# **************************************************************************
# *
# * Authors:    Scipion Team (scipion@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from dynamo.protocols.protocol_extraction import OTHER
from dynamo.tests.test_dynamo_sta_base import TestDynamoStaBase
from pyworkflow.tests import DataSet
from pyworkflow.utils import magentaStr
from tomo.constants import TR_DYNAMO
from tomo.protocols import ProtTomoExtractCoords
from tomo.protocols.protocol_extract_coordinates import Output3dCoordExtraction
from tomo.tests import EMD_10439, DataSetEmd10439


class TestDynamoSubtomoExtraction(TestDynamoStaBase):

    ds = None
    tomoImported = None
    coordsImported = None
    tomosBinned = None
    subtomosSameAsPicking = None
    subtomosAnotherTomo = None
    bin2BoxSize = None
    unbinnedBoxSize = None
    bin2SRate = None
    nParticles = None
    unbinnedSRate = None

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.ds = DataSet.getDataSet(EMD_10439)
        cls.bin2BoxSize = DataSetEmd10439.bin2BoxSize.value
        cls.unbinnedBoxSize = DataSetEmd10439.unbinnedBoxSize.value
        cls.bin2SRate = DataSetEmd10439.bin2SRate.value
        cls.nParticles = DataSetEmd10439.nParticles.value
        cls.unbinnedSRate = DataSetEmd10439.unbinnedSRate.value
        cls.runPreviousProtocols()

    @classmethod
    def runPreviousProtocols(cls):
        # Import the tomogram
        cls.tomoImported = super().runImportTomograms(tomoFiles=cls.ds.getFile(DataSetEmd10439.tomoEmd10439.name),
                                                      sRate=DataSetEmd10439.unbinnedSRate.value)
        # Bin the tomogram to make it smaller
        cls.tomosBinned = super().runBinTomograms(inTomos=cls.tomoImported,
                                                  binning=DataSetEmd10439.binFactor.value)
        # Import the coordinates from the binned tomogram
        cls.coordsImported = super().runImport3dCoords(sqliteFile=cls.ds.getFile(DataSetEmd10439.coords39Sqlite.name),
                                                       inTomos=cls.tomosBinned,
                                                       boxSize=cls.bin2BoxSize)
        cls.subtomosSameAsPicking = super().runExtractSubtomograms(cls.coordsImported,
                                                                   boxSize=cls.bin2BoxSize)
        cls.subtomosAnotherTomo = super().runExtractSubtomograms(cls.coordsImported,
                                                                 tomoSource=OTHER,
                                                                 tomograms=cls.tomoImported,
                                                                 boxSize=cls.unbinnedBoxSize)

    @classmethod
    def runExtract3dCoords(cls, inputSubTomos=None, inputTomos=None, boxSize=None):
        print(magentaStr("\n==> Extracting the 3D coordinates:"))
        protExtract3dCoords = cls.newProtocol(ProtTomoExtractCoords,
                                              inputSubTomos=inputSubTomos,
                                              inputTomos=inputTomos,
                                              boxSize=boxSize)
        cls.launchProtocol(protExtract3dCoords)
        coordsExtracted = getattr(protExtract3dCoords, Output3dCoordExtraction.coordinates3d.name, None)
        cls.assertIsNotNone(coordsExtracted, "There was a problem with the 3d coordinates extraction")
        return coordsExtracted

    def test_extractParticlesSameAsPicking(self):
        # The imported 3d coordinates were picked from the binned tomogram
        super().checkExtractedSubtomos(self.coordsImported,
                                       self.subtomosSameAsPicking,
                                       expectedSetSize=self.nParticles,
                                       expectedSRate=self.bin2SRate,
                                       expectedBoxSize=self.bin2BoxSize,
                                       convention=TR_DYNAMO,
                                       orientedParticles=True)  # Picked with PySeg

    def test_extractParticlesOtherTomoSource(self):
        super().checkExtractedSubtomos(self.coordsImported,
                                       self.subtomosAnotherTomo,
                                       expectedSetSize=self.nParticles,
                                       expectedSRate=self.unbinnedSRate,
                                       expectedBoxSize=self.unbinnedBoxSize,
                                       convention=TR_DYNAMO,
                                       orientedParticles=True)  # Picked with PySeg

    # __________________________________________________________________________________________________________________
    # NOTE:
    # Although the coordinates extraction is not a part of the plugin emantomo, a part of its functionality
    # will be tested here as it has a very direct relation with the particle extraction. In the coordinates
    # extraction a new set of coordinates is generated, so in this case not only the shifts but the coordinates
    # are expected to be scaled properly according to the sampling rate of the introduced tomograms.
    # __________________________________________________________________________________________________________________

    def test_extract3dCoordsToBiggerTomo(self):
        """Subtomos extracted from the same tomo used for the picking, which was at bin 2. Coordinates will
        be extracted to the original size (unbinned)."""
        boxSize = self.unbinnedBoxSize
        extractedCoords = self.runExtract3dCoords(inputSubTomos=self.subtomosSameAsPicking,
                                                  inputTomos=self.tomoImported,
                                                  boxSize=boxSize)
        super().checkExtracted3dCoordinates(self.subtomosSameAsPicking,
                                            extractedCoords,
                                            expectedSetSize=self.nParticles,
                                            expectedBoxSize=boxSize,
                                            expectedSRate=self.unbinnedSRate,
                                            convention=TR_DYNAMO,
                                            orientedParticles=True)  # Picked with PySeg

    def test_extract3dCoordsToSmallerTomo(self):
        """Subtomos extracted from the another tomo source, which was unbinned. Coordinates will
        be extracted to bin 2."""
        boxSize = self.bin2BoxSize
        extractedCoords = self.runExtract3dCoords(inputSubTomos=self.subtomosAnotherTomo,
                                                  inputTomos=self.tomosBinned,
                                                  boxSize=boxSize)
        self.checkExtracted3dCoordinates(self.subtomosAnotherTomo,
                                         extractedCoords,
                                         expectedSetSize=self.nParticles,
                                         expectedBoxSize=boxSize,
                                         expectedSRate=self.bin2SRate,
                                         convention=TR_DYNAMO,
                                         orientedParticles=True)  # Picked with PySeg

    def test_extract3dCoordsToTheSameTomo(self):
        """Subtomos extracted from the same tomo used for the picking, which was at bin 2. Coordinates will
        be extracted to the same tomogram."""
        boxSize = self.bin2BoxSize
        extractedCoords = self.runExtract3dCoords(inputSubTomos=self.subtomosSameAsPicking,
                                                  inputTomos=self.tomosBinned,
                                                  boxSize=boxSize)
        super().checkExtracted3dCoordinates(self.subtomosSameAsPicking,
                                            extractedCoords,
                                            expectedSetSize=self.nParticles,
                                            expectedBoxSize=boxSize,
                                            expectedSRate=self.bin2SRate,
                                            convention=TR_DYNAMO,
                                            orientedParticles=True)  # Picked with PySeg
