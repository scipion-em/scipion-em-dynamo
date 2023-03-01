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
from dynamo.protocols import DynamoBinTomograms
from dynamo.protocols.protocol_extraction import SAME_AS_PICKING, OTHER, DynamoExtraction
from pyworkflow.tests import DataSet, setupTestProject
from pyworkflow.utils import magentaStr
from tomo.protocols import ProtImportTomograms, ProtImportCoordinates3DFromScipion
from tomo.protocols.protocol_import_coordinates_from_scipion import outputObjs
from tomo.tests import EMD_10439, DataSetEmd10439
from tomo.tests.test_base_coordinates_and_subtomos import TestUtilsAverageOfSubtomos, TestUtilsExtractSubtomos, \
    TestUtilsExtractCoords


class TestDynamoStaBase(TestUtilsExtractSubtomos, TestUtilsExtractCoords, TestUtilsAverageOfSubtomos):

    ds = DataSet.getDataSet(EMD_10439)
    nParticles = 39
    bin2BoxSize = 44
    unbinnedBoxSize = 88
    unbinnedSRate = 13.68
    bin2SRate = 27.36

    @classmethod
    def runImportTomograms(cls):
        # Import tomograms
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=cls.ds.getFile(DataSetEmd10439.tomoEmd10439.name),
                                             samplingRate=cls.unbinnedSRate, )

        cls.launchProtocol(protImportTomogram)
        tomoImported = protImportTomogram.outputTomograms
        cls.assertIsNotNone(tomoImported, "There was a problem with tomogram output")
        return tomoImported

    @classmethod
    def runImport3dCoords(cls, tomoImported):
        # Import coordinates
        print(magentaStr("\n==> Importing the 3D coordinates:"))
        protImportCoordinates3d = cls.newProtocol(ProtImportCoordinates3DFromScipion,
                                                  sqliteFile=cls.ds.getFile(DataSetEmd10439.coords39Sqlite.name),
                                                  importTomograms=tomoImported,
                                                  boxSize=cls.bin2BoxSize)

        cls.launchProtocol(protImportCoordinates3d)
        coordsImported = getattr(protImportCoordinates3d, outputObjs.coordinates.name, None)
        cls.assertIsNotNone(coordsImported, "There was a problem with the 3D coordinates output")
        return coordsImported

    @classmethod
    def runBinTomograms(cls, tomoImported, binning=2):
        # Bin the tomogram to make it smaller
        print(magentaStr("\n==> Tomogram binning:"))
        protBinTomos = cls.newProtocol(DynamoBinTomograms,
                                       inputTomos=tomoImported,
                                       binning=binning)

        cls.launchProtocol(protBinTomos)
        tomosBinned = getattr(protBinTomos, protBinTomos._possibleOutputs.tomograms.name, None)
        cls.assertIsNotNone(tomosBinned, 'No tomograms were binned.')
        return tomosBinned

    @classmethod
    def runExtractSubtomograms(cls, coordsImported, tomoSource=SAME_AS_PICKING, tomograms=None, boxSize=None):
        # Extract subtomograms
        print(magentaStr("\n==> Extracting the subtomograms:"))
        protLabel = 'Extraction - same as picking'
        argsDict = {'inputCoordinates': coordsImported,
                    'tomoSource': tomoSource,
                    'boxSize': boxSize,
                    'doInvert': True}
        if tomoSource != SAME_AS_PICKING:
            argsDict['tomoSource'] = OTHER
            argsDict['inputTomograms'] = tomograms
            protLabel = 'Extraction - another tomo source'

        protTomoExtraction = cls.newProtocol(DynamoExtraction, **argsDict)
        protTomoExtraction.setObjLabel(protLabel)
        cls.launchProtocol(protTomoExtraction)
        subtomosExtracted = getattr(protTomoExtraction, protTomoExtraction._possibleOutputs.subtomograms.name, None)
        cls.assertIsNotNone(subtomosExtracted, "There was a problem with the subtomogram extraction")
        return subtomosExtracted




