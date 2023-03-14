# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *              Scipion Team (scipion@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
from dynamo.protocols import DynamoImportSubtomos
from dynamo.protocols.protocol_import_subtomos import DynImportSubtomosOuts
from dynamo.tests.test_dynamo_base import TestDynamoStaBase
from pyworkflow.tests import DataSet
from tomo.objects import SetOfSubTomograms
from tomo.tests import EMD_10439, DataSetEmd10439


class TestDynamoImportSubTomograms(TestDynamoStaBase):
    """ This class check if the protocol to import subtomograms from Dynamo works
     properly."""

    bin2SRate = None
    bin2BoxSize = None
    tomosBinned = None
    extractProtocol = None
    ds = None

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
        tomoImported = super().runImportTomograms(tomoFiles=cls.ds.getFile(DataSetEmd10439.tomoEmd10439.name),
                                                      sRate=DataSetEmd10439.unbinnedSRate.value)
        # Bin the tomogram to make it smaller
        cls.tomosBinned = super().runBinTomograms(inTomos=tomoImported,
                                                  binning=DataSetEmd10439.binFactor.value)
        # Import the coordinates from the binned tomogram
        coordsImported = super().runImport3dCoordsSqlite(
            sqliteFile=cls.ds.getFile(DataSetEmd10439.coords39Sqlite.name),
            inTomos=cls.tomosBinned,
            boxSize=cls.bin2BoxSize)
        cls.extractProtocol = super().runExtractSubtomograms(coordsImported,
                                                             boxSize=cls.bin2BoxSize,
                                                             returnProtocol=True)

    @classmethod
    def _runImportDynSubTomograms(cls, tomos=None, sRate=None) -> SetOfSubTomograms:
        tblFile = cls.extractProtocol._getExtraPath('Crop1', 'crop.tbl')
        protImportSubtomos = cls.newProtocol(DynamoImportSubtomos,
                                             filesPath=tblFile,
                                             tomoSet=tomos,
                                             samplingRate=sRate)
        cls.launchProtocol(protImportSubtomos)
        outSubtomos = getattr(protImportSubtomos, DynImportSubtomosOuts.subtomograms.name, None)
        cls.assertIsNotNone(outSubtomos, 'There was a problem importing the subtomograms')
        return outSubtomos

    def test_import_without_associated_tomos(self):
        sRate = self.bin2SRate
        outsubtomos = self._runImportDynSubTomograms(sRate=sRate)
        super().checkSubtomograms(outsubtomos,
                                  expectedSetSize=self.nParticles,
                                  expectedSRate=sRate,
                                  expectedBoxSize=self.bin2BoxSize)

    def test_import_with_associated_tomos(self):
        outsubtomos = self._runImportDynSubTomograms(tomos=self.tomosBinned)
        super().checkSubtomograms(outsubtomos,
                                  expectedSetSize=self.nParticles,
                                  expectedSRate=self.bin2SRate,
                                  expectedBoxSize=self.bin2BoxSize,
                                  tomograms=self.tomosBinned)