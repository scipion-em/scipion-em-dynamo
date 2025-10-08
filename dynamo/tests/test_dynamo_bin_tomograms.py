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
import numpy as np
from dynamo.tests.test_dynamo_base import TestDynamoStaBase
from pyworkflow.tests import DataSet
from tomo.tests import RE4_STA_TUTO, DataSetRe4STATuto, TS_03, TS_54


class TestDynamoBinTomograms(TestDynamoStaBase):
    tomoFiles = None
    binningFactor = 2
    origSize = np.array([1024, 1024, 512])
    bin4sRate = float(DataSetRe4STATuto.sRateBin4.value)

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.ds = DataSet.getDataSet(RE4_STA_TUTO)
        cls.bin8SRate = cls.bin4sRate * cls.binningFactor
        cls.testAcqDict, cls.testDimsDict = DataSetRe4STATuto.genTestTomoDicts(
            tsIdList=(TS_03, TS_54),
            binning=cls.binningFactor)
        cls.importedTomos = super().runImportTomograms(
            tomoFiles=cls.ds.getFile(DataSetRe4STATuto.tomosPath.value),
            filesPattern="*.mrc",
            sRate=cls.bin4sRate,
            exclusionWords=str(DataSetRe4STATuto.exclusionWordsTs03ts54.value))

    def testBinTomograms(self):
        binnedTomograms = super().runBinTomograms(self.importedTomos, binning=self.binningFactor)
        # Check the results
        self.checkTomograms(binnedTomograms,
                            expectedSetSize=2,
                            expectedSRate=self.bin8SRate,
                            expectedDimensions=self.testDimsDict,
                            testAcqObj=self.testAcqDict)
