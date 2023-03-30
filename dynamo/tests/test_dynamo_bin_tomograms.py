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


class TestDynamoBinTomograms(TestDynamoStaBase):
    tomoFiles = None
    binningFactor = 2
    origSize = np.array([1024, 1024, 512])
    sRate = 5
    bin2SRate = 10
    importedTomos = None

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        ds = DataSet.getDataSet('tomo-em')
        cls.tomoFiles = ds.getFile('tomo1')
        cls.importedTomos = super().runImportTomograms(tomoFiles=cls.tomoFiles, sRate=cls.sRate)

    def testBinTomograms(self):
        binnedTomograms = super().runBinTomograms(self.importedTomos, binning=self.binningFactor)
        # Check the results
        self.assertSetSize(binnedTomograms, 1)
        self.assertEqual(binnedTomograms.getSamplingRate(), self.bin2SRate)
        self.assertTrue(np.array_equal(np.array(binnedTomograms.getDimensions()), self.origSize / self.binningFactor))