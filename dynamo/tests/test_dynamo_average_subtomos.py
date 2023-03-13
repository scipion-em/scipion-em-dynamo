# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team  (scipion@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnología (CSIC), Madrid, Spain
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
from dynamo.protocols import DynamoProtAvgSubtomograms
from dynamo.protocols.protocol_extraction import SAME_AS_PICKING
from dynamo.tests.test_dynamo_base import TestDynamoStaBase
from pyworkflow.tests import DataSet
from pyworkflow.utils import magentaStr
from tomo.tests import EMD_10439, DataSetEmd10439


class TestDynamoAverageSubtomograms(TestDynamoStaBase):

    ds = None
    tomoFiles = None
    sqliteFile = None
    binFactor = None
    bin2BoxSize = None
    subtomosExtracted = None

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.ds = DataSet.getDataSet(EMD_10439)
        cls.bin2BoxSize = DataSetEmd10439.bin2BoxSize.value
        cls.subtomosExtracted = cls.runPreviousProtocols()

    @classmethod
    def runPreviousProtocols(cls):
        # Import the tomogram
        tomoImported = super().runImportTomograms(tomoFiles=cls.ds.getFile(DataSetEmd10439.tomoEmd10439.name),
                                                  sRate=DataSetEmd10439.unbinnedSRate.value)
        # Bin the tomogram to make it smaller
        tomosBinned = super().runBinTomograms(inTomos=tomoImported,
                                              binning=DataSetEmd10439.binFactor.value)
        # Import the coordinates from the binned tomogram
        coordsImported = super().runImport3dCoordsSqlite(sqliteFile=cls.ds.getFile(DataSetEmd10439.coords39Sqlite.name),
                                                         inTomos=tomosBinned,
                                                         boxSize=cls.bin2BoxSize)
        # Extract subtomograms
        return super().runExtractSubtomograms(inCoords=coordsImported,
                                              tomoSource=SAME_AS_PICKING,
                                              boxSize=cls.bin2BoxSize)

    @classmethod
    def runAverageSubtomograms(cls):
        print(magentaStr("\n==> Averaging the subtomograms:"))
        protAvgSubtomo = cls.newProtocol(DynamoProtAvgSubtomograms, inSubtomos=cls.subtomosExtracted)
        cls.launchProtocol(protAvgSubtomo)
        return getattr(protAvgSubtomo, protAvgSubtomo._possibleOutputs.average.name, None)

    def test_average(self):
        avg = self.runAverageSubtomograms()  # The imported coordinates correspond to a binned 2 tomogram
        super().checkAverage(avg,
                             expectedSRate=DataSetEmd10439.bin2SRate.value,
                             expectedBoxSize=self.bin2BoxSize,
                             hasHalves=False)  # Dynamo average protocol doesn't generate halves
