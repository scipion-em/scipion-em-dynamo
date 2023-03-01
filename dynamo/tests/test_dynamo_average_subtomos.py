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
from dynamo.tests.test_dynamo_sta_base import TestDynamoStaBase
from pyworkflow.tests import setupTestProject
from pyworkflow.utils import magentaStr


class TestDynamoAverageSubtomograms(TestDynamoStaBase):

    subtomosExtracted = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.subtomosExtracted = cls.runPreviousProtocols()

    @classmethod
    def runPreviousProtocols(cls):
        tomoImported = super().runImportTomograms()  # Import tomograms
        tomosBinned = super().runBinTomograms(tomoImported)  # Bin the tomogram to make it smaller
        coordsImported = super().runImport3dCoords(tomosBinned)  # Import the coordinates from the binned tomogram
        # Extract subtomograms
        return super().runExtractSubtomograms(coordsImported,
                                              tomoSource=SAME_AS_PICKING,
                                              boxSize=super().bin2BoxSize)

    @classmethod
    def runAverageSubtomograms(cls):
        print(magentaStr("\n==> Averaging the subtomograms:"))
        protAvgSubtomo = cls.newProtocol(DynamoProtAvgSubtomograms, inSubtomos=cls.subtomosExtracted)
        cls.launchProtocol(protAvgSubtomo)
        return getattr(protAvgSubtomo, protAvgSubtomo._possibleOutputs.average.name, None)

    def test_average(self):
        avg = self.runAverageSubtomograms()  # The imported coordinates correspond to a binned 2 tomogram
        super().checkAverage(avg,
                             expectedSRate=super().bin2SRate,
                             expectedBoxSize=super().bin2BoxSize,
                             hasHalves=False)  # Dynamo average protocol doesn't generate halves
