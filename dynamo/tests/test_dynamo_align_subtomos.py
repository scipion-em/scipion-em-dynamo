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
from dynamo.protocols import DynamoSubTomoMRA
from dynamo.protocols.protocol_extraction import SAME_AS_PICKING
from dynamo.protocols.protocol_subtomo_MRA import FROM_PREVIOUS_ESTIMATION, NO_THRESHOLD
from dynamo.tests.test_dynamo_base import TestDynamoStaBase
from pyworkflow.tests import DataSet
from pyworkflow.utils import magentaStr
from tomo.constants import TR_DYNAMO
from tomo.tests import EMD_10439, DataSetEmd10439


class TestDynamoAlignSubtomograms(TestDynamoStaBase):
    ds = None
    tomoFiles = None
    sqliteFile = None
    binFactor = None
    bin2BoxSize = None
    bin2SRate = None
    nParticles = None
    subtomosExtracted = None
    avg = None

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.ds = DataSet.getDataSet(EMD_10439)
        cls.bin2BoxSize = DataSetEmd10439.bin2BoxSize.value
        cls.bin2SRate = DataSetEmd10439.bin2SRate.value
        cls.nParticles = DataSetEmd10439.nParticles.value
        cls.subtomosExtracted, cls.avg = cls.runPreviousProtocols()

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
        subtomos = super().runExtractSubtomograms(inCoords=coordsImported,
                                                  tomoSource=SAME_AS_PICKING,
                                                  boxSize=cls.bin2BoxSize)
        # Average the extracted subtomograms
        avg = super().runAverageSubtomograms(subtomos)

        return subtomos, avg

    @classmethod
    def runAlignSubtomos(cls, nIters='3', dims='0', thMode=str(NO_THRESHOLD),
                         areaSearchMode=str(FROM_PREVIOUS_ESTIMATION), protLabel=None):
        protAlign = cls.newProtocol(DynamoSubTomoMRA,
                                    inputVolumes=cls.subtomosExtracted,
                                    templateRef=cls.avg,
                                    numberOfIters=nIters,
                                    dim=dims,
                                    thresholdMode=thMode,
                                    limm=areaSearchMode,
                                    useGpu=True)
        if protLabel:
            protAlign.setObjLabel(protLabel)
        cls.launchProtocol(protAlign)
        subtomos = getattr(protAlign, protAlign._possibleOutputs.subtomograms.name, None)
        avg = getattr(protAlign, protAlign._possibleOutputs.average.name, None)
        return subtomos, avg

    def test_alignSubtomos_oneRound(self):
        print(magentaStr("\n==> aligning the subtomograms, 1 round:"))
        subtomos, avg = self.runAlignSubtomos(protLabel='Subtomo align, 1 round')
        self.checkResults(avg, subtomos)

    def test_alignSubtomos_twoRounds(self):
        print(magentaStr("\n==> aligning the subtomograms, 2 rounds:"))
        subtomos, avg = self.runAlignSubtomos(nIters='2 2', dims='0', protLabel='Subtomo align, 2 rounds, dim=0')
        self.checkResults(avg, subtomos)

    def test_alignSubtomos_threeRounds(self):
        print(magentaStr("\n==> aligning the subtomograms, 3 rounds:"))
        subtomos, avg = self.runAlignSubtomos(nIters='1 1 1', dims='22 44',
                                              protLabel='Subtomo align, 3 rounds, dimSize=2')
        self.checkResults(avg, subtomos)

    def test_alignSubtomos_3Rounds_3Dims(self):
        print(magentaStr("\n==> aligning the subtomograms, 3 rounds, 3 dims:"))
        subtomos, avg = self.runAlignSubtomos(nIters='1 1 1', dims='22 44 44',
                                              protLabel='Subtomo align, 3 rounds, dimSize=3')
        self.checkResults(avg, subtomos)

    def checkResults(self, avg, subtomos):
        # Check the average
        super().checkAverage(avg,
                             expectedSRate=self.bin2SRate,
                             expectedBoxSize=self.bin2BoxSize,
                             hasHalves=False)  # Dynamo average protocol doesn't generate halves
        # Check the subtomograms
        super().checkRefinedSubtomograms(self.subtomosExtracted, subtomos,
                                         expectedSetSize=self.nParticles,
                                         expectedBoxSize=self.bin2BoxSize,
                                         expectedSRate=self.bin2SRate,
                                         convention=TR_DYNAMO,
                                         orientedParticles=True)  # Picked with PySeg
