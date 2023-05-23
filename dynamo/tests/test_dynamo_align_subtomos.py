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
from xmipp3.constants import MASK3D_CYLINDER
from xmipp3.protocols import XmippProtCreateMask3D
from xmipp3.protocols.protocol_preprocess.protocol_create_mask3d import SOURCE_GEOMETRY

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
                         areaSearchMode=str(FROM_PREVIOUS_ESTIMATION), alignMask=None, protLabel=None):
        argsDict = {'inputVolumes': cls.subtomosExtracted,
                    'templateRef': cls.avg,
                    'numberOfIters': nIters,
                    'dim': dims,
                    'thresholdMode': thMode,
                    'limm': areaSearchMode,
                    'useGpu': True}
        if alignMask:
            argsDict['alignMask'] = alignMask
        protAlign = cls.newProtocol(DynamoSubTomoMRA, **argsDict)
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

    def test_alignSubtomos_oneRound_With_AlignMask(self):
        print(magentaStr("\n==> aligning the subtomograms with alignment mask, 1 round:"))
        # Generate the mask
        protMask3D = self.newProtocol(XmippProtCreateMask3D,
                                      source=SOURCE_GEOMETRY,
                                      samplingRate=self.bin2SRate,
                                      size=self.bin2BoxSize,
                                      geo=MASK3D_CYLINDER,
                                      radius=15,
                                      shiftCenter=True,
                                      centerZ=6,
                                      height=20,
                                      doSmooth=True)
        self.launchProtocol(protMask3D)
        alignMask = getattr(protMask3D, 'outputMask', None)
        # Align the subtomograms
        subtomos, avg = self.runAlignSubtomos(protLabel='Subtomo align, 1 round', alignMask=alignMask)
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
