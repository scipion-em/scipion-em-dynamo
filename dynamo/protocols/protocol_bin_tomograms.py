# **************************************************************************
# *
# * Authors:    David Herreros Calero (dherreros@cnb.csic.es)
# *             Scipion Team (scipion@cnb.csic.es)
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
from enum import Enum
from os.path import abspath

from dynamo.protocols.protocol_base_dynamo import DynamoProtocolBase
from pwem.emlib.image import ImageHandler
from pyworkflow.object import Set
from pyworkflow.protocol import params, GT, STEPS_PARALLEL
from pyworkflow.utils import removeBaseExt, getExt
from tomo.objects import Tomogram, SetOfTomograms
from dynamo import Plugin


class DynamoBinOuts(Enum):
    tomograms = SetOfTomograms


class DynamoBinTomograms(DynamoProtocolBase):
    """Reduce the size of a SetOfTomograms by a binning factor"""

    _label = 'bin tomograms'
    _possibleOutputs = DynamoBinOuts
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.finalTomoNamesDict = {}
        self.ih = None
        self.sRate = None
        self.doConvertFiles = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputTomos', params.PointerParam,
                      pointerClass='SetOfTomograms',
                      label="Input Tomograms",
                      important=True)
        form.addParam('binning', params.IntParam,
                      default=2,
                      validators=[GT(0)],
                      label="Binning Factor",
                      help="A Binning Factor of 1 means that no binning will be carried out.")
        form.addParam('zChunk', params.IntParam,
                      default=300,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Number of slices kept in memory",
                      help="Maximum number of Z slices that are kept simultaneously in the memory during the "
                           "binning process. This parameter might be important for larger size tomograms, making "
                           "possible to process them in vertical slabs of thickness = value introduced in the  "
                           "current parameter. This procedure can be accelerated using the multiple threads to engage "
                           "several cores in parallel. However, this will only make sense if the total memory occupied "
                           "by all the slabs simultaneously in memory in a given time fits in the RAM of the machine.")
        form.addParam('binThreads', params.IntParam,
                      label='Dynamo threads',
                      default=3,
                      help='Number of threads used by Dynamo each time it is called in the protocol execution. For '
                           'example, if 2 Scipion threads and 3 Dynamo threads are set, the tomograms will be '
                           'processed in groups of 2 at the same time with a call of tomo3d with 3 threads each, so '
                           '6 threads will be used at the same time. Beware the memory of your machine has '
                           'memory enough to load together the number of tomograms specified by Scipion threads.')
        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        stepIds = []
        for tsId in self.tomoDict.keys():
            cInPid = self._insertFunctionStep(self.convertInputStep, tsId,
                                              prerequisites=[],
                                              needsGPU=False)
            binId = self._insertFunctionStep(self.binTomosStep, tsId,
                                             prerequisites=cInPid,
                                             needsGPU=False)
            cOutId = self._insertFunctionStep(self.createOutputStep, tsId,
                                              prerequisites=binId,
                                              needsGPU=False)
            stepIds.append(cOutId)
        self._insertFunctionStep(self._closeOutputSet,
                                 prerequisites=stepIds,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.ih = ImageHandler()
        self.tomoDict = {tomo.getTsId(): tomo.clone() for tomo in self.inputTomos.get()}
        self.doConvertFiles = not self.isCompatibleFileFormat()
        self.sRate = self.inputTomos.get().getSamplingRate() * self.getBinningFactor(fromDynamo=False)

    def convertInputStep(self, tsId: str):
        if self.doConvertFiles:
            tomo = self.tomoDict[tsId]
            origName = tomo.getFileName()
            finalName = self.getConvertedOrLinkedTsFn(tsId)
            self.ih.convert(origName, finalName)

    def binTomosStep(self, tsId: str):
        mFile = self.createMCodeFile(tsId)
        args = ' %s' % mFile
        self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    def createOutputStep(self, tsId: str):
        with self._lock:
            outTomos = self.getOutputSetOfTomograms()
            inTomo = self.tomoDict[tsId]
            tomo = Tomogram()
            tomo.copyInfo(inTomo)
            tomo.setSamplingRate(self.sRate)
            tomo.setFileName(self.getOutTsFn(tsId))
            outTomos.append(tomo)
            outTomos.update(tomo)
            outTomos.write()
            self._store(outTomos)

    # --------------------------- DEFINE utils functions ----------------------
    def isCompatibleFileFormat(self):
        """Compatible with MRC and em (MRC with that extension)"""
        compatibleExts = ['.em', '.mrc']
        return True if (getExt(self.inputTomos.get().getFirstItem().getFileName())
                        not in compatibleExts) else False

    def getConvertedOrLinkedTsFn(self, tsId: str):
        return self._getExtraPath(f'in_{tsId}.mrc')

    def getOutTsFn(self, tsId: str):
        return self._getExtraPath(f'{tsId}.mrc')

    def createMCodeFile(self, tsId: str):
        # FROM DYNAMO:
        # ______________________________________________________________________________________________
        # bin(fileIn,fileOut,binFactor,varargin)
        # Create a Volume template object
        # p.addParamValue('slabSize',[],'short','ss');
        # p.addParamValue('matlabWorkers',0,'short','mw');
        # p.addParamValue('maximumMegaBytes',[]);
        # p.addParamValue('showStatistics',false,'short','sst');
        # ______________________________________________________________________________________________
        origName = abspath(self.getConvertedOrLinkedTsFn(tsId))
        finalName = abspath(self.getOutTsFn(tsId))
        mFile = self._getExtraPath('binTomograms.m')
        with open(mFile, 'w') as codeFile:
            content = ("dpktomo.tools.bin('%s', '%s', %i, 'slabSize', %i, 'matlabWorkers', %i, "
                       "'showStatistics', true)\n") % (origName, finalName, super().getBinningFactor(),
                                                       self.zChunk.get(), self.binThreads.get())
            codeFile.write(content)
        return mFile

    def getOutputSetOfTomograms(self) -> SetOfTomograms:
        outTomograms = getattr(self, self._possibleOutputs.tomograms.name, None)
        if outTomograms:
            outTomograms.enableAppend()
            tomograms = outTomograms
        else:
            tomograms = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
            tomograms.copyInfo(self.inputTomos.get())
            tomograms.setSamplingRate(self.sRate)
            tomograms.setStreamState(Set.STREAM_OPEN)
            setattr(self, self._possibleOutputs.tomograms.name, tomograms)
            self._defineOutputs(**{self._possibleOutputs.tomograms.name: tomograms})
            self._defineSourceRelation(self.inputTomos, tomograms)

        return tomograms


    # --------------------------- DEFINE info functions ----------------------
    def _methods(self):
        methodsMsgs = ["*Binning Factor*: %s" % self.binning.get()]
        return methodsMsgs

    def _summary(self):
        summary = []
        if self.getOutputsSize() >= 1:
            for _, outTomos in self.iterOutputAttributes():
                summary.append("Output *%s*:" % outTomos.getNameId().split('.')[1])
                summary.append("    * Binning Factor: *%s*" % self.binning.get())
                summary.append("    * Number of Tomograms Binned: *%s*" %
                               outTomos.getSize())
        else:
            summary.append("Output tomograms not ready yet.")
        return summary
