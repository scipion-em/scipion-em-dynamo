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
from dynamo.protocols.protocol_base_dynamo import DynamoProtocolBase
from pwem.emlib.image import ImageHandler
from pyworkflow.protocol import params, GT
from pyworkflow.utils import removeBaseExt, getExt
from tomo.objects import Tomogram, SetOfTomograms
from dynamo import Plugin


class DynamoBinOuts(Enum):
    tomograms = SetOfTomograms


class DynamoBinTomograms(DynamoProtocolBase):
    """Reduce the size of a SetOfTomograms by a binning factor"""

    _label = 'bin tomograms'
    _possibleOutputs = DynamoBinOuts

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.content = ''
        self.finalTomoNamesDict = {}

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputTomos', params.PointerParam,
                      pointerClass='SetOfTomograms',
                      label="Tomograms",
                      important=True)
        self._addBinningParam(form)
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
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        for tomo in self.inputTomos:
            origTomoName = tomo.getFileName()
            finalTomoName = self.getFinalTomoName(tomo)
            self.finalTomoNamesDict[finalTomoName] = tomo
            if self.doConvertFiles:
                self._insertFunctionStep(self.convertInputStep, origTomoName, finalTomoName)
            # Generate one unique file with all the tomograms to be binned
            self._insertFunctionStep(self.createMCodeStep, origTomoName, finalTomoName)
        # That way, we can carry out the binning of all the tomograms provided with only one call to MATLAB, improving
        # the efficiency, as this call is very slow
        self._insertFunctionStep(self.binTomosStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.inputTomos = self.inputTomos.get()
        self.doConvertFiles = not self.isCompatibleFileFormat()

    @staticmethod
    def convertInputStep(origName, finalName):
        ih = ImageHandler()
        ih.convert(origName, finalName)

    def createMCodeStep(self, origTomoName, finalTomoName):
        # FROM DYNAMO:
        # ______________________________________________________________________________________________
        # bin(fileIn,fileOut,binFactor,varargin)
        # Create a Volume template object
        # p.addParamValue('slabSize',[],'short','ss');
        # p.addParamValue('matlabWorkers',0,'short','mw');
        # p.addParamValue('maximumMegaBytes',[]);
        # p.addParamValue('showStatistics',false,'short','sst');
        # ______________________________________________________________________________________________
        self.content += "dpktomo.tools.bin('%s', '%s', %i, 'slabSize', %i, 'matlabWorkers', %i, " \
                        "'showStatistics', true)\n" % \
                        (origTomoName, finalTomoName, super().getBinningFactor(), self.zChunk.get(),
                         self.numberOfThreads.get())

    def binTomosStep(self):
        mFile = self._getExtraPath('binTomograms.m')
        with open(mFile, 'w') as codeFile:
            codeFile.write(self.content)
        args = ' %s' % mFile
        self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    def createOutputStep(self):
        sr = self.inputTomos.getSamplingRate() * super().getBinningFactor(forDynamo=False)  # Not for, but from
        outTomos = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
        outTomos.setSamplingRate(sr)
        for tomoName, inTomo in self.finalTomoNamesDict.items():
            tomo = Tomogram()
            tomo.copyInfo(inTomo)
            tomo.setSamplingRate(sr)
            tomo.setFileName(tomoName)
            outTomos.append(tomo)

        self._defineOutputs(**{DynamoBinOuts.tomograms.name: outTomos})
        self._defineSourceRelation(self.inputTomos, outTomos)

    # --------------------------- DEFINE utils functions ----------------------
    def getConvertedOutFileName(self, inFileName):
        return self._getExtraPath(removeBaseExt(inFileName) + '.mrc')

    def getFinalTomoName(self, tomo):
        return self.getConvertedOutFileName(tomo.getFileName())

    def isCompatibleFileFormat(self):
        """Compatible with MRC and em (MRC with that extension)"""
        compatibleExts = ['.em', '.mrc']
        return True if getExt(self.inputTomos.getFirstItem().getFileName()) not in compatibleExts else False

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
