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
from pwem.emlib.image import ImageHandler
from pyworkflow import BETA
from pyworkflow.protocol import params, GT
from pwem.protocols import EMProtocol
from pyworkflow.utils import getExt, removeBaseExt

from tomo.protocols import ProtTomoBase
from tomo.objects import Tomogram, SetOfTomograms

from dynamo import Plugin


class DynBinTomosOutputs(Enum):
    tomograms = SetOfTomograms


class DynamoBinTomograms(EMProtocol, ProtTomoBase):
    """Reduce the size of a SetOfTomograms by a binning factor"""

    _label = 'bin tomograms'
    _devStatus = BETA
    _possibleOutputs = DynBinTomosOutputs
    inTomosDir = 'inTomograms'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.content = ''
        self.finalTomoNamesDict = {}

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
                           "binning process. This parameter might be important for larger size tomograms.")
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
        self.content += "dpktomo.tools.bin('%s', '%s', %i, 'ss', %i, 'mw', %i, 'sst', true)\n" % \
                        (origTomoName, finalTomoName, self.getBinningFactor(), self.zChunk.get(),
                         self.numberOfThreads.get())

    def binTomosStep(self):
        mFile = self._getExtraPath('binTomograms.m')
        with open(mFile, 'w') as codeFile:
            codeFile.write(self.content)
        args = ' %s' % mFile
        self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    def createOutputStep(self):
        sr = self.inputTomos.getSamplingRate() * self.getBinningFactor(forDynamo=False)  # Not for, but from
        outTomos = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
        outTomos.setSamplingRate(sr)
        for tomoName, inTomo in self.finalTomoNamesDict.items():
            tomo = Tomogram()
            tomo.copyInfo(inTomo)
            tomo.setSamplingRate(sr)
            tomo.setFileName(tomoName)
            outTomos.append(tomo)

        self._defineOutputs(**{DynBinTomosOutputs.tomograms.name: outTomos})
        self._defineSourceRelation(self.inputTomos, outTomos)

    # --------------------------- DEFINE utils functions ----------------------
    def isCompatibleFileFormat(self):
        """Compatible with MRC and em"""
        compatibleExts = ['.em', '.mrc']
        return True if getExt(self.inputTomos.getFirstItem().getFileName()) not in compatibleExts else False

    def getConvertedOutFileName(self, inFileName):
        return self._getExtraPath(removeBaseExt(inFileName) + '.mrc')

    def getBinningFactor(self, forDynamo=True):
        """From Dynamo: a binning Factor of 1 will decrease the size of the Tomograms by 2,
        a Binning Factor of 2 by 4... So Dynamo interprets the binning factor as 2**binFactor, while IMOD
        interprets it literally. Thus, this method will convert the binning introduced by the user in the
        Dynamo convention"""
        return (2**(self.binning.get() - 1)) / 2 if forDynamo else self.binning.get()

    def getFinalTomoName(self, tomo):
        return self.getConvertedOutFileName(tomo.getFileName())

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
