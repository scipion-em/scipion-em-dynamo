# **************************************************************************
# *
# * Authors:    Scipion Team (scipion@cnb.csic.es)
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
from collections import namedtuple
from enum import Enum
from os.path import join, basename, abspath

from dynamo.protocols.protocol_base_dynamo import DynamoProtocolBase
from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import params, GT, GE
from pyworkflow.utils import removeBaseExt, getExt, Message, createLink, makePath
from tomo.objects import Tomogram, SetOfTomograms, SetOfTiltSeries, TiltSeriesBase, TiltImageBase
from dynamo import Plugin


class DynamoTsAlignOuts(Enum):
    tiltSeries = SetOfTiltSeries


class DynamoTsAlign(DynamoProtocolBase):
    """Tilt  series alignment"""

    _label = 'tilt series alignment'
    _possibleOutputs = DynamoTsAlignOuts
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._acq = None
        # self.content = ''
        # self.finalTomoNamesDict = {}

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputTs', params.PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label="Tilt series",
                      important=True)

        form.addSection('Detection settings')
        form.addParam('binning', params.IntParam,
                      default=1,
                      validators=[GE(1)],
                      label='Detection binning factor',
                      help='Binning is used internally only for quick and robust bead detection.')
        form.addParam('beadRadius', params.IntParam,
                      default=-1,
                      validators=[GE(1)],
                      label='Unbinned TS gold bead radius [pix.]')
        form.addParam('maskRadius', params.IntParam,
                      default=-1,
                      validators=[GE(1)],
                      label='Radius of the mask around gold bead [pix.]',
                      help='The general police is to choose a radius that covers the white "halo" around the bead and '
                           'a couple of additional pixels')
        form.addParam('templateSideLength', params.IntParam,
                      default=-1,
                      label='Template side length [pix.]',
                      help='It should be at least twice the current value of the mask. It is used for detection, '
                           'alignment and depiction of sets of markers.\nIf set to -1, it will be considered as '
                           'twice of the current value of the mask.')

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        mdObjDict = self._initialize()
        for mdObj in mdObjDict.values():
            self._insertFunctionStep(self.convertInputStep, mdObj)
            self._insertFunctionStep(self.runTsAlignStep, mdObj)
        #     origTomoName = tomo.getFileName()
        #     finalTomoName = self.getFinalTomoName(tomo)
        #     self.finalTomoNamesDict[finalTomoName] = tomo
        #     if self.doConvertFiles:
        #         self._insertFunctionStep(self.convertInputStep, origTomoName, finalTomoName)
        #     # Generate one unique file with all the tomograms to be binned
        #     self._insertFunctionStep(self.createMCodeStep, origTomoName, finalTomoName)
        # # That way, we can carry out the binning of all the tomograms provided with only one call to MATLAB, improving
        # # the efficiency, as this call is very slow
        # self._insertFunctionStep(self.binTomosStep)
        # self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        inTsSet = self.inputTs.get()
        self._acq = inTsSet.getAcquisition()
        mdObjDict = {}
        for ts in inTsSet:
            ts = ts.clone(ignoreAttrs=[])
            tsId = ts.getTsId()
            tsDir = self._getExtraPath(tsId)
            mdObjDict[tsId] = DynTsAliMdObj(ts=ts,
                                            tsDir=tsDir,
                                            workflowDir=join(tsDir, 'workflow'),
                                            tltFile=join(tsDir, tsId + '.tlt'),
                                            matlabFile=join(tsDir, 'alignTs.m'))
        return mdObjDict

    @staticmethod
    def convertInputStep(mdObj):
        # Create a directory for the current TS under extra dir
        makePath(mdObj.tsDir)
        # Create the tlt file that corresponds to the current ts
        mdObj.ts.generateTltFile(mdObj.tltFile)

    def runTsAlignStep(self, mdObj):
        mCode = self._genMatlabCode(mdObj)
        with open(mdObj.matlabFile, 'w') as codeFile:
            codeFile.write(mCode)
        args = ' %s' % mdObj.matlabFile
        self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    # def _initialize(self):
    #     self.inputTomos = self.inputTomos.get()
    #     self.doConvertFiles = not self.isCompatibleFileFormat()
    #
    # @staticmethod
    # def convertInputStep(origName, finalName):
    #     ih = ImageHandler()
    #     ih.convert(origName, finalName)
    #
    # def createMCodeStep(self, origTomoName, finalTomoName):
    #     # FROM DYNAMO:
    #     # ______________________________________________________________________________________________
    #     # bin(fileIn,fileOut,binFactor,varargin)
    #     # Create a Volume template object
    #     # p.addParamValue('slabSize',[],'short','ss');
    #     # p.addParamValue('matlabWorkers',0,'short','mw');
    #     # p.addParamValue('maximumMegaBytes',[]);
    #     # p.addParamValue('showStatistics',false,'short','sst');
    #     # ______________________________________________________________________________________________
    #     self.content += "dpktomo.tools.bin('%s', '%s', %i, 'slabSize', %i, 'matlabWorkers', %i, " \
    #                     "'showStatistics', true)\n" % \
    #                     (origTomoName, finalTomoName, super().getBinningFactor(), self.zChunk.get(),
    #                      self.numberOfThreads.get())
    #
    # def binTomosStep(self):
    #     mFile = self._getExtraPath('binTomograms.m')
    #     with open(mFile, 'w') as codeFile:
    #         codeFile.write(self.content)
    #     args = ' %s' % mFile
    #     self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())
    #
    # def createOutputStep(self):
    #     sr = self.inputTomos.getSamplingRate() * super().getBinningFactor(forDynamo=False)  # Not for, but from
    #     outTomos = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
    #     outTomos.setSamplingRate(sr)
    #     for tomoName, inTomo in self.finalTomoNamesDict.items():
    #         tomo = Tomogram()
    #         tomo.copyInfo(inTomo)
    #         tomo.setSamplingRate(sr)
    #         tomo.setFileName(tomoName)
    #         outTomos.append(tomo)
    #
    #     self._defineOutputs(**{DynamoTsAlignOuts.tomograms.name: outTomos})
    #     self._defineSourceRelation(self.inputTomos, outTomos)
    #
    # --------------------------- DEFINE utils functions ----------------------
    def _genMatlabCode(self, mdObj):
        # Create the workflow
        cmd = "name = 'scipionDynamoTsAlign';\n"
        cmd += "folder = '%s';\n" % mdObj.workflowDir
        cmd += "u = dtsa(name,'--nogui','-path',folder, 'fp',1);\n"
        # Entering the data ######################################
        # Basic data
        cmd += "u.enter.tiltSeries('%s');\n" % abspath(mdObj.ts.getFirstItem().getFileName())
        cmd += "u.enter.tiltAngles('%s');\n" % mdObj.tltFile
        # Acquisition data
        cmd += "u.enter.settingAcquisition.apix(%f);\n" % self.inputTs.get().getSamplingRate()
        cmd += "u.enter.settingAcquisition.sphericalAberration(%f);\n" % self._acq.getSphericalAberration()
        cmd += "u.enter.settingAcquisition.amplitudeContrast(%f);\n" % self._acq.getAmplitudeContrast()
        cmd += "u.enter.settingAcquisition.voltage(%f);\n" % self._acq.getVoltage()
        # Detection settings
        # cmd += "u.enter.settingDetection.detectionBinningFactor(%i);\n" % self.getBinningFactor(
        #     self.binning.get())
        cmd += "u.enter.settingDetection.detectionBinningFactor(%i);\n" % self.binning.get()
        cmd += "u.enter.settingDetection.beadRadius(%i);\n" % self.beadRadius.get()
        cmd += "u.enter.settingDetection.maskRadius(%i);\n" % self.maskRadius.get()
        cmd += "u.enter.templateSidelength(%i);\n" % self._getTemplateSideLength()
        # Computing settings
        cmd += "u.enter.settingComputing.parallelCPUUse(0);\n"  # enable the use of parallel cores
        cmd += "u.enter.settingComputing.cpus('*');\n"
        # cmd += "u.enter.settingComputing.cpus(%i);\n" % self.numberOfThreads.get()
        # Run the workflow

        # cmd += "fid = fopen('/home/jjimenez/test_JJ.txt', 'wt')\n"
        # cmd += "fprintf(fid, pwd);\n"
        # cmd += "fclose(fid);\n"

        cmd += "u.area.indexing.step.tiltGapFiller.parameterSet.residualsThreshold(8);\n"
        # cmd += "workflow.area.refinement.step.trimMarkers.parameterSet.maximalResidualObservation(5);"

        cmd += "u.run.all('noctf', 1);\n"
        # JORGE_END
        return cmd

    def _getTemplateSideLength(self):
        tsl = self.templateSideLength.get()
        return tsl if tsl != -1 else 2 * self.maskRadius.get()

    # def getConvertedOutFileName(self, inFileName):
    #     return self._getExtraPath(removeBaseExt(inFileName) + '.mrc')
    #
    # def getFinalTomoName(self, tomo):
    #     return self.getConvertedOutFileName(tomo.getFileName())
    #
    # def isCompatibleFileFormat(self):
    #     """Compatible with MRC and em (MRC with that extension)"""
    #     compatibleExts = ['.em', '.mrc']
    #     return True if getExt(self.inputTomos.getFirstItem().getFileName()) not in compatibleExts else False
    #
    # # --------------------------- DEFINE info functions ----------------------
    # def _methods(self):
    #     methodsMsgs = ["*Binning Factor*: %s" % self.binning.get()]
    #     return methodsMsgs
    #
    # def _summary(self):
    #     summary = []
    #     if self.getOutputsSize() >= 1:
    #         for _, outTomos in self.iterOutputAttributes():
    #             summary.append("Output *%s*:" % outTomos.getNameId().split('.')[1])
    #             summary.append("    * Binning Factor: *%s*" % self.binning.get())
    #             summary.append("    * Number of Tomograms Binned: *%s*" %
    #                            outTomos.getSize())
    #     else:
    #         summary.append("Output tomograms not ready yet.")
    #     return summary


class DynTsAliMdObj:

    def __init__(self, ts=None, tsDir=None, workflowDir=None, tltFile=None, matlabFile=None):
        self.ts = ts
        self.tsDir = tsDir
        self.workflowDir = workflowDir
        self.tltFile = tltFile
        self.matlabFile = matlabFile
