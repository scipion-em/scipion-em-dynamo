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
from enum import Enum
from os.path import join, abspath
import numpy as np
from emtable import Table
from dynamo.protocols.protocol_base_dynamo import DynamoProtocolBase
from pwem import ALIGN_NONE, ALIGN_2D
from pwem.emlib.image import ImageHandler
from pwem.objects import Transform
from pyworkflow import BETA
from pyworkflow.object import Set, String
from pyworkflow.protocol import GE, LEVEL_ADVANCED, IntParam, PointerParam, BooleanParam, EnumParam
from pyworkflow.utils import Message, makePath
from tomo.objects import SetOfTiltSeries, TiltSeries, SetOfTomograms, Tomogram
from dynamo import Plugin

# Fields in Dynamo align star files
ALI_TABLE = 'align'
INDEX = 'index'
USED = 'used'
TILT_ANGLE = 'thetas'
TILT_AXIS_ANGLE = 'psis'
SHIFT_X = 'x'
SHIFT_Y = 'y'
TYPE_DICT_FOR_PARSING = {
    INDEX: int,
    USED: int,
    TILT_ANGLE: float,
    TILT_AXIS_ANGLE: float,
    SHIFT_X: float,
    SHIFT_Y: float
}


class DynRecTomoChoices(Enum):
    SIRT = 0
    WBP = 1


class DynamoTsAlignOuts(Enum):
    tiltSeries = SetOfTiltSeries()
    tiltSeriesInterpolated = SetOfTiltSeries()
    tomograms = SetOfTomograms()


class DynamoTsAlign(DynamoProtocolBase):
    """Tilt  series alignment"""

    _label = 'tilt series alignment'
    _possibleOutputs = DynamoTsAlignOuts
    _devStatus = BETA
    tsAliPrjName = 'scipionDynamoTsAlign.AWF'
    tsOutFileName = 'alignedBinnedStack.mrc'
    tsOutAliFileName = 'stackAlignerImod.star'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._inTsSet = None
        self._sRate = -1
        self._beadRadiusPx = -1
        self._maskRadiusPx = -1
        self.excludedViewsMsg = String()

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputTs', PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label="Tilt series",
                      important=True)
        form.addParam('genInterpolated', BooleanParam,
                      default=True,
                      label='Generate the interpolated aligned TS?')
        group = form.addGroup('Detection settings')
        group.addParam('beadDiamNm', IntParam,
                       default=-1,
                       validators=[GE(1)],
                       label='Gold bead diameter nm]')
        group.addParam('maskDiamNm', IntParam,
                       default=-1,
                       expertLevel=LEVEL_ADVANCED,
                       label='Mask around gold bead diameter [Å]',
                       help='The general police is to choose a radius that diameter the white "halo" around the '
                            'bead and a couple of additional pixels.|nIf set to -1, it will be assumed as 1.5 * '
                            'gold bead diameter')
        group.addParam('templateSideLengthPix', IntParam,
                       default=-1,
                       expertLevel=LEVEL_ADVANCED,
                       label='Template side length [px]',
                       help='It should be at least twice the current value of the mask. It is used for detection, '
                            'alignment and depiction of sets of markers.\nIf set to -1, it will be considered as '
                            'twice of the current value of the mask.')
        form.addParam('recTomo', BooleanParam,
                      default=False,
                      label='Reconstruct the tomogram?')
        group = form.addGroup('Reconstruction settings', condition='recTomo')
        group.addParam('binning', IntParam,
                       default=4,
                       label='binning factor',
                       help='The tomogram reconstruction implies the generation of the interpolated tilt series. The '
                            'binning factor introduced will be applied to both.')
        group.addParam('recMethod', EnumParam,
                       choices=[DynRecTomoChoices.SIRT.name, DynRecTomoChoices.WBP.name],
                       display=EnumParam.DISPLAY_HLIST,
                       label="Reconstruction method",
                       default=DynRecTomoChoices.SIRT.value,
                       help="Choose either SART or weighted back "
                            "projection (WBP).")
        group.addParam('tomoThk', IntParam,
                       default=-1,
                       label='Thickness (binned) [px]',
                       help='If set to -1, the thickness will be internally calculated as 1/3 of the X '
                            'dimension of the binned tilt series.')
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        mdObjDict = self._initialize()
        for mdObj in mdObjDict.values():
            self._insertFunctionStep(self.convertInputStep, mdObj)
            self._insertFunctionStep(self.runTsAlignStep, mdObj)
            self._insertFunctionStep(self.createOutputStep, mdObj)
        self._insertFunctionStep(self._closeOutputSet)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self._inTsSet = self.inputTs.get()
        self._sRate = self._inTsSet.getSamplingRate()
        self._binnedSRate = self._sRate * 2 ** self.getBinningFactor()
        self._beadRadiusPx = self._getRadiusInPix(self.beadDiamNm.get())
        self._maskRadiusPx = self._getMaskRadiusPix()
        mdObjDict = {}
        for ts in self._inTsSet:
            ts = ts.clone(ignoreAttrs=[])
            tsId = ts.getTsId()
            tsDir = self._getExtraPath(tsId)
            outAliDir = join(tsDir, self.tsAliPrjName, 'align')
            mdObjDict[tsId] = DynTsAliMdObj(ts=ts,
                                            tsDir=tsDir,
                                            tltFile=join(tsDir, tsId + '.tlt'),
                                            matlabFile=join(tsDir, 'alignTs.m'),
                                            tsInterpFileName=join(outAliDir, self.tsOutFileName),
                                            outAliFile=join(outAliDir, self.tsOutAliFileName),
                                            tomoFileName=join(tsDir, self.tsAliPrjName, 'reconstruction',
                                                              self._getTomoFileName())
                                            )
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

    def createOutputStep(self, mdObj):
        # Read the alignment file
        aliData = self._readAliFile(mdObj.outAliFile)
        # Tilt series
        outTsSet = self._getOutputSetOfTs()
        self._fillTsAndUpdateTsSet(mdObj, outTsSet, aliData)
        # Interpolated tilt series
        if self.genInterpolated.get():
            outTsSetInterp = self._getOutputSetOfTs(interpolated=True)
            self._fillTsAndUpdateTsSet(mdObj, outTsSetInterp, aliData, interpolated=True)
        # Tomograms
        if self.recTomo.get():
            outTomos = self._getOutputSetOfTomos()
            self._createTomoAndUpdateTomoSet(mdObj, outTomos)

        self._store()

    # --------------------------- DEFINE utils functions ----------------------
    def _genMatlabCode(self, mdObj):

        acq = self._inTsSet.getAcquisition()
        # Create the workflow
        cmd = "name = '%s';\n" % self.tsAliPrjName
        cmd += "folder = '%s';\n" % mdObj.tsDir
        cmd += "u = dtsa(name,'--nogui','-path',folder, 'fp',1);\n"
        # Entering the data ######################################
        # Basic data
        cmd += "u.enter.tiltSeries('%s');\n" % abspath(mdObj.ts.getFirstItem().getFileName())
        cmd += "u.enter.tiltAngles('%s');\n" % mdObj.tltFile
        # Acquisition data
        cmd += "u.enter.settingAcquisition.apix(%f);\n" % self._inTsSet.getSamplingRate()
        cmd += "u.enter.settingAcquisition.sphericalAberration(%f);\n" % acq.getSphericalAberration()
        cmd += "u.enter.settingAcquisition.amplitudeContrast(%f);\n" % acq.getAmplitudeContrast()
        cmd += "u.enter.settingAcquisition.voltage(%f);\n" % acq.getVoltage()
        # Detection settings
        cmd += "u.enter.settingDetection.beadRadius(%i);\n" % self._beadRadiusPx
        cmd += "u.enter.settingDetection.maskRadius(%i);\n" % self._maskRadiusPx
        cmd += "u.enter.templateSidelength(%i);\n" % self._getTemplateSideLengthPix()
        # Computing settings
        cmd += "u.enter.settingComputing.parallelCPUUse(1);\n"  # enable the use of parallel cores
        cmd += "u.enter.settingComputing.cpus(%i);\n" % self.numberOfThreads.get()
        # Alignment
        cmd += "u.run.area.uptoRefinement();\n"
        cmd += "u.area.alignment.step.fixAlignmentMarkers.run();\n"
        if self.genInterpolated.get() or self.recTomo.get():  # The reconstruction needs the interpolated
            cmd += ("u.area.alignment.step.alignWorkingStack.parameterSet.alignmentBinLevel(%i);\n" %
                    self.getBinningFactor())
            cmd += "u.area.alignment.step.alignWorkingStack.run()\n"

        # Reconstruction
        # -------------------------------------------------------------------------------------------------------------------------------------------------------------
        # {'reconstructionFullSize'      }    {'400  400  400'}    {'400  400  400'}    {'reconstruction size'                      }    {'pixels'                   }
        # {'reconstructionShiftCenter'   }    {'0  0  0'      }    {'0  0  0'      }    {'shift tomogram from center'               }    {'pixels'                   }
        # {'useCenterOnbinnedCoordinates'}    {'0'            }    {'0'            }    {'use center on binned coordinates'         }    {'pixels in binned tomogram'}
        # {'centerBinnedCoordinatesValue'}    {'0  0  0'      }    {'0  0  0'      }    {'coordinates in binned tomogram'           }    {'pixels in binned tomogram'}
        # {'reconstructFullSIRT'         }    {'1'            }    {'1'            }    {'non-binned sized SIRT-like reconstruction'}    {0x0 char                   }
        # {'reconstructFullWBP'          }    {'1'            }    {'1'            }    {'WBP reconstruction'                       }    {0x0 char                   }
        # {'reconstructFullWBPCTF'       }    {'0'            }    {'0'            }    {'CTF-corrected WBP reconstruction'         }    {0x0 char                   }
        # -------------------------------------------------------------------------------------------------------------------------------------------------------------
        if self.recTomo.get():
            # Disable the full size reconstruction
            cmd += "u.area.reconstruction.step.fullReconstruction.parameterSet.reconstructFullWBP(0);\n"
            cmd += "u.area.reconstruction.step.fullReconstruction.parameterSet.reconstructFullSIRT(0);\n"
            sirtValue = 1 if self.recMethod.get() == DynRecTomoChoices.SIRT.value else 0
            wbpValue = 1 if self.recMethod.get() == DynRecTomoChoices.WBP.value else 0
            cmd += ("u.area.reconstruction.step.binnedReconstruction.parameterSet.reconstructBinnedWBP(%i);\n"
                    % wbpValue)
            cmd += ("u.area.reconstruction.step.binnedReconstruction.parameterSet.reconstructBinnedSIRT(%i);\n"
                    % sirtValue)
            cmd += ("u.area.reconstruction.step.binnedReconstruction.parameterSet.reconstructionBinnedHeight(%i);\n"
                    % self._getThickness())
            cmd += "u.area.reconstruction.run();\n"

        return cmd

    def _getOutputSetOfTs(self, interpolated=False) -> SetOfTiltSeries:

        if interpolated:
            outSetName = self._possibleOutputs.tiltSeriesInterpolated.name
            suffix = 'interpolated'
            sRate = self._binnedSRate
            ali = ALIGN_NONE
        else:
            outSetName = self._possibleOutputs.tiltSeries.name
            suffix = ''
            sRate = self._sRate
            ali = ALIGN_2D

        outTsSet = getattr(self, outSetName, None)
        if outTsSet:
            outTsSet.enableAppend()
        else:
            outTsSet = SetOfTiltSeries.create(self._getPath(), template='tiltseries', suffix=suffix)
            outTsSet.copyInfo(self._inTsSet)
            outTsSet.setSamplingRate(sRate)
            outTsSet.setAlignment(ali)
            outTsSet.setStreamState(Set.STREAM_OPEN)
            # Create the outputs and define the relations
            self._defineOutputs(**{outSetName: outTsSet})
            self._defineSourceRelation(self._inTsSet, outTsSet)

        return outTsSet

    def _getOutputSetOfTomos(self) -> SetOfTomograms:

        outSetName = self._possibleOutputs.tomograms.name
        outputSetOfTomograms = getattr(self, outSetName, None)
        if outputSetOfTomograms:
            outputSetOfTomograms.enableAppend()
        else:
            outputSetOfTomograms = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
            outputSetOfTomograms.copyInfo(self._inTsSet)
            outputSetOfTomograms.setSamplingRate(self._binnedSRate)
            outputSetOfTomograms.setStreamState(Set.STREAM_OPEN)
            self._defineOutputs(**{outSetName: outputSetOfTomograms})
            self._defineSourceRelation(self._inTsSet, outputSetOfTomograms)

        return outputSetOfTomograms

    def _readAliFile(self, aliFile):
        # The alignment file is generated as MATLAB files that can be parsed calling Dynamo, and as star files, that
        # can be read using the emtable package. The second one is chosen because of the efficiency of a call to
        # emtable compared to a call to Dynamo

        # First of all, the alignment file needs to be copied and edited, so they have table names that are accepted by
        # emtable. This is what we have:
        # data_
        # _centerX 1855
        # _centerY 1919.5
        # _centerZ 0
        # loop_
        # _index
        # _used
        # _thetas
        # _psis
        # _x
        # _y
        # 1     1    -57    85.2    -235.8922       -399.7366
        # 2     1    -54    85.2    -90.04349       -378.4401
        # [...]
        #
        # The expected names must be data_[name], and they must appear before the loop_ line, that indicates the
        # beginning of the contents of a new table
        fixedAliFile = self._fixStarTableName(aliFile)
        dataTable = Table()
        dataTable.read(fixedAliFile, tableName=ALI_TABLE, types=TYPE_DICT_FOR_PARSING)
        return dataTable

    @staticmethod
    def _fixStarTableName(inFile):
        newFileName = inFile.replace('.star', '_fixed.star')
        # Read the file
        with open(inFile, "r") as f:
            lines = f.readlines()
        # Locate the line with the loop_ statement
        ind = lines.index('loop_\n')
        # Add the table name in that position, so emtable can read it correctly
        lines.insert(ind, 'data_align\n')
        # Create the new file with the updated contents (table names)
        with open(newFileName, "w") as f:
            f.writelines(lines)
        return newFileName

    def _createTomoAndUpdateTomoSet(self, mdObj, outTomos):
        ts = mdObj.ts
        newTomogram = Tomogram()
        newTomogram.setLocation(mdObj.tomoFileName)
        newTomogram.setSamplingRate(self._binnedSRate)
        newTomogram.setOrigin()
        newTomogram.setAcquisition(ts.getAcquisition())
        newTomogram.setTsId(ts.getTsId())

        outTomos.append(newTomogram)
        outTomos.update(newTomogram)
        outTomos.write()

    def _fillTsAndUpdateTsSet(self, mdObj, outTsSet, aliData, interpolated=False):
        ts = mdObj.ts
        tsId = ts.getTsId()
        # The interpolated TS may contain fewer images if there have been image exclusion, so it's created from scratch
        # instead of cloning it
        outTs = TiltSeries(tsId=tsId) if interpolated else ts.clone()
        outTs.copyInfo(ts)
        outTsSet.append(outTs)
        identityMatrix = np.eye(3)
        outTiltAxisAngle = aliData[0].get(TILT_AXIS_ANGLE)  # It's the same for all the tilt images
        acq = outTs.getAcquisition()
        acq.setTiltAxisAngle(outTiltAxisAngle)
        outTs.setAcquisition(acq)
        outTs.setInterpolated(interpolated)
        excludedViewsList = []
        # Tilt series
        for i, ti in enumerate(ts.iterItems()):
            aliRow = aliData[i]
            enabled = bool(aliRow.get(USED))
            if interpolated and not enabled:
                excludedViewsList.append(i)
                continue
            else:
                transform = Transform()
                outTi = ti.clone()
                outTi.copyInfo(ti, copyId=True)
                trMatrix = identityMatrix
                if interpolated:
                    outFileName = mdObj.tsInterpFileName
                else:
                    outFileName = ti.getFileName()
                    if enabled:
                        acq = outTi.getAcquisition()
                        acq.setTiltAxisAngle(outTiltAxisAngle)
                        # Alignment data
                        tilt = aliRow.get(TILT_ANGLE)
                        sx = aliRow.get(SHIFT_X)
                        sy = aliRow.get(SHIFT_Y)
                        outTi.setTiltAngle(tilt)
                        outTi.setAcquisition(acq)
                        trMatrix = np.eye(3)
                        trMatrix[0, 0] = trMatrix[1, 1] = np.cos(np.deg2rad(outTiltAxisAngle))
                        trMatrix[0, 1] = np.sin(np.deg2rad(outTiltAxisAngle))
                        trMatrix[1, 0] = -trMatrix[0, 1]
                        trMatrix[0, 2] = sx
                        trMatrix[1, 2] = sy

            outTi.setFileName(outFileName)
            transform.setMatrix(trMatrix)
            outTi.setTransform(transform)
            outTi.setEnabled(enabled)
            outTs.append(outTi)

        outTsSet.update(outTs)
        outTsSet.write()

        if excludedViewsList:
            prevMsg = self.excludedViewsMsg.get() if self.excludedViewsMsg.get() else ''
            self.excludedViewsMsg.set(prevMsg + f'\n{tsId}: {excludedViewsList}')
            self._store(self.excludedViewsMsg)

    def _getRadiusInPix(self, inDiameterInNm):
        return 10 * inDiameterInNm / (2 * self._sRate)

    def _getMaskRadiusPix(self):
        maskDiamNm = self.maskDiamNm.get()
        return 1.5 * self._beadRadiusPx if maskDiamNm == -1 else self._getRadiusInPix(maskDiamNm)

    def _getTemplateSideLengthPix(self):
        tslPix = self.templateSideLengthPix.get()
        return 4 * self._maskRadiusPx if tslPix == -1 else tslPix

    def _getThickness(self):
        tomoThk = self.tomoThk.get()
        if tomoThk != -1:
            return tomoThk
        else:
            ih = ImageHandler()
            x, _, _, _ = ih.getDimensions(self._inTsSet.getFirstItem().getFirstItem().getFileName())
            return round(x / (2 ** self.getBinningFactor() * 3))

    def _getTomoFileName(self):
        if self.recMethod.get() == DynRecTomoChoices.SIRT.value:
            return 'binnedReconstructionSIRT.mrc'
        else:
            return 'binnedReconstructionWBP.mrc'

    # --------------------------- DEFINE info functions ----------------------
    def _validate(self):
        errorMsg = []
        if self.inputTs.get().interpolated():
            errorMsg.append("The introduced tilt series are interpolated. Please introduce non-interpolated.")
        return errorMsg

    def _summary(self):
        msg = []
        exludedViewsMsg = self.excludedViewsMsg.get()
        if exludedViewsMsg:
            msg.append("*Interpolated TS stacks have a few tilt images removed.*\n" +
                       self.excludedViewsMsg.get())
        return msg


class DynTsAliMdObj:

    def __init__(self, ts=None, tsDir=None, tltFile=None, matlabFile=None, tsInterpFileName=None, tomoFileName=None,
                 outAliFile=None):
        self.ts = ts
        self.tsDir = tsDir
        self.tltFile = tltFile
        self.matlabFile = matlabFile
        self.tsInterpFileName = tsInterpFileName
        self.tomoFileName = tomoFileName
        self.outAliFile = outAliFile
