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
from pwem.objects import Transform
from pyworkflow import BETA
from pyworkflow.object import Set
from pyworkflow.protocol import GE, LEVEL_ADVANCED, IntParam, PointerParam
from pyworkflow.utils import Message, makePath
from tomo.objects import SetOfTiltSeries, TiltSeries, TiltImage
from dynamo import Plugin


class DynamoTsAlignOuts(Enum):
    tiltSeries = SetOfTiltSeries()
    tiltSeriesInterpolated = SetOfTiltSeries()


class DynamoTsAlign(DynamoProtocolBase):
    """Tilt  series alignment"""

    _label = 'tilt series alignment'
    _possibleOutputs = DynamoTsAlignOuts
    _devStatus = BETA
    tsAliPrjName = 'scipionDynamoTsAlign.AWF'
    tsOutFileName = 'alignedFullStack.mrc'
    tsOutAliFileName = 'stackAlignerImod.star'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._acq = None
        self._sRate = -1
        self._beadRadiusPx = -1
        self._maskRadiusPx = -1
        # self.content = ''
        # self.finalTomoNamesDict = {}

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputTs', PointerParam,
                      pointerClass='SetOfTiltSeries',
                      label="Tilt series",
                      important=True)
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

        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        mdObjDict = self._initialize()
        for mdObj in mdObjDict.values():
            self._insertFunctionStep(self.convertInputStep, mdObj)
            self._insertFunctionStep(self.runTsAlignStep, mdObj)
            self._insertFunctionStep(self.createOutputStep, mdObj)
            self._insertFunctionStep(self.closeOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        inTsSet = self.inputTs.get()
        self._acq = inTsSet.getAcquisition()
        self._sRate = inTsSet.getSamplingRate()
        self._beadRadiusPx = self._getRadiusInPix(self.beadDiamNm.get())
        self._maskRadiusPx = self._getMaskRadiusPix()
        mdObjDict = {}
        for ts in inTsSet:
            ts = ts.clone(ignoreAttrs=[])
            tsId = ts.getTsId()
            tsDir = self._getExtraPath(tsId)
            outAliDir = join(tsDir, self.tsAliPrjName, 'align')
            mdObjDict[tsId] = DynTsAliMdObj(ts=ts,
                                            tsDir=tsDir,
                                            tltFile=join(tsDir, tsId + '.tlt'),
                                            matlabFile=join(tsDir, 'alignTs.m'),
                                            tsInterpFileName=join(outAliDir, self.tsOutFileName),
                                            outAliFile=join(outAliDir, self.tsOutAliFileName))
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
        # Note: it always generates the interpolated
        # Tilt series
        outTsSet = self._getOutputSetOfTs()
        self._fillTsAndUpdateTsSet(mdObj, outTsSet)
        # Interpolated tilt series
        outTsSetInterp = self._getOutputSetOfTs(interpolated=True)
        self._fillTsAndUpdateTsSet(mdObj, outTsSetInterp, interpolated=True)

    def closeOutputStep(self):
        inTsSet = self.inputTs.get()
        # Tilt series -> close
        outTsSet = self._getOutputSetOfTs()
        outTsSet.setStreamState(Set.STREAM_CLOSED)
        # Interpolated tilt series -> close
        outTsSetInterp = self._getOutputSetOfTs(interpolated=True)
        outTsSetInterp.setStreamState(Set.STREAM_CLOSED)
        # Create the outputs and define the relations
        self._defineOutputs(**{self._possibleOutputs.tiltSeries.name: outTsSet,
                               self._possibleOutputs.tiltSeriesInterpolated.name: outTsSetInterp})
        self._defineSourceRelation(inTsSet, outTsSet)
        self._defineSourceRelation(inTsSet, outTsSetInterp)

    # --------------------------- DEFINE utils functions ----------------------
    def _genMatlabCode(self, mdObj):
        # Create the workflow
        cmd = "name = '%s';\n" % self.tsAliPrjName
        cmd += "folder = '%s';\n" % mdObj.tsDir  #workflowDir
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
        # cmd += "u.enter.settingDetection.detectionBinningFactor(%i);\n" % self.binning.get()
        cmd += "u.enter.settingDetection.beadRadius(%i);\n" % self._beadRadiusPx
        cmd += "u.enter.settingDetection.maskRadius(%i);\n" % self._maskRadiusPx
        cmd += "u.enter.templateSidelength(%i);\n" % self._getTemplateSideLengthPix()
        # Computing settings
        cmd += "u.enter.settingComputing.parallelCPUUse(0);\n"  # enable the use of parallel cores
        # cmd += "u.enter.settingComputing.cpus('*');\n"
        cmd += "u.enter.settingComputing.cpus(%i);\n" % self.numberOfThreads.get()
        # Run the workflow

        # cmd += "fid = fopen('/home/jjimenez/test_JJ.txt', 'wt')\n"
        # cmd += "fprintf(fid, pwd);\n"
        # cmd += "fclose(fid);\n"

        # cmd += "u.area.indexing.step.tiltGapFiller.parameterSet.residualsThreshold(8);\n"
        # cmd += "workflow.area.refinement.step.trimMarkers.parameterSet.maximalResidualObservation(5);"

        # cmd += "u.run.all('noctf', 1);\n"

        # Alignment
        cmd += "u.run.area.uptoAlignment();\n"
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

        # cmd += "u.area.reconstruction.step.fullReconstruction.parameterSet.reconstructionShiftCenter([1855, 1919, 150]);\n"
        # cmd += "u.area.reconstruction.step.fullReconstruction.parameterSet.reconstructionFullSize([928, 960, 300]);\n"
        # cmd += "u.area.reconstruction.step.fullReconstruction.parameterSet.reconstructFullWBP(0);\n"
        # cmd += "u.area.reconstruction.step.fullReconstruction.parameterSet.reconstructFullSIRT(1);\n"
        # cmd += "u.area.reconstruction.step.binnedReconstruction.parameterSet.reconstructBinnedWBP(0);\n"
        # cmd += "u.area.reconstruction.step.binnedReconstruction.parameterSet.reconstructBinnedSIRT(0);\n"
        # cmd += "u.area.reconstruction.step.binnedReconstruction.parameterSet.reconstructionBinnedHeight(300);\n"
        # cmd += "u.area.reconstruction.run();\n"

        return cmd

    def _getOutputSetOfTs(self, interpolated=False):
        inTsSet = self.inputTs.get()
        if interpolated:
            outSetName = self._possibleOutputs.tiltSeriesInterpolated.name
            suffix = 'interpolated'
        else:
            outSetName = self._possibleOutputs.tiltSeries.name
            suffix = ''
        outTsSet = getattr(self, outSetName, None)
        if outTsSet:
            outTsSet.enableAppend()
        else:
            outTsSet = SetOfTiltSeries.create(self._getPath(), template='tiltseries', suffix=suffix)
            outTsSet.copyInfo(inTsSet)
            outTsSet.setSamplingRate(self._sRate)
            outTsSet.setStreamState(Set.STREAM_OPEN)
            setattr(self, outSetName, inTsSet)

        return outTsSet

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
        dataTable.read(fixedAliFile, tableName='align')
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

    def _fillTsAndUpdateTsSet(self, mdObj, outTsSet, interpolated=False):
        aliData = self._readAliFile(mdObj.outAliFile)
        tiltSeries = TiltSeries()
        tiltSeries.copyInfo(mdObj.ts)
        outTi = TiltImage()
        transform = Transform()
        # Tilt series
        for aliRow, ti in zip(aliData, mdObj.ts.iterItems()):
            trMatrix = np.eye(3)
            outTi.copyInfo(ti, copyId=True)
            enabled = bool(aliRow.get('used'))
            if interpolated:
                outFileName = mdObj.tsOutFileName
                if not enabled:
                    continue
            else:
                outFileName = ti.getFileName()
                if enabled:
                    # Alignment data
                    rot = aliRow.get('thetas')
                    tilt = aliRow.get('psis')
                    sx = aliRow.get('x')
                    sy = aliRow.get('y')
                    outTi.setTiltAngle(tilt)
                    outTi.getAcquisition().setTiltAxisAngle(rot)
                    trMatrix[0, 0] = trMatrix[1, 1] = np.cos(np.deg2rad(rot))
                    trMatrix[0, 1] = np.sin(np.deg2rad(rot))
                    trMatrix[1, 0] = -trMatrix[0, 1]
                    trMatrix[0, 2] = sx
                    trMatrix[1, 2] = sy

            outTi.setFileName(outFileName)
            transform.setMatrix(trMatrix)
            outTi.setTransform(transform)
            outTi.setEnabled(enabled)
            tiltSeries.append(outTi)

        outTsSet.update(tiltSeries)

    def _getRadiusInPix(self, inDiameterInNm):
        return 10 * inDiameterInNm / (2 * self._sRate)

    def _getMaskRadiusPix(self):
        maskDiamNm = self.maskDiamNm.get()
        return 1.5 * self._beadRadiusPx if maskDiamNm == -1 else self._getRadiusInPix(maskDiamNm)

    def _getTemplateSideLengthPix(self):
        tslPix = self.templateSideLengthPix.get()
        return 4 * self._maskRadiusPx if tslPix == -1 else tslPix

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

    def __init__(self, ts=None, tsDir=None, tltFile=None, matlabFile=None, tsInterpFileName=None, outAliFile=None):
        self.ts = ts
        self.tsDir = tsDir
        self.tltFile = tltFile
        self.matlabFile = matlabFile
        self.tsInterpFileName = tsInterpFileName
        self.outAliFile = outAliFile
