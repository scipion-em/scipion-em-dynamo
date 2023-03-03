# **************************************************************************
# *
# * Authors:    David Herreros Calero (dherreros@cnb.csic.es)
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

import os
import glob
import numpy as np

from pyworkflow import BETA
import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pwem.protocols import EMProtocol

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfCoordinates3D
import tomo.constants as const

from dynamo import Plugin
from dynamo.convert import matrix2eulerAngles, readDynCoord


class DynamoSubBoxing(EMProtocol, ProtTomoBase):
    """SubBoxing using Dynamo"""

    _label = 'subBoxing'
    _devStatus = BETA
    OUTPUT_PREFIX = ""

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam, label="Input 3D Coordinates", important=True,
                      pointerClass='SetOfCoordinates3D', help='Select the SetOfCoordinates3D.')
        form.addParam('inputAverage', params.PointerParam, label="Input Average", important=True,
                      pointerClass='AverageSubTomogram', help='Select the AverageSubTomogram.')
        form.addParam('symmetry', params.StringParam, label='Point Group Symmetry', important=True,
                      default='C2')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('writeSetOfCoordinates3D')
        self._insertFunctionStep('launchDynamoSubBoxingStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def writeSetOfCoordinates3D(self):
        coords = self.inputCoordinates.get()
        tomos = coords.getPrecedents()
        samplingRateCoord = coords.getSamplingRate()
        samplingRateTomo = tomos.getFirstItem().getSamplingRate()
        self.factor = float(samplingRateCoord / samplingRateTomo)
        self.lines = []
        inputSet = tomos.getFiles()
        self.coordsFileName = self._getExtraPath('coords.txt')
        self.anglesFileName = self._getExtraPath('angles.txt')
        outC = open(self.coordsFileName, "w")
        outA = open(self.anglesFileName, 'w')
        for idt, item in enumerate(inputSet):
            coordDict = []
            for coord3DSet in coords.iterCoordinates():
                if pwutils.removeBaseExt(item) == pwutils.removeBaseExt(coord3DSet.getVolName()):
                    angles_coord = matrix2eulerAngles(coord3DSet.getMatrix())
                    x = round(self.factor * coord3DSet.getX(const.BOTTOM_LEFT_CORNER))
                    y = round(self.factor * coord3DSet.getY(const.BOTTOM_LEFT_CORNER))
                    z = round(self.factor * coord3DSet.getZ(const.BOTTOM_LEFT_CORNER))
                    outC.write("%d\t%d\t%d\t%d\n" % (x, y, z, idt + 1))
                    outA.write("%f\t%f\t%f\n" % (angles_coord[0], angles_coord[1], angles_coord[2]))
                    coordDict.append(coord3DSet.getObjId())
            if coordDict:
                self.lines.append(coordDict)
        outC.close()
        outA.close()

    def launchDynamoSubBoxingStep(self):
        matlabFile = self.writeMatlabCode()

        args = ' %s' % matlabFile
        self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

        pwutils.cleanPattern('*.m')

    def createOutputStep(self):
        inputCoords = self.inputCoordinates.get()
        tomos = inputCoords.getPrecedents()
        suffix = self._getOutputSuffix(SetOfCoordinates3D)
        outSet = self._createSetOfCoordinates3D(tomos, suffix)
        outSet.setBoxSize(np.loadtxt(self._getExtraPath('recBoxSize.txt')))
        outSet.setSamplingRate(inputCoords.getSamplingRate())
        outSet.setPrecedents(tomos)
        for coordFile in glob.glob(self._getExtraPath('*.tbl')):
            idt = int(''.join(filter(str.isdigit, pwutils.removeBaseExt(coordFile))))
            readDynCoord(coordFile, outSet, tomos[idt].clone())
        args = {}
        args['outputCoordinates3D'] = outSet
        self._defineOutputs(**args)
        self._defineSourceRelation(inputCoords, outSet)
        self._defineSourceRelation(self.inputAverage, outSet)

    # --------------------------- DEFINE utils functions ----------------------
    def writeMatlabCode(self):
        # Initialization params
        average = self.inputAverage.get()
        codeFilePath = self._getExtraPath("DynamoSubBoxing.m")

        # Write code to Matlab code file
        codeFid = open(codeFilePath, 'w')
        content = "box=%d\n" \
                  "aux=readmatrix('%s')\n" \
                  "angles=readmatrix('%s')\n" \
                  "coords=aux(:,1:3)\n" \
                  "tags=aux(:,4)'\n" \
                  "dynamo_mapview('%s')\n" \
                  "hGUI = findobj(0)\n" \
                  "hGUI.ShowHiddenHandles=true\n" \
                  "h=findobj(0,'type','figure')\n" \
                  "click_panel=findobj(h.Children,'tag','uipanel_clicks')\n" \
                  "north_ui=findobj(click_panel.Children,'tag','edit_anchor_north')\n" \
                  "distance_ui=findobj(click_panel.Children,'tag','text_anchor_distance')\n" \
                  "while ishandle(h)\n" \
                  "north=str2num(north_ui.String)\n" \
                  "distance=2*floor(str2num(distance_ui.String))\n" \
                  "pause(0.01)\n" \
                  "end\n" \
                  "rSubunitFromCenter=north-box/2\n" \
                  "writematrix(distance,'%s')\n" \
                  "for tag=unique(tags)\n" \
                  "tomoCoords=coords(tags==tag,:)\n" \
                  "tomoAngles=angles(tags==tag,:)\n" \
                  "t=dynamo_table_blank(size(tomoCoords,1),'r',tomoCoords,'angles',tomoAngles)\n" \
                  "ts=dynamo_subboxing_table(t,rSubunitFromCenter,'sym','%s')\n" \
                  "dynamo_write(ts,['%s','_ID',num2str(tag),'.tbl'])\n" \
                  "end\n" \
                  % (average.getXDim(), os.path.abspath(self.coordsFileName),
                     os.path.abspath(self.anglesFileName), os.path.abspath(average.getFileName()),
                     os.path.abspath(self._getExtraPath('recBoxSize.txt')),
                     self.symmetry.get(), os.path.abspath(self._getExtraPath('subBoxTable')))

        codeFid.write(content)
        codeFid.close()
        return codeFilePath

    # --------------------------- DEFINE info functions ----------------------
    def _validate(self):
        pass

    def _methods(self):
        methodsMsgs = []
        boxSize = np.loadtxt(self._getExtraPath('recBoxSize.txt'))
        if self.getOutputsSize() >= 1:
            methodsMsgs.append("Particle box size: *%s*" % boxSize)
            methodsMsgs.append("Total particles found: *%s*" %
                           self.outputCoordinates3D.getSize())
            methodsMsgs.append("Symmetry: *%s*" % self.symmetry.get())
        else:
            methodsMsgs.append("Output coordinates not ready yet.")
        return methodsMsgs

    def _summary(self):
        summary = []
        boxSize = np.loadtxt(self._getExtraPath('recBoxSize.txt'))
        if self.getOutputsSize() >= 1:
            summary.append("Particle box size: *%s*" % boxSize)
            summary.append("Total particles found: *%s*" %
                           self.outputCoordinates3D.getSize())
            summary.append("Symmetry: *%s*" % self.symmetry.get())
        else:
            summary.append("Output coordinates not ready yet.")
        return summary