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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

import os

import pyworkflow.utils as pwutils
from pyworkflow.protocol.params import PointerParam, IntParam
from pyworkflow.utils.properties import Message
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.em.metadata import (MetaData, MDL_XCOOR, MDL_YCOOR, MDL_ZCOOR,
                                    MDL_PICKING_PARTICLE_SIZE)

from tomo.protocols import ProtTomoPicking
from tomo.objects import SetOfCoordinates3D, Coordinate3D

from dynamo import Plugin


class DynamoBoxing(ProtTomoPicking):
    """Manual picker from Dynamo"""

    _label = 'dynamo boxer'

    def __init__(self, **kwargs):
        ProtTomoPicking.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtTomoPicking._defineParams(self, form)

        form.addParam('boxSize', IntParam, label="Box Size")

        # form.addParam('inputCoordinates', PointerParam, label="Input Coordinates",
        #               allowsNull=True, pointerClass='SetOfCoordinates3D',
        #               help='Select the SetOfCoordinates3D.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        # # Copy input coordinates to Extra Path
        # if self.inputCoordinates.get():
        #     self._insertFunctionStep('copyInputCoords')

        # Launch Boxing GUI
        self._insertFunctionStep('launchDynamoBoxingStep', interactive=True)

    # --------------------------- STEPS functions -----------------------------
    # def copyInputCoords(self):
    #     setCoords2Jsons(self.inputTomograms.get(), self.inputCoordinates.get(), self._getExtraPath())

    def launchDynamoBoxingStep(self):
        codeFilePath = self.writeMatlabCode()

        args = ' %s' % codeFilePath
        self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

        # Open dialog to request confirmation to create output
        if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, None):
            self._createOutput()

    def _createOutput(self):
        coord3DSetDict = {}
        coord3DMap = {}
        setTomograms = self.inputTomograms.get()
        suffix = self._getOutputSuffix(SetOfCoordinates3D)
        coord3DSet = self._createSetOfCoordinates3D(setTomograms, suffix)
        coord3DSet.setName("tomoCoord")
        coord3DSet.setPrecedents(setTomograms)
        coord3DSet.setSamplingRate(setTomograms.getSamplingRate())
        coord3DSet.setBoxSize(self.boxSize.get())
        for tomo in setTomograms.iterItems():
            outFile = pwutils.join(self._getExtraPath(), pwutils.removeBaseExt(tomo.getFileName()) + '.txt')
            if not os.path.isfile(outFile):
                continue

            # Populate Set of 3D Coordinates with 3D Coordinates
            md = MetaData()
            md.readPlain(outFile, "xcoor ycoor zcoor")
            for objId in md:
                x = md.getValue(MDL_XCOOR, objId)
                y = md.getValue(MDL_YCOOR, objId)
                z = md.getValue(MDL_ZCOOR, objId)
                coord = Coordinate3D()
                coord.setPosition(x, y, z)
                coord.setVolume(tomo)
                coord3DSet.append(coord)

            coord3DSetDict['00'] = coord3DSet

        name = self.OUTPUT_PREFIX + suffix
        args = {}
        args[name] = coord3DSet
        self._defineOutputs(**args)
        self._defineSourceRelation(setTomograms, coord3DSet)

    # --------------------------- DEFINE utils functions ----------------------
    def writeMatlabCode(self):
        # Initialization params
        # codeFilePath = os.path.join(os.environ.get("SCIPION_HOME"), "software", "tmp", "DynamoPicker.m")
        codeFilePath = os.path.join(os.getcwd(), "DynamoPicker.m")
        listTomosFile = os.path.join(os.environ.get("SCIPION_HOME"), "software", "tmp", "tomos.vll")
        catalogue = os.path.abspath(self._getExtraPath("tomos"))

        # Create list of tomos file
        tomoFid = open(listTomosFile, 'w')
        for tomo in self.inputTomograms.get().iterItems():
            tomoPath = os.path.abspath(tomo.getFileName())
            tomoFid.write(tomoPath + '\n')
        tomoFid.close()

        # Write code to Matlab code file
        codeFid = open(codeFilePath, 'w')
        content = "path='%s'\n" \
                  "savePath = '%s'\n" \
                  "catalogue_name='%s'\n" \
                  "dcm -create %s -fromvll %s\n" \
                  "hGUI=findobj(0)\n" \
                  "hGUI.ShowHiddenHandles=true\n" \
                  "dcm -c %s\n" \
                  "dcm_handles=findobj(0,'name','catalogue manager: %s')\n" \
                  "l=dynamo_read(fullfile(path,'dcmData.m'))\n" \
                  "l=[l{:}]\n" \
                  "eval(l)\n" \
                  "items = []\n" \
                  "for c=dcm_data\n" \
                  "items(end+1)=str2num(strtrim(c{1}{1}))\n" \
                  "end\n" \
                  "models=cell(1,max(items))\n" \
                  "dcmOpen=true\n" \
                  "while dcmOpen\n" \
                  "handles=dpkslicer.getHandles()\n" \
                  "l=dynamo_read(fullfile(path,'waitForPicker.m'))\n" \
                  "l=[l{:}]\n" \
                  "eval(l)\n" \
                  "l=dynamo_read(fullfile(path,'checkDCM.m'))\n" \
                  "l=[l{:}]\n" \
                  "eval(l)\n" \
                  "end\n" \
                  "numTomos=length(dcm_data)\n" \
                  "for idt=1:numTomos\n" \
                  "[~,outFile,~]=fileparts(dcm_data{idt}{end})\n" \
                  "outFile=[outFile '.txt']\n" \
                  "crop_points=[]\n" \
                  "for model=models{idt}\n" \
                  "if iscell(model)\n" \
                  "crop_points=[crop_points; model{end}.crop_points]\n" \
                  "else\n" \
                  "crop_points=[crop_points; model.crop_points]\n" \
                  "end\n" \
                  "end\n" \
                  "if ~isempty(crop_points)\n" \
                  "writematrix(crop_points,fullfile(savePath,outFile),'Delimiter',' ')\n" \
                  "end\n" \
                  "end\n" % (os.path.abspath(os.getcwd()), self._getExtraPath(),
                             catalogue, catalogue, listTomosFile, catalogue, catalogue)

        codeFid.write(content)
        codeFid.close()

        # Create function files
        functions = ['checkDCM.m', 'waitForPicker.m', 'dcmData.m', 'extractModelsCatalogue.m']
        contents = ["try;\ndcm_handles.Name;\ndcmOpen=true;\ncatch;\ndcmOpen=false;\n"
                    "end;\npause(0.25);\n",
                    "if ~isempty(handles);\nfastslicer=handles.figure_fastslicer;\n"
                    "items=split(fastslicer.Name,',');\nidx=sscanf(strtrim(items{2}),'volume index %d');\n"
                    "uiwait(handles.figure_fastslicer);\nm=dynamo_read(fullfile(path,'extractModelsCatalogue.m'));\n"
                    "m=[m{:}];\neval(m);\nend;\n",
                    "aux=findobj(dcm_handles,'tag','uitable_main');\naux=aux.DisplayData();\n"
                    "items=size(aux,1);\ndcm_data=cell(1,items);\nfor idx=1:items;\ndcm_data{idx}={aux{idx,3} aux{idx,1} aux{idx,2}};\n"
                    "end;\n",
                    "model=dread(dcmodels(catalogue_name,'i',idx));\n"
                    "models{idx}=model;\n"]
        for function, content in zip(functions, contents):
            fid = open(os.path.join(os.getcwd(), function), 'w')
            fid.write(content)
            fid.close()

        return codeFilePath