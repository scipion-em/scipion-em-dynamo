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
import numpy as np

from pyworkflow import BETA
import pyworkflow.utils as pwutils
from pyworkflow.protocol.params import PointerParam, IntParam, EnumParam
from pyworkflow.utils.properties import Message
from pyworkflow.gui.dialog import askYesNo

from tomo.protocols import ProtTomoPicking
from tomo.viewers.views_tkinter_tree import TomogramsTreeProvider
import tomo.constants as const

from dynamo import Plugin
from dynamo.viewers.views_tkinter_tree import DynamoTomoDialog
from dynamo.convert import textFile2Coords

class DynamoBoxing(ProtTomoPicking):
    """Manual vectorial picker from Dynamo. After choosing the Tomogram to be picked, the tomo slicer from Dynamo will be
    direclty loaded with all the models previously saved in the disk (if any).
    This picking will only save the "user points" defined in a set of models. It is possible to
    create several models at once in a given tomogram. Once the coordinates are defined,
    the models are automatically saved in the catalogue and registered.

    Currently the following Dynamo models are supported:
        - Ellipsoidal Vesicle"""

    _label = 'vectorial picking'
    _devStatus = BETA

    modelChoices = ["Ellipsoidal Vesicle", "Surface", "General"]
    modelNames = {modelChoices[0]: "ellipsoidalVesicle", modelChoices[1]: "surface", modelChoices[2]: "general"}
    OUTPUT_PREFIX = 'outputMeshes'

    def __init__(self, **kwargs):
        ProtTomoPicking.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtTomoPicking._defineParams(self, form)

        form.addParam('boxSize', IntParam, label="Box Size")
        form.addParam('selection', EnumParam, choices=['Yes', 'No'], default=1,
                      label='Modify previous meshes?', display=EnumParam.DISPLAY_HLIST,
                      help='This option allows to add and/or remove coordinates to a previous SetOfMeshes')
        form.addParam('inputMeshes', PointerParam, label="Input Meshes", condition='selection == 0',
                      allowsNull=True, pointerClass='SetOfMeshes',
                      help='Select the previous SetOfMeshes you want to modify')
        form.addParam('modelType', EnumParam,
                      choices=self.modelChoices, default=0,
                      label='Model type',
                      help='Select the type of model defined in the Tomograms.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        # Copy input coordinates to Extra Path
        self._insertFunctionStep('copyInputMeshes')

        # Launch Boxing GUI
        self._insertFunctionStep('launchDynamoBoxingStep', interactive=True)

    # --------------------------- STEPS functions -----------------------------

    def copyInputMeshes(self):
        # Initialize the catalogue
        listTomosFile = self._getExtraPath("tomos.vll")
        catalogue = os.path.abspath(self._getExtraPath("tomos"))

        # Create list of tomos file
        tomoFid = open(listTomosFile, 'w')
        for tomo in self.inputTomograms.get().iterItems():
            tomoPath = os.path.abspath(tomo.getFileName())
            tomoFid.write(tomoPath + '\n')
        tomoFid.close()

        # Save coordinates into .txt file for each tomogram
        codeFile = self._getExtraPath('coords2model.m')
        if self.selection.get() == 0:
            inputMeshes = self.inputMeshes.get()
            inputTomograms = self.inputTomograms.get()
            for tomo in inputTomograms.iterItems(iterate=False):
                outFileCoord = self._getExtraPath(pwutils.removeBaseExt(tomo.getFileName())) + ".txt"
                coords_tomo = []
                for coord in inputMeshes.iterCoordinates(tomo.getObjId()):
                    coords_tomo.append(list(coord.getPosition(const.BOTTOM_LEFT_CORNER)) + [coord.getGroupId()])
                if coords_tomo:
                    np.savetxt(outFileCoord, np.asarray(coords_tomo), delimiter=' ')

            # Create small program to tell Dynamo to save the inputMeshes in a Ellipsoidal Vesicle Model
            contents = "dcm -create %s -fromvll %s\n" \
                       "path='%s'\n" \
                       "catalogue=dread(['%s' '.ctlg'])\n" \
                       "nVolumes=length(catalogue.volumes)\n" \
                       "for idv=1:nVolumes\n" \
                       "tomoPath=catalogue.volumes{idv}.fullFileName()\n" \
                       "tomoIndex=catalogue.volumes{idv}.index\n" \
                       "[~,tomoName,~]=fileparts(tomoPath)\n" \
                       "coordFile=[path '/' tomoName '.txt']\n" \
                       "if ~isfile(coordFile)\n" \
                       "continue\n" \
                       "end\n" \
                       "coords_ids=readmatrix(coordFile,'Delimiter',' ')\n" \
                       "idm_vec=unique(coords_ids(:,4))'\n" \
                       "for idm=idm_vec\n" \
                       "model_name=['model_',num2str(idm)]\n" \
                       "coords=coords_ids(coords_ids(:,4)==idm,1:3)\n" \
                       "vesicle=dmodels.%s()\n" \
                       "vesicle.name=model_name\n" \
                       "addPoint(vesicle,coords(:,1:3),coords(:,3))\n" \
                       "vesicle.linkCatalogue('%s','i',tomoIndex,'s',1)\n" \
                       "vesicle.saveInCatalogue()\n" \
                       "end\n" \
                       "end\n" \
                       "exit\n" % (catalogue, listTomosFile, os.path.abspath(self._getExtraPath()), catalogue,
                                   self.modelNames[self.modelChoices[self.modelType.get()]], catalogue)
        else:
            contents = "dcm -create %s -fromvll %s\n" % (catalogue, listTomosFile)

        codeFid = open(codeFile, 'w')
        codeFid.write(contents)
        codeFid.close()

        # Tell Dynamo to create the catalogue with the models
        args = ' %s' % codeFile
        self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    def launchDynamoBoxingStep(self):

        tomoList = []
        for tomo in self.inputTomograms.get().iterItems():
            tomogram = tomo.clone()
            tomoName = pwutils.removeBaseExt(tomo.getFileName())
            outFile = self._getExtraPath(tomoName + '.txt')
            if os.path.isfile(outFile):
                tomogram.count = np.loadtxt(outFile, delimiter=' ').shape[0]
            else:
                tomogram.count = 0
            tomoList.append(tomogram)

        tomoProvider = TomogramsTreeProvider(tomoList, self._getExtraPath(), "txt")

        self.dlg = DynamoTomoDialog(None, self._getExtraPath(), provider=tomoProvider)

        # Open dialog to request confirmation to create output
        import tkinter as tk
        if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, tk.Frame()):
            self._createOutput()

        pwutils.cleanPattern(self._getExtraPath('*.m'))

    def _createOutput(self):
        textFile2Coords(self, self.inputTomograms.get(), self._getExtraPath(), False, True)

    # --------------------------- DEFINE info functions ----------------------
    def getMethods(self, output):
        msg = 'User picked %d particles ' % output.getSize()
        msg += 'with a particle size of %s.' % output.getBoxSize()
        return msg

    def _methods(self):
        methodsMsgs = []
        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes():
                msg = self.getMethods(output)
                methodsMsgs.append("%s: %s" % (self.getObjectTag(output), msg))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs

