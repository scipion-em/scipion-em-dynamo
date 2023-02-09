# **************************************************************************
# *
# * Authors:    David Herreros Calero (dherreros@cnb.csic.es)
# *             Scipion Team (scipion@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import glob
from enum import Enum
from os.path import exists, join, abspath

import numpy as np

from dynamo.utils import getCatalogFile, genMCode4ReadDynModel, genMCode4ReadAndSaveData, \
    getCurrentTomoPointsFile, getCurrentTomoCroppedFile
from dynamo.viewers.DynamoTomoProvider import DynamoTomogramProvider
from pyworkflow import BETA
import pyworkflow.utils as pwutils
from pyworkflow.object import String
from pyworkflow.protocol import LEVEL_ADVANCED
from pyworkflow.protocol.params import PointerParam, IntParam, BooleanParam
from pyworkflow.utils import removeBaseExt
from pyworkflow.utils.properties import Message
from pyworkflow.gui.dialog import askYesNo
from tomo.objects import SetOfMeshes, Coordinate3D

from tomo.protocols import ProtTomoPicking
import tomo.constants as const

from dynamo import Plugin, VLL_FILE, CATALOG_BASENAME
from dynamo.viewers.views_tkinter_tree import DynamoTomoDialog


class OutputDynPicking(Enum):
    meshes = SetOfMeshes


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
    _possibleOutputs = OutputDynPicking

    def __init__(self, **kwargs):
        ProtTomoPicking.__init__(self, **kwargs)
        self.dlg = None
        self.dynModelsPathDict = {}  # Used to store the path where the corresponding models to a tomo are stored

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtTomoPicking._defineParams(self, form)

        form.addParam('boxSize', IntParam, label="Box Size")
        form.addParam('modPrevMeshes', BooleanParam,
                      default=False,
                      label='Modify previous meshes?',
                      help='This option allows to add and/or remove coordinates to a previous SetOfMeshes')
        form.addParam('inputMeshes', PointerParam,
                      label="Input Meshes",
                      condition='modPrevMeshes',
                      allowsNull=True,
                      pointerClass='SetOfMeshes',
                      help='Select the previous SetOfMeshes you want to modify')
        form.addParam('deleteGenMFiles', BooleanParam,
                      default=True,
                      label='Remove the .m files generated after the execution?',
                      expertLevel=LEVEL_ADVANCED,
                      help='It can be useful for developers to check exactly what was .m files were generated by '
                           'Scipion and executed by Dynamo.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.launchDynamoBoxingStep, interactive=True)

    # --------------------------- STEPS functions -----------------------------
    def convertInputStep(self):
        """Initialize the catalogue"""
        # Create the vll (list of tomos) file
        vllFile = self._getExtraPath(VLL_FILE)
        tomoCounter = 1  # Matlab begins counting in 1
        with open(vllFile, 'w') as tomoFid:
            for tomo in self.inputTomograms.get().iterItems():
                tomoPath = abspath(tomo.getFileName())
                tomoFid.write(tomoPath + '\n')
                self.dynModelsPathDict[tomo.getTsId()] = self._getExtraPath(CATALOG_BASENAME, 'tomograms',
                                                                            'volume_%i' % tomoCounter, 'models')
                tomoCounter += 1

        catalogFile = getCatalogFile(self._getExtraPath(), withExt=False)
        if self.modPrevMeshes.get():
            # Save coordinates into .txt file for each tomogram and pass them to dynamo
            inputMeshes = self.inputMeshes.get()
            inputTomograms = self.inputTomograms.get()
            for tomo in inputTomograms:
                outFileCoord = self._getExtraPath(pwutils.removeBaseExt(tomo.getFileName())) + ".txt"
                coordsInCurrentTomo = []
                for coord in inputMeshes.iterCoordinates(tomo.getObjId()):
                    coordsInCurrentTomo.append(list(coord.getPosition(const.BOTTOM_LEFT_CORNER)) + [coord.getGroupId()])
                if coordsInCurrentTomo:
                    np.savetxt(outFileCoord, np.asarray(coordsInCurrentTomo), delimiter=' ')

            # Create small program to tell Dynamo to save the inputMeshes in a Ellipsoidal Vesicle Model
            contents = "dcm -create %s -fromvll %s\n" % (catalogFile, vllFile)
            contents += "catalogue=dread('%s')\n" % getCatalogFile(self._getExtraPath())
            contents += "nVolumes=length(catalogue.volumes)\n"
            contents += "for idv=1:nVolumes\n"
            contents += "tomoPath=catalogue.volumes{idv}.fullFileName()\n"
            contents += "tomoIndex=catalogue.volumes{idv}.index\n"
            contents += "[~,tomoName,~]=fileparts(tomoPath)\n"
            contents += "coordFile=fullfile('%s', tomoName '.txt']\n" % self._getExtraPath()
            contents += "if ~isfile(coordFile)\n"
            contents += "continue\n"
            contents += "end\n"
            contents += "coords_ids=readmatrix(coordFile,'Delimiter',' ')\n"
            contents += "idm_vec=unique(coords_ids(:,4))'\n"
            contents += "end\n"
            contents += "exit\n"
        else:
            contents = "dcm -create %s -vll %s\n" % (catalogFile, vllFile)

        codeFile = self._getExtraPath('coords2model.m')
        with open(codeFile, 'w') as codeFid:
            codeFid.write(contents)

        # Tell Dynamo to create the catalogue with the models
        args = ' %s' % codeFile
        self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    def launchDynamoBoxingStep(self):
        tomoList = []
        for tomo in self.inputTomograms.get().iterItems():
            tomogram = tomo.clone()
            tomoList.append(tomogram)

        tomoProvider = DynamoTomogramProvider(tomoList, self._getExtraPath(), "txt")
        self.dlg = DynamoTomoDialog(None, self._getExtraPath(), provider=tomoProvider)

        # Open dialog to request confirmation to create output
        import tkinter as tk
        if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, tk.Frame()):
            self._createOutput()
            # Delete the .m generated files if requested
        if self.deleteGenMFiles.get():
            pwutils.cleanPattern(self._getExtraPath('*.m'))

    def _createOutput(self):
        precedents = self.inputTomograms.get()
        meshes = SetOfMeshes.create(self._getPath(), template='meshes%s.sqlite')
        meshes.setPrecedents(precedents)
        meshes.setSamplingRate(precedents.getSamplingRate())
        meshes.setBoxSize(self.boxSize.get())
        meshes._dynCatalogue = String(getCatalogFile(self._getExtraPath()))  # Extended attribute
        tomoIdDict = {tomo.getTsId(): tomo for tomo in precedents}
        for tomoId, modelsDir in self.dynModelsPathDict.items():
            tomo = tomoIdDict[tomoId]
            pointsFile = getCurrentTomoPointsFile(self._getExtraPath(), tomo)  # For the clicked points
            croppedFile = getCurrentTomoCroppedFile(self._getExtraPath(), tomo)  # For the cropped points and angles in case of mesh calculation from the GUI
            modelsDir = self.readModels(modelsDir, tomoId, pointsFile=pointsFile, croppedFile=croppedFile)
            if modelsDir:
                with open(pointsFile, 'r') as coordFile:
                    for line in coordFile:
                        coord = Coordinate3D()
                        values = line.replace('\n', '').split('\t')
                        coord.setVolume(tomo)
                        coord.setPosition(float(values[0]), float(values[1]), float(values[2]),
                                          const.BOTTOM_LEFT_CORNER)
                        coord.setGroupId(int(values[3]))
                        # Extended attributes
                        coord._dynModelName = String(values[4])
                        coord._dynModelFile = String(values[5])
                        meshes.append(coord)

                # Generate a set of coordinates in case the user generated the meshes from the Dynamo GUI
                if exists(croppedFile):
                    pass

        self._defineOutputs(**{OutputDynPicking.meshes.name: meshes})
        self._defineSourceRelation(self.inputTomograms, meshes)
        self._updateOutputSet(OutputDynPicking.meshes.name, meshes, state=meshes.STREAM_CLOSED)

    def readModels(self, modelsDir, tomoId, pointsFile=None, croppedFile=None):
        """Read the models generated for each tomograms and write the info to a the corresponding file,
        depending if there was only a picking or a picking and a mesh calculation"""
        modelsInDir = False
        modelFilesInDir = glob.glob(join(modelsDir, '*.omd'))
        if modelFilesInDir:  # The models directories are created before the annotation step, so they can be empty
            modelsInDir = True
            vesicleInd = 1
            for modelFile in modelFilesInDir:  # There's only one model per vesicle, a file is generated for each model
                contents = genMCode4ReadDynModel(modelFile)  # Read current model
                contents += genMCode4ReadAndSaveData(vesicleInd, modelFile, pointsFile=pointsFile, croppedFile=croppedFile)
                vesicleInd += 1

                codeFile = self._getExtraPath('parseModels_%s_%s.m' % (tomoId, removeBaseExt(modelFile)))
                with open(codeFile, 'w') as codeFid:
                    codeFid.write(contents)
                args = ' %s' % codeFile
                self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())
        return modelsInDir

    # --------------------------- DEFINE info functions ----------------------
    @staticmethod
    def getMethods(output):
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

