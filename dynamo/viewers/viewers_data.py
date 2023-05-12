# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
# *              Scipion Team (scipion@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import datetime
import glob
from os import remove
from os.path import join, abspath, exists
from dynamo.utils import createBoxingOutputObjects, getNewestModelModDate, getDynamoModels
from dynamo.viewers.DynamoTomoProvider import DynamoTomogramProvider
from pwem.viewers import ObjectView
import pyworkflow.viewer as pwviewer
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.utils import makePath, removeBaseExt
from pyworkflow.utils.properties import Message
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.objects import SetOfCoordinates3D, SetOfMeshes
from dynamo.viewers.views_tkinter_tree import DynamoTomoDialog
from dynamo import VLL_FILE, MB_BY_LEVELS, MB_ELLIPSOIDAL, MB_VESICLE, MB_ELLIPSOIDAL_MARKED, MB_GENERAL, \
    MB_GENERAL_BOXES, PROJECT_DIR, DATA_MODIFIED_FROM_VIEWER

# Dynamo model names mapping
DYN_GEN_MODEL_OBJ_NAME = 'general'
DYN_EV_MODEL_OBJ_NAME = 'ellipsoidalVesicle'
dynModelsDict = {
    MB_BY_LEVELS: 'membraneByLevels',
    MB_VESICLE: 'vesicle',
    MB_ELLIPSOIDAL: DYN_EV_MODEL_OBJ_NAME,
    MB_ELLIPSOIDAL_MARKED: DYN_EV_MODEL_OBJ_NAME,
    MB_GENERAL: DYN_GEN_MODEL_OBJ_NAME,
    MB_GENERAL_BOXES: DYN_GEN_MODEL_OBJ_NAME,
}


class DynamoDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [SetOfMeshes, SetOfCoordinates3D]
    _name = 'Dynamo'

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return ObjectView(
            self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, coord=None, **kwargs):
        views = []
        path = self.protocol._getExtraPath()

        if isinstance(obj, SetOfCoordinates3D):  # SetofMeshes will also be True as it inherits from SetOfCoordinates3D
            outputCoords = obj
            precedentsPointer = outputCoords._precedentsPointer
            precedents = outputCoords.getPrecedents()
            tomoIdDict = {tomo.getTsId(): tomo.clone() for tomo in precedents}
            tomoList = list(tomoIdDict.values())
            nCoordsDict = {}
            # Write the coordinates to text files, as expected by Dynamo
            for tomo in tomoList:
                coordCounter = 0
                tomoId = tomo.getTsId()
                with open(join(path, removeBaseExt(tomo.getFileName()) + '.txt'), 'w') as fid:
                    for coord in outputCoords.iterCoordinates(tomo):
                        coords = list(coord.getPosition(BOTTOM_LEFT_CORNER))
                        modelType = self._getModelType(coord)
                        # angles = list(matrix2eulerAngles(coord.getMatrix()))[:3]  # Keep the angles, not the shifts (Dynamo boxing GUI is not expecting them at this point)
                        fid.write(",".join(map(str, coords + [coord.getGroupId(), modelType])) + "\n")
                        coordCounter += 1
                nCoordsDict[tomoId] = coordCounter

            # Remove the possible particle count files from previous executions that may cause incorrect values in the
            # tomo provider
            # [remove(partFile) for partFile in glob.glob(join(path, '*_count.txt'))]
            tomoProvider = DynamoTomogramProvider(tomoList, path, nParticlesDict=nCoordsDict)
            listTomosFile = join(path, VLL_FILE)

            # Create list of tomos file (VLL file), only required if using the Dynamo viewer in a protocol from any
            # other plugin
            if not exists(listTomosFile):
                with open(listTomosFile, 'w') as tomoFid:
                    for tomo in tomoList:
                        tomoPath = abspath(tomo.getFileName())
                        tomoFid.write(tomoPath + '\n')

            dynamoDialogCallingTime = datetime.datetime.now()
            self.dlg = DynamoTomoDialog(self._tkRoot, path,
                                        provider=tomoProvider,
                                        calledFromViewer=True,
                                        suffix='fixed')

            # Check if the modification file of the newest model file is higher than the time capture right before
            # calling the Dynamo dialog. In that case, it means that some modification was carried out by the user from
            # it and we have to ask if the changes should be saved
            if exists(self.protocol._getExtraPath(PROJECT_DIR, DATA_MODIFIED_FROM_VIEWER)):
                import tkinter as tk
                from dynamo.protocols.protocol_boxing import DynPickingOuts
                from dynamo.protocols.protocol_model_workflow import DynModelWfOuts
                from dynamo.protocols import DynamoBoxing

                frame = tk.Frame()
                # Because the coordinates are written as general models, they'll have cropped points and angles defined
                # (particularity of Dynamo for this kind of models). Hence, coordinates will be saved as if the model
                # workflow had been carried out
                if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, frame):
                    tmpPath = self.protocol._getTmpPath()
                    if not exists(tmpPath):  # Some files will be created there, and it may not exist
                        makePath(tmpPath)
                    outMeshes, outCoords = createBoxingOutputObjects(self.protocol, precedentsPointer,
                                                                     boxSize=outputCoords.getBoxSize(),
                                                                     savePicked=True,
                                                                     suffix='fixed')
                    # Preserve previous outputs and generate new ones if necessary. The possible outputs are
                    # different depending on if the viewer was called from the boxing protocol or from the model
                    # workflow protocol results
                    coordsName = DynPickingOuts.coordinates.name
                    prevCoords = getattr(self.protocol, coordsName, None)
                    if isinstance(self.protocol, DynamoBoxing):
                        meshesName = DynPickingOuts.meshes.name
                        prevMeshes = getattr(self.protocol, meshesName, None)
                        prevOutputs = [prevCoords, prevMeshes]
                        prevOutputNames = [coordsName, meshesName]
                    else:  # Model workflow outputs
                        # failedMeshesName = DynModelWfOuts.failedMeshes.name
                        # prevFailedMeshes = getattr(self.protocol, failedMeshesName, None)
                        prevOutputs = [prevCoords]
                        prevOutputNames = [coordsName]

                    outputsDict = {}
                    if outCoords:
                        outputsDict[DynModelWfOuts.fixedCoordinates.name] = outCoords
                    if outMeshes:
                        outputsDict[DynModelWfOuts.fixedMeshes.name] = outMeshes
                    for key, val in zip(prevOutputNames, prevOutputs):
                        if val:
                            outputsDict[key] = val

                    self.protocol._defineOutputs(**outputsDict)
                    for key, val in outputsDict.items():
                        if val:
                            self.protocol._defineSourceRelation(precedentsPointer, val)
            return views

    @staticmethod
    def _getModelType(coord):
        """Map the Dynamo model names into the protocol encoding model value"""
        # If more than one model of the same type, they're stored as modelName_num
        modelName = getattr(coord, '_dynModelName', None)
        return dynModelsDict.get(modelName.get().split('_')[0] if modelName is not None else MB_GENERAL)
