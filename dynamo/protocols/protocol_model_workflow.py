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
from os import remove
from os.path import abspath, exists
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import IntParam, PointerParam, BooleanParam, LEVEL_ADVANCED
from tomo.objects import SetOfCoordinates3D, SetOfMeshes, Coordinate3D
from tomo.protocols import ProtTomoBase
from dynamo import Plugin, M_GENERAL_DES, M_GENERAL_WITH_BOXES_DES, M_GENERAL_NAME, M_SURFACE_NAME, \
    M_ELLIPSOIDAL_VESICLE_NAME, M_GENERAL_WITH_BOXES_NAME, M_SURFACE_DES, \
    M_SPH_VESICLE_NAME, M_VESICLE_DES, M_ELLIPSOIDAL_VESICLE_DES, M_MARKED_ELLIP_VESICLE_NAME, \
    M_MARKED_ELLIP_VESICLE_DES, \
    MB_BY_LEVELS, MB_ELLIPSOIDAL, MB_GENERAL, MB_GENERAL_BOXES, MB_VESICLE, MB_ELLIPSOIDAL_MARKED, \
    MODELS_NOT_PROCESSED_IN_MW, MODELS_ALLOWED_IN_MW_NAMES, M_VESICLE_NAME
from ..utils import genMCode4ReadAndSaveData, dynamoCroppingResults2Scipion, createSetOfOutputCoords, getCroppedFile

# Model types mapping
MODEL_CHOICES = [M_ELLIPSOIDAL_VESICLE_NAME, M_SURFACE_NAME, M_GENERAL_NAME]

# Model types encoding
M_VESICLE = 0
M_SURFACE = 1
M_GENERAL = 2

# Dynamo model names mapping
dynModelsDict = {
    MB_BY_LEVELS: M_SURFACE,
    MB_VESICLE: M_VESICLE,
    MB_ELLIPSOIDAL: M_VESICLE,
    MB_ELLIPSOIDAL_MARKED: M_VESICLE,
    MB_GENERAL: M_GENERAL,
    MB_GENERAL_BOXES: M_GENERAL,
}

# Attribute names of extended attributes that are specific of Dynamo
DYN_MODEL_NAME = '_dynModelName'
DYN_MODEL_FILE = '_dynModelFile'

# Keys for the failed models
TOMO_ID = 'tomoId'
MODEL_NAME = 'modelName'
MODEL_FILE = 'modelFile'
FAILED_MODEL_KEYS = [TOMO_ID, MODEL_NAME, MODEL_FILE]


class DynModelWfOuts(Enum):
    # Instantiation needed in case the of multiple outputs of the same type (overridden if not)
    coordinates = SetOfCoordinates3D()
    fixedCoordinates = SetOfCoordinates3D()
    failedMeshes = SetOfMeshes()
    fixedMeshes = SetOfMeshes()


class DynamoModelWorkflow(EMProtocol, ProtTomoBase):
    """
    Apply a model workflow to a SetOfMeshes generated by Dynamo Boxing protocol.
    This workflow will use the models created by the user to create the
    corresponding cropping meshes needed to extract the crop points.
    Considerations:
        1. The geometry will be automatically approximated from the clicked points.
        2. The meshes for the general models will be calculated treating them as surface models.
    """

    _label = 'model workflow'
    _devStatus = BETA
    _possibleOutputs = DynModelWfOuts

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.failedList = []

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMeshes', PointerParam,
                      pointerClass='SetOfMeshes',
                      label="Input Meshes",
                      important=True,
                      help="Input Meshes that will be used to create the cropping geometry and "
                           "to extract the crop points")
        form.addParam('boxSize', IntParam,
                      default=0,
                      label='Box size',
                      important=True)

        form.addSection('Model parameters')
        form.addLine('Accepted Dynamo models: description in help --> ', help=self._genModelsDescriptionMsg())
        group = form.addGroup('Mesh creation for depiction')
        group.addParam('meshParameter', IntParam,
                       default=5,
                       label='Mesh parameter',
                       help='Intended mesh parameter for the "mesh" that supports the depiction of the model. It '
                            'governs the number or triangles.')
        group.addParam('maxTr', IntParam,
                       default=100000,
                       label="Maximun number of triangles",
                       help='Maximum number of triangles allowed during generation of a depiction mesh')
        group.addParam('doRefineMesh', BooleanParam,
                       default=False,
                       label='Refine mesh?',
                       expertLevel=LEVEL_ADVANCED,
                       help='If set to Yes, it will refine both the mesh and the cropped mesh, which means that the '
                            'depiction grid will be subdivided (each triangle will generate four children). Hence, a '
                            'higher number of points will be generated, but the computation time will be significantly '
                            'increased.')
        group.addParam('subDivision', IntParam,
                       default=2,
                       expertLevel=LEVEL_ADVANCED,
                       condition='doRefineMesh',
                       label="Subdivision iterations",
                       help="Specify the number of times the Mesh geometry will be subdivided. This will increase the "
                            "number of triangles in the mesh, making it smoother. However, it will also increase the "
                            "number of cropping points")
        group = form.addGroup('Particle cropping')
        group.addParam('cropping', IntParam,
                       default=10,
                       label="Cropping parameter",
                       help='Intended mesh parameter for the "crop_mesh" that defined a cropping geometry on a surface')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        modelsDict = self.inputMeshes.get().getUniqueValues([Coordinate3D.TOMO_ID_ATTR,
                                                             DYN_MODEL_NAME,
                                                             DYN_MODEL_FILE,
                                                             '_groupId'])
        for tomoId, modelName, modelFile in zip(modelsDict[Coordinate3D.TOMO_ID_ATTR],
                                                modelsDict[DYN_MODEL_NAME],
                                                modelsDict[DYN_MODEL_FILE]):
            self._insertFunctionStep(self.applyWorkflowStep, tomoId, modelName, modelFile)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def applyWorkflowStep(self, tomoId, modelName, modelFile):
        commandsFile = self.writeMatlabFile(tomoId, modelName, modelFile)
        args = ' %s' % commandsFile
        try:
            self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())
        except:
            self.failedList.append(dict(zip(FAILED_MODEL_KEYS, [tomoId, modelName, modelFile])))

    def createOutputStep(self):
        croppedFile = getCroppedFile(self._getTmpPath())  # All the calculated mesh points will be in this file
        inputMeshes = self.inputMeshes.get()
        outCoords = None
        failedMeshes = None
        outputsDict = {}
        if exists(croppedFile):  # If all the models failed in the model workflow
            precedentsPointer = inputMeshes._precedentsPointer
            precedents = precedentsPointer.get()
            tomoList = [tomo.clone() for tomo in precedents]
            tomoFileDict = {abspath(tomo.getFileName()): tomo for tomo in tomoList}
            outCoords = createSetOfOutputCoords(self._getPath(), self._getExtraPath(), precedentsPointer,
                                                boxSize=self.boxSize.get())
            dynamoCroppingResults2Scipion(outCoords, croppedFile, tomoFileDict)
            # Remove previous file to avoid data repetition because of the append mode
            remove(croppedFile)
            outputsDict[self._possibleOutputs.coordinates.name] = outCoords

        # Create a set with the failed meshes if there are any, so the user can correct them
        if self.failedList:
            failedMeshes = SetOfMeshes.create(self._getPath(), template='meshes%s.sqlite', suffix='failed')
            failedMeshes.copyInfo(inputMeshes)
            failedMeshes._dynCatalogue = inputMeshes._dynCatalogue
            for d in self.failedList:
                whereCond = "_tomoId=='%s' AND _dynModelName='%s' AND _dynModelFile=='%s'" % \
                            (d[TOMO_ID], d[MODEL_NAME], d[MODEL_FILE])
                for clickedPoint in inputMeshes.iterItems(where=whereCond):
                    failedMeshes.append(clickedPoint)
            outputsDict[self._possibleOutputs.failedMeshes.name] = failedMeshes

        # Define outputs and relations
        self._defineOutputs(**outputsDict)
        if outCoords:
            self._defineSourceRelation(self.inputMeshes, outCoords)
            self._updateOutputSet(self._possibleOutputs.coordinates.name, outCoords, state=outCoords.STREAM_CLOSED)
        if failedMeshes:
            self._defineSourceRelation(self.inputMeshes, failedMeshes)
            self._updateOutputSet(self._possibleOutputs.failedMeshes.name, failedMeshes,
                                  state=failedMeshes.STREAM_CLOSED)

    # --------------------------- DEFINE utils functions ----------------------
    def writeMatlabFile(self, tomoId, modelName, modelFile):
        content = ''
        codeFilePath = self._getExtraPath('modelWf_%s_%s.m' % (tomoId, modelName))
        modelType = self._getModelType(modelName)
        if modelType == M_VESICLE:
            content = self.genVesicleCmdFileContents()
        elif modelType == M_SURFACE:
            content = self.genSCmdFileContents()
        elif modelType == M_GENERAL:
            # Change its type to surface model and process it
            content = self.genGen2SurfCmdFileContents()
        content = genMCode4ReadAndSaveData(self._getTmpPath(), modelFile, savePicked=False, saveCropped=True,
                                           modelWfCode=content)
        with open(codeFilePath, 'w') as codeFid:
            codeFid.write(content)
        return codeFilePath

    def _genCommonModelWfSteps(self):
        """See
        https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Example_of_membrane_model_workflow_through_the_command_line"""
        contentMWf = "m.mesh_parameter = %i\n" % self.meshParameter.get()
        contentMWf += "m.crop_mesh_parameter = %i\n" % self.cropping.get()
        contentMWf += "m.mesh_maximum_triangles = %i\n" % self.maxTr.get()
        contentMWf = "m.subdivision_iterations=%i\n" % self.subDivision.get()
        contentMWf += "m.controlUpdate()\n"
        contentMWf += "m.createMesh()\n"
        if self.doRefineMesh.get():
            contentMWf += "m.refineMesh()\n"
        contentMWf += "m.createCropMesh()\n"
        if self.doRefineMesh.get():
            contentMWf += "m.refineCropMesh()\n"
        contentMWf += "m.updateCrop()\n"
        contentMWf += "m.grepTable()\n"
        return contentMWf

    def genVesicleCmdFileContents(self):
        # Let Dynamo approximate the geometry based on the points annotated in the boxing protocol
        contentMWf = "m.approximateGeometryFromPoints()\n"
        contentMWf += self._genCommonModelWfSteps()
        return contentMWf

    def genSCmdFileContents(self):
        return self._genCommonModelWfSteps()

    def genGen2SurfCmdFileContents(self):
        # Change model type from general to surface
        contentMWf = "m = model.changeType(m, '%s')\n" % 'membraneByLevels'
        contentMWf += "zCoords = m.points(:, 3)\n"
        contentMWf += "zUVals = unique(zCoords)\n"
        contentMWf += "nParticles = length(zCoords)\n"
        contentMWf += "groupLabels = zeros(1, nParticles)\n"
        contentMWf += "for j=1:length(zUVals)\n"
        contentMWf += "currentZ = zUVals(i)\n"
        contentMWf += "groupLabels(zCoords == currentZ) = j\n"
        contentMWf += "end\n"
        contentMWf += "m.group_labels = groupLabels\n"
        contentMWf += "m.last_group_label = i\n"
        # Mesh creation steps
        contentMWf += self.genSCmdFileContents()
        # Format and write the output data in a text file that will be read in the step create output
        return contentMWf

    @staticmethod
    def _genModelsNotationMsg():
        modelsHelp = '[G]: %s,  ' % M_GENERAL_NAME
        modelsHelp += '[S]: %s,  ' % M_SURFACE_NAME
        modelsHelp += '[V]: %s' % M_VESICLE_NAME
        return modelsHelp

    @staticmethod
    def _genModelsDescriptionMsg():
        modelsHelp = '*%s* (NOTE: they will be converted into surface model to generate the mesh):\n\n' % M_GENERAL_NAME.upper()
        modelsHelp += '\t1) *%s*: %s\n\n' % (M_GENERAL_NAME, M_GENERAL_DES)
        modelsHelp += '\t2) *%s*: %s\n\n\n' % (M_GENERAL_WITH_BOXES_NAME, M_GENERAL_WITH_BOXES_DES)
        modelsHelp += '*%s*:\n%s\n\n\n' % (M_SURFACE_NAME.upper(), M_SURFACE_DES)
        modelsHelp += '*%s*:\n\n' % M_VESICLE_NAME.upper()
        modelsHelp += '\t1) *%s*: %s\n\n' % (M_SPH_VESICLE_NAME, M_VESICLE_DES)
        modelsHelp += '\t2) *%s*: %s\n\n' % (M_ELLIPSOIDAL_VESICLE_NAME, M_ELLIPSOIDAL_VESICLE_DES)
        modelsHelp += '\t3) *%s*: %s\n\n\n' % (M_MARKED_ELLIP_VESICLE_NAME, M_MARKED_ELLIP_VESICLE_DES)
        return modelsHelp

    @staticmethod
    def _getModelType(modelName):
        # Map the Dynamo model names into the protocol encoding model values
        return dynModelsDict[modelName.split('_')[0]]  # If more than one model of the same type, they're stored as modelName_num

    def getMeshResultFile(self, tomoId):
        return abspath(self._getExtraPath('%s.txt' % tomoId))

    # --------------------------- DEFINE INFO functions ----------------------
    def _summary(self):
        summary = []
        if self.getOutputsSize() >= 1:
            for _, outCoords in self.iterOutputAttributes():
                summary.append("Output *%s*:" % outCoords.getNameId().split('.')[1])
                summary.append("    * Particle box size: *%s*" % self.boxSize.get())
                summary.append("    * Coordinates defined by geometry: *%s*" %
                               outCoords.getSize())
        else:
            summary.append("Output coordinates not ready yet.")
        return summary

    def _validate(self):
        errorMsg = []
        # Only sets of meshes generated using the Dynamo picking protocol are accepted (they must contain
        # an attribute named '_dynCatalogue')
        if not getattr(self.inputMeshes.get(), '_dynCatalogue', None):
            errorMsg.append('Only sets of meshes generated using the Dynamo picking protocol are accepted')
        return errorMsg

    def _warnings(self):
        warnMsg = []
        presentModelList = self.inputMeshes.get().getUniqueValues([DYN_MODEL_NAME])
        for presentModel in presentModelList:
            if presentModel in MODELS_NOT_PROCESSED_IN_MW:
                warnMsg.append('Some of the models provided are not allowed in this protocol. Allowed models are:\n%s'
                               % '\n - '.join(MODELS_ALLOWED_IN_MW_NAMES.insert(0, '')))
                break
        return warnMsg
