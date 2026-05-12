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
import logging
from enum import Enum
from os import remove
from os.path import abspath, exists
from pwem.protocols import EMProtocol
from pyworkflow.protocol import IntParam, PointerParam, BooleanParam, LEVEL_ADVANCED, STEPS_PARALLEL
from pyworkflow.utils import Message, cyanStr
from tomo.objects import SetOfCoordinates3D, SetOfMeshes, Coordinate3D
from tomo.protocols import ProtTomoBase
from dynamo import Plugin, M_GENERAL_DES, M_GENERAL_WITH_BOXES_DES, M_GENERAL_NAME, M_SURFACE_NAME, \
    M_ELLIPSOIDAL_VESICLE_NAME, M_GENERAL_WITH_BOXES_NAME, M_SURFACE_DES, \
    M_SPH_VESICLE_NAME, M_VESICLE_DES, M_ELLIPSOIDAL_VESICLE_DES, M_MARKED_ELLIP_VESICLE_NAME, \
    M_MARKED_ELLIP_VESICLE_DES, \
    MB_BY_LEVELS, MB_ELLIPSOIDAL, MB_GENERAL, MB_GENERAL_BOXES, MB_VESICLE, MB_ELLIPSOIDAL_MARKED, \
    MODELS_NOT_PROCESSED_IN_MW, M_VESICLE_NAME
from ..utils import genMCode4ReadAndSaveData, dynamoCroppingResults2Scipion, createSetOfOutputCoords, getCroppedFile

logger = logging.getLogger(__name__)

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
    Applies Dynamo model workflows to tomography-derived meshes in order to
    generate cropping geometries and extract particle coordinates suitable
    for downstream subtomogram analysis and structural interpretation.

    AI Generated:

    Dynamo Model Workflow (DynamoModelWorkflow) - User Manual
        Overview

        The Dynamo Model Workflow protocol processes collections of meshes
        generated from Dynamo boxing procedures and converts them into
        biologically meaningful cropping geometries and particle coordinates.
        Its primary purpose is to transform user-defined structural models
        into organized sampling regions that can later be used for particle
        extraction, subtomogram averaging, or geometric analysis within
        tomographic datasets.

        In practical cryo-electron tomography workflows, this protocol is
        especially useful when particles are distributed along membranes,
        vesicles, curved surfaces, or manually annotated structural features.
        Instead of relying on isolated coordinate picking, the protocol
        interprets the biological geometry represented by the meshes and
        generates coherent cropping positions that follow the modeled
        structure.

        Biological Context and Intended Usage

        Many biological assemblies in tomography are not isolated particles
        but organized systems embedded within membranes or extended cellular
        environments. Membrane proteins, vesicular systems, organelle
        surfaces, and filament-associated complexes often require geometric
        interpretation before particles can be extracted consistently.

        This protocol addresses that need by allowing the user to define
        structural models through Dynamo-based annotations and then generate
        refined geometrical representations suitable for systematic sampling.
        The resulting coordinates follow the topology of the modeled object,
        enabling extraction strategies that preserve biological organization
        and orientation consistency.

        The workflow is particularly valuable for membrane-associated studies,
        where particle positions are expected to lie on curved surfaces rather
        than in arbitrary volumetric regions.

        Accepted Model Types

        The protocol supports several categories of Dynamo models, including
        general models, surface models, and vesicle-oriented models. Each
        model type represents a different biological interpretation of the
        annotated points.

        Surface models are intended for structures that naturally form
        continuous biological boundaries, such as membranes or layered
        assemblies. Vesicle models are optimized for approximately spherical
        or ellipsoidal biological objects, allowing the workflow to infer a
        smooth geometry from sparse annotations. General models are accepted
        as well and are internally interpreted as surface-like geometries so
        they can participate in mesh generation and coordinate extraction.

        From a biological perspective, choosing the correct model type is
        important because it determines how the geometry is reconstructed and
        how particle positions are distributed along the structure.

        Geometry Approximation and Mesh Generation

        One of the central concepts of the protocol is geometric
        approximation. The workflow interprets the annotated points as
        evidence of an underlying biological surface and reconstructs a mesh
        representation that approximates that structure.

        This approach is especially useful when the biological object cannot
        be represented by a simple regular geometry. Curved membranes,
        irregular vesicles, and partially segmented cellular structures can
        all be approximated through mesh reconstruction.

        The user can control the density and smoothness of the generated
        meshes through mesh-related parameters. Lower mesh densities generally
        produce faster computations and simpler geometries, while denser
        meshes provide smoother and more accurate representations of complex
        biological surfaces.

        Increasing mesh refinement improves geometric detail but also
        generates substantially larger numbers of cropping points. In
        practical workflows, users should balance geometric precision against
        computational cost and downstream dataset size.

        Refinement and Surface Sampling

        Optional mesh refinement allows the protocol to subdivide the
        geometry into finer sampling regions. This operation increases the
        number of surface elements and therefore increases the spatial density
        of extracted coordinates.

        Biologically, refinement becomes important when studying densely
        distributed particles or highly curved membranes where coarse
        approximations may fail to capture local geometry accurately.
        However, excessive refinement may create unnecessarily large particle
        datasets and increase redundancy in regions with limited structural
        variability.

        In many practical cryo-ET studies, moderate refinement is sufficient
        to preserve biological continuity while keeping extraction manageable.

        Cropping Geometry and Particle Extraction

        After generating the mesh representation, the workflow creates
        cropping geometries that define where particle coordinates will be
        extracted. These cropping geometries follow the reconstructed surface
        and provide systematic sampling positions across the modeled
        structure.

        The resulting coordinates can later be used for subtomogram
        extraction, alignment, classification, or averaging workflows. Since
        the coordinates originate from biologically informed geometries, the
        extracted particles often preserve meaningful spatial organization.

        The protocol also stores the particle box size associated with the
        extraction process, ensuring consistency between coordinate generation
        and subsequent subtomogram processing.

        Handling Failed Models

        Some meshes may fail during workflow execution because of incomplete
        annotations, incompatible geometries, or unstable model definitions.
        Instead of discarding these cases silently, the protocol preserves
        failed meshes in a dedicated output collection.

        This behavior is particularly important in biological projects where
        only a subset of tomograms or structures may require manual
        correction. By isolating problematic meshes, the workflow allows the
        user to inspect and refine annotations without interrupting successful
        processing of the remaining dataset.

        Validation and Data Consistency

        The protocol expects meshes generated through Dynamo-compatible
        boxing workflows. This restriction ensures that all required metadata,
        model annotations, and geometric descriptors are available during
        processing.

        Consistency between tomograms, mesh annotations, and model
        definitions is critical for obtaining biologically meaningful
        coordinates. Incorrect or incomplete metadata may lead to invalid
        geometry reconstruction or extraction errors.

        Practical Recommendations

        For routine membrane-associated workflows, users should begin with
        moderate mesh density and limited refinement in order to evaluate the
        resulting geometry visually before generating very large coordinate
        datasets.

        Vesicle-oriented workflows generally benefit from well-distributed
        annotations that cover the entire biological structure. Sparse or
        uneven annotations may produce unstable surface approximations.

        When processing heterogeneous cellular environments, it is advisable
        to inspect failed meshes separately and verify whether the biological
        structure is adequately represented by the chosen model type.

        In downstream subtomogram averaging projects, maintaining consistent
        extraction geometry across tomograms is essential for preserving
        biological interpretability and improving alignment stability.

        Final Perspective

        The Dynamo Model Workflow protocol bridges manual geometric
        annotation and automated coordinate generation within cryo-electron
        tomography. Rather than treating particles as isolated points, it
        interprets them as components of larger biological surfaces and
        spatial organizations.

        For many tomographic studies involving membranes, vesicles, or
        extended cellular architectures, this geometry-aware extraction
        strategy provides a more biologically meaningful foundation for
        subtomogram analysis and structural interpretation.
    """

    _label = 'model workflow'
    _possibleOutputs = DynModelWfOuts
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.failedList = []

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
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
        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        modelsDict = self.inputMeshes.get().getUniqueValues([Coordinate3D.TOMO_ID_ATTR,
                                                             DYN_MODEL_NAME,
                                                             DYN_MODEL_FILE,
                                                             '_groupId'])
        pIdList = []
        for tomoId, modelName, modelFile in zip(modelsDict[Coordinate3D.TOMO_ID_ATTR],
                                                modelsDict[DYN_MODEL_NAME],
                                                modelsDict[DYN_MODEL_FILE]):
            wfId = self._insertFunctionStep(self.applyWorkflowStep, tomoId, modelName, modelFile,
                                     prerequisites=[],
                                     needsGPU=False)
            pIdList.append(wfId)
        self._insertFunctionStep(self.createOutputStep,
                                 prerequisites=pIdList,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def applyWorkflowStep(self, tomoId, modelName, modelFile):
        logger.info(cyanStr(f'===> {tomoId}: Running the model workflow:'))
        logger.info(cyanStr(f'======> Model name = {modelName}'))
        logger.info(cyanStr(f'======> Model file = {modelFile}'))
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

    def _genCommonModelWfSteps(self, isVesicle=False):
        """See
        https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Example_of_membrane_model_workflow_through_the_command_line"""
        contentMWf = "m.mesh_parameter = %i\n" % self.meshParameter.get()
        contentMWf += "m.crop_mesh_parameter = %i\n" % self.cropping.get()
        contentMWf += "m.mesh_maximum_triangles = %i\n" % self.maxTr.get()
        contentMWf = "m.subdivision_iterations=%i\n" % self.subDivision.get()
        if not isVesicle:
            contentMWf += "m.controlUpdate()\n"  # This step fails for vesicle model on the Dynamo side
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
        contentMWf += self._genCommonModelWfSteps(isVesicle=True)
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
        preMsg = 'Some of the models provided are not allowed in this protocol. Allowed models are:\n'
        for presentModel in presentModelList:
            if presentModel in MODELS_NOT_PROCESSED_IN_MW:
                warnMsg.append(f'{preMsg}{MODELS_NOT_PROCESSED_IN_MW}')
                break
        return warnMsg
