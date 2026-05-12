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
import datetime
from enum import Enum
from os.path import abspath
from dynamo.protocols.protocol_base_dynamo import DynamoProtocolBase
from dynamo.utils import createBoxingOutputObjects, getDynamoModels, getNewestModelModDate
from dynamo.viewers.DynamoTomoProvider import DynamoTomogramProvider
import pyworkflow.utils as pwutils
from pyworkflow.protocol import LEVEL_ADVANCED
from pyworkflow.protocol.params import IntParam, BooleanParam
from pyworkflow.utils.properties import Message
from pyworkflow.gui.dialog import askYesNo
from tomo.objects import SetOfCoordinates3D, SetOfMeshes
from tomo.protocols import ProtTomoPicking
from dynamo import VLL_FILE, CATALOG_BASENAME
from dynamo.viewers.views_tkinter_tree import DynamoTomoDialog


class DynPickingOuts(Enum):
    coordinates = SetOfCoordinates3D()
    meshes = SetOfMeshes()


class DynamoBoxing(ProtTomoPicking, DynamoProtocolBase):
    """
    Provides an interactive environment for manual vectorial picking of
    structures inside tomograms using Dynamo visualization and annotation
    tools. The protocol enables users to define trajectories, surfaces, or
    particle distributions directly within tomographic volumes and transform
    these annotations into biologically meaningful coordinate datasets.

    AI Generated:

    Vectorial Picking (DynamoBoxing) - User Manual
        Overview

        The Vectorial Picking protocol enables interactive manual annotation
        of tomographic data using the Dynamo graphical environment. Its main
        purpose is to assist users in identifying biological structures,
        defining trajectories, and generating coordinate information directly
        from cryo-electron tomography datasets. This workflow is especially
        useful when studying elongated assemblies, membrane-associated
        complexes, cytoskeletal filaments, or curved biological surfaces that
        cannot be easily described through isolated particle picking alone.

        In practical cryo-ET studies, many biological structures follow
        continuous spatial organizations rather than appearing as discrete
        independent particles. Examples include actin filaments,
        microtubules, membrane tubules, viral lattices, or repeating
        macromolecular assemblies distributed along curved geometries. This
        protocol allows researchers to manually define these spatial patterns
        and transform them into coordinate systems suitable for downstream
        subtomogram extraction and averaging workflows.

        Interactive Picking Workflow

        The protocol operates through an interactive Dynamo visualization
        session in which users manually inspect tomograms and define models
        directly within the volumetric data. Previously created annotations
        can be reloaded automatically, allowing iterative refinement over
        multiple sessions and facilitating long-term curation of complex
        biological datasets.

        During the annotation process, users can generate one or more models
        for each tomogram. These models may represent biological trajectories,
        surfaces, tubular arrangements, filament traces, or other geometric
        organizations relevant to the experiment. The workflow is highly
        flexible and supports exploratory interpretation of crowded cellular
        environments where automated picking may be unreliable.

        Biological Importance of Vectorial Picking

        Manual vectorial picking is particularly important in situations where
        structural organization follows continuous or spatially constrained
        patterns. Standard particle picking approaches are often insufficient
        for describing filamentous systems, membrane-bound assemblies, or
        repeating structures distributed along curved cellular geometries.

        By defining biological trajectories directly within tomograms,
        researchers can preserve contextual information about spatial
        organization, curvature, polarity, and relative orientation. This
        information becomes especially valuable in studies involving cellular
        ultrastructure, cytoskeletal dynamics, membrane remodeling, or
        supramolecular assemblies embedded within native environments.

        In many biological workflows, vectorial picking serves as the bridge
        between visual interpretation and quantitative structural analysis.
        The resulting coordinates can later support subtomogram averaging,
        classification, alignment, or geometric measurements.

        Meshes and Coordinate Generation

        The protocol produces geometric representations derived from the user
        annotations together with coordinate datasets suitable for downstream
        processing. These outputs preserve the biological organization
        identified during manual inspection and provide structured spatial
        information that can be reused throughout the tomography workflow.

        Depending on the annotation strategy, the generated coordinates may
        represent interpolated particle positions distributed along curves,
        surfaces, or trajectories. This is particularly useful for elongated
        assemblies where biological particles follow repeating arrangements
        rather than isolated spatial positions.

        The resulting meshes can also assist in visualization and contextual
        interpretation of the tomographic scene. In cellular cryo-ET studies,
        these geometric outputs frequently help researchers understand spatial
        organization within crowded intracellular environments.

        Particle Size and Sampling Considerations

        The protocol allows users to define an approximate particle size
        associated with the generated coordinates. This parameter influences
        how extracted regions are interpreted during downstream subtomogram
        processing and should reflect the expected dimensions of the
        biological target.

        Choosing an appropriate particle size is biologically important
        because values that are too small may truncate relevant structural
        features, while excessively large values may include neighboring
        densities or unrelated cellular material. For filamentous or
        membrane-associated systems, the selected size should adequately cover
        the repeating structural motif of interest.

        Iterative Annotation and Curation

        One of the strengths of the protocol is its support for iterative
        annotation workflows. Researchers can reopen datasets, refine
        trajectories, correct annotations, and progressively improve the
        biological interpretation of difficult tomograms. This iterative
        strategy is especially valuable in challenging cellular datasets where
        structural identification evolves during analysis.

        In collaborative environments, the ability to preserve and revisit
        annotation models also facilitates quality control and reproducibility.
        Expert users may curate biologically meaningful trajectories while
        less experienced users contribute preliminary annotations for later
        refinement.

        Outputs and Their Interpretation

        The protocol produces annotated geometric models together with
        coordinate datasets representing biologically relevant positions
        inside the tomograms. These outputs are designed to integrate
        naturally with downstream subtomogram extraction and averaging
        workflows.

        The coordinate outputs preserve both spatial organization and
        orientation information derived from the manual annotations. This is
        particularly important for studies involving directional assemblies,
        filament polarity, membrane curvature, or ordered macromolecular
        arrays.

        The generated meshes additionally provide a convenient visual summary
        of the annotated biological structures and may support interpretation,
        presentation, or further geometric analysis.

        Practical Recommendations

        In routine cryo-ET workflows, it is generally advisable to begin with
        exploratory annotation sessions in order to understand the structural
        organization of the tomograms before attempting large-scale particle
        extraction. Careful visual inspection often reveals biological
        heterogeneity, curvature, or spatial constraints that automated
        methods may overlook.

        Users working with filamentous systems should maintain smooth and
        biologically consistent trajectories during annotation. For membrane
        systems or curved assemblies, preserving continuity and contextual
        interpretation is often more important than maximizing coordinate
        density.

        Iterative refinement of annotations is strongly recommended,
        particularly in noisy tomograms or crowded cellular environments.
        Revisiting annotations after preliminary averaging or classification
        can substantially improve biological accuracy.

        Final Perspective

        Manual vectorial picking remains one of the most biologically
        informative annotation strategies in cryo-electron tomography because
        it combines expert structural interpretation with quantitative spatial
        analysis. By enabling researchers to define trajectories, surfaces,
        and ordered spatial arrangements directly inside tomograms, the
        protocol supports biologically meaningful extraction of structural
        information from highly complex cellular datasets.
    """

    _label = 'vectorial picking'
    _possibleOutputs = DynPickingOuts

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.dlg = None
        self.dynModelsPathDict = {}  # Used to store the path where the corresponding models to a tomo are stored

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        ProtTomoPicking._defineParams(self, form)

        form.addParam('boxSize', IntParam, expertLevel=LEVEL_ADVANCED,
                      default=32, label="Box Size")
        form.addParam('deleteGenMFiles', BooleanParam,
                      default=True,
                      label='Remove the .m files generated after the execution?',
                      expertLevel=LEVEL_ADVANCED,
                      help='It can be useful for developers to check exactly what was .m files were generated by '
                           'Scipion and executed by Dynamo.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep,
                                 needsGPU=False)
        self._insertFunctionStep(self.launchDynamoBoxingStep,
                                 needsGPU = False,
                                 interactive=True)

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

    def launchDynamoBoxingStep(self):
        tomoList = []
        for tomo in self.inputTomograms.get().iterItems():
            tomogram = tomo.clone()
            tomoList.append(tomogram)

        tomoProvider = DynamoTomogramProvider(tomoList, self._getExtraPath(), "txt")
        dynamoDialogCallingTime = datetime.datetime.now()
        self.dlg = DynamoTomoDialog(None, self._getExtraPath(),
                                    provider=tomoProvider,
                                    calledFromViewer=False)

        modelList = getDynamoModels(self._getExtraPath())
        if modelList:
            # Check if the modification file of the newest model file is higher than the time capture right before
            # calling the Dynamo dialog. In that case, it means that some modification was carried out by the user from
            # it and we have to ask if the changes should be saved
            if dynamoDialogCallingTime < getNewestModelModDate(modelList):
                # Open dialog to request confirmation to create output
                import tkinter as tk
                if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, tk.Frame()):
                    self._createOutput()
                    # Delete the .m generated files if requested
                if self.deleteGenMFiles.get():
                    pwutils.cleanPattern(self._getExtraPath('*.m'))

    def _createOutput(self):
        precedentsPointer = self.inputTomograms
        meshes, outCoords = createBoxingOutputObjects(self, precedentsPointer, boxSize=self.boxSize.get())
        # Define outputs and relations
        outputDict = {self._possibleOutputs.meshes.name: meshes}
        if outCoords:
            outputDict[self._possibleOutputs.coordinates.name] = outCoords
        self._defineOutputs(**outputDict)
        self._defineSourceRelation(precedentsPointer, meshes)
        if outCoords:
            self._defineSourceRelation(precedentsPointer, outCoords)
        self._updateOutputSet(self._possibleOutputs.meshes.name, meshes, state=meshes.STREAM_CLOSED)

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

