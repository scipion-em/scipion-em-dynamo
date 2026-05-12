# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *              Scipion Team (scipion@cnb.csic.es)
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
import glob
from enum import Enum
from os.path import basename, dirname, join, normpath

from pwem.emlib.image.image_readers import EmImageReader
from pyworkflow.object import Float
from pyworkflow.protocol.params import PointerParam, FloatParam
from pyworkflow.utils.path import createLink, makePath
from pwem.emlib.image import ImageHandler
from tomo.protocols.protocol_base import ProtTomoImportFiles
from tomo.objects import SubTomogram, SetOfSubTomograms, SetOfCoordinates3D
from .protocol_base_dynamo import DynamoProtocolBase, IN_TOMOS
from ..convert import dynTableLine2Subtomo


class DynImportSubtomosOuts(Enum):
    coordinates = SetOfCoordinates3D
    subtomograms = SetOfSubTomograms


class DynamoImportSubtomos(ProtTomoImportFiles, DynamoProtocolBase):
    """
    Imports subtomograms and their associated metadata from Dynamo table files
    into a Scipion tomography workflow environment. The protocol is intended to
    organize particle-centered subtomographic datasets generated during Dynamo
    processing and transform them into structured objects that can be directly
    used for downstream subtomogram averaging, classification, alignment, or
    visualization tasks.

    AI Generated:

    Import Subtomograms From Dynamo Tables (DynamoImportSubtomos) — User Manual
        Overview

        The Import Subtomograms From Dynamo Tables protocol is designed to
        incorporate subtomogram datasets generated in Dynamo into Scipion while
        preserving their spatial and experimental context. In subtomogram
        averaging workflows, particles are often extracted from tomograms and
        stored together with metadata describing their coordinates, orientations,
        and relationships to the parent tomograms. This protocol provides a
        bridge between Dynamo-based particle extraction workflows and Scipion-
        based processing pipelines.

        From a biological perspective, this import step is essential because it
        preserves the correspondence between each subtomogram and the original
        cellular or structural environment from which it was extracted. This
        contextual information becomes especially important in studies involving
        membrane-associated complexes, viral assemblies, cytoskeletal networks,
        ribosomes in situ, or any heterogeneous macromolecular population
        distributed throughout a tomographic volume.

        Inputs and Dataset Organization

        The protocol expects Dynamo table files together with the associated
        subtomogram volumes generated during particle extraction. Each table
        typically represents the metadata associated with a tomogram and contains
        the positional and orientational information required to reconstruct the
        relationship between extracted particles and their parent volume.

        Optionally, users may provide the original tomograms used during
        extraction. When these tomograms are available, the protocol establishes
        explicit spatial relationships between the imported subtomograms and
        their source tomograms. This improves interoperability with downstream
        Scipion tools that rely on coordinate information, visualization of
        particle distributions, or contextual biological interpretation.

        If tomograms are not available, the protocol can still import the
        subtomograms independently as long as an appropriate sampling rate is
        provided. This mode is useful when working with archived particle
        datasets, collaborative exchanges, or partial processing workflows where
        the original tomograms are unavailable.

        Sampling Rate and Physical Interpretation

        Correct sampling rate assignment is biologically critical because it
        defines the physical scale of the imported subtomograms. Accurate voxel
        size information ensures that structural dimensions, alignment searches,
        masking procedures, and resolution estimations remain physically
        meaningful throughout downstream processing.

        When tomograms are supplied, the sampling information is inherited
        automatically, reducing the risk of inconsistencies between extracted
        particles and their parent volumes. When importing standalone
        subtomograms, users should carefully verify the voxel size used during
        extraction to avoid scaling errors that could compromise alignment,
        averaging, or interpretation.

        Coordinate Preservation and Tomographic Context

        One of the major advantages of this protocol is the preservation of
        three-dimensional particle coordinates. Maintaining this information
        allows researchers to revisit the spatial organization of macromolecular
        complexes inside the native cellular environment after subtomogram
        processing has been completed.

        In biological studies, coordinate preservation enables analyses such as
        mapping ribosome distributions, studying membrane curvature dependencies,
        identifying lattice organization in viral capsids, or correlating
        molecular states with local ultrastructural environments. This spatial
        continuity is particularly important in cryo-electron tomography because
        biological interpretation often depends not only on particle structure
        but also on particle localization within the cell or specimen.

        File Conversion and Compatibility

        The protocol supports heterogeneous subtomogram file formats and prepares
        them for standardized processing within Scipion. This compatibility layer
        simplifies integration of datasets generated across different Dynamo
        workflows or legacy processing environments.

        From a practical perspective, this flexibility is useful when combining
        datasets originating from different laboratories, microscope facilities,
        or historical projects. Standardized import reduces downstream
        compatibility issues and facilitates reproducible processing pipelines.

        Outputs and Their Interpretation

        After execution, the protocol generates a structured set of imported
        subtomograms together with their associated metadata. When tomograms are
        provided, an additional set of three-dimensional coordinates is produced,
        linking every particle to its original location in the tomographic
        volume.

        The resulting dataset can be directly used for subtomogram alignment,
        classification, averaging, focused refinement, or visualization. Because
        particle identities and spatial relationships are preserved, downstream
        analyses remain biologically interpretable and traceable to the original
        tomographic experiment.

        Practical Recommendations

        In routine cryo-electron tomography workflows, it is highly recommended
        to import subtomograms together with their parent tomograms whenever
        possible. This preserves the complete experimental context and maximizes
        compatibility with downstream protocols requiring coordinate information.

        Users should verify that all subtomograms belonging to a dataset share a
        consistent voxel size and extraction convention before import. Mixing
        particles generated under incompatible extraction parameters may lead to
        unreliable alignments or biologically inconsistent averages.

        When working with collaborative datasets or archived material, careful
        validation of metadata consistency is especially important. Incorrect
        coordinate systems, mismatched sampling rates, or incomplete metadata can
        significantly affect downstream structural interpretation.

        Final Perspective

        Importing subtomograms is more than a technical conversion step. It is
        the process that reconnects extracted particle data with the biological
        specimen from which it originated. Proper preservation of spatial
        metadata, voxel size information, and tomographic relationships is
        essential for reliable subtomogram averaging and for meaningful
        interpretation of molecular organization inside native cellular
        environments.
    """

    _label = 'import subtomograms from tbl files'
    _possibleOutputs = DynImportSubtomosOuts

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.matchingFiles = None
        self.tblTomoDict = None
        self.tblSubtomoFilesDict = {}
        self.sRate = Float()

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        super()._defineImportParams(form)
        form.addParam(IN_TOMOS, PointerParam,
                      pointerClass='SetOfTomograms',
                      allowsNull=True,
                      label="Tomograms (opt.)",
                      help="If not provided, the subtomograms won't be referred to any tomogram. "
                           "If provided, the sampling rate value will be read from them.")
        form.addParam('samplingRate', FloatParam,
                      condition='not %s' % IN_TOMOS,
                      allowsNull=True,
                      label='Sampling rate [Å/px] (opt.)')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep,
                                 needsGPU=False)
        self._insertFunctionStep(self.importSubTomogramsStep,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        tomograms = self.getInTomos()
        self.matchingTblFiles = self.getMatchFiles()
        # If tomograms provided, the sampling rate value will be read from them. If not, the user
        # will have to introduce it manually
        if tomograms:
            self.tblTomoDict = {tblFile: tomo.clone() for tblFile, tomo in zip(self.matchingTblFiles, tomograms)}
            self.sRate.set(tomograms.getSamplingRate())
        else:
            self.tblTomoDict = {tblFile: None for tblFile in self.matchingTblFiles}
            self.sRate.set(self.samplingRate.get())

    def convertInputStep(self):
        # Check the subtomograms file extension in the first directory expected to contain some (the directory in which
        # the corresponding tbl file is contained)
        for tblFile in self.matchingTblFiles:
            filesDir = dirname(tblFile)
            dirBaseName = self._getExtraPath(basename(normpath(filesDir)))
            makePath(dirBaseName)
            # Get the subtomo from the current tbl corresponding directory excluding the tbl file
            subtomoFiles = list(set(glob.glob(join(filesDir, '*'))) - set(glob.glob(join(filesDir, '*.tbl'))))
            files2store = []
            if subtomoFiles[0].endswith('.mrc'):
                # Create links
                for subtomoFile in sorted(subtomoFiles):
                    newFileName = join(dirBaseName, basename(subtomoFile))
                    createLink(subtomoFile, newFileName)
                    files2store.append(newFileName)
            else:
                # Convert the files
                for subtomoFile in subtomoFiles:
                    mrcTomoFile = join(dirBaseName, basename(subtomoFile) + '.mrc')
                    emFileHandler = EmImageReader()
                    emFileHandler.emToMrc(subtomoFile, mrcTomoFile)
                    files2store.append(mrcTomoFile)
            self.tblSubtomoFilesDict[tblFile] = files2store

    def importSubTomogramsStep(self):
        imgh = ImageHandler()
        samplingRate = self.sRate.get()
        coordSet = None
        subtomoSet = SetOfSubTomograms.create(self._getPath(), template='subtomograms%s.sqlite')
        subtomoSet.setSamplingRate(samplingRate)
        subtomo = SubTomogram()
        subtomo.setSamplingRate(samplingRate)
        tomograms = self.getInTomos()
        if tomograms:
            coordSet = SetOfCoordinates3D.create(self._getPath(), template='coordinates3d%s.sqlite')
            coordSet.setPrecedents(tomograms)
            coordSet.setSamplingRate(samplingRate)
            coordSet.setBoxSize(20)
        for tblFile, tomo in self.tblTomoDict.items():
            with open(tblFile, 'r') as fhTable:
                lines = fhTable.readlines()
                # The subtomograms files are generated in the same directory as the .tbl file, one for each tomogram
                for line, fileName in zip(lines, self.tblSubtomoFilesDict[tblFile]):
                    _, _, _, n = imgh.getDimensions(fileName)
                    self._fillSubtomogram(line, subtomo, subtomoSet, fileName, tomo=tomo, coordSet=coordSet)
        # Needs to be registered before assigning it to the set of subtomograms
        self._defineOutputs(**{self._possibleOutputs.coordinates.name: coordSet})
        if tomograms:
            subtomoSet.setCoordinates3D(coordSet)
        self._defineOutputs(**{self._possibleOutputs.subtomograms.name: subtomoSet})
        if tomograms:
            self._defineSourceRelation(tomograms, subtomoSet)

    @staticmethod
    def _fillSubtomogram(line, subtomo, subtomoSet, newFileName, tomo=None, coordSet=None):
        """ adds a subtomogram to a set """
        subtomo.cleanObjId()
        subtomo.setFileName(newFileName)
        dynTableLine2Subtomo(line, subtomo, subtomoSet=subtomoSet, tomo=tomo, coordSet=coordSet)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        matchingFiles = self.getMatchFiles()
        tblFilesNotFoundMsg = 'No Dynamo .tbl files were detected in the introduced path/s.'
        if matchingFiles:
            if not matchingFiles[0].endswith('.tbl'):
                errors.append(tblFilesNotFoundMsg)
        else:
            errors.append(tblFilesNotFoundMsg)
        if not self.getInTomos() and not self.samplingRate.get():
            errors.append('If tomograms are not provided, a sampling rate value is required.')
        return errors
