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
import copy
import glob
import logging
from enum import Enum
from os.path import abspath, join
from typing import List
import mrcfile
from dynamo.protocols.protocol_base_dynamo import IN_TOMOS, IN_COORDS, DynamoProtocolBase
from dynamo.utils import getCatalogFile
from pwem.objects import Transform
from pyworkflow.object import Boolean
from pyworkflow.protocol import PointerParam, EnumParam, IntParam, BooleanParam, STEPS_PARALLEL
from pyworkflow.utils import removeExt, Message, makePath, cyanStr, redStr
from tomo.constants import BOTTOM_LEFT_CORNER, TR_DYNAMO
from tomo.objects import SetOfSubTomograms, SubTomogram
from dynamo import Plugin, VLL_FILE
from dynamo.convert import matrix2eulerAngles
from tomo.utils import scaleTrMatrixShifts

logger = logging.getLogger(__name__)
CROP_DIR = 'Crop'
LOG_FILE_NAME = 'log.txt'

# Tomogram type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1


class DynSubtomoExtractOuts(Enum):
    subtomograms = SetOfSubTomograms


class DynamoExtraction(DynamoProtocolBase):
    """
    Extracts subtomograms from tomographic datasets using coordinate
    information and Dynamo-based particle cropping workflows. The protocol
    generates individual 3D particles that can later be used for alignment,
    classification, averaging, or structural interpretation in subtomogram
    averaging pipelines.

    AI Generated:

    Subtomogram Extraction (DynamoExtraction) — User Manual
        Overview

        The Subtomogram Extraction protocol is designed to isolate local
        three-dimensional regions from tomograms using previously defined
        particle coordinates. In cryo-electron tomography workflows, this
        operation represents one of the central transitions between particle
        detection and downstream structural analysis. The extracted
        subtomograms become the working units for alignment, averaging,
        classification, and refinement.

        From a biological perspective, the protocol allows users to recover
        molecular complexes directly from crowded cellular environments or
        from purified in situ preparations. Each extracted subtomogram
        preserves the structural context and orientation associated with its
        original tomographic location, making this step essential for
        studying macromolecular organization inside native environments.

        Inputs and Coordinate Interpretation

        The protocol requires a set of three-dimensional coordinates that
        identify the positions of the particles of interest. These
        coordinates are typically generated during manual or automated
        particle picking workflows and are associated with one or more
        tomograms.

        By default, the extraction is performed from the same tomograms used
        during the coordinate generation stage. This is generally the safest
        and most biologically consistent option because it guarantees that
        coordinates and tomograms share the same spatial reference system.
        However, the protocol also supports extraction from alternative
        tomograms when users wish to use differently processed datasets,
        such as denoised tomograms, reconstructed versions with different
        filtering strategies, or tomograms generated at another binning
        level.

        When alternative tomograms are used, users should verify that the
        spatial correspondence between coordinates and tomograms remains
        accurate. Small mismatches in sampling rate, alignment, or geometry
        can produce incorrect particle localization and compromise the
        biological validity of the extracted subtomograms.

        Box Size and Biological Context

        One of the most important parameters is the extraction box size.
        This parameter determines the cubic region surrounding each particle
        that will be isolated from the tomogram. The chosen size should be
        large enough to fully contain the macromolecular complex and a small
        amount of surrounding contextual density, but not so large that it
        introduces excessive background noise or neighboring structures.

        In practical cryo-ET workflows, selecting an appropriate box size
        depends on the biological target. Small soluble proteins may require
        relatively compact boxes, whereas membrane-associated assemblies,
        ribosomes, viral particles, or large molecular machines may need
        substantially larger regions. If the extraction region extends
        beyond the tomogram boundaries, some particles may become invalid
        and excluded from the final output.

        Biological users should therefore inspect particle positions near
        tomogram edges carefully. Excessively large boxes can reduce the
        number of valid particles and increase storage and computational
        costs during downstream averaging procedures.

        Preservation of Orientations

        The protocol preserves the spatial orientation associated with each
        coordinate. This is especially important for workflows involving
        vectorial picking, membrane geometry estimation, or template-guided
        localization strategies where particle orientation carries
        biologically meaningful information.

        Maintaining consistent orientation metadata enables subsequent
        subtomogram alignment procedures to converge more efficiently and
        allows users to perform analyses related to membrane topology,
        filament directionality, or molecular organization inside cells.

        Contrast Inversion

        The protocol optionally supports contrast inversion after particle
        extraction. This option is important because different tomography
        reconstruction pipelines may represent particles either as dark
        densities on bright backgrounds or the opposite.

        Many downstream subtomogram averaging and single-particle software
        packages assume that biological densities appear bright over a dark
        background. Inverting the contrast ensures compatibility with these
        conventions and can improve visualization and interpretation during
        later processing stages.

        From a biological interpretation standpoint, contrast inversion does
        not alter structural information. It simply standardizes the density
        representation expected by subsequent software environments.

        Parallel Processing and Large Datasets

        Modern cryo-electron tomography projects frequently involve hundreds
        or thousands of particles distributed across many tomograms. The
        protocol is therefore designed to process multiple tomograms
        efficiently using parallel execution strategies.

        This capability becomes especially important in cellular tomography
        studies where large datasets are common and extraction can otherwise
        become computationally demanding. Efficient parallelization reduces
        processing time while maintaining reproducibility across all
        tomograms in the experiment.

        Outputs and Their Interpretation

        The primary output is a set of extracted subtomograms associated
        with the original particle coordinates and tomographic references.
        Each subtomogram preserves the spatial metadata needed for further
        alignment and averaging procedures.

        Biologically, the extracted particles should be interpreted as local
        volumetric observations of the targeted molecular complexes. Their
        quality depends strongly on the original tomogram reconstruction,
        coordinate precision, box size selection, and local structural
        heterogeneity.

        In some cases, particles located near tomogram boundaries or inside
        damaged regions may be excluded automatically. This behavior is
        generally beneficial because incomplete particles can negatively
        affect downstream averaging quality and introduce structural
        artifacts.

        Practical Recommendations

        In routine biological workflows, users should begin by selecting a
        conservative box size that comfortably contains the particle while
        minimizing surrounding clutter. Visual inspection of several
        extracted subtomograms is strongly recommended before proceeding to
        large-scale averaging or classification.

        When using coordinates derived from one tomogram set together with
        another tomogram source, users should carefully verify sampling rate
        consistency and spatial registration. Small geometric mismatches may
        not be immediately obvious but can significantly degrade alignment
        quality later in the workflow.

        Contrast inversion should only be enabled when necessary for
        compatibility with downstream software or visualization conventions.
        If uncertainty exists, inspecting a subset of extracted particles
        before large-scale processing is advisable.

        Final Perspective

        Subtomogram extraction is a foundational step in cryo-electron
        tomography because it transforms coordinate-based particle
        localization into analyzable structural data. The biological quality
        of all downstream subtomogram averaging results depends strongly on
        the accuracy and consistency achieved during this stage. Careful
        selection of extraction parameters, validation of coordinate
        integrity, and thoughtful interpretation of particle context are
        essential for obtaining reliable structural information from
        tomographic datasets.
    """

    _label = 'subtomogram extraction'
    _possibleOutputs = DynSubtomoExtractOuts
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.dynamoTomoIdDict = {}
        self.scaleFactor = None
        self.tomoTsIdDict = None
        self.coordsRemoved = Boolean()
        self.failedItems = []

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_COORDS, PointerParam,
                      label="Coordinates",
                      important=True,
                      pointerClass='SetOfCoordinates3D')
        form.addParam('tomoSource', EnumParam,
                      choices=['same as picking', 'other'],
                      default=SAME_AS_PICKING,
                      display=EnumParam.DISPLAY_HLIST,
                      label='Tomogram source',
                      help='By default the subtomograms will be extracted from the tomogram used in the picking '
                           'step ( _same as picking_ option ).\nIf you select _other_ option, you must provide '
                           'a different tomogram to extract from.\n*Note*: In the _other_ case, ensure that provided '
                           'tomogram and coordinates are related ')
        form.addParam(IN_TOMOS, PointerParam,
                      pointerClass='SetOfTomograms',
                      condition='tomoSource != %s' % SAME_AS_PICKING,
                      label='Tomograms',
                      help='Select the tomogram from which to extract.')
        form.addParam('boxSize', IntParam,
                      label='Box size',
                      help='The subtomograms are extracted as a cubic box of this size. '
                           'The wizard will select the box size considering the sampling rate ratio between the '
                           'introduced coordinates and the tomograms that will br used for the extraction.')
        form.addSection(label='Postprocess')
        form.addParam('doInvert', BooleanParam,
                      default=True,
                      label='Invert contrast?',
                      help='Invert the contrast if your tomogram is black '
                           'over a white background.  Xmipp, Spider, Relion '
                           'and Eman require white particles over a black '
                           'background.')
        form.addParallelSection(threads=1, mpi=0, binThreads=3)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        closeSetStepDeps = []
        for tsId in self.tomoTsIdDict.keys():
            pId = self._insertFunctionStep(self.writeSetOfCoordinates3D,
                                           tsId,
                                           prerequisites=[],
                                           needsGPU=False)
            pId = self._insertFunctionStep(self.launchDynamoExtractStep,
                                           tsId,
                                           prerequisites=pId,
                                           needsGPU=False)
            if self.doInvert.get():
                pId = self._insertFunctionStep(self.invertContrastStep,
                                               tsId,
                                               prerequisites=pId,
                                               needsGPU=False)
            closeSetStepDeps.append(pId)
        self._insertFunctionStep(self.createOutputStep,
                                 prerequisites=closeSetStepDeps,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        # Get the intersection between the tomograms and coordinates provided (this covers possible subsets made)
        inTomos = self.getInputTomograms()
        inCoords = self.getInCoords()
        coordsPresentTomoIds = inCoords.getTSIds()
        tomosPresentTsIds = inTomos.getTSIds()
        commonTomoIds = set(coordsPresentTomoIds) & set(tomosPresentTsIds)
        self.tomoTsIdDict = {tomo.getTsId(): tomo.clone() for tomo in inTomos if tomo.getTsId() in commonTomoIds}
        # Calculate the scale factor
        samplingRateCoord = inCoords.getSamplingRate()
        samplingRateTomo = inTomos.getFirstItem().getSamplingRate()
        self.scaleFactor = float(samplingRateCoord / samplingRateTomo)

    def writeSetOfCoordinates3D(self, tsId: str) -> None:
        logger.info(cyanStr(f"tsId = {tsId} - Writing the coordinates of tomogram into Dynamo format..."))
        try:
            makePath(self._getTomoResultsDir(tsId))
            tomoFile = self._getVllFileName(tsId)
            tomo = self.tomoTsIdDict[tsId]
            # Write the VLL file, the coords file and the angles file as expected by Dynamo
            with open(self._getCoordsFileName(tsId), "w") as outC, \
                    open(self._getAnglesFileName(tsId), 'w') as outA, \
                    open(tomoFile, 'w') as tomoFid:
                tomoFid.write(f'{abspath(tomo.getFileName())}\n')
                with self._lock:
                    for coord in self.getInCoords().iterCoordinates(tomo):
                        angles_coord = matrix2eulerAngles(coord.getMatrix())
                        x = self.scaleFactor * coord.getX(BOTTOM_LEFT_CORNER)
                        y = self.scaleFactor * coord.getY(BOTTOM_LEFT_CORNER)
                        z = self.scaleFactor * coord.getZ(BOTTOM_LEFT_CORNER)
                        outC.write("%.2f\t%.2f\t%.2f\t1\n" % (x, y, z))
                        outA.write("%.2f\t%.2f\t%.2f\n" % (angles_coord[0], angles_coord[1], angles_coord[2]))
        except Exception as e:
            self.failedItems.append(tsId)
            logger.error(redStr(f'tsId = {tsId} -> input conversion failed with the exception -> {e}'))

    def launchDynamoExtractStep(self, tsId: str):
        logger.info(cyanStr(f"tsId = {tsId} - Extracting the particles from tomogram..."))
        if tsId not in self.failedItems:
            try:
                codeFilePath = self.writeMatlabCode(tsId)
                args = ' %s > %s' % (codeFilePath, self._getLogFileName(tsId))
                self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())
            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> Dynamo extraction failed with the exception -> {e}'))

    def invertContrastStep(self, tsId: str):
        logger.info(cyanStr(f"tsId = {tsId} - Inverting the contrast of the particles extracted..."))
        if tsId not in self.failedItems:
            try:
                for subTomoFile in self._getSubtomoFileNames(tsId):
                    # Read the subtomo file
                    with mrcfile.mmap(subTomoFile, mode='r', permissive=True) as mrc:
                        invertedData = -1 * mrc.data
                    # Write the result
                    with mrcfile.new_mmap(subTomoFile, overwrite=True, shape=invertedData.shape,
                                          mrc_mode=2) as mrc:  # Mode 2 is float32 (see new_mmap)
                        for i in range(len(invertedData)):
                            mrc.data[i, :, :] = invertedData[i, :, :]
                        mrc.update_header_from_data()
                        mrc.voxel_size = self.getInputTomograms().getSamplingRate()
            except Exception as e:
                self.failedItems.append(tsId)
                logger.error(redStr(f'tsId = {tsId} -> Invert contrast failed with the exception -> {e}'))

    def createOutputStep(self):
        logger.info(cyanStr("Registering the results..."))
        outSubtomos = self.getOutSetOfSubtomos()
        for tsId in self.tomoTsIdDict.keys():
            if tsId in self.failedItems:
                continue
            try:
                tomo = self.tomoTsIdDict[tsId]
                tomoFileName = tomo.getFileName()
                sRate = tomo.getSamplingRate()
                currentSubtomoFiles = sorted(self._getSubtomoFileNames(tsId))
                excludedIndices = self._getDynamoExcludedPartInds(tsId)
                if excludedIndices:
                    logger.info(cyanStr(f"===> tsId = {tsId} - Excluded indices [{len(excludedIndices)}] by "
                                        f"Dynamo {excludedIndices}"))
                else:
                    logger.info(cyanStr(f"tsId = {tsId} - No indices were excluded by Dynamo..."))
                coordCounter = 0
                for i, inCoord in enumerate(self.getInCoords().iterCoordinates(volume=tomo)):
                    if i in excludedIndices:
                        continue
                    subtomogram = SubTomogram()
                    transform = Transform()
                    subtomoFile = currentSubtomoFiles[coordCounter]
                    subtomogram.setSamplingRate(sRate)
                    subtomogram.setFileName(subtomoFile)
                    subtomogram.setVolName(tomoFileName)
                    subtomogram.setCoordinate3D(inCoord)
                    trMatrix = copy.copy(inCoord.getMatrix())
                    transform.setMatrix(scaleTrMatrixShifts(trMatrix, self.scaleFactor))
                    subtomogram.setTransform(transform, convention=TR_DYNAMO)
                    outSubtomos.append(subtomogram)
                    outSubtomos.update(subtomogram)
                    coordCounter += 1
            except Exception as e:
                logger.error(redStr(f'tsId = {tsId} -> Unable to register the output with '
                                    f'exception {e}. Skipping... '))
                continue
        self._defineOutputs(**{self._possibleOutputs.subtomograms.name: outSubtomos})
        self._defineSourceRelation(self.getInCoords(isPointer=True), outSubtomos)

    # --------------------------- DEFINE utils functions ----------------------
    def getInputTomograms(self):
        """ Return the tomogram associated to the 'SetOfCoordinates3D' or 'Other' tomograms. """
        if self.tomoSource.get() == SAME_AS_PICKING:
            return self.getInCoords().getPrecedents()
        else:
            return self.getInTomos()

    def writeMatlabCode(self, tsId: str) -> str:
        codeFilePath = self._getMatlabFileCode(tsId)
        # Write code to Matlab code file
        with open(codeFilePath, 'w') as codeFid:
            catalogue = getCatalogFile(self._getTomoResultsDir(tsId))
            content = "savePath = '%s'\n" % self._getCroppedParticlesDir(tsId)
            content += "box = %i\n" % self.boxSize.get()
            content += "dcm -create '%s' -fromvll '%s'\n" % (removeExt(catalogue), self._getVllFileName(tsId))
            content += "c = dread('%s')\n" % catalogue
            content += "coordsData = readmatrix('%s')\n" % self._getCoordsFileName(tsId)
            content += "angles = readmatrix('%s')\n" % self._getAnglesFileName(tsId)
            content += "coords = coordsData(:,1:3)\n"
            content += "tags = coordsData(:,4)'\n"
            content += "parfor(tag=unique(tags), %i)\n" % self.binThreads.get()
            content += "tomoCoords = coords(tags == tag, :)\n"
            content += "tomoAngles = angles(tags == tag, :)\n"
            content += "t = dynamo_table_blank(size(tomoCoords, 1), 'r', tomoCoords, 'angles', tomoAngles)\n"
            content += "dtcrop(c.volumes{tag}.fullFileName, t, strcat(savePath, num2str(tag)), box, 'ext', 'mrc')\n"
            content += "end\n"
            codeFid.write(content)

        return codeFilePath

    def _getTomoResultsDir(self, tsId: str) -> str:
        return self._getExtraPath(tsId)

    def _getCoordsFileName(self, tsId: str) -> str:
        return join(self._getTomoResultsDir(tsId), 'coords.txt')

    def _getAnglesFileName(self, tsId: str) -> str:
        return join(self._getTomoResultsDir(tsId), 'angles.txt')

    def _getVllFileName(self, tsId: str) -> str:
        return join(self._getTomoResultsDir(tsId), VLL_FILE)

    def _getMatlabFileCode(self, tsId: str) -> str:
        return join(self._getTomoResultsDir(tsId), 'DynamoExtraction.m')

    def _getCroppedParticlesDir(self, tsId: str) -> str:
        return join(self._getTomoResultsDir(tsId), CROP_DIR)

    def _getLogFileName(self, tsId: str) -> str:
        return join(self._getTomoResultsDir(tsId), LOG_FILE_NAME)

    def getOutSetOfSubtomos(self) -> SetOfSubTomograms:
        outSubtomos = SetOfSubTomograms.create(self._getPath(), template='submograms%s.sqlite')
        inTomos = self.getInputTomograms()
        finalSRate = inTomos.getSamplingRate()
        outSubtomos.setSamplingRate(finalSRate)
        outSubtomos.setCoordinates3D(self.getInCoords())
        inTomosAcq = inTomos.getAcquisition()
        if inTomosAcq:
            outSubtomos.setAcquisition(inTomosAcq)
        return outSubtomos

    def _getSubtomoFileNames(self, tsId: str) -> List[str]:
        return glob.glob(join(f'{self._getCroppedParticlesDir(tsId)}*', '*.mrc'))

    def _getDynamoExcludedPartInds(self, tsId: str) -> List[int]:
        indices = []
        logFile = self._getLogFileName(tsId)
        with open(logFile, 'r') as file:
            for line in file:
                if line.startswith("ATTENTION: cannot crop particle"):
                    parts = line.strip().split()
                    # The particle index is expected to be the last element in the line
                    index = int(parts[-1])
                    indices.append(index)
        return indices

        # --------------------------- DEFINE info functions ----------------------
    def _methods(self):
        methodsMsgs = []
        if self.getOutputsSize() >= 1:
            msg = ("A total of %i subtomograms of size %i were extracted"
                   % (self.getInCoords().getSize(), self.boxSize.get()))
            if self.tomoSource.get() == OTHER:
                msg += (" from another set of tomograms: %s"
                        % self.getObjectTag('inputTomogram'))
            msg += " using coordinates %s" % self.getObjectTag('inputCoordinates')
            msg += self.methodsVar.get('')
            methodsMsgs.append(msg)
        else:
            methodsMsgs.append("Set of Subtomograms not ready yet")
        return methodsMsgs

    def _summary(self):
        summary = []
        if self.isFinished():
            summary.append("Tomogram source: *%s*" % self.getEnumText("tomoSource"))
            if self.coordsRemoved.get():
                summary.append('*Some coordinates were removed* (This is because the box associated to those '
                               'coordinates partially or totally lay out of the corresponding tomogram. A good way to '
                               'avoid this is to decrease the box size).')

            if self.doInvert:
                summary.append('*Contrast was inverted.*')
        return summary
