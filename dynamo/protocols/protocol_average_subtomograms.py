# ding=utf-8
# **************************************************************************
# *
# * Authors:     Scipion Team  (scipion@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnología (CSIC), Madrid, Spain
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
from os.path import join
from dynamo import Plugin
from dynamo.convert import writeDynTable, writeSetOfVolumes
from dynamo.protocols.protocol_base_dynamo import DynamoProtocolBase
from pwem.convert.headers import setMRCSamplingRate
from pwem.emlib.image.image_readers import EmImageReader
from pyworkflow.protocol import PointerParam, BooleanParam
from pyworkflow.utils import Message, makePath
from tomo.objects import AverageSubTomogram


class DynAvgOuts(Enum):
    average = AverageSubTomogram


class DynamoProtAvgSubtomograms(DynamoProtocolBase):
    """
    Generates an averaged 3D subtomogram map from a collection of aligned
    subtomograms. The protocol combines multiple tomographic particles into a
    consensus reconstruction that enhances common structural features while
    reducing noise and missing information.

    AI Generated:

    Average Subtomograms (DynamoProtAvgSubtomograms) - User Manual
        Overview

        The Average Subtomograms protocol computes a consensus 3D average from
        a collection of subtomograms using Dynamo subtomogram averaging tools.
        Its main objective is to improve the signal-to-noise ratio of cryo-electron
        tomography data by reinforcing structural information shared across many
        particles. This type of averaging is one of the central steps in
        subtomogram analysis workflows because individual subtomograms are often
        extremely noisy and affected by incomplete angular sampling.

        In biological practice, subtomogram averaging is commonly used to study
        macromolecular complexes directly inside cells or native environments.
        Typical applications include membrane protein analysis, ribosome
        organization, viral assembly studies, and in situ structural biology.
        By combining many particles representing the same biological structure,
        the protocol produces a clearer and more interpretable density map that
        can later be refined, classified, or compared with known atomic models.

        Inputs and Data Preparation

        The protocol requires a set of subtomograms that already contain
        positional and orientational information. These particles are interpreted
        as multiple observations of the same biological object and are merged into
        a common structural reference frame before averaging.

        From a practical perspective, the quality of the final average strongly
        depends on the consistency of the input particles. Accurate particle
        picking, reliable alignment, and homogeneous biological content are
        essential for obtaining meaningful results. Mixing particles from
        different conformational states or different assemblies may reduce the
        interpretability of the final reconstruction by blurring structural
        details.

        Orientation Handling and Missing Wedge Compensation

        One of the most important biological considerations in cryo-electron
        tomography is the missing wedge effect, which arises from the limited
        angular range accessible during tomographic acquisition. This effect
        introduces anisotropic resolution and directional artifacts in the final
        averages.

        The protocol provides an option to randomize subtomogram orientations
        before averaging. This strategy is particularly useful when particles are
        distributed with preferred orientations because randomization helps
        compensate for directional bias and produces a more isotropic average.
        In many experimental datasets, especially membrane-associated systems,
        this approach can significantly improve the visual appearance and overall
        interpretability of the reconstruction.

        When orientation randomization is disabled, the original particle
        orientations are preserved. This is generally preferable when the
        orientations already reflect biologically meaningful alignment results or
        when directional organization within the sample must remain intact.

        Implicit Rotation Masking

        The protocol also supports implicit rotation masking during averaging.
        This option restricts the effective contribution of each subtomogram to a
        spherical region, reducing the influence of undefined regions generated
        during rotational transformations.

        Biologically, implicit masking is often beneficial because it suppresses
        edge artifacts and minimizes the propagation of empty regions into the
        final average. This is particularly important for particles extracted from
        crowded cellular environments or for datasets where particles occupy only
        a limited fraction of the extraction box.

        In many routine subtomogram averaging workflows, enabling implicit
        masking improves average stability and leads to cleaner density maps.
        However, users should remain aware that aggressive masking may also
        suppress peripheral structural features if the extraction region is too
        small.

        Averaging Workflow

        During execution, the protocol prepares the subtomograms in Dynamo
        format, organizes the associated metadata, and launches the averaging
        process using the Dynamo computational framework. The averaging combines
        all particles into a single consensus volume representing the common
        structural information across the dataset.

        The procedure is designed to integrate naturally into iterative
        subtomogram refinement workflows. Users often alternate averaging with
        alignment refinement, classification, or particle cleaning in order to
        progressively improve structural resolution and biological homogeneity.

        Outputs and Interpretation

        The protocol produces an averaged subtomogram volume that can be directly
        visualized in standard cryo-EM software packages. The resulting map
        preserves the sampling information of the original subtomograms and is
        suitable for downstream processing steps such as refinement, focused
        classification, segmentation, or atomic fitting.

        From a biological perspective, the averaged map should be interpreted as
        a consensus representation of the input population. Well-defined regions
        generally correspond to structurally conserved features, whereas blurred
        or weak regions may indicate conformational variability, compositional
        heterogeneity, flexibility, or insufficient alignment accuracy.

        Practical Recommendations

        In routine subtomogram averaging workflows, it is generally advisable to
        begin with carefully curated particles and verify alignment quality before
        averaging large datasets. Orientation randomization is often useful for
        strongly preferred orientations, while preserving orientations is more
        appropriate for biologically ordered assemblies.

        Implicit rotation masking is commonly recommended because it improves the
        robustness of the averaging process and reduces rotational artifacts.
        Nevertheless, users should visually inspect the final averages to ensure
        that masking does not remove biologically relevant peripheral densities.

        For highly heterogeneous datasets, performing classification before
        averaging usually leads to more interpretable reconstructions. Averaging
        structurally mixed populations may artificially smooth important
        conformational differences and reduce the achievable resolution.

        Final Perspective

        Subtomogram averaging is one of the foundational methodologies in
        cryo-electron tomography because it transforms noisy individual particles
        into biologically meaningful 3D reconstructions. Careful particle
        selection, thoughtful handling of orientation distributions, and proper
        masking strategies are critical for obtaining reliable structural
        information from complex cellular environments.
    """

    _label = 'Average subtomograms'
    _possibleOutputs = DynAvgOuts
    tableName = 'initial.tbl'
    dataDirName = 'data'
    averageDirName = 'average'

    # --------------- DEFINE param functions ---------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inSubtomos', PointerParam,
                      pointerClass='SetOfSubTomograms',
                      strict=True,
                      important=True,
                      label='Subtomograms')
        form.addParam('impRotMasking', BooleanParam,
                      default=True,
                      label='Do implicit rotation masking? (opt.)',
                      help='If set to Yes, in the rotated particles, the material outside a spherical mask will not '
                           'be computed. The particles will de facto appear with a spherical mask.')
        form.addParam('randomizeOrientation', BooleanParam,
                      default=False,
                      label='Randomize the subtomos orientation?',
                      help='If set to Yes, the orientation of the picked subtomograms will be randimized. This ensures'
                           'to fill the missing wedge obtaining a ball in the average. If set to No, then the orientation'
                           'of the subtomos will be preserve in the average.')
        self.insertBinThreads(form)

    # --------------- INSERT steps functions ----------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep, needsGPU=False)
        self._insertFunctionStep(self.avgStep, needsGPU=False)
        self._insertFunctionStep(self.convertOutputStep, needsGPU=False)
        self._insertFunctionStep(self.createOutputStep, needsGPU=False)

    # --------------- STEPS functions -----------------------
    def convertInputStep(self):
        inSubtomos = self.inSubtomos.get()
        tableName = self._getExtraPath(self.tableName)
        dataDir = self._getExtraPath(self.dataDirName)
        avgDir = self._getExtraPath(self.dataDirName)
        makePath(*[dataDir, avgDir])
        # Generate the data folder formatted along the Dynamo convention
        writeSetOfVolumes(inSubtomos, join(dataDir, 'particle_'), 'id')
        # Generate the Dynamo data table
        with open(tableName, 'w') as fhTable:
            writeDynTable(fhTable, inSubtomos, randomizeOrientation=self.randomizeOrientation)

    def avgStep(self):
        codeFileName = self._getExtraPath('inCode.doc')
        cmd = self.genAvgCmd()
        with open(codeFileName, 'w') as cFile:
            cFile.write(cmd)

        Plugin.runDynamo(self, codeFileName)

    def convertOutputStep(self):
        # Replacing directly the .em to .mrc generates headers with non-valid dimensions (1 x 1 x boxSize),
        # So the average is converted explicitly using the Image Handler
        emFileHandler = EmImageReader()
        genFile = self.getOutputFile()
        outputFile = self.getOutputFile('mrc')
        emFileHandler.emToMrc(genFile, outputFile)

    def createOutputStep(self):
        inSubtomos = self.inSubtomos.get()
        avg = AverageSubTomogram()
        outFn = self.getOutputFile('mrc')
        outSr = inSubtomos.getSamplingRate()
        setMRCSamplingRate(outFn, outSr)  # Update the apix value in file header
        avg.setFileName(outFn)
        avg.setSamplingRate(outSr)
        self._defineOutputs(**{DynAvgOuts.average.name: avg})
        self._defineSourceRelation(inSubtomos, avg)

    # --------------------------- INFO functions ------------------------------

    # --------------------------- UTILS functions ----------------------------
    def getOutputFile(self, ext='em'):
        return self._getExtraPath(self.averageDirName, self.averageDirName + '.' + ext)

    def genAvgCmd(self):
        cmd = "daverage('%s', " % self._getExtraPath(self.dataDirName)
        cmd += "'table', '%s', " % self._getExtraPath(self.tableName)
        cmd += "'o', '%s', " % self.getOutputFile()  # output file
        if self.impRotMasking.get():
            cmd += "'implicitRotationMasking', 1, "
        cmd += "'extension', 'mrc', "  # informs Dynamo that the data folder uses an ext different to the default .em
        cmd += "'matlab_workers', %i, " % self.binThreads.get()
        cmd += "'v', 1"  # Verbose
        cmd += ")"
        return cmd
