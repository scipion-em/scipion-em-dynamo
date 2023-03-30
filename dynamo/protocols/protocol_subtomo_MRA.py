# -*- coding: utf-8 -*-
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
import os
from enum import Enum
from os.path import join, abspath
import numpy as np
from pwem.objects.data import SetOfVolumes
from pyworkflow import BETA
from pyworkflow.object import Set, String
from pyworkflow.protocol import GPU_LIST, USE_GPU
from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam, StringParam, LEVEL_ADVANCED, \
    NumericListParam, Form
from pyworkflow.utils import Message
from pyworkflow.utils.path import makePath
from dynamo import Plugin
from dynamo.convert import writeSetOfVolumes, writeDynTable, readDynTable
from tomo.protocols.protocol_base import ProtTomoSubtomogramAveraging
from tomo.objects import AverageSubTomogram, SetOfSubTomograms

IMPORT_CMD_FILE = 'prepareProject.m'
SHOW_PROJECT_CMD_FILE = "showProject.m"
ALIGNMENT_CMD_FILE = "runAlign.m"
DEFAULT_DIM = "0"
DYNAMO_ALIGNMENT_PROJECT = 'dynamoAlignmentProject'
INI_TABLE = "initialTable.tbl"
DATADIR_NAME = "data"
MASKSDIR_NAME = "masks"
TEMPLATESDIR_NAME = 'templates'

# Cone flip modes
NO_INVERSION = 0
INVERTED_COARSEST = 1
INVERTED_ALL_REF_LEVELS = 2
FLIP_MODES = [NO_INVERSION, INVERTED_COARSEST, INVERTED_ALL_REF_LEVELS]

# Particles for averaging thresholding modes
NO_THRESHOLD = 0
ABS_THRESHOLD = 1
EFF_THRESHOLD_1 = 2
EFF_THRESHOLD_2 = 3
TOTAL_NO_PARTICLES_BY_CC = 4
FRACTION_NO_PARTICLES = 5
AVG_THRESHOLD_MODES = [NO_THRESHOLD, ABS_THRESHOLD, EFF_THRESHOLD_1, EFF_THRESHOLD_2, TOTAL_NO_PARTICLES_BY_CC,
                       FRACTION_NO_PARTICLES]

# Area search modes
NO_LIMITS = 0
CENTER_OF_THE_BOX = 1
FROM_PREVIOUS_ESTIMATION = 2
FROM_FIRST_ITER_ROUND = 3
FROM_FIRST_ITER_PRJ = 4
AREA_SEARCH_MODES = [NO_LIMITS, CENTER_OF_THE_BOX, FROM_PREVIOUS_ESTIMATION, FROM_FIRST_ITER_ROUND, FROM_FIRST_ITER_PRJ]


class DynRefineOuts(Enum):
    subtomograms = SetOfSubTomograms
    average = AverageSubTomogram


class DynamoSubTomoMRA(ProtTomoSubtomogramAveraging):
    """This protocol will align subtomograms using Dynamo"""  # MRA Subtomogram Averaging"""

    _label = 'Subtomogram alignment'
    _devStatus = BETA
    _possibleOutputs = DynRefineOuts

    @classmethod
    def getUrl(cls):
        return "https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Alignment_project"

    def __init__(self, **args):
        ProtTomoSubtomogramAveraging.__init__(self, **args)
        self.dimRounds = String()
        self.fhTable = None
        self.masksDir = None
        self.doMra = None

    # --------------------------- DEFINE param functions ------------------------

    def _defineParams(self, form: Form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addHidden(USE_GPU, BooleanParam,
                       default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation. Select the one you want to use.")
        form.addHidden(GPU_LIST, StringParam, default='0',
                       expertLevel=LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used")
        form.addParam('inputVolumes', PointerParam,
                      pointerClass="SetOfSubTomograms",
                      important=True,
                      label='Set of subtomograms',
                      help="Set of subtomograms to align with dynamo")
        form.addParam('templateRef', PointerParam,
                      pointerClass='Volume, SubTomogram',
                      important=True,
                      label="Template",
                      help='The size of the template should be equal or smaller than the size of the particles.')  # If you '
        # 'pass a single file in multi-reference modus (MRA), Dynamo will just made copies of it.')
        form.addParam('numberOfIters', NumericListParam,
                      label='Iterations (R)',
                      important=True,
                      default=5,
                      help="Number of iterations per round (R). Used to identify the number of rounds desired to be "
                           "carried out.")
        form.addParam('dim', NumericListParam,
                      label='Particle dimensions (R)',
                      default=DEFAULT_DIM,
                      help="If only one round, leave 0 to use the size of your particle. If working with multiple "
                           "rounds, the size of the particles for each round is expected to be explicitly specified. "
                           "If not, the size of the input particles will be used for all the rounds. This option can "
                           "be used, for example, to reduce the particles size for a particular round and increase the "
                           "speed. E.g.: 64 128 128.")
        form.addBooleanParam('useDynamoGui', 'Launch dynamo GUI',
                             default=False,
                             expertLevel=LEVEL_ADVANCED,
                             help="Launches Dynamo's alignment project GUI. Do not 'Run' the project,"
                                  " Scipion will do it for you.")

        form.addSection(label='Masks')
        form.addParam('alignMask', PointerParam,
                      pointerClass='Volume',  # SetOfVolumes',
                      label="Alignment mask (opt)",
                      allowsNull=True,
                      help='This is the MOST important mask from the different types that can be used by Dynamo. '
                           'It is used for the local correlation. The template is rotated and shifted (virtually, '
                           'through fourier acceleration). At each posible combination of rotation angles and shift, '
                           'the mask is also rotated and shifted, defining a moving region inside the template. '
                           'The rotated and shifted template is compared to the data particle only inside this moving '
                           'region. See more details in '
                           'https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Alignment_mask.')
        # TODO: only for MRA --> uncomment when it's offered again
        # form.addParam('cmask', PointerParam,
        #               pointerClass='Volume, SetOfVolumes',
        #               label="Classification mask (opt)",
        #               allowsNull=True,
        #               condition='Only for multirreference')
        form.addParam('fmask', PointerParam,
                      pointerClass='Volume',  # SetOfVolumes',
                      label="Fourier mask on template (opt)",
                      allowsNull=True,
                      help='Used in very few special cases. The Fourier Mask that you define on a template during an '
                           'alignment project describes the Fourier content of the template, not the one of the data '
                           'particles. It does not reflect directly the missing wedge of the tomogram. See more '
                           'details in '
                           'https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Fourier_mask_on_template.')
        form.addParam('smask', PointerParam,
                      pointerClass='Volume',  # SetOfVolumes',
                      label="Fourier Shell Correlation mask (opt)",
                      allowsNull=True,
                      help='Used in the context of adaptive bandpass filtering. This procedure needs an automatic '
                           'evaluation of the attained resolution at each iteration step. This is performed through '
                           'and FSC computation of the averages computed independently in the different channels of '
                           'the odd/even computation. See more details in '
                           'https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Smoothing_mask')

        form.addSection(label='Angular scanning')
        form.addParam('cr', NumericListParam,
                      label='Cone aperture (R)',
                      default=360,
                      help="The ‘cone_’ parameters define a set of orientations that will be sampled around a "
                           "previously determined old orientation. Here we speak of orientations of the vertical axis "
                           "of the template, not full unconstrained rotations. Allowing this axis to move inside a "
                           "cone involves only two Euler angles (tdrot and tilt).\n\n"
                           "‘cone_range’, is the extent of this cone in degrees (360 being the full possible range of "
                           "axis orientations).")
        form.addParam('cs', NumericListParam,
                      label='Cone sampling (R)',
                      default=45,
                      help="The ‘cone_’ parameters define a set of orientations that will be sampled around a "
                           "previously determined old orientation. Here we speak of orientations of the vertical axis "
                           "of the template, not full unconstrained rotations. Allowing this axis to move inside a "
                           "cone involves only two Euler angles (tdrot and tilt).\n\n"
                           "’cone_sampling’ is the step inside the cone defined by ´cone_range´, also in degrees. "
                           "The orientations are generated so as to provide an uniform coverage.")
        form.addParam('cf', NumericListParam,
                      label='Cone flip (R)',
                      default=0,
                      expertLevel=LEVEL_ADVANCED,
                      help='Generates a mirrored scanning geometry: the "cone" of directions is complemented with the '
                           'diametrally oposed cone.\n'
                           'This is useful when averaging elongated particles in the case in which the direction of '
                           'each one is not certain, i.e., the initial table catches the overall orientation of each '
                           'particle, but it is not certain on which end is the "head" and which is the "tail", so '
                           'that the refinement should allow "flippling" the particles (but still produce a scanning '
                           'set of angles concentrated along the axis of the particle). Values:\n'
                           '\n\t* 0:  No inversion of the cone (default!)'
                           '\n\t* 1:  The cone is inverted only for the coarsest level of the multigrid refinement'
                           '\n\t* 2:  The cone is inverted for all refinement levels')
        form.addParam('inplane_range', NumericListParam,
                      label='Azimuth rotation range (R)',
                      default=360,
                      help='The ‘inplane_’ parameters complete the set of scanned Euler triplets. After each of the '
                           'axis reorentations defined by the ‘cone_’ parameters, the template will be rotated about '
                           'the new orientation of its axis. This involves only the ‘narot’ angle.\n\n'
                           '‘inplane_range’ defines the angular interval to be scanned around the old value of narot.')
        form.addParam('inplane_sampling',
                      NumericListParam,
                      label='Azimuth rotation sampling (R)',
                      default=45,
                      help='The ‘inplane_’ parameters complete the set of scanned Euler triplets. After each of the '
                           'axis reorentations defined by the ‘cone_’ parameters, the template will be rotated about '
                           'the new orientation of its axis. This involves only the ‘narot’ angle.\n\n'
                           'The project parameter ‘inplane_range’ defines the angular interval to be scanned around '
                           'the old value of narot, and ‘inplane_sampling’ defines the interval.')
        form.addParam('inplane_flip', NumericListParam,
                      label='Azimuth flip (R)',
                      default=0,
                      expertLevel=LEVEL_ADVANCED,
                      help='Flips the set of inplane rotations.\n'
                           'The set of inplane rotations to scan will be the original set plus the flipped '
                           'orientations.\n'
                           'This is useful when the particles have a directionality, but it is not very well defined. '
                           'For instance, if you have decorations on a microtubule, you might expect that most of them '
                           'will have similar orientations along the direction of the tube. However, this direction '
                           'might not be obvious: if you have particles coming from two different tubes, reducing the '
                           'span of azimuthal rotations to a narrow set along the estimated direction of the tube '
                           'might be dangerous, as the particles from the two tubes might be oriented in different '
                           'directions.Values:\n'
                           '\n\t* 0  :  no flip (default).'
                           '\n\t* 1  :  flips the coarsest level  in the multilevel grid.'
                           '\n\t* 2  :  flips the full set (all levels).')
        form.addParam('rf', NumericListParam,
                      default=5,
                      label='Refine iterations per particle (R)',
                      help="How many refinement iterations are carried out on each single particle. This refinement "
                           "when comparing rotations of the reference against the data, takes the best orientation and "
                           "looks again with a finer sampling. The sampling in the refined search will be half of the "
                           "sampling used in the original one.  The range of the refined search encompasses all the "
                           "orientations that neighobur the best orientation found in the original search.")
        form.addParam('rff', NumericListParam,
                      label='Refine factor (R)',
                      default=2,
                      help="Controls the size of the angular neighborhood during the local refinement of the angular "
                           "grid.")

        form.addSection('Shift limits')
        form.addParam('lim', NumericListParam,
                      label='Shift limits (R)',
                      default='4',
                      help='Restricts the search area to an sphere of the given radius centered and oriented in the '
                           'last found position (Ellipsoid will be implemented soon).')
        # TODO: manage this properly or directly offer only the sphere
        # 'Restricts the search area to an ellipsoid centered and oriented in the last found '
        #      'position. The three parameters are the semi-axes of the ellipsoid.')
        form.addParam('limm', NumericListParam,
                      default=NO_LIMITS,
                      label='Shift limiting way (R)',
                      help='States how exactly the shifts (area search) will be interpreted:\n'
                           '\n\t* 0:  no limitations (can easily produce artifacts if the initial reference is bad).'
                           '\n\t* 1:  limits are understood from the center of the particle cube.'
                           '\n\t* 2:  limits are understood from the previous estimation on the particle position '
                           '(i.e., the shifts available) With this option, the origin of the shifts changes at every '
                           'iteration.'
                           '\n\t* 3:  limis are understood from the estimation provided for the first iteration of the '
                           'round. The origin of the shifts will change at each round.'
                           '\n\t* 4:  limis are understood from the estimation provided for the first iteration of the '
                           'project. The origin of the shifts is thus defined for the full project, and stays static '
                           'all during the full computation.\n\n'
                           'Note that options 3 and 4 are useful to avoid particles gradually shifting away from '
                           'the initially user-entered locations.')
        form.addParam('separation', IntParam,
                      label='Separation in tomogram [pix.] (R)',
                      default=0,
                      help='When tuned to  positive number, it will check the relative positions (positions in the '
                           'tomogram+shifts) of all the particles in each tomogram separately. Whenever two particles '
                           'are closer together than "separation_in_tomogram", only the particle with the higher '
                           'correlation will stay.')

        form.addSection(label='Thresholding')
        form.addParam('threshold', NumericListParam,
                      label='Threshold parameter (R)',
                      default=0.2,
                      help='Different thresholding policies can be used in order to select which particles are averaged'
                           ' in view of their CC (cross correlation value) . The value of the thresholding parameter '
                           'defined here  will be interpreted differently depending on the "threshold_modus"')
        form.addParam('thresholdMode', NumericListParam,
                      default=NO_THRESHOLD,
                      label='Threshold modus (R)',
                      help='Specify which particles contribute to the average at the end of each iteration. '
                           'Different thresholding policies can be used to select particles according to their CC '
                           'value. Thus value of the "threshold" parameter you input  (denoted as THRESHOLD below) '
                           'will be interpreted differently depending on the "threshold_modus" defined here.\n\n'
                           'Possible values of the thresholding policy "threshold_modus":\n'
                           '\n\t* 0: no thresholding policy'
                           '\n\t* 1: THRESHOLD is an absolute threshold (only particles with CC above this value are '
                           'selected).'
                           '\n\t* 2: effective threshold = mean(CC) * THRESHOLD.'
                           '\n\t* 3: effective threshold = mean(CC) + std(CC) * THRESHOLD.'
                           '\n\t* 4: THRESHOLD is the total number of particles (ordered by CC ).'
                           '\n\t* 5: THRESHOLD ranges between 0 and 1  and sets the fraction of particles.')
        # '\n\t* 11,21,31,34,41,51: select the same particles as 1,2,3,4 or 5, BUT non selected '
        # 'particles will be excluded:'
        # '\n\t\t- from averaging in the present iteration, and'
        # '\n\t\t- ALSO from alignment in the next iteration (unlike 1,2,3,4,5).')
        form.addParam('threshold2', NumericListParam,
                      label='Second threshold parameter (R)',
                      default=0.2,
                      expertLevel=LEVEL_ADVANCED,
                      help="Thresholding II is operated against the average produced by the particles that survived"
                           "the first thresholding.")
        form.addParam('thresholdMode2', NumericListParam,
                      default=NO_THRESHOLD,
                      label='Second threshold modus (R)',
                      expertLevel=LEVEL_ADVANCED,
                      help="Thresholding II is operated against the average produced by the particles that survived "
                           "the first thresholding. It uses the same syntax as Threshold I")

        form.addSection(label='Cross correlation')
        form.addParam('high', NumericListParam,
                      label='High pass (R)',
                      default=2,
                      help='Cut off frequency for high pass filtering. The units are pixels in the Fourier space of '
                           'the particle (if your template is smaller than the data, the bandpass parameters will be '
                           'rescaled). The bandpass works with a soft mask, allowing a smoothing window of two pixels '
                           'around the cut frequency.')
        form.addParam('low', NumericListParam,
                      label='Low frequency (R)',
                      default=32,
                      help='Cut off frequency for low pass filtering. The units are pixels in the Fourier space of the '
                           'particle (if your template is smaller than the data, the bandpass parameters will be '
                           'rescaled). The bandpass works with a soft mask, allowing a smothing window of two pixels '
                           'around the cut frequency.')
        form.addParam('sym', StringParam,
                      default='c1',
                      label='Symmetry group (R)',
                      help="Symmetrization is applied at the beginning of the round to the input reference. It is also "
                           "used during the computation of the [eo_fsc] or the iteration.\n"
                           "First chars in string indicate thetype of symmetry operator:\n"
                           "\n\t* 'c'  rotational symmetry around z."
                           "\n\t* 'h'  helical symmetry around z  (parameters: dpsi; dz)."
                           "\n\t* 'ico'  icosahedral symmetry."
                           "\n\t* 'cbo' for cuboctahedral symmetry.\n"
                           "The rest of the string contains the respective parameters. Examples: 'c1', 'h60,5' ,"
                           "'h[60,5]','h[-60,5]'")
        form.addParam('ccmatrixBatch', IntParam,
                      label='Cross-correlation matrix batch',
                      default=128,
                      expertLevel=LEVEL_ADVANCED,
                      help="Number of particles to be kept in memory simultaneously during the computation of the "
                           "ccmatrix. The larger this number, the more efficient the algorithm performance, as more "
                           "computations can be kept for reuse.However, trying to keep all the particles in memory can "
                           "lead to saturate it,blocking the CPU. Additionally, a small batch allows to divide the "
                           "matrix in more blocks. This might be useful in parallel computations.")
        form.addParallelSection(threads=4, mpi=0)

        # form.addParam('pca', BooleanParam,
        #               label='Perform PCA',
        #               default=False,
        #               help="If selected, principal component analysis (PCA) is carried out.")
        # form.addSection(label='Initial density')
        # form.addParam('compensateMissingWedge', BooleanParam,
        #               label='Compensate for missing wedge',
        #               default=False,
        #               condition="refMode == %i" % MODE_GEN_REF,
        #               help="If 'Yes' is selected, this operation will be performed using the method "
        #                    "dynamo_table_randomize_azimuth.")
        # form.addParam('ccmatrix', BooleanParam,
        #               label='Compute  cross-correlation matrix',
        #               default=False,
        #               help="Computation of a Cross-Correlation matrix among the aligned particles.")
        # form.addParam('ccmatrixType', StringParam,
        #               label='Cross-correlation matrix type',
        #               default='align',
        #               condition="ccmatrix",
        #               expertLevel=LEVEL_ADVANCED,
        #               help="string with three characters, each position controling a different aspect: thresholding, "
        #                    "symmetrization, compensation")

    # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self.initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.alignStep)
        self._insertFunctionStep(self.createOutputStep)
        # if self.doMra:
        #     self._insertFunctionStep(self.closeSetsStep)

    # --------------------------- STEPS functions -------------------------------
    def initialize(self):
        # self.doMra = isinstance(self.templateRef.get(), SetOfVolumes)
        nRounds = len(self.numberOfIters.getListFromValues())
        self.dimRounds.set(self.getDimRounds(nRounds))

    def convertInputStep(self):
        dataDir = self._getExtraPath(DATADIR_NAME)
        self.masksDir = self._getExtraPath(MASKSDIR_NAME)
        makePath(*[dataDir, self.masksDir])
        inputVols = self.inputVolumes.get()
        # Convert the input particles into .em if necessary
        areInEmFormat = inputVols.getFirstItem().getFileName().endswith('.em')
        dataDirName = join(dataDir, "particle_")
        writeSetOfVolumes(inputVols, dataDirName, 'id')

        # Write the tbl file with the data read from the introduced particles
        fnTable = self._getExtraPath(INI_TABLE)
        with open(fnTable, 'w') as fhTable:
            writeDynTable(fhTable, inputVols)

        # NOTE for rounds: There are up to 8 rounds. round's params can be specified like:
        # dvput('dynamoAlignmentProject', 'cr_r2', '360');  --> note "_r2" for round 2
        with open(self._getExtraPath(IMPORT_CMD_FILE), 'w') as fhCommands:
            content = "dcp.new('%s', 'data', '%s', 'gui', 0)\n" % (DYNAMO_ALIGNMENT_PROJECT, DATADIR_NAME)
            if not areInEmFormat:
                content += "dynamo_data_format('%s/particle_*.mrc', 'data', 'modus', 'convert', 'extension', '.em')\n" \
                           % DATADIR_NAME
            template = self.templateRef.get()
            if self.doMra:
                templatesDir = self._getExtraPath(TEMPLATESDIR_NAME)
                makePath(templatesDir)
                # Reference management
                refFileNames = [abspath(ref.getFileName()) for ref in self.templateRef.get()]
                content += "dynamo_write_multireference(%s, 'template', '%s')\n" % \
                           (str(refFileNames).replace('[', '{').replace(']', '}'), TEMPLATESDIR_NAME)
                content += self.get_dvput('table', INI_TABLE)
            else:
                content += self.get_dvput('table', INI_TABLE)
                if template:
                    # The template will be a volume (validation method ensures it)
                    referenceName = 'template.em'
                    templateOrigFName = template.getFileName()
                    if not templateOrigFName.endswith('.em'):
                        content += "object = dynamo_read('%s')\n" % abspath(templateOrigFName)
                        content += "dynamo_write(object, '%s')\n" % abspath(self._getExtraPath(referenceName))
                    # convertOrLinkVolume(template, self._getExtraPath(referenceName))
                    content += self.get_dvput('template', referenceName)

            # Masks management
            masks = [self.alignMask.get(), self.fmask.get(), self.smask.get()]  # self.cmask.get()
            for mask in masks:
                content += self.prepareMask(mask)

            # Write the file that will be passed to Dynamo
            fhCommands.write(content)

        Plugin.runDynamo(self, IMPORT_CMD_FILE, cwd=self._getExtraPath())

    def showDynamoGUI(self):
        fhCommands2 = open(self._getExtraPath(SHOW_PROJECT_CMD_FILE), 'w')
        content2 = "dcp '%s';" % DYNAMO_ALIGNMENT_PROJECT
        fhCommands2.write(content2)
        fhCommands2.close()
        Plugin.runDynamo(self, SHOW_PROJECT_CMD_FILE, cwd=self._getExtraPath())

    def alignStep(self):
        with open(self._getExtraPath(ALIGNMENT_CMD_FILE), 'w') as fhCommands2:
            alignmentCommands = self.get_computing_command()
            if not self.useDynamoGui:
                # alignmentCommands += "dynamo_vpr_run('%s','check',true,'unfold',true)" % DYNAMO_ALIGNMENT_PROJECT
                alignmentCommands += "dvcheck('%s')\n" % DYNAMO_ALIGNMENT_PROJECT
                alignmentCommands += "dvunfold('%s')\n" % DYNAMO_ALIGNMENT_PROJECT
            fhCommands2.write(alignmentCommands)

        Plugin.runDynamo(self, ALIGNMENT_CMD_FILE, cwd=self._getExtraPath())
        if self.useDynamoGui:
            self.showDynamoGUI()

        # This way shows output more or less on the fly.
        self.runJob("./%s.exe" % DYNAMO_ALIGNMENT_PROJECT, [], env=Plugin.getEnviron(gpuId=self.getGpuList()[0]),
                    cwd=self._getExtraPath())

        resultsDir = self.getLastIterResultsDir()
        if not os.path.exists(resultsDir):
            raise RuntimeError("Results folder (%s) no generated. "
                               "Probably there has been an error while running the alignment in Dynamo. "
                               "Please, see run.stdout log for more details." % resultsDir)

    def getTotalIterations(self):
        iters = self.numberOfIters.getListFromValues()
        return sum(iters)

    def createOutputStep(self):

        niters = self.getTotalIterations()

        if self.doMra:
            fhSurvivRefs = open(self._getExtraPath('%s/results/ite_%04d/currently_surviving_references_ite_%04d.txt')
                                % (DYNAMO_ALIGNMENT_PROJECT, niters, niters), 'r')
            nline = next(fhSurvivRefs).rstrip()
            for nref in range(len(nline.split())):
                ref = int(nline.split()[nref])
                subtomoSet = self._createSetOfSubTomograms()
                inputSet = self.inputVolumes.get()
                subtomoSet.copyInfo(inputSet)
                self.fhTable = open(
                    self._getExtraPath('%s/results/ite_%04d/averages/refined_table_ref_%03d_ite_%04d.tbl') %
                    (DYNAMO_ALIGNMENT_PROJECT, niters, ref, niters), 'r')
                subtomoSet.copyItems(inputSet, updateItemCallback=self._updateItem)
                self.fhTable.close()
                averageSubTomogram = AverageSubTomogram()
                averageSubTomogram.setFileName(
                    self._getExtraPath('%s/results/ite_%04d/averages/average_ref_%03d_ite_%04d.em')
                    % (DYNAMO_ALIGNMENT_PROJECT, niters, ref, niters))
                averageSubTomogram.setSamplingRate(inputSet.getSamplingRate())

                name = 'outputSubtomogramsRef%s' % str(ref)
                args = {name: subtomoSet}
                subtomoSet.setStreamState(Set.STREAM_OPEN)
                self._defineOutputs(**args)
                self._defineSourceRelation(inputSet, subtomoSet)

                name2 = 'AverageRef%s' % str(ref)
                args2 = {name2: averageSubTomogram}
                self._defineOutputs(**args2)
                self._defineSourceRelation(inputSet, averageSubTomogram)
        else:
            averageSubTomogram = AverageSubTomogram()
            outSubtomos = SetOfSubTomograms.create(self._getPath(), template='subtomograms%s.sqlite')
            inputSet = self.inputVolumes.get()
            outSubtomos.copyInfo(inputSet)
            resTbl = join(self.getLastIterAvgsDir(), 'refined_table_ref_001_ite_%04d.tbl' % niters)
            avgFile = join(self.getLastIterAvgsDir(), 'average_symmetrized_ref_001_ite_%04d.em' % niters)
            # Open the final particles table, that will be used in the updateItemCallback
            self.fhTable = open(resTbl, 'r')
            outSubtomos.copyItems(inputSet, updateItemCallback=self._updateItem)
            self.fhTable.close()
            # Fill the resulting average object
            averageSubTomogram.setFileName(self.convertToMrc(avgFile))
            averageSubTomogram.setSamplingRate(inputSet.getSamplingRate())
            # Define outputs and relations
            outsDict = {self._possibleOutputs.subtomograms.name: outSubtomos,
                        self._possibleOutputs.average.name: averageSubTomogram}
            self._defineOutputs(**outsDict)
            self._defineSourceRelation(self.inputVolumes, outSubtomos)
            self._defineSourceRelation(self.inputVolumes, averageSubTomogram)

    def closeSetsStep(self):
        for outputset in self._iterOutputsNew():
            if isinstance(outputset, SetOfSubTomograms):
                outputset[1].setStreamState(Set.STREAM_CLOSED)
        self._store()

    # --------------------------- UTILS functions --------------------------------
    def getDimRounds(self, nRounds):
        """The number of rounds will be determined the same as Dynamo, which is through the number of elements
        introduced in label 'Iterations'. The parameter 'Particle dimensions' present a specific behavior: if one round,
        a zero value will indicate that the size of the particles is the same as the size of the input subtomogrmas.
        But, if there are more than one rounds, the specific size considered for the particles at that round must be
        specified. This method manages this scenario to avoid the user be forced to introduce the dimensions everytime.
        """
        if nRounds == 1:
            if self.dim.get() != DEFAULT_DIM:
                return self.dim.get()
            else:
                dim, _, _ = self.inputVolumes.get().getDimensions()
                return dim
        else:
            dimVals = self.dim.getListFromValues()
            nDims = len(dimVals)
            if nRounds == nDims:
                return self.dim.get()
            elif nRounds > nDims:
                inParticleSize, _, _ = self.inputVolumes.get().getDimensions()
                dimPattern = str(inParticleSize) + ' '
                if nDims == 1 and dimVals[0] == 0:
                    return dimPattern * nRounds
                else:
                    return self.dim.get() + ' ' + dimPattern * (nRounds - nDims)

    def getRoundParams(self, dynamoParamName, param: String, projectName=DYNAMO_ALIGNMENT_PROJECT, caster=int):
        """ Returns the dynamo command for any of the params that can be specified in the rounds.
        See --> https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Starters_guide#Alignment_projects

        :param dynamoParamName: Dynamo's parameter name
        :param param: Scipion's param containing the values for the rounds
        :param projectName: Optional, defaults to DYNAMO_ALIGNMENT_PROJECT. Name of the dynamo alignment project

        :return the dvput commands as a string
        """
        # Get the values as list
        values = param.getListFromValues(caster=caster)
        command = ""
        for index, value in enumerate(values):
            finalParamName = dynamoParamName
            if index != 0:
                finalParamName += '_r' + str(index + 1)
            command += self.get_dvput(finalParamName, value, projectName=projectName)
        return command

    @staticmethod
    def get_dvput(paramName, value, projectName=DYNAMO_ALIGNMENT_PROJECT):
        """From Dynamo: Changes parameters of a project residing in disk.
        INPUT:
            project  : an existing project
            option   :  a string defining what to do with the modified project:
                       'disk','d'    :  save modified project
                       'unfold','u'  :  unfold modified project
                       'cd'          :  save only if passes check
                       'cu'          :  unfold only if passes check

            parameter/value couples:
                'inround' : integer that identifies a round. Restricts all the subsequent round parameter modifications
                            to one round.
                You can input any series of  parameter/value couples where the parameter is a legal Dynamo project
                parameters (type dynamo_vpr_help for a complete list).
        """
        return "dvput('%s', 'disk', '%s', '%s')\n" % (projectName, paramName, value)

    def get_computing_command(self):
        """ Returns the dynamo commands related to the angular search, threashold, GPu, ..."""
        command = self.getRoundParams("dim", self.dimRounds)
        command += self.get_dvput("apix", self.inputVolumes.get().getSamplingRate())
        command += self.getRoundParams('sym', self.sym, caster=str)
        command += self.getRoundParams("ite", self.numberOfIters)
        # command += self.get_dvput('mra', int(self.doMra))
        # command += self.get_dvput('pcas', int(self.pca.get()))
        # --- Angular scanning ---
        command += self.getRoundParams("cr", self.cr)
        command += self.getRoundParams("cs", self.cs)
        command += self.getRoundParams("cf", self.cf)
        # command += self.get_dvput('ccp', self.ccp)
        command += self.getRoundParams('rf', self.rf)
        command += self.getRoundParams('rff', self.rff)
        command += self.getRoundParams("ir", self.inplane_range)
        command += self.getRoundParams("is", self.inplane_sampling)
        command += self.getRoundParams("if", self.inplane_flip)
        # command += self.get_dvput('icp', self.inplane_check_peak)
        # --- Thresholding ---
        command += self.get_dvput('stm', self.separation)
        if not self.anyValActiveInNumListParam(self.thresholdMode) \
                and not self.anyValActiveInNumListParam(self.thresholdMode2) \
                and not self.anyValActiveInNumListParam(self.limm):
            # Don't compute the CC matrix
            command += self.get_dvput('ccms', 0)
        else:
            # CC matrix stuff
            command += self.get_dvput('ccms', 1)
            command += self.get_dvput('ccmt', 'align')
            command += self.get_dvput('batch', self.ccmatrixBatch)
            # Thresholding stuff
            command += self.getRoundParams('thrm', self.thresholdMode)
            command += self.getRoundParams('thr', self.threshold, caster=float)
            command += self.getRoundParams('thr2m', self.thresholdMode2)
            command += self.getRoundParams('thr2', self.threshold2, caster=float)
        # --- Area search ---
        command += self.getRoundParams('limm', self.limm)
        command += self.getRoundParams('lim', self.lim)
        # --- Filtering ---
        command += self.getRoundParams('low', self.low)
        command += self.getRoundParams('high', self.high)

        # --- Processing software + hardware resources ---
        command += self.get_dvput('mwa', self.numberOfThreads.get())  # Cores used to calculate the average in each iter
        if self.useGpu.get():
            # Param 'cores' is used to specify the number of CPUs involved in the alignment. If GPU is used, Dynamo
            # only works well setting it to 1.
            command += self.get_dvput('cores', 1)  # Not working with more than 1 CPU when using GPU
            command += self.get_dvput('destination', 'standalone_gpu')
            command += self.get_dvput('gpu_motor', 'spp')
            # command += self.get_dvput('gpu_identifier_set', self.getGpuList()[0])
        else:
            command += self.get_dvput('cores', self.numberOfThreads.get())
            command += self.get_dvput('destination', 'standalone')

        return command

    def _updateItem(self, item, row):
        readDynTable(self, item)

    def prepareMask(self, maskObj):
        if maskObj:
            if isinstance(maskObj, SetOfVolumes):
                writeSetOfVolumes(maskObj, join(self._getExtraPath(), 'fmasks/fmask_initial_ref_'), 'ix')
                return self.get_dvput('fmask', self.masksDir)
            else:
                return self.get_dvput('fmask', maskObj.getFileName())
        else:
            return ''

    def getLastIterResultsDir(self):
        return self._getExtraPath('%s', 'results', 'ite_%04d') % \
               (DYNAMO_ALIGNMENT_PROJECT, self.getTotalIterations())

    def getLastIterAvgsDir(self):
        return join(self.getLastIterResultsDir(), 'averages')

    def convertToMrc(self, inFileName):
        import xmipp3
        program = 'xmipp_image_convert'
        outFName = inFileName.replace('.em', '.mrc')
        args = '-i %s ' % inFileName
        args += '-o %s ' % outFName
        args += '-t vol'
        self.runJob(program, args, env=xmipp3.Plugin.getEnviron())
        return outFName

    def sizesOk(self, inVolume, checkLE=True):
        """Method to check the size conditions from Dynamo:
        - The size of the template should be lower or equal than the size of the particles.
        - The size of the masks must be equal to the size of the template."""
        refDims = self.templateRef.get().getDimensions()
        testDims = inVolume.getDimensions()
        if checkLE:
            return testDims <= refDims
        else:
            return testDims == refDims

    @staticmethod
    def anyValActiveInNumListParam(iParam) -> bool:
        """Checks if any of the values of a NumericListParam is greater than 0, which means to
        be active at least in one of the rounds"""
        return np.any(np.array(iParam.getListFromValues()) > 0)

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        # doMra = isinstance(self.templateRef.get(), SetOfVolumes)
        # ref = self.templateRef.get()
        subtomo = self.inputVolumes.get().getFirstItem()
        masks = [self.alignMask.get(), self.fmask.get(), self.smask.get()]  # self.cmask.get()
        introducedMasks = any(masks)
        # if doMra:
        #
        #     def getVolumesSetSize(iVol):
        #         return 1 if isinstance(iVol, Volume) else iVol.getSize()
        #
        #     refSize = getVolumesSetSize(ref)
        #     masksSetSizes = [getVolumesSetSize(mask) for mask in masks if mask]
        #     if introducedMasks:
        #         if any(masksSetSizes) != refSize:
        #             validateMsgs.append('All the optional introduced masks must be sets of volumes of the same '
        #                                 'size as the set of references')

        # Check the reference
        if not self.sizesOk(subtomo):
            validateMsgs.append('The size of the template should be equal or smaller than the size of the particles.')
        # Check the masks
        if introducedMasks:
            for mask in masks:
                if not self.sizesOk(mask):
                    validateMsgs.append('The introduced masks must be of the same size as the template.')
                    break
        # Check the dims values
        dimValues = self.dim.getListFromValues()
        nDims = len(dimValues)
        nIters = len(self.numberOfIters.getListFromValues())
        if nIters > 1 and nDims > 1 and np.any(np.array(dimValues) == 0):
            validateMsgs.append('If working with multiple rounds, the *particle dimensions* for each round must be '
                                'greater than 0.')
        # Check the flip modes
        coneFlipModes = self.cf.getListFromValues()
        aziFlipModes = self.inplane_flip.getListFromValues()
        for roundVal in coneFlipModes + aziFlipModes:
            if roundVal not in FLIP_MODES:
                validateMsgs.append('Non-valid value detected for one of the *cone flip*. Please check the help to '
                                    'see the admitted values.')

        # Check the thresholds modes values
        th1Modes = self.thresholdMode.getListFromValues()
        th2Modes = self.thresholdMode2.getListFromValues()
        areaSearchValues = self.limm.getListFromValues()
        for roundVal in th1Modes + th2Modes:
            if roundVal not in AVG_THRESHOLD_MODES:
                validateMsgs.append('Non-valid value detected for one of the *thresholds*. Please check the help to '
                                    'see the admitted values.')
                break
        for roundVal in areaSearchValues:
            if roundVal not in AREA_SEARCH_MODES:
                validateMsgs.append('Non-valid value detected for the *area search mode*. Please check the help to see '
                                    'the admitted values.')
                break
        return validateMsgs
