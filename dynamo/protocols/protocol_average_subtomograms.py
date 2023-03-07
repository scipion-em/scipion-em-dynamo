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
from os.path import join
from dynamo import Plugin
from dynamo.convert import writeDynTable, writeSetOfVolumes
from dynamo.protocols.protocol_base_dynamo import DynamoProtocolBase
from pwem.emlib.image import ImageHandler
from pyworkflow.protocol import PointerParam, BooleanParam, IntParam, GE, FloatParam, EnumParam, LE
from pyworkflow.utils import Message, makePath
from tomo.objects import AverageSubTomogram


# Values for enum param fcCompensateOpt
PART_NUMBER = 0
PART_FRACTION = 1


class DynamoProtAvgSubtomograms(DynamoProtocolBase):
    _label = 'Average subtomograms'
    tableName = 'initial.tbl'
    dataDirName = 'data'
    averageDirName = 'average'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    # --------------- DEFINE param functions ---------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inSubtomos', PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True,
                      label='Subtomograms')
        form.addSection(label='Parameters')
        form.addParam('fcompensate', BooleanParam,
                      default=True,
                      label='Do Fourier compensation?',
                      important=True,
                      help='If set to Yes, it will divide each component of the "raw average" of the particles by '
                           'the number of particles actually contributing to this particular component')
        group = form.addGroup('Fourier compensation',
                              condition='fcompensate')
        group.addParam('fcCompensateOpt', EnumParam,
                       display=EnumParam.DISPLAY_HLIST,
                       label='Min. particles for a Fourier coefficient expressed as',
                       choices=['Number', 'Fraction'],
                       default=PART_NUMBER)
        group.addParam('fmin', IntParam,
                       default=1,
                       validators=[GE(1)],
                       condition='fcCompensateOpt == %i' % PART_NUMBER,
                       label='Min. no. of particles',
                       help='Used to manage the elimination of underrepresented Fourier components. Minimum number of '
                            'particles that must contribute into a Fourier coefficient in order to accept it')
        group.addParam('fminFraction', FloatParam,
                       default=0,
                       validators=[GE(0), LE(1)],
                       condition='fcCompensateOpt == %i' % PART_FRACTION,
                       label='Min. fraction of particles [0, 1]',
                       help='Used to manage the elimination of underrepresented Fourier components. Minimum fraction '
                            'of particles that must contribute into a Fourier coefficient in order to accept it')
        group.addParam('fCompensationSmoothingMask', PointerParam,
                       pointerClass='VolumeMask',
                       allowsNull=True,
                       label='Mask for Fourier compensation (opt.)',
                       help='If a mask is provided, the raw average will be multiplied with this mask before '
                            'the Fourier compensation')
        form.addParam('nmask', PointerParam,
                      pointerClass='VolumeMask',
                      label='Normalization mask (opt.)',
                      allowsNull=True,
                      help='If a mask is provided, each particle will be normalized to have zero mean and standard '
                           'deviation 1 inside this mask before adding it to the average.')
        form.addParam('impRotMasking', BooleanParam,
                      default=True,
                      label='Do implicit rotation masking? (opt.)',
                      help='If set to Yes, in the rotated particles, the material outside a spherical mask will not '
                           'be computed. The particles will de facto appear with a spherical mask.')
        form.addParallelSection(threads=4, mpi=0)

    # --------------- INSERT steps functions ----------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.avgStep)
        self._insertFunctionStep(self.convertOutputStep)
        self._insertFunctionStep(self.createOutputStep)

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
            writeDynTable(fhTable, inSubtomos)

    def avgStep(self):
        codeFileName = self._getExtraPath('inCode.doc')
        cmd = self.genAvgCmd()
        with open(codeFileName, 'w') as cFile:
            cFile.write(cmd)

        Plugin.runDynamo(self, codeFileName)

    def convertOutputStep(self):
        # Replacing directly the .em to .mrc generates headers with non-valid dimensions (1 x 1 x boxSize),
        # So the average is converted explicitly using the Image Handler
        ih = ImageHandler()
        genFile = self.getOutputFile()
        outputFile = self.getOutputFile('mrc')
        ih.convert(genFile, outputFile)

    def createOutputStep(self):
        inSubtomos = self.inSubtomos.get()
        avg = AverageSubTomogram()
        avg.setFileName(self.getOutputFile('mrc'))
        avg.setSamplingRate(inSubtomos.getSamplingRate())
        self._defineOutputs(**{super()._possibleOutputs.average.name: avg})
        self._defineSourceRelation(inSubtomos, avg)

    # --------------------------- INFO functions ------------------------------

    # --------------------------- UTILS functions ----------------------------
    def getOutputFile(self, ext='em'):
        return self._getExtraPath(self.averageDirName, self.averageDirName + '.' + ext)

    def genAvgCmd(self):
        fminFraction = self.fminFraction.get()
        fCompSmoothMask = self.fCompensationSmoothingMask.get()
        nmask = self.nmask.get()
        cmd = "daverage('%s', " % self._getExtraPath(self.dataDirName)
        cmd += "'table', '%s', " % self._getExtraPath(self.tableName)
        cmd += "'o', '%s', " % self.getOutputFile()  # output file
        if self.fcompensate.get():
            cmd += "'fcompensate', 1, "
            if fminFraction > 0:
                # fminFraction overrides fmin and sets the fraction of particles that are needed to accept
                # a Fourier coeffient.
                cmd += "'fminFraction', %.2f, " % fminFraction
            else:
                cmd += "'fmin', %i, " % self.fmin.get()
            if fCompSmoothMask:
                cmd += "'fCompensationSmoothingMask', '%s', " % fCompSmoothMask.getFileName()
        if nmask:
            cmd += "'nmask', '%s', " % nmask.getFileName()
        if self.impRotMasking.get():
            cmd += "'implicitRotationMasking', 1, "
        cmd += "'matlab_workers', %i, " % self.numberOfThreads.get()

        cmd += "'v', 1"  # Verbose
        cmd += ")"
        return cmd
