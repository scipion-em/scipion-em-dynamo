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
from pwem.emlib.image import ImageHandler
from pyworkflow.protocol import PointerParam, BooleanParam
from pyworkflow.utils import Message, makePath
from tomo.objects import AverageSubTomogram


class DynAvgOuts(Enum):
    average = AverageSubTomogram


class DynamoProtAvgSubtomograms(DynamoProtocolBase):

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
        ih = ImageHandler()
        genFile = self.getOutputFile()
        outputFile = self.getOutputFile('mrc')
        ih.convert(genFile, outputFile)

    def createOutputStep(self):
        inSubtomos = self.inSubtomos.get()
        avg = AverageSubTomogram()
        avg.setFileName(self.getOutputFile('mrc'))
        avg.setSamplingRate(inSubtomos.getSamplingRate())
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
