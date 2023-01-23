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
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam
from pyworkflow.utils import Message, makePath, createLink
from tomo.objects import AverageSubTomogram
from tomo.protocols import ProtTomoBase


class dynAvgOutputs(Enum):
    average = AverageSubTomogram


class DynamoProtAvgSubtomograms(EMProtocol, ProtTomoBase):

    _label = 'Average subtomograms'
    _devStatus = BETA
    _possibleOutputs = dynAvgOutputs
    tableName = 'initial.tbl'
    dataDirName = 'data'
    averageDirName = 'average'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------- DEFINE param functions ---------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inSubtomos', PointerParam,
                      pointerClass='SetOfSubTomograms',
                      important=True,
                      label='Subtomograms')

    # --------------- INSERT steps functions ----------------
    def _insertAllSteps(self):
        # self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.avgStep)
        # self._insertFunctionStep(self.convertOutputStep)
        # self._insertFunctionStep(self.createOutputStep)

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
        # matlabFile = join(Plugin.getHome(), 'matlab', 'src', 'daverage.m')
        codeFileName = 'inCode.doc'
        cmd = "daverage('%s', " % self.dataDirName
        cmd += "'table', '%s', " % self.tableName
        cmd += "'o', '%s', " % self.getOutputFile()  # output file
        cmd += "'extension', 'mrc', "  # Dynamo that the data folder uses an extension different to the default .em
        cmd += "'v', true"  # Verbose
        cmd += ")"
        with open(self._getExtraPath(codeFileName), 'w') as cFile:
            cFile.write(cmd)

        Plugin.runDynamo(self, codeFileName, cwd=self._getExtraPath())

    # --------------------------- INFO functions ------------------------------

    # --------------------------- UTILS functions ----------------------------
    def getOutputFile(self):
        return join(self.averageDirName, self.averageDirName + '.mrc')
