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
from pwem.emlib.image import ImageHandler
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import PointerParam, BooleanParam, IntParam, FloatParam, StringParam
from pyworkflow.utils import Message, makePath
from tomo.objects import AverageSubTomogram
from tomo.protocols import ProtTomoBase


class DynTempMatchOuts(Enum):
    average = AverageSubTomogram


class DynamoProtTemplateMatching(EMProtocol, ProtTomoBase):
    _label = 'Template matching'
    _devStatus = BETA
    _possibleOutputs = DynTempMatchOuts

    # --------------- DEFINE param functions ---------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inTomos', PointerParam,
                      pointerClass='SetOfTomograms',
                      important=True,
                      label='Tomograms')
        form.addParam('template', PointerParam,
                      pointerClass='Volume',
                      important=True,
                      label='Template')
        form.addParam('mask', PointerParam,
                      pointerClass='VolumeMask',
                      label='Mask (opt.)',
                      allowsNull=True,
                      help='If no mask is given, a spherical mask with the size of the template will be used.')
        form.addParam('yTiltRange', StringParam,
                      default='-60 60',
                      label='Y tilt range')
        form.addParam('sizeChunk', IntParam,
                      default=0,
                      label='Chunk size',
                      help='If a chunk size is provided, the volume will be split internally in parts of the '
                           'specified size. The dimensions of the chunk are referred to the unbinned volume.')
        form.addParam('coneRange', IntParam,
                      default=0,
                      label='Orientation search cone aperture (deg.)',
                      help='Orientations will be looked for inside a cone. In this context, the most usual values '
                           'are 360 (sample the full sphere) or 0 (do not scan orientations).')
        form.addParam('coneSampling', IntParam,
                      default=0,
                      label='Angular interval for the orientation search (deg.)',
                      condition='coneRange > 0',
                      help='If set to 0, no search will be performed.')
        form.addParam('inPlaneRange', IntParam,
                      default=0,
                      label='In plane rotation range (deg.)')
        form.addParam('InPlaneSampling', IntParam,
                      default=0,
                      label='Angular interval for the in plane rotation search (deg.)',
                      condition='inPlaneRange > 0')
        form.addParallelSection(threads=4, mpi=0)

    # --------------- INSERT steps functions ----------------
    def _insertAllSteps(self):
        # self._insertFunctionStep(self.convertInputStep)
        tomoList = [tomo.clone() for tomo in self.inTomos.get()]
        for tomo in tomoList:
            self._insertFunctionStep(self.tempMatchStep, tomo)
        # self._insertFunctionStep(self.convertOutputStep)
        # self._insertFunctionStep(self.createOutputStep)

    # --------------- STEPS functions -----------------------
    # def convertInputStep(self):
    #     inSubtomos = self.inSubtomos.get()
    #     tableName = self._getExtraPath(self.tableName)
    #     dataDir = self._getExtraPath(self.dataDirName)
    #     avgDir = self._getExtraPath(self.dataDirName)
    #     makePath(*[dataDir, avgDir])
    #     # Generate the data folder formatted along the Dynamo convention
    #     writeSetOfVolumes(inSubtomos, join(dataDir, 'particle_'), 'id')
    #     # Generate the Dynamo data table
    #     with open(tableName, 'w') as fhTable:
    #         writeDynTable(fhTable, inSubtomos)

    def tempMatchStep(self, tomo):
        codeFileName = self._getExtraPath('tempMatch.m')
        cmd = self.genTempMatchCmd(tomo.getFileName())
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
        self._defineOutputs(**{DynTempMatchOuts.average.name: avg})
        self._defineSourceRelation(inSubtomos, avg)

    # --------------------------- INFO functions -----------------------------
    def _warnings(self):
        warnMsgs = []
        coneRange = self.coneRange.get()
        coneSampling = self.coneSampling.get()
        inPlaneRange = self.inPlaneRange.get()
        inPlaneSampling = self.InPlaneSampling.get()
        if coneRange > 0 and coneSampling == 0:
            warnMsgs.append('The cone range is greater than 0, but not the cone sampling. A cone sampling value '
                            'greater than 0 must be provided to perform the orientation search.')
        if inPlaneRange > 0 and inPlaneSampling == 0:
            warnMsgs.append('The in plane range is greater than 0, but not the in plane sampling. An in plane sampling '
                            'value greater than 0 must be provided to perform the plane rotation search.')
        return warnMsgs

    # --------------------------- UTILS functions ----------------------------
    def getOutputFile(self, ext='em'):
        return self._getExtraPath(self.averageDirName, self.averageDirName + '.' + ext)

    def genTempMatchCmd(self, tomoFileName):
        """Syntax example of the native code:

        pts = dynamo_match('t20s.mrc','average32.em','mask','maskTight32.em',...
        'outputFolder','cs30', 'ytilt',[-39,36],'sc',[256,256,256],'cr',360,'cs',30,'bin',1);

        Where:

        - 'sc' : size of chunk
          This is the maximum extent of a block in the unbinned tomogram that will be kept in memory by each one of
          the processors working in parallel.
        - 'cr' : cone range
          Orientations will be looked for inside a cone. In this context, the most usual values are 360 (sample the
          full sphere) or 0 (do not scan orientations)
        - 'cs' : cone sampling
          Determines the scanning density inside the sphere.
        - 'ir' : in plane rotation range
        - 'is' : in plane sampling. It determines the scanning density for in plane rotations
        - 'mw' : number of matlab workers for parallel computation
        """
        sizeChunk = self.sizeChunk.get()
        coneRange = self.coneRange.get()
        coneSampling = self.coneSampling.get()
        inPlaneRange = self.inPlaneRange.get()
        inPlaneSampling = self.InPlaneSampling.get()
        cmd = [f"dynamo_match('{tomoFileName}', '{self.template.get().getFileName()}', "
               f"'outputFolder', '{self._getExtraPath()}', " 
               f"'ytilt', {[int(i) for i in self.yTiltRange.get().split(' ')]}, "
               f"'mw', {self.numberOfThreads.get()}"]
        if sizeChunk:
            cmd.append(f", 'sc', [{sizeChunk}, {sizeChunk}, {sizeChunk}]")
        if coneRange > 0 and coneSampling > 0:
            cmd.append(f", 'cr', {coneRange}, 'cs', {coneSampling}")
        if inPlaneRange > 0 and inPlaneSampling > 0:
            cmd.append(f", 'ir', {inPlaneRange}, 'is', {inPlaneSampling}")
        return ' '.join(cmd) + ')'

