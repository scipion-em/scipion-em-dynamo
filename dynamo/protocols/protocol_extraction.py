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
import copy
import glob
from enum import Enum
from os.path import abspath, join, dirname, basename

from dynamo.protocols.protocol_base_dynamo import IN_TOMOS, IN_COORDS, DynamoProtocolBase
from dynamo.utils import getCatalogFile
from pwem.objects import Transform
from pwem.protocols import EMProtocol
from pyworkflow.object import String
from pyworkflow.protocol import PointerParam, EnumParam, IntParam, BooleanParam
from pyworkflow.utils import removeExt, Message
from tomo.constants import BOTTOM_LEFT_CORNER, TR_DYNAMO
from tomo.objects import SetOfSubTomograms, SubTomogram, Coordinate3D, Tomogram
from dynamo import Plugin, VLL_FILE
from dynamo.convert import matrix2eulerAngles
from tomo.protocols import ProtTomoBase
from tomo.utils import scaleTrMatrixShifts


CROP_DIR = 'Crop'

# Tomogram type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1


class DynSubtomoExtractOuts(Enum):
    subtomograms = SetOfSubTomograms


class DynamoExtraction(DynamoProtocolBase):
    """Extraction of subtomograms using Dynamo"""

    _label = 'subtomogram extraction'
    _possibleOutputs = DynSubtomoExtractOuts

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.dynamoTomoIdDict = {}
        self.scaleFactor = None
        self.tomoTsIdDict = None
        self.coordsFileName = None
        self.anglesFileName = None
        self.cropDirName = None
        self.coordList = []
        self.removedCoordsIndices = String()

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_COORDS, PointerParam,
                      label="Input Coordinates",
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
                      label='Input tomogram',
                      help='Select the tomogram from which to extract.')
        form.addParam('boxSize', IntParam,
                      label='Box size',
                      help='The subtomograms are extracted as a cubic box of this size. '
                           'The wizard will select the box size considering the sampling rate ratio between the '
                           'introduced coordinates and the tomograms that will br used for the extraction.')
        form.addSection(label='Postprocess')
        form.addParam('doInvert', BooleanParam,
                      default=False,
                      label='Invert contrast?',
                      help='Invert the contrast if your tomogram is black '
                           'over a white background.  Xmipp, Spider, Relion '
                           'and Eman require white particles over a black '
                           'background.')
        self.insertBinThreads(form)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.writeSetOfCoordinates3D,
                                 needsGPU=False)
        self._insertFunctionStep(self.launchDynamoExtractStep,
                                 needsGPU=False)
        if self.doInvert.get():
            self._insertFunctionStep(self.invertContrastStep,
                                     needsGPU=False)
        self._insertFunctionStep(self.createOutputStep,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.cropDirName = self._getExtraPath(CROP_DIR)
        # Get the intersection between the tomograms and coordinates provided (this covers possible subsets made)
        inTomos = self.getInputTomograms()
        inCoords = self.getInCoords()
        coordsPresentTomoIds = inCoords.getTSIds()
        tomosPresentTsIds = inTomos.getTSIds()
        commonTomoIds = set(coordsPresentTomoIds).intersection(set(tomosPresentTsIds))
        self.tomoTsIdDict = {tomo.getTsId(): tomo.clone() for tomo in inTomos if tomo.getTsId() in commonTomoIds}
        self.coordsFileName = self._getExtraPath('coords.txt')
        self.anglesFileName = self._getExtraPath('angles.txt')
        # Calculate the scale factor
        samplingRateCoord = inCoords.getSamplingRate()
        samplingRateTomo = inTomos.getFirstItem().getSamplingRate()
        self.scaleFactor = float(samplingRateCoord / samplingRateTomo)

    def writeSetOfCoordinates3D(self):
        # Write the VLL file, the coords file and the angles file as expected by Dynamo
        listTomoFile = self._getExtraPath(VLL_FILE)
        tomoDims = self.getInputTomograms().getDimensions()
        removedCoordsIndices = []
        with open(self.coordsFileName, "w") as outC, open(self.anglesFileName, 'w') as outA, \
                open(listTomoFile, 'w') as tomoFid:
            for idt, tomo in enumerate(self.tomoTsIdDict.values()):
                self.dynamoTomoIdDict[idt + 1] = tomo.getFileName()
                tomoPath = abspath(tomo.getFileName())
                tomoFid.write(tomoPath + '\n')
                for coord in self.getInCoords().iterCoordinates(tomo):
                    angles_coord = matrix2eulerAngles(coord.getMatrix())
                    x = self.scaleFactor * coord.getX(BOTTOM_LEFT_CORNER)
                    y = self.scaleFactor * coord.getY(BOTTOM_LEFT_CORNER)
                    z = self.scaleFactor * coord.getZ(BOTTOM_LEFT_CORNER)
                    if self.isParticleOutOfTomo((x, y, z), tomoDims):
                        removedCoordsIndices.append(coord.getObjId())
                        continue
                    outC.write("%.2f\t%.2f\t%.2f\t%i\n" % (x, y, z, idt + 1))
                    outA.write("%.2f\t%.2f\t%.2f\n" % (angles_coord[0], angles_coord[1], angles_coord[2]))
                    # Add the coordinates of each in-tomo particle to a list that will be used in the CreateOutputStep
                    # to avoid index mismatching by preventing the particles out of the tomograms to be sent to be
                    # processed by Dynamo
                    self.coordList.append(coord.clone())
        if not self.coordList:
            raise Exception("All the coordinates were removed (This is because the box associated to those "
                            "coordinates partially or totally lay out of the corresponding tomogram. A good way to "
                            "avoid this is to decreae the box size)")
        if self.removedCoordsIndices:
            self.removedCoordsIndices.set(" ".join(map(str, removedCoordsIndices)))
            self._store(self.removedCoordsIndices)

    def launchDynamoExtractStep(self):
        codeFilePath = self.writeMatlabCode()
        args = ' %s' % codeFilePath
        self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    def invertContrastStep(self):
        import xmipp3
        program = 'xmipp_image_operate'
        for subTomoFile in glob.glob(join(self.cropDirName + '*/**.mrc')):
            args = "-i %s --mult -1" % subTomoFile
            self.runJob(program, args, env=xmipp3.Plugin.getEnviron())

    def createOutputStep(self):
        outSubtomos = SetOfSubTomograms.create(self._getPath(), template='submograms%s.sqlite')
        finalSRate = self.getInputTomograms().getSamplingRate()
        outSubtomos.setSamplingRate(finalSRate)
        outSubtomos._coordsPointer = self.getInCoords(isPointer=True)
        if self.getInputTomograms().getFirstItem().getAcquisition():
            outSubtomos.setAcquisition(self.getInputTomograms().getFirstItem().getAcquisition())
        currentSubtomoFiles = sorted(glob.glob(self._getExtraPath('**', '*.mrc')))
        for currentCoord, subtomoFile in zip(self.coordList, currentSubtomoFiles):
            subtomogram = SubTomogram()
            transform = Transform()
            subtomogram.setSamplingRate(finalSRate)
            subtomogram.setFileName(subtomoFile)
            subtomogram.setVolName(self.getTomogramFileFromSubtomoFile(subtomoFile))
            subtomogram.setCoordinate3D(currentCoord)
            trMatrix = copy.copy(currentCoord.getMatrix())
            transform.setMatrix(scaleTrMatrixShifts(trMatrix, self.scaleFactor))
            subtomogram.setTransform(transform, convention=TR_DYNAMO)
            outSubtomos.append(subtomogram)

        if len(outSubtomos) == 0:
            raise "No particles were generated. Please check if Crop directories exist in this protocol's extra " \
                  "folder. If it's the case, check that the tomograms from were the particles are desired to be " \
                  "extracted are the ones in which the picking was performed, or a binned version of them."

        self._defineOutputs(**{self._possibleOutputs.subtomograms.name: outSubtomos})
        self._defineSourceRelation(self.getInCoords(isPointer=True), outSubtomos)

    # --------------------------- DEFINE utils functions ----------------------
    def getInputTomograms(self):
        """ Return the tomogram associated to the 'SetOfCoordinates3D' or 'Other' tomograms. """
        if self.tomoSource.get() == SAME_AS_PICKING:
            return self.getInCoords().getPrecedents()
        else:
            return self.getInTomos()

    def writeMatlabCode(self):
        codeFilePath = self._getExtraPath("DynamoExtraction.m")
        # Write code to Matlab code file
        with open(codeFilePath, 'w') as codeFid:
            catalogue = getCatalogFile(self._getExtraPath())
            content = "savePath = '%s'\n" % self.cropDirName
            content += "box = %i\n" % self.boxSize.get()
            content += "dcm -create '%s' -fromvll '%s'\n" % (removeExt(catalogue), self._getExtraPath(VLL_FILE))
            content += "c = dread('%s')\n" % catalogue
            content += "coordsData = readmatrix('%s')\n" % self.coordsFileName
            content += "angles = readmatrix('%s')\n" % self.anglesFileName
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

    def isParticleOutOfTomo(self, scaledXYZ: tuple, tomoDims: tuple) -> bool:
        """Dynamo skips the particles whose box size or part of it is out of the tomogram. For example if a particle is
        located 10 voxels from one of the tomogram edges and the box size introduced is 30 voxels (15 up and down, left
        and right, from the coordinate) the part of the box corresponding to the 5 voxels out of the tomogram would
        cut out of the tomogram and Dynamo would skip it. Coordinates are assumed to be at the same size scale as the
        tomograms from which they are going to be extracted

        :param scaledXYZ: a tuple containing the properly scaled x, y, and z coordinates of a coordinate.
        :param tomoDims: a tuple with the width, height and thickness of the tomograms from which the particles will be
        extracted."""
        halfBoxSize = self.boxSize.get() / 2
        coordX, coordY, coordZ = scaledXYZ[:]
        tomoWidth, tomoHeight, tomoThickness = tomoDims[:]
        outCheckList = [abs(coordX) + halfBoxSize > tomoWidth,
                        abs(coordY) + halfBoxSize > tomoHeight,
                        abs(coordZ) + halfBoxSize > tomoThickness]
        return True if any(outCheckList) else False

    def getTomogramFileFromSubtomoFile(self, subtomoFile):
        """Dynamo generates, for each tomogram, a directory named CropN, where N is the number of the tomogram (stored
         by Scipion as keys in self.dynamoTomoIdDict). This method gets N from the directory name in which the given
         subtomogram filename is stored and uses it to get the tomgram value corresponding to that key from the
         mentioned dictionary."""
        return self.dynamoTomoIdDict[int(basename(dirname(subtomoFile)).replace(CROP_DIR, ''))]

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
            inCoordsSize = self.getInCoords().getSize()
            outCoordsSize = getattr(self, self._possibleOutputs.subtomograms.name).getSize()
            summary.append("Tomogram source: *%s*" % self.getEnumText("tomoSource"))
            if inCoordsSize > outCoordsSize:
                summary.append('*%i coordinates were removed* (This is because the box associated to those '
                               'coordinates partially or totally lay out of the corresponding tomogram. A good way to '
                               'avoid this is to decreae the box size).\nRemoved coordinates ids: %s.' %
                               (inCoordsSize - outCoordsSize, self.removedCoordsIndices.get()))

            if self.doInvert:
                summary.append('*Contrast was inverted.*')
        return summary
