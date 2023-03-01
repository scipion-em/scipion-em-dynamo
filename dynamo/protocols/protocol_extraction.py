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

import glob
from os.path import abspath, join, basename
from dynamo.protocols.protocol_base_dynamo import DynamoProtocolBase
from dynamo.utils import getCatalogFile
from pwem.emlib.image import ImageHandler
from pwem.objects import Transform
from pyworkflow.protocol import PointerParam, EnumParam, FloatParam, IntParam, BooleanParam, LEVEL_ADVANCED
from pyworkflow.utils import removeExt
from tomo.objects import SetOfSubTomograms, SubTomogram
import tomo.constants as const
from dynamo import Plugin, VLL_FILE
from dynamo.convert import matrix2eulerAngles

# Tomogram type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1


class DynamoExtraction(DynamoProtocolBase):
    """Extraction of subtomograms using Dynamo"""

    _label = 'vectorial extraction'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.dynamoTomoIdDict = None
        self.scaleFactor = None
        self.tomoTsIdDict = None
        self.coordsFileName = None
        self.anglesFileName = None
        self.cropDirName = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', PointerParam,
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
        form.addParam('inputTomograms', PointerParam,
                      pointerClass='SetOfTomograms',
                      condition='tomoSource != %s' % SAME_AS_PICKING,
                      label='Input tomogram',
                      help='Select the tomogram from which to extract.')
        form.addParam('boxSize', IntParam,
                      label='Box size',
                      help='The subtomograms are extracted as a cubic box of this size. '
                           'The wizard will select the box size considering the sampling rate ratio between the '
                           'introduced coordinates and the tomograms that will br used for the extraction.')
        form.addParam('downFactor', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      default=1.0,
                      label='Downsampling factor',
                      help='It can be used to directly interpolate the extracted subtomograms to '
                           'another size that may correspond to the size of a tomogram that has not been '
                           'reconstructed and this way it is not necessary to reconstruct it. '
                           'If 1.0 is used, no downsampling is applied. Non-integer downsampling factors are '
                           'allowed.\n\nIMPORTANT: the box size of the generated particles will be the result of'
                           'multiplying the introduced box size by this factor.')
        form.addSection(label='Postprocess')
        form.addParam('doInvert', BooleanParam,
                      default=False,
                      label='Invert contrast?',
                      help='Invert the contrast if your tomogram is black '
                           'over a white background.  Xmipp, Spider, Relion '
                           'and Eman require white particles over a black '
                           'background.')
        form.addParallelSection(threads=4, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.writeSetOfCoordinates3D)
        self._insertFunctionStep(self.launchDynamoExtractStep)
        self._insertFunctionStep(self.invertContrastStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.cropDirName = self._getExtraPath('Crop')
        # Get the intersection between the tomograms and coordinates provided (this covers possible subsets made)
        inTomos = self.getInputTomograms()
        inCoords = self.inputCoordinates.get()
        coordsPresentTomoIds = inCoords.getUniqueValues(['_tomoId'])
        tomosPresentTsIds = inTomos.getUniqueValues(['_tsId'])
        commonTomoIds = list(set(coordsPresentTomoIds).intersection(set(tomosPresentTsIds)))
        self.tomoTsIdDict = {tomo.getTsId(): tomo for tomo in inTomos if tomo.getTsId() in commonTomoIds}
        self.coordsFileName = self._getExtraPath('coords.txt')
        self.anglesFileName = self._getExtraPath('angles.txt')
        # Calculate the scale factor
        samplingRateCoord = inCoords.getSamplingRate()
        samplingRateTomo = inTomos.getFirstItem().getSamplingRate()
        self.scaleFactor = float(samplingRateCoord / samplingRateTomo)

    def writeSetOfCoordinates3D(self):
        # Write the VLL file, the coords file and the angles file as expected by Dynamo
        listTomoFile = self._getExtraPath(VLL_FILE)
        self.dynamoTomoIdDict = {}
        with open(self.coordsFileName, "w") as outC, open(self.anglesFileName, 'w') as outA, \
                open(listTomoFile, 'w') as tomoFid:
            for idt, tomo in enumerate(self.tomoTsIdDict.values()):
                self.dynamoTomoIdDict[idt + 1] = tomo
                tomoPath = abspath(tomo.getFileName())
                tomoFid.write(tomoPath + '\n')
                for coord3DSet in self.inputCoordinates.get().iterCoordinates(tomo):
                    angles_coord = matrix2eulerAngles(coord3DSet.getMatrix())
                    x = self.scaleFactor * coord3DSet.getX(const.BOTTOM_LEFT_CORNER)
                    y = self.scaleFactor * coord3DSet.getY(const.BOTTOM_LEFT_CORNER)
                    z = self.scaleFactor * coord3DSet.getZ(const.BOTTOM_LEFT_CORNER)
                    outC.write("%.2f\t%.2f\t%.2f\t%i\n" % (x, y, z, idt + 1))
                    outA.write("%.2f\t%.2f\t%.2f\n" % (angles_coord[0], angles_coord[1], angles_coord[2]))

    def launchDynamoExtractStep(self):
        codeFilePath = self.writeMatlabCode()
        args = ' %s' % codeFilePath
        self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())
        # pwutils.cleanPattern('*.m')

    def invertContrastStep(self):
        import xmipp3
        program = 'xmipp_image_operate'
        for subTomoFile in glob.glob(join(self.cropDirName + '*/**.mrc')):
            args = "-i %s --mult -1" % subTomoFile
            self.runJob(program, args, env=xmipp3.Plugin.getEnviron())

    def createOutputStep(self):
        outSubtomos = SetOfSubTomograms.create(self._getPath(), template='submograms%s.sqlite')
        finalSRate = self.getInputTomograms().getSamplingRate() / self.downFactor.get()
        outSubtomos.setSamplingRate(finalSRate)
        outSubtomos._coordsPointer = self.inputCoordinates
        if self.getInputTomograms().getFirstItem().getAcquisition():
            outSubtomos.setAcquisition(self.getInputTomograms().getFirstItem().getAcquisition())

        for ind in range(len(self.dynamoTomoIdDict.values())):
            currentParticlesDir = join(self.cropDirName + str(ind + 1))
            currentTomo = self.dynamoTomoIdDict[ind + 1]
            currentTomoFName = currentTomo.getFileName()
            currentSubtomoFiles = glob.glob(join(currentParticlesDir, '*.mrc'))
            currentCoords = [coord.clone() for coord in self.inputCoordinates.get().iterCoordinates(currentTomo)]
            for i, subtomoFile in enumerate(currentSubtomoFiles):
                subtomogram = SubTomogram()
                transform = Transform()
                currentCoord = currentCoords[i]
                subtomogram.setSamplingRate(finalSRate)
                dfactor = self.downFactor.get()
                if dfactor != 1:
                    ImageHandler.scaleSplines(subtomoFile + ':mrc', subtomoFile, dfactor)
                subtomogram.setFileName(subtomoFile)
                subtomogram.setVolName(currentTomoFName)
                subtomogram.setCoordinate3D(currentCoord)
                transform.setMatrix(currentCoord.getMatrix())
                subtomogram.setTransform(transform)
                outSubtomos.append(subtomogram)

        if len(outSubtomos) == 0:
            raise "No particles were generated. Please check if Crop directories exist in this protocol's extra " \
                  "folder. If it's the case, check that the tomograms from were the particles are desired to be " \
                  "extracted are the ones in which the picking was performed, or a binned version of them."

        self._defineOutputs(**{super()._possibleOutputs.subtomograms.name: outSubtomos})
        self._defineSourceRelation(self.inputCoordinates, outSubtomos)

    # --------------------------- DEFINE utils functions ----------------------
    def getInputTomograms(self):
        """ Return the tomogram associated to the 'SetOfCoordinates3D' or 'Other' tomograms. """
        if self.tomoSource.get() == SAME_AS_PICKING:
            return self.inputCoordinates.get().getPrecedents()
        else:
            return self.inputTomograms.get()

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
            content += "parfor(tag=unique(tags), %i)\n" % self.numberOfThreads.get()
            content += "tomoCoords = coords(tags == tag, :)\n"
            content += "tomoAngles = angles(tags == tag, :)\n"
            content += "t = dynamo_table_blank(size(tomoCoords, 1), 'r', tomoCoords, 'angles', tomoAngles)\n"
            content += "dtcrop(c.volumes{tag}.fullFileName, t, strcat(savePath, num2str(tag)), box, 'ext', 'mrc')\n"
            content += "end\n"
            codeFid.write(content)

        return codeFilePath

    def readSetOfSubTomograms(self, workDir, outputSubTomogramsSet, coordSet):
        coords = self.inputCoordinates.get()
        for ids, subTomoFile in enumerate(sorted(glob.glob(join(workDir, '*' + '.mrc')))):
            subtomogram = SubTomogram()
            subtomogram.cleanObjId()
            subtomogram.setLocation(subTomoFile)
            dfactor = self.downFactor.get()
            if dfactor != 1:
                fnSubtomo = self._getExtraPath(basename(workDir.strip("/")) + "_downsampled_subtomo%d.mrc" % (ids + 1))
                ImageHandler.scaleSplines(subtomogram.getFileName() + ':mrc', fnSubtomo, dfactor)
                subtomogram.setLocation(fnSubtomo)
            subtomogram.setCoordinate3D(coords[coordSet[ids]])
            subtomogram.setVolName(coords[coordSet[ids]].getVolName())
            transform = Transform()
            transform.setMatrix(coords[coordSet[ids]].getMatrix())
            subtomogram.setTransform(transform)
            outputSubTomogramsSet.append(subtomogram)
        return outputSubTomogramsSet

    # --------------------------- DEFINE info functions ----------------------
    def _methods(self):
        methodsMsgs = []
        if self.getOutputsSize() >= 1:
            msg = ("A total of %i subtomograms of size %i were extracted"
                   % (self.inputCoordinates.get().getSize(), self.boxSize.get()))
            if self.tomoSource.get() == OTHER:
                msg += (" from another set of tomograms: %s"
                        % self.getObjectTag('inputTomogram'))
            msg += " using coordinates %s" % self.getObjectTag('inputCoordinates')
            msg += self.methodsVar.get('')
            methodsMsgs.append(msg)
        else:
            methodsMsgs.append("Set of Subtomograms not ready yet")
        if self.downFactor.get() != 1:
            methodsMsgs.append("Subtomograms downsample by factor %d."
                               % self.downFactor.get())
        return methodsMsgs

    def _summary(self):
        summary = ["Tomogram source: *%s*" % self.getEnumText("tomoSource")]
        # if self.getOutputsSize() >= 1:
        #     summary.append("Particle box size: *%s*" % self.boxSize.get())
        #     summary.append("Subtomogram extracted: *%s*" % self.inputCoordinates.get().getSize())
        # else:
        #     summary.append("Output subtomograms not ready yet.")
        if self.doInvert:
            summary.append('*Contrast was inverted.*')
        return summary
