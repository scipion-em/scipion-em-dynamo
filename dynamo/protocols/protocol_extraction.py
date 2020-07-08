# **************************************************************************
# *
# * Authors:    David Herreros Calero (dherreros@cnb.csic.es)
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
import glob

import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
from pwem.emlib.image import ImageHandler
from pwem.objects import Transform
from pwem.protocols import EMProtocol

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfSubTomograms, SubTomogram, TomoAcquisition

from dynamo import Plugin
from dynamo.convert import matrix2eulerAngles


# Tomogram type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1


class DynamoExtraction(EMProtocol, ProtTomoBase):
    """Extraction of subtomograms using Dynamo"""

    _label = 'vectorial extraction'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoordinates', params.PointerParam, label="Input Coordinates", important=True,
                      pointerClass='SetOfCoordinates3D', help='Select the SetOfCoordinates3D.')
        form.addParam('tomoSource', params.EnumParam,
                      choices=['same as picking', 'other'],
                      default=0,
                      display=params.EnumParam.DISPLAY_HLIST,
                      label='Tomogram source',
                      help='By default the subtomograms will be extracted '
                           'from the tomogram used in the picking '
                           'step ( _same as picking_ option ). \n'
                           'If you select _other_ option, you must provide '
                           'a different tomogram to extract from. \n'
                           '*Note*: In the _other_ case, ensure that provided '
                           'tomogram and coordinates are related ')

        form.addParam('inputTomograms', params.PointerParam,
                      pointerClass='SetOfTomograms',
                      condition='tomoSource != %s' % SAME_AS_PICKING,
                      label='Input tomogram',
                      help='Select the tomogram from which to extract.')

        form.addParam('boxSize', params.FloatParam,
                      label='Box size',
                      help='The subtomograms are extracted as a cubic box of this size.\n'
                           'The wizard selects same box size as picking.\n'
                           'Dynamo only accepts even numbers as box sizes. Otherwise, an exception will '
                           'be returned.')

        form.addParam('downFactor', params.FloatParam, default=1.0,
                      label='Downsampling factor',
                      help='If 1.0 is used, no downsample is applied. '
                           'Non-integer downsample factors are possible. ')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('writeSetOfCoordinates3D')
        self._insertFunctionStep('launchDynamoExtractStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def writeSetOfCoordinates3D(self):
        samplingRateCoord = self.inputCoordinates.get().getSamplingRate()
        samplingRateTomo = self.getInputTomograms().getFirstItem().getSamplingRate()
        self.factor = float(samplingRateCoord / samplingRateTomo)
        self.lines = []
        inputSet = self.getInputTomograms().getFiles()
        self.coordsFileName = self._getExtraPath('coords.txt')
        self.anglesFileName = self._getExtraPath('angles.txt')
        outC = open(self.coordsFileName, "w")
        outA = open(self.anglesFileName, 'w')
        for idt, item in enumerate(inputSet):
            coordDict = []
            for coord3DSet in self.inputCoordinates.get().iterCoordinates():
                if pwutils.removeBaseExt(item) == pwutils.removeBaseExt(coord3DSet.getVolName()):
                    angles_coord = matrix2eulerAngles(coord3DSet.getMatrix())
                    x = round(self.factor * coord3DSet.getX())
                    y = round(self.factor * coord3DSet.getY())
                    z = round(self.factor * coord3DSet.getZ())
                    outC.write("%d\t%d\t%d\t%d\n" % (x, y, z, idt+1))
                    outA.write("%f\t%f\t%f\n" % (angles_coord[0], angles_coord[1], angles_coord[2]))
                    coordDict.append(coord3DSet.getObjId())
            if coordDict:
                self.lines.append(coordDict)
        outC.close()
        outA.close()

    def launchDynamoExtractStep(self):
        codeFilePath = self.writeMatlabCode()

        args = ' %s' % codeFilePath
        self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

        pwutils.cleanPattern('*.m')

    def createOutputStep(self):
        self.outputSubTomogramsSet = self._createSetOfSubTomograms(self._getOutputSuffix(SetOfSubTomograms))
        self.outputSubTomogramsSet.setSamplingRate(self.getInputTomograms().getSamplingRate() * self.downFactor.get())
        self.outputSubTomogramsSet.setCoordinates3D(self.inputCoordinates)
        acquisition = TomoAcquisition()
        acquisition.setAngleMin(self.getInputTomograms().getFirstItem().getAcquisition().getAngleMin())
        acquisition.setAngleMax(self.getInputTomograms().getFirstItem().getAcquisition().getAngleMax())
        acquisition.setStep(self.getInputTomograms().getFirstItem().getAcquisition().getStep())
        self.outputSubTomogramsSet.setAcquisition(acquisition)
        for ind, tomoFile in enumerate(self.tomoFiles):
            cropPath = os.path.join(self._getExtraPath('Crop%d' % (ind+1)), '')
            if os.path.isdir(cropPath):
                coordSet = self.lines[ind]
                outputSet = self.readSetOfSubTomograms(cropPath, self.outputSubTomogramsSet, coordSet)

        self._defineOutputs(outputSetOfSubtomogram=outputSet)
        self._defineSourceRelation(self.inputCoordinates, outputSet)

    # --------------------------- DEFINE utils functions ----------------------
    def getInputTomograms(self):
        """ Return the tomogram associated to the 'SetOfCoordinates3D' or 'Other' tomograms. """
        if self.tomoSource.get() == SAME_AS_PICKING:
            return self.inputCoordinates.get().getPrecedents()
        else:
            return self.inputTomograms.get()

    def writeMatlabCode(self):
        # Initialization params
        codeFilePath = os.path.join(os.getcwd(), "DynamoExtraction.m")
        listTomosFile = self._getTmpPath("tomos.vll")
        catalogue = os.path.abspath(self._getExtraPath("tomos"))
        self.tomoFiles = sorted(self.getInputTomograms().getFiles())

        # Create list of tomos file
        tomoFid = open(listTomosFile, 'w')
        for tomoFile in self.tomoFiles:
            tomoPath = os.path.abspath(tomoFile)
            tomoFid.write(tomoPath + '\n')
        tomoFid.close()

        # Check of box size is an even number (as required by Dynamo)
        boxSize = round(self.factor * self.boxSize.get())
        if boxSize % 2 != 0:
            boxSize += 1

        # Write code to Matlab code file
        codeFid = open(codeFilePath, 'w')
        content = "path='%s'\n" \
                  "savePath = '%s'\n" \
                  "catalogue_name='%s'\n" \
                  "box=%d\n" \
                  "dcm -create %s -fromvll %s\n" \
                  "c=dread(strcat(catalogue_name,'.ctlg'))\n" \
                  "aux=readmatrix('%s')\n" \
                  "angles=readmatrix('%s')\n" \
                  "coords=aux(:,1:3)\n" \
                  "tags=aux(:,4)'\n" \
                  "for tag=unique(tags)\n" \
                  "tomoCoords=coords(tags==tag,:)\n" \
                  "tomoAngles=angles(tags==tag,:)\n" \
                  "t=dynamo_table_blank(size(tomoCoords,1),'r',tomoCoords,'angles',rad2deg(tomoAngles))\n" \
                  "dtcrop(c.volumes{tag}.fullFileName,t,strcat(savePath,num2str(tag)),box)\n" \
                  "end\n" \
                   % (os.path.abspath(os.getcwd()), self._getExtraPath('Crop'),
                      catalogue,boxSize, catalogue, listTomosFile, os.path.abspath(self.coordsFileName),
                      os.path.abspath(self.anglesFileName))

        codeFid.write(content)
        codeFid.close()

        return codeFilePath

    def readSetOfSubTomograms(self, workDir, outputSubTomogramsSet, coordSet):
        coords = self.inputCoordinates.get()
        for ids, subTomoFile in enumerate(sorted(glob.glob(os.path.join(workDir, '*.em')))):
            subtomogram = SubTomogram()
            subtomogram.cleanObjId()
            subtomogram.setLocation(subTomoFile)
            dfactor = self.downFactor.get()
            if dfactor != 1:
                fnSubtomo = self._getExtraPath(os.path.basename(workDir.strip("/")) + "_downsampled_subtomo%d.mrc" % (ids+1))
                ImageHandler.scaleSplines(subtomogram.getLocation(), fnSubtomo, dfactor)
                subtomogram.setLocation(fnSubtomo)
            subtomogram.setCoordinate3D(coords[coordSet[ids]])
            subtomogram.setVolName(coords[coordSet[ids]].getVolName())
            transform = Transform()
            transform.setMatrix(coords[coordSet[ids]].getMatrix())
            subtomogram.setTransform(transform)
            outputSubTomogramsSet.append(subtomogram)
        return outputSubTomogramsSet

    # --------------------------- DEFINE info functions ----------------------
    def _validate(self):
        errors = []
        if not self.boxSize.get() % 2 == 0:
            errors.append('The box size provided to Dynamo is an odd number. '
                          'This will rise an execption when Dynamo is invoked.\n'
                          'Please, provide an even number to continue with the extraction.')
        return errors

    def _tomosOther(self):
        """ Return True if other tomograms are used for extract. """
        return self.tomoSource == OTHER

    def _methods(self):
        methodsMsgs = []
        if self.getOutputsSize() >= 1:
            msg = ("A total of %s subtomograms of size %s were extracted"
                   % (str(self.inputCoordinates.get().getSize()), self.boxSize.get()))
            if self._tomosOther():
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
        summary = []
        summary.append("Tomogram source: %s"
                       % self.getEnumText("tomoSource"))
        if self.getOutputsSize() >= 1:
            summary.append("Particle box size: %s" % self.boxSize.get())
            summary.append("Subtomogram extracted: %s" %
                           self.inputCoordinates.get().getSize())
        else:
            summary.append("Output subtomograms not ready yet.")
        return summary
