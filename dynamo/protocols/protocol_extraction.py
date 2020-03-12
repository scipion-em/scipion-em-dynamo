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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

import os
import numpy as np

import pyworkflow.utils as pwutils
import pyworkflow.em as pwem
import pyworkflow.protocol.params as params

from tomo.protocols import ProtTomoBase
from tomo.objects import SetOfCoordinates3D, Coordinate3D

from dynamo import Plugin


# Tomogram type constants for particle extraction
SAME_AS_PICKING = 0
OTHER = 1


class DynamoExtraction(pwem.EMProtocol, ProtTomoBase):
    """Extraction of subtomograms using Dynamo"""

    _label = 'dynamo extract'

    def __init__(self, **kwargs):
        pwem.EMProtocol.__init__(self, **kwargs)

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
                           'In Dynamo, this parameter must be an even number.')

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
        self.lines = []
        self.tomoFiles = []
        inputSet = self.getInputTomograms()
        self.coordsFileName = self._getExtraPath('coords.txt')
        self.anglesFileName = self._getExtraPath('angles.txt')
        outC = file(self.coordsFileName, "w")
        outA = file(self.anglesFileName, 'w')
        for idt, item in enumerate(inputSet):
            coordDict = []
            tomo = item.clone()
            for coord3DSet in self.inputCoordinates.get().iterCoordinates():
                if os.path.basename(tomo.getFileName()) == os.path.basename(coord3DSet.getVolName()):
                    angles_coord = coord3DSet.eulerAngles()
                    outC.write("%d\t%d\t%d\t%d\n" % (coord3DSet.getX(), coord3DSet.getY(), coord3DSet.getZ(), idt+1))
                    outA.write("%f\t%f\t%f\n" % (angles_coord[0], angles_coord[1], angles_coord[2]))
                    coordDict.append(coord3DSet.clone())
            if coordDict:
                self.lines.append(coordDict)
                self.tomoFiles.append(tomo.getFileName())
                self.samplingRateTomo = tomo.getSamplingRate()
        outC.close()
        outA.close()

    def launchDynamoExtractStep(self):
        codeFilePath = self.writeMatlabCode()

        args = ' %s' % codeFilePath
        self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

        pwutils.cleanPattern('*.m')

    def createOutputStep(self):
        pass

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
        listTomosFile = os.path.join(os.environ.get("SCIPION_HOME"), "software", "tmp", "tomos.vll")
        catalogue = os.path.abspath(self._getExtraPath("tomos"))

        # Create list of tomos file
        tomoFid = open(listTomosFile, 'w')
        for tomoFile in self.tomoFiles:
            tomoPath = os.path.abspath(tomoFile)
            tomoFid.write(tomoPath + '\n')
        tomoFid.close()

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
                      catalogue, self.boxSize.get(), catalogue, listTomosFile, os.path.abspath(self.coordsFileName),
                      os.path.abspath(self.anglesFileName))

        codeFid.write(content)
        codeFid.close()

        return codeFilePath

    # --------------------------- DEFINE info functions ----------------------
    def getMethods(self, output):
        msg = 'User picked %d particles ' % output.getSize()
        msg += 'with a particle size of %s.' % output.getBoxSize()
        return msg

    def _methods(self):
        methodsMsgs = []
        if self.inputTomograms is None:
            return ['Input tomogram not available yet.']

        methodsMsgs.append("Input tomograms imported of dims %s." %(
                              str(self.inputTomograms.get().getDim())))

        if self.getOutputsSize() >= 1:
            for key, output in self.iterOutputAttributes():
                msg = self.getMethods(output)
                methodsMsgs.append("%s: %s" % (self.getObjectTag(output), msg))
        else:
            methodsMsgs.append(Message.TEXT_NO_OUTPUT_CO)

        return methodsMsgs