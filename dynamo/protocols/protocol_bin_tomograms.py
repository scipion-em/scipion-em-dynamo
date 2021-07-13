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

from pyworkflow import BETA
from pyworkflow.protocol import params, STEPS_PARALLEL
import pyworkflow.utils as pwutils

from pwem.protocols import EMProtocol

from tomo.protocols import ProtTomoBase
from tomo.objects import Tomogram

from dynamo import Plugin


class DynamoBinTomograms(EMProtocol, ProtTomoBase):
    """Reduce the size of a SetOfTomograms by a binning factor"""

    _label = 'bin tomograms'
    _devStatus = BETA

    OUTPUT_PREFIX = 'resizedTomograms'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

   # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputTomos', params.PointerParam, pointerClass='SetOfTomograms',
                      label="Input Tomograms", important=True,
                      help="Select the input tomograms to be resized.")
        form.addParam('binning', params.IntParam, default=2,
                      label="Binning Factor",
                      help="Binning Factor to be applied to the Tomograms. A Binning Factor of 1 will decrease the "
                           "size of the Tomograms by 2, a Binning Factor of 2 by 4...")
        form.addParam('zChunk', params.IntParam, default=300, expertLevel=params.LEVEL_ADVANCED,
                      label="Number of slices kept in memory",
                      help="Maximum number of Z slices that are kept simultaneously in the memory during the "
                           "binning process. This parameter might be important for larger size tomograms.")
        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        binning_steps = []
        self.inputTomos = self.inputTomos.get()
        for tomo in self.inputTomos:
            binning_steps.append(self._insertFunctionStep('binTomosStep',
                                                          tomo.getFileName(), tomo.getObjId()))
        self._insertFunctionStep('createOutputStep', prerequisites=binning_steps)

    # --------------------------- STEPS functions -----------------------------
    def binTomosStep(self, tomoName, tomoId):
        commandsFile = self.writeMatlabFile(os.path.abspath(tomoName), tomoId)
        args = ' %s' % commandsFile
        self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    def createOutputStep(self):
        # Create a Volume template object
        binning = self.binning.get()
        sr = self.inputTomos.getSamplingRate() * 2 ** binning
        tomo = Tomogram()
        tomo.setSamplingRate(sr)
        binned_tomos = self._createSetOfTomograms()
        binned_tomos.setSamplingRate(sr)

        for inTomo in self.inputTomos:
            inTomo_file = inTomo.getFileName()
            binned_scipion_name = pwutils.removeBaseExt(inTomo_file) + '_binned{}'.format(binning) + \
                                  pwutils.getExt(inTomo_file)
            binned_scipion_path = self._getExtraPath(binned_scipion_name)
            tomo.cleanObjId()
            tomo.setLocation(binned_scipion_path)
            tomo.setOrigin(inTomo.getOrigin())
            binned_tomos.append(tomo)

        args = {self.OUTPUT_PREFIX: binned_tomos}
        self._defineOutputs(**args)
        self._defineSourceRelation(self.inputTomos, binned_tomos)

    # --------------------------- DEFINE utils functions ----------------------
    def writeMatlabFile(self, tomoPath, tomoId):
        codeFilePath = self._getExtraPath('bin_Tomo_%d.m' % tomoId)
        binned_scipion_name = pwutils.removeBaseExt(tomoPath) + '_binned{}'.format(self.binning.get()) + \
                              pwutils.getExt(tomoPath)
        binned_scipion_path = os.path.abspath(self._getExtraPath(binned_scipion_name))
        content = "dpktomo.tools.bin('%s','%s',%d,'ss',%d)\n" \
                  "exit\n" % (tomoPath, binned_scipion_path, self.binning.get(), self.zChunk.get())
        codeFid = open(codeFilePath, 'w')
        codeFid.write(content)
        codeFid.close()
        return codeFilePath

    # --------------------------- DEFINE info functions ----------------------
    def _methods(self):
        methodsMsgs = []
        methodsMsgs.append("*Binning Factor*: %s" % self.binning.get())
        return methodsMsgs

    def _summary(self):
        summary = []
        if self.getOutputsSize() >= 1:
            for _, outTomos in self.iterOutputAttributes():
                summary.append("Output *%s*:" % outTomos.getNameId().split('.')[1])
                summary.append("    * Binning Factor: *%s*" % self.binning.get())
                summary.append("    * Number of Tomograms Binned: *%s*" %
                               outTomos.getSize())
        else:
            summary.append("Output tomograms not ready yet.")
        return summary