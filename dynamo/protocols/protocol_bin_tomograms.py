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
from enum import Enum
from os.path import join

from pwem.emlib.image import ImageHandler
from pyworkflow import BETA
from pyworkflow.protocol import params, STEPS_PARALLEL
import pyworkflow.utils as pwutils

from pwem.protocols import EMProtocol
from pyworkflow.utils import getExt, makePath, removeBaseExt

from tomo.protocols import ProtTomoBase
from tomo.objects import Tomogram, SetOfTomograms

from dynamo import Plugin


class DynBinTomosOutputs(Enum):
    tomograms = SetOfTomograms


class DynamoBinTomograms(EMProtocol, ProtTomoBase):
    """Reduce the size of a SetOfTomograms by a binning factor"""

    _label = 'bin tomograms'
    _devStatus = BETA
    _possibleOutputs = DynBinTomosOutputs
    inTomosDir = 'inTomograms'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

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
        self._initialize()
        if self.doConvertFiles:
            makePath(self.getInTomosDir())
            for tomo in self.inputTomos:
                fName = tomo.getFileName()
                self._insertFunctionStep(self.convertInputStep, fName)
                self._insertFunctionStep(self.binTomosStep, fName, tomo.getObjId())
        else:
            for tomo in self.inputTomos:
                self._insertFunctionStep(self.binTomosStep, tomo.getFileName(), tomo.getObjId())

        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        self.inputTomos = self.inputTomos.get()
        self.doConvertFiles = not self.isCompatibleFileFormat()

    def convertInputStep(self, inFileName):
        ih = ImageHandler()
        ih.convert(inFileName, self.getConvertedOutFileName(inFileName))

    def binTomosStep(self, tomoName, tomoId):
        if self.doConvertFiles:
            tomoName = self.getConvertedOutFileName(tomoName)
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
            tomo.setAcquisition(inTomo.getAcquisition())
            binned_tomos.append(tomo)

        self._defineOutputs(**{DynBinTomosOutputs.tomograms.name: binned_tomos})
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

    def isCompatibleFileFormat(self):
        """Compatible with MRC and em"""
        compatibleExts = ['.em', '.mrc']
        return True if getExt(self.inputTomos.getFirstItem().getFileName()) not in compatibleExts else False

    def getInTomosDir(self, *pathList):
        return self._getExtraPath(self.inTomosDir, *pathList)

    def getConvertedOutFileName(self, inFileName):
        outBaseName = removeBaseExt(inFileName)
        return self.getInTomosDir(outBaseName + '.mrc')

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
