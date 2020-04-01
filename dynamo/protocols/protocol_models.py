# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
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

import os, glob
import numpy as np

from pwem.protocols import ProtAnalysis3D

from pyworkflow.protocol.params import PointerParam
import pyworkflow.utils as pwutils
from pyworkflow.utils.properties import Message

from dynamo.viewers.views_tkinter_tree import DynamoDialog

from tomo.viewers.views_tkinter_tree import TomogramsTreeProvider
from tomo.objects import Mesh
from tomo.protocols import ProtTomoBase


"""
Protocols to create models in Dynamo
"""

class DynamoModels(ProtAnalysis3D, ProtTomoBase):
    """Model manger from Dynamo for Mesh creation.
    Once the desired Tomogram is opened, a surface model will be created and loaded in Dynamo automatically. Once the points of
    the Mesh have been defined, close the window to save automatically your data."""
    _label = 'model manager'

    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)
        self.OUTPUT_PREFIX = 'outputMeshes'

        # --------------------------- DEFINE param functions ------------------------
    def _defineParams(self, form):
        form.addSection(label='Input Parameters')
        form.addParam('inputTomograms', PointerParam,
                      pointerClass="SetOfTomograms",
                      label='Set of tomograms',
                      help="Set of tomograms to create a model")

        # --------------------------- INSERT steps functions --------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('launchDynamoGUIStep', interactive=True)

        # --------------------------- STEPS functions -------------------------------
    def launchDynamoGUIStep(self):
        tomoList = [tomo.clone() for tomo in self.inputTomograms.get().iterItems()]

        tomoProvider = TomogramsTreeProvider(tomoList, self._getExtraPath(), "txt")

        self.dlg = DynamoDialog(None, self._getExtraPath(), False, provider=tomoProvider,)

        if glob.glob(self._getExtraPath("*.txt")):
            self._createOutput()

    def _createOutput(self):
        outSet = self._createSetOfMeshes()
        for tomo in self.inputTomograms.get().iterItems():
            tomoName = pwutils.removeBaseExt(tomo.getFileName())
            outFile = self._getExtraPath(tomoName + '.txt')
            if os.path.isfile(outFile):
                data = np.loadtxt(outFile, delimiter=',')
                groups = np.unique(data[:, 3]).astype(int)
                for group in groups:
                    mesh = Mesh(group=group, path=outFile)
                    mesh.setVolume(tomo.clone())
                    outSet.append(mesh)
        outSet.setVolumes(self.inputTomograms.get())
        name = self.OUTPUT_PREFIX
        args = {name: outSet}
        self._defineOutputs(**args)
        self._defineSourceRelation(self.inputTomograms.get(), outSet)

        # Update Outputs
        outSet.setObjComment(self.getSummary())
        self._updateOutputSet(name, outSet, state=outSet.STREAM_CLOSED)

        # --------------------------- INFO functions --------------------------------
    def _summary(self):
        if self.getOutputsSize() >= 1:
            summary = []
            for key, output in self.iterOutputAttributes():
                summary.append("*%s:* \n %s \n" % (key, output.getObjComment()))
        else:
            summary = []
            summary.append(Message.TEXT_NO_OUTPUT_CO)

        return summary

    def getSummary(self):
        summary = []
        count = 0
        for file in os.listdir(self._getExtraPath()):
            if file.endswith(".txt"):
                count += 1
        summary.append("Meshes defined for %d/%d files have been saved in Scipion (%s)" % (
            count/2, self.inputTomograms.get().getSize(), self._getExtraPath()))

        return "\n".join(summary)

    def _methods(self):
        tomos = self.inputTomograms.get()
        return [
            "Model creation using Dynamo",
            "A total of %d tomograms of dimensions %s were used"
            % (tomos.getSize(), tomos.getDimensions()),
        ]

    def _citations(self):
        return ['CASTANODIEZ2012139']
