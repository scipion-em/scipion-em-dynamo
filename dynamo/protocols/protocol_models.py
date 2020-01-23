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

import os

from pyworkflow.em import ProtAnalysis3D
from pyworkflow.protocol.params import PointerParam
from pyworkflow.utils import importFromPlugin
import pyworkflow.utils as pwutils

from dynamo.viewers.views_tkinter_tree import DynamoDialog

from tomo.viewers.views_tkinter_tree import TomogramsTreeProvider

Mesh = importFromPlugin("tomo.objects", "Mesh")
SetOfMeshes = importFromPlugin("tomo.objects", "SetOfMeshes")
ProtTomoBase = importFromPlugin("tomo.protocols", "ProtTomoBase")


"""
Protocols to create models in Dynamo
"""

class DynamoModels(ProtAnalysis3D, ProtTomoBase):
    """ It will align subtomograms using Dynamo"""
    _label = 'model manager'

    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)

        # --------------------------- DEFINE param functions ------------------------

    def _defineParams(self, form):
        form.addSection(label='Input Parameters')
        form.addParam('inputTomograms', PointerParam,
                      pointerClass="SetOfTomograms",
                      label='Set of tomograms',
                      help="Set of tomograms to create a model")


        # form.addParallelSection(threads=0, mpi=8)

        # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('launchDynamoGUIStep', interactive=True)

        # --------------------------- STEPS functions -------------------------------
    def launchDynamoGUIStep(self):
        tomoList = [tomo.clone() for tomo in self.inputTomograms.get().iterItems()]

        tomoProvider = TomogramsTreeProvider(tomoList, self._getExtraPath(), "txt")

        self.dlg = DynamoDialog(None, self._getExtraPath(), provider=tomoProvider,)

        self._createOutput()

    def _createOutput(self):
        outSet = self._createSetOfMeshes()
        for tomo in self.inputTomograms.get().iterItems():
                tomoName = pwutils.removeBaseExt(tomo.getFileName())
                outFile = self._getExtraPath(tomoName + '.txt')
                roi = Mesh(outFile)
                roi.setVolume(tomo)
        outSet.append(roi)
        outSet.setVolumes(self.inputTomograms.get())
        self._defineOutputs(outputMeshes=outSet)
        self._defineSourceRelation(self.inputTomograms.get(), outSet)


        # --------------------------- INFO functions --------------------------------

    def _summary(self):
        pass

    def _methods(self):
        pass

    def _citations(self):
        pass
