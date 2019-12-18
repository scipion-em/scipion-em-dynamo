# -*- coding: utf-8 -*-
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
from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam
"""
Protocols to run Dynamo methods
"""

class DynamoModels(ProtAnalysis3D):
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
        self._insertFunctionStep('writeInputFile')
        self._insertFunctionStep('modelStep')
        self._insertFunctionStep('createOutput')

        # --------------------------- STEPS functions -------------------------------

    def writeInputFile(self):
        self.inputFilePath = os.path.join(os.environ.get("SCIPION_HOME"), "software", "tmp", "commands.doc")
        tomoFile = os.path.join(os.environ.get("SCIPION_HOME"), "software", "tmp", "tomos.vll")
        catalogue = self._getExtraPath('tomos')
        tomoFid = open(tomoFile, 'w')
        for tomo in self.inputTomograms.get():
            tomoFile = tomo.getFileName()
            tomoName = os.path.basename(tomoFile)
            tomoFid.write(tomoName + '\n')
        tomoFid.close()
        inputFid = open(self.inputFilePath, 'w')
        content = 'dcm -create %s -fromvll %s \n' \
                  'dtmslice %s -c %s \n' \
                  'uiwait(msgbox(\'Click Ok when done\')) \n' \
                  'm = dread(dcmodels(\'%s\',\'i\',1)); \n' \
                  'writematrix(m.points, \'%s\'); \n' \
                  'exit' % (catalogue, tomoFile, tomoFile, catalogue, catalogue, self._getExtraPath('extra.txt'))
        inputFid.write(content)
        inputFid.close()


    def modelStep(self):
        program = '/home/davidh/dynamo_v1.146/matlab/bin/dynamo'
        args = ' %s' % self.inputFilePath
        dynamo = self.runJob(program, args)

    def createOutput(self):
        pass
        # self.subtomoSet = self._createSetOfSubTomograms()
        # inputSet = self.inputVolumes.get()
        # self.subtomoSet.copyInfo(inputSet)
        # self.fnDoc = '%s/mltomo_it00000%d.doc' % (
        # self._getExtraPath(), self.numberOfIters)
        # self.docFile = open(self.fnDoc)
        # self.subtomoSet.copyItems(inputSet, updateItemCallback=self._updateItem)
        # self.docFile.close()
        # classesSubtomoSet = self._createSetOfClassesSubTomograms(
        #     self.subtomoSet)
        # classesSubtomoSet.classifyItems(updateClassCallback=self._updateClass)
        # self._defineOutputs(outputSubtomograms=self.subtomoSet)
        # self._defineSourceRelation(self.inputVolumes, self.subtomoSet)
        # self._defineOutputs(outputClassesSubtomo=classesSubtomoSet)
        # self._defineSourceRelation(self.inputVolumes, classesSubtomoSet)
        # self._cleanFiles()

        # --------------------------- INFO functions --------------------------------

    def _summary(self):
        pass
        # summary = []
        # if hasattr(self, 'outputClassesSubtomo'):
        #     summary.append(
        #         "Input subtomograms: *%d* \nRequested classes: *%d*\nGenerated classes: *%d* in *%d* iterations\n"
        #         % (self.inputVolumes.get().getSize(), self.numberOfReferences,
        #            self.outputClassesSubtomo.getSize(), self.numberOfIters))
        # else:
        #     summary.append("Output classes not ready yet.")
        # return summary

    def _methods(self):
        pass
        # methods = []
        # if hasattr(self, 'outputClassesSubtomo'):
        #     methods.append(
        #         'We classified %d subtomograms from %s into %d classes %s using *MLTomo*.'
        #         % (self.inputVolumes.get().getSize(),
        #            self.getObjectTag('inputVolumes'),
        #            self.outputClassesSubtomo.getSize(),
        #            self.getObjectTag('outputClassesSubtomo')))
        # else:
        #     methods.append("Output classes not ready yet.")
        # return methods

    def _citations(self):
        pass
        # return ['Scheres2009c']
