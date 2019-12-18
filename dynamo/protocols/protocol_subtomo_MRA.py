# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
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

from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam
from tomo.protocols.protocol_base import ProtTomoSubtomogramAveraging
"""
Protocols to run Dynamo methods
"""

class DynamoSubTomoMRA(ProtTomoSubtomogramAveraging):
    """ It will align subtomograms using Dynamo"""
    _label = 'MRA alignment'

    def __init__(self, **args):
        ProtTomoSubtomogramAveraging.__init__(self, **args)

        # --------------------------- DEFINE param functions ------------------------

    def _defineParams(self, form):
        form.addSection(label='Input subtomograms')
        form.addParam('inputVolumes', PointerParam,
                      pointerClass="SetOfSubTomograms",
                      label='Set of volumes',
                      help="Set of subtomograms to align with dynamo")

        form.addSection(label='Templates')
        form.addParam('generateTemplate', BooleanParam, default=True,
                      label='Generate a template reference:', help="Generate a reference template based on parameters")

        # Pointer to a reference
        form.addParam('templateRef', PointerParam, label="Template",
                      condition="not generateTemplate", pointerClass='Volume',
                      help='Template to be used')

        # Generation template options
        form.addParam('useChosenSoP', BooleanParam, label='Use random chosen set of particles', default=False,
                      condition="generateTemplate", help="Use a random set of particles")

        form.addParam('numberOfParticles', IntParam, label='Number of references', default=50,
                      condition="generateTemplate and useChosenSoP", help="Number of references to generate automatically")

        form.addParam('useRandomTable', BooleanParam, label='Use a random table (rotate particles randomly)', default=False,
                      condition="generateTemplate", help="")

        form.addParam('compensateMissingWedge', BooleanParam, label='Compensate for missing wedge', default=False,
                      condition="generateTemplate", help="")

        # Masks
        form.addSection(label='Masks')

        form.addParam('alignmentMask', PointerParam, label="Alignment mask",
                      pointerClass='VolumeMask',
                      help='Mask for the alignment')

        form.addParam('classificationMask', PointerParam, label="Classification mask",
                      pointerClass='VolumeMask', allowsNull=True,
                      help='Mask for the classification steps')

        form.addParam('fourierMask', PointerParam, label="Fourier mask on template",
                      pointerClass='VolumeMask', allowsNull=True,
                      help='A binary mask describing which fourier components of'
                           ' the initial average are known.')

        form.addParam('fscMask', PointerParam, label="FSC mask",
                      pointerClass='VolumeMask', allowsNull=True,
                      help='a direct space mask that will be imposed onto any couple of volumes when computing their FSC.')


        # form.addParallelSection(threads=0, mpi=8)

        # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('alignStep')
        self._insertFunctionStep('createOutput')

        # --------------------------- STEPS functions -------------------------------

    def convertInputStep(self):
        pass
        # TODO: Add convertion of image files to a compatible with Dynamo
        # fnDir = self._getExtraPath("inputVolumes")
        # makePath(fnDir)
        # fnRoot = os.path.join(fnDir, "subtomo")
        # writeSetOfVolumes(self.inputVolumes.get(), fnRoot)
        # self.fnSel = self._getExtraPath("subtomograms.sel")
        # self.runJob("xmipp_selfile_create",
        #             '"%s*.vol">%s' % (fnRoot, self.fnSel), numberOfMpi=1)

    def alignStep(self):
        pass
        # fnIn = self._getExtraPath("subtomograms.sel")
        # fnOut = self._getExtraPath("mltomo")
        # self._createFilesForMLTomo()
        # fhDoc = self._getExtraPath("subtomograms.doc")
        # args = ' -i ' + fnIn + \
        #        ' -o ' + fnOut + \
        #        ' -nref ' + str(self.numberOfReferences.get()) + \
        #        ' -doc ' + fhDoc + \
        #        ' -iter ' + str(self.numberOfIters.get()) + \
        #        ' -ang ' + str(self.angularSampling.get()) + \
        #        ' -dim ' + str(self.downscDim.get()) + \
        #        ' ' + self.extraParams.get()
        # fhWedge = self._getExtraPath("wedge.doc")
        # if exists(fhWedge):
        #     args = args + ' -missing ' + fhWedge
        # self.runJob("xmipp_ml_tomo", args, numberOfMpi=self.numberOfMpi.get())

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
