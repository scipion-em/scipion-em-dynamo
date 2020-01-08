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

from os.path import join, basename
from os import environ
from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam, StringParam, FloatParam, LEVEL_ADVANCED
from pyworkflow.utils import importFromPlugin
from pyworkflow.utils.path import makePath
from dynamo import Plugin
from dynamo.convert import writeSetOfVolumes, writeDynTable, readDynTable
ProtTomoSubtomogramAveraging = importFromPlugin("tomo.protocols.protocol_base", "ProtTomoSubtomogramAveraging")


class DynamoSubTomoMRA(ProtTomoSubtomogramAveraging):
    """ This protocol will align subtomograms using Dynamo MRA Subtomogram Averaging"""

    _label = 'MRA alignment'

    def __init__(self, **args):
        ProtTomoSubtomogramAveraging.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------

    def _defineParams(self, form):
        form.addSection(label='Input subtomograms')
        form.addParam('projName', StringParam, label='Project Name',
                      help='Name for the dynamo align project that will be generated')
        form.addParam('inputVolumes', PointerParam, pointerClass="SetOfSubTomograms", label='Set of volumes',
                      help="Set of subtomograms to align with dynamo")
        form.addParam('numberOfRounds', IntParam, label='Rounds', default=1, expertLevel=LEVEL_ADVANCED,
                      help="Number of rounds (from 1 to 8). Each round consists in X iterations with the same "
                           "parameters, but parameters could vary in different rounds")
        form.addParam('numberOfIters', IntParam, label='Iterations', default=5, help="Number of iterations per round")
        form.addParam('numberOfRefs', IntParam, label='Number of references for MRA', default=1,
                      help="Number of references for multi-reference alignment(MRA)")
        form.addParam('pca', BooleanParam, label='Perform PCA', default=False, expertLevel=LEVEL_ADVANCED,
                      help="If selected, principal component analysis alignment is performed")

        form.addSection(label='Templates')
        form.addParam('generateTemplate', BooleanParam, default=False,
                      label='Generate reference templates:', help="Generate a reference template based on parameters")
        form.addParam('templateRef', PointerParam, label="Template",
                      condition="not generateTemplate", pointerClass='Volume',
                      help='Template to be used')
        form.addParam('useChosenSoP', BooleanParam, label='Use random chosen set of particles', default=False,
                      condition="generateTemplate", help="Use a random set of particles")
        form.addParam('numberOfParticles', IntParam, label='Number of references', default=50,
                      condition="generateTemplate and useChosenSoP", help="Number of references to generate automatically")
        form.addParam('useRandomTable', BooleanParam, label='Use a random table (rotate particles randomly)', default=False,
                      condition="generateTemplate", help="")
        form.addParam('compensateMissingWedge', BooleanParam, label='Compensate for missing wedge', default=False,
                      condition="generateTemplate", help="")

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
                      pointerClass='VolumeMask', allowsNull=True, help='A direct space mask that will be imposed onto '
                                                                       'any couple of volumes when computing their FSC.')

        form.addSection(label='Angular scanning')
        form.addParam('coneAperture', IntParam, label='Cone aperture', default=360, help=" ")
        form.addParam('coneSampling', IntParam, label='Cone sampling', default=60, help=" ")
        form.addParam('coneFlip', IntParam, label='Cone flip', default=0, expertLevel=LEVEL_ADVANCED, help=" ")
        form.addParam('freezeRef', IntParam, label='Freeze reference', default=0, expertLevel=LEVEL_ADVANCED, help=" ")
        form.addParam('azymuthRange', IntParam, label='Azimuth rotation range', default=0, help=" ")
        form.addParam('azymuthSampling', IntParam, label='Azimuth rotation sampling', default=1, help=" ")
        form.addParam('azymuthFlip', IntParam, label='Azymuth flip', default=0, expertLevel=LEVEL_ADVANCED, help=" ")
        form.addParam('azymuthFreeze', IntParam, label='Freeze azymuth reference', default=0, expertLevel=LEVEL_ADVANCED, help=" ")
        form.addParam('refine', IntParam, label='Refine', default=5, help=" ")
        form.addParam('refineFactor', IntParam, label='Refine factor', default=2, help=" ")

        form.addSection(label='Shift/Threshold')
        form.addParam('shiftLims', StringParam, label='Shift limits', default='1 1 1', help=" ")
        form.addParam('shiftWay', IntParam, label='Shift limiting way', default=0, help=" ")
        form.addParam('separation', IntParam, label='Separation in tomogram', default=0, help=" ")
        form.addParam('threshold', FloatParam, label='Threshold', default=0.2, help=" ")
        form.addParam('thresholdMode', IntParam, label='Threshold modus', default=0, help=" ")
        form.addParam('threshold2', FloatParam, label='Second threshold', default=0.2, expertLevel=LEVEL_ADVANCED, help=" ")
        form.addParam('thresholdMode2', IntParam, label='Second threshold modus', default=0, expertLevel=LEVEL_ADVANCED, help=" ")

        form.addSection(label='Classification')
        form.addParam('ccmatrix', BooleanParam, label='Compute  cross-correlation matrix', default=False, help=" ")
        form.addParam('ccmatrixType', StringParam, label='Cross-correlation matrix type', default='bin 1; sym c1;', help=" ")
        form.addParam('ccmatrixBatch', IntParam, label='Cross-correlation matrix batch', default=10, expertLevel=LEVEL_ADVANCED, help=" ")

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('alignStep')
        self._insertFunctionStep('createOutput')

    # --------------------------- STEPS functions -------------------------------

    def convertInputStep(self):
        inputVols = self.inputVolumes.get()
        fnDir = self._getExtraPath("data")
        fnTable = self._getExtraPath("initial.tbl")
        makePath(fnDir)
        fnRoot = join(fnDir, "particle_")
        writeSetOfVolumes(inputVols, fnRoot)
        fhTable = open(fnTable, 'w')
        writeDynTable(fhTable, inputVols)
        fhTable.close()
        pcaInt = int(self.pca.get())
        fhCommands = open(self._getExtraPath("commands.doc"), 'w')
        content = "dcp.new('%s','table','initial.tbl','data','data','masks','default','gui',0);" % self.projName + \
                  "dvput('%s', 'ite_r1', '%s');" % (self.projName, self.numberOfIters) + \
                  "dvput('%s', 'mra', %s);" % (self.projName, self.numberOfRefs) + \
                  "dvput('%s', 'pcas', %d);" % (self.projName, pcaInt) + \
                  "dvput('%s', 'template', '%s');" % (self.projName, self.templateRef) + \
                  "dvcheck('%s');" % self.projName
                  # "dvput('%s', 'template', '%s');" % (self.projName, self.templateRef) + \
                  # "dvput('%s', 'nref_r1', '%s');" % (self.projName, self.numberOfRefs) + \
                  # "dvput('%s', 'mask', 1);" % self.projName + \
                  # "dvput('%s', 'cr', '%s');" % (self.projName, self.coneAperture) + \
                  # "dvput('%s', 'cs', '%s');" % (self.projName, self.coneSampling) + \
                  # "dvput('%s', 'inplane_range', 0);" % self.projName + \
                  # "dvput('%s', 'inplane_sampling', 1);" % vprojName + \
                  # "dvput('%s', 'refine', '%s');" % (self.projName, self.refine) + \
                  # "dvput('%s', 'low', 10);" % self.projName + \
                  # "dvput('%s', 'sym', 'c57');" % self.projName + \
                  # "dvput('%s', 'dim', '%s');" % (self.projName, self.inputVolumes.get().getDimensions()) + \
                  # "dvput('%s', 'area_search', 10);" % self.projName + \
                  # "dvput('%s', 'area_search_modus', 1);" % self.projName + \
                  # "dvcheck('%s');dvunfold('%s');('%s')" % (self.projName, projName, projName)
        fhCommands.write(content)
        fhCommands.close()
        # # dvunfold fails because template and data dimensions and more things detected by dvcheck fails
        # => it should work when all the correct parameters are set

        # Change default template and mask:
        # 'mask', 'seed1/settings/editorMaskEllipsoid.em', 'cmask', 'my_mask.em', 'fmask', 'seed1_c57_Zfilt_eo/settings', 'smask', 'seed1/settings/smoothingMaskOnes.em'
        # 'mask', 'seed1/settings/editorMaskEllipsoid.em', 'cmask', 'default', 'fmask', 'default', 'smask', 'default'

        # # Alignment params (from form) + boxSize (see dvhelp in Dynamo) (all of them should go before dvcheck)
        # # System Parameters
        # dvput('%s', 'destination', 'matlab_gpu');
        # dvput('%s', 'cores', 1);
        # dvput('%s', 'matlab_workers_average', 6);

    def alignStep(self):
        Plugin.runDynamo(self, 'commands.doc', cwd=self._getExtraPath())

    def createOutput(self):
        self.fhTable = open(self._getExtraPath("initial.tbl"), 'r')  # Change to "real.tbl"
        self.subtomoSet = self._createSetOfSubTomograms()
        inputSet = self.inputVolumes.get()
        self.subtomoSet.copyInfo(inputSet)
        self.subtomoSet.copyItems(inputSet, updateItemCallback=self._updateItem)
        classesSubtomoSet = self._createSetOfClassesSubTomograms(self.subtomoSet)
        classesSubtomoSet.classifyItems(updateClassCallback=self._updateClass)
        self.fhTable.close()
        self._defineOutputs(outputSubtomograms=self.subtomoSet)
        self._defineSourceRelation(self.inputVolumes, self.subtomoSet)
        self._defineOutputs(outputClassesSubtomo=classesSubtomoSet)
        self._defineSourceRelation(self.inputVolumes, classesSubtomoSet)

    # --------------------------- INFO functions --------------------------------

    def _summary(self):
        summary = []
        if hasattr(self, 'outputClassesSubtomo'):
            summary.append(
                "Input subtomograms: *%d* \nGenerated classes: *%d*"
                % (self.inputVolumes.get().getSize(), self.outputClassesSubtomo.getSize()))
        else:
            summary.append("Output classes not ready yet.")
        return summary

    def _methods(self):
        methods = []
        if hasattr(self, 'outputClassesSubtomo'):
            methods.append(
                'We classified %d subtomograms from %s into %d classes %s using Dynamo *MRA* Subtomogram averaging.'
                % (self.inputVolumes.get().getSize(),
                   self.getObjectTag('inputVolumes'),
                   self.outputClassesSubtomo.getSize(),
                   self.getObjectTag('outputClassesSubtomo')))
        else:
            methods.append("Output classes not ready yet.")
        return methods

    def _citations(self):
        return ['CASTANODIEZ2012139']

    # --------------------------- UTILS functions ----------------------------------

    def _updateItem(self, item, row):
        readDynTable(self, item)

    def _updateClass(self, item):
        pass
        # update class info (setRepresentative), see mltomo
