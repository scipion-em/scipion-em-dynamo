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

from os.path import join
from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam, StringParam, FloatParam, LEVEL_ADVANCED
from pyworkflow.utils import importFromPlugin
from pyworkflow.utils.path import makePath
from dynamo import Plugin
from dynamo.convert import writeVolume, writeSetOfVolumes, writeDynTable, readDynTable
ProtTomoSubtomogramAveraging = importFromPlugin("tomo.protocols.protocol_base", "ProtTomoSubtomogramAveraging")
AverageSubTomogram = importFromPlugin("tomo.objects", "AverageSubTomogram")
SetOfAverageSubTomograms = importFromPlugin("tomo.objects", "SetOfAverageSubTomograms")


class DynamoSubTomoMRA(ProtTomoSubtomogramAveraging):
    """ This protocol will align subtomograms using Dynamo MRA Subtomogram Averaging"""

    _label = 'Subtomogram alignment'

    def __init__(self, **args):
        ProtTomoSubtomogramAveraging.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------

    def _defineParams(self, form):
        form.addSection(label='Input subtomograms')
        form.addParam('inputVolumes', PointerParam, pointerClass="SetOfSubTomograms", label='Set of volumes',
                      help="Set of subtomograms to align with dynamo")
        form.addParam('sym', StringParam, default='c1', label='Symmetry group',
                      help="Specify the article's symmetry. Symmetrization is applied at the beginning of the round to "
                           "the input reference.")
        form.addParam('numberOfRounds', IntParam, label='Rounds', default=1, expertLevel=LEVEL_ADVANCED,
                      help="Number of rounds (from 1 to 8). Each round consists in X iterations with the same "
                           "parameters, but parameters could vary in different rounds")
        form.addParam('numberOfIters', IntParam, label='Iterations', default=5, help="Number of iterations per round")
        form.addParam('pca', BooleanParam, label='Perform PCA', default=False,
                      help="If selected, principal component analysis alignment is performed")

        form.addSection(label='Templates')
        form.addParam('mra', BooleanParam, label='Perform MRA', default=False,
                      help="If selected, multi-reference alignment (MRA) is performed")
        form.addParam('nref', IntParam, label='Number of references for MRA', default=1, condition="mra",
                      help="Number of references for multi-reference alignment(MRA)")
        form.addParam('generateTemplate', BooleanParam, default=False, label='Generate reference templates:',
                      help="Generate a reference template based on parameters")
        form.addParam('templateRef', PointerParam, label="Template", condition="not generateTemplate",
                      pointerClass='Volume', allowsNull=True,
                      help='The size of the template should be equal or smaller than the size of the particles. If you '
                           'pass a single file in multireference modus (MRA), Dynamo will just made copies of it.')
        # form.addParam('templateSetRef', PointerParam, label="Templates", condition="not generateTemplate and mra",
        #               pointerClass='SetOfVolumes', allowsNull=True,
        #               help='The size of the template should be equal or smaller than the size of the particles. If you '
        #                    'pass a single file in multireference modus (MRA), Dynamo will just made copies of it.')
        form.addParam('useChosenSoP', BooleanParam, label='Use random chosen set of particles', default=False,
                      condition="generateTemplate", help="Use a random set of particles")
        form.addParam('numberOfParticles', IntParam, label='Number of references', default=50,
                      condition="generateTemplate and useChosenSoP", help="Number of references to generate "
                                                                          "automatically")
        form.addParam('useRandomTable', BooleanParam, label='Use a random table (rotate particles randomly)',
                      default=False, condition="generateTemplate", help="")
        form.addParam('compensateMissingWedge', BooleanParam, label='Compensate for missing wedge', default=False,
                      condition="generateTemplate", help="")

        form.addSection(label='Masks')
        form.addParam('mask', PointerParam, label="Alignment mask", pointerClass='VolumeMask', allowsNull=True,
                      help='Needs to have the same dimmensionality as the template. It does NOT need to be smoothened. '
                           'A hard mask (only 0 and 1 values) is expected for the Roseman"s-like normalization')
        form.addParam('cmask', PointerParam, label="Classification mask", pointerClass='VolumeMask', allowsNull=True,
                      help='Needs to have the same dimmensionality as the template. It does NOT need to be smoothened. '
                           'A hard mask (only 0 and 1 values) is expected for the Roseman"s-like normalization')
        form.addParam('fmask', PointerParam, label="Fourier mask on template", pointerClass='VolumeMask',
                      allowsNull=True,
                      help='The fmask indicates which fourier coefficients are present at the starting reference volume'
                           '. The file should contain only ones or zeros.')
        # form.addParam('setfmask', PointerParam, label="Fourier mask on template", pointerClass='SetOfVolumes',
        #               allowsNull=True, condition="mra",
        #               help='The fmask indicates which fourier coefficients are present at the starting reference volume'
        #                    '. The file should contain only ones or zeros.')
        form.addParam('smask', PointerParam, label="FSC mask", pointerClass='VolumeMask', allowsNull=True,
                      help='A direct space mask that will be imposed onto any couple of volumes when computing their '
                           'FSC (smask)')

        form.addSection(label='Angular scanning')
        form.addParam('cr', IntParam, label='Cone range', default=360,
                      help="The first two Euler angles are used to define the orientation of the vertical axis of the "
                           "protein. First Euler angle (tdrot) rotates the template around its z axis. Second Euler "
                           "angle (tilt) rotates the template around its x axis.Dynamo scans for this axis inside a "
                           "cone: The 'cone_range' parameter defines the angular aperture of this cone.360 degrees is "
                           "thus the value for a global scan. To skip the part of the angular search that looks for "
                           "orientations, you have to set 'cone range' to 0, and 'cone_sampling' to 1")
        form.addParam('cs', IntParam, label='Cone sampling', default=60,
                      help="This parameter expresses the discretization inside this cone. The sampling is given in "
                           "degrees, and corresponds to a representative angular distance between two neighboring "
                           "orientations inside the cone.")
        form.addParam('cf', IntParam, label='Cone flip', default=0, expertLevel=LEVEL_ADVANCED,
                      help="Generates a mirrored scanning geometry: the 'cone' of directions is complemented with the "
                           "diametrally oposed cone. This is useful when averaging elongated particles in the case in "
                           "which the direction of each one is not certain, i.e., the initial table catches the overall"
                           " orientation of each particle, but it is not certain on which end is the 'head' and which "
                           "is the 'tail', so that the refinement should allow 'flippling' the particles (but still "
                           "produce a scanning set of angles concentrated along the axis of the particle). 0:  No "
                           "inversion of the cone (default!) 1:  The cone is inverted only for the coarsest level of "
                           "the multigrid refinement 2:  The cone is inverted for all refinement levels")
        form.addParam('ccp', IntParam, label='Cone check peak', default=0, expertLevel=LEVEL_ADVANCED,
                      help="Controls peak quality inside the scanned directions. Useful when the particles are in close"
                           " contact with some feature of high intensity. Ensures that only angles yielding a real peak"
                           " will be considered. If the angle that yields the maximum of CC is in the boundary, the "
                           "particle will get aligned with its original alignment parameters. Values: 0: normal "
                           "behaviour (no peak quality control, default) any integer  : degrees that define the peak. A"
                           " maximum occurring within this distance to the boundary of the defined cone will be "
                           "discarded.")
        form.addParam('inplane_range', IntParam, label='Inplane rotation range', default=360,
                      help='The third Euler angle ("narot") defines rotations about the new axis of the reference. 360 '
                           'degrees is the value for a global scan. To skip the part of the angular search that looks '
                           'for azimuthal rotations , you have to set 1)  "inplane range" to zero, and 2)  '
                           '"inplane_sampling" to 1. Likewise,  to skip the part of the angular search that looks for '
                           'orientations, you have to manipulate the "cone" parameters: 1)  "cone range" to zero, and 2'
                           ') "cone_sampling" to 1.')
        form.addParam('inplane_sampling', IntParam, label='Inplane rotation sampling', default=45,
                      help='Parameter "inplane_sampling" is associated with the "narot" angle (New Axis ROTation). 1) '
                           'the axis of the template is rotated to a new orientation (defined by "tdrot" and "tilt"). 2'
                           ') the rotated template rotates again on its new axis (an "inplane" rotation). "inplane_'
                           'sampling" defines the angular interval (in degrees) between two of these inplane rotations')
        form.addParam('inplane_flip', IntParam, label='Inplane flip', default=0, expertLevel=LEVEL_ADVANCED,
                      help='Flips the set of inplane rotations. The set of inplane rotations to scan will be the '
                           'original set plus the flipped orientations.This is useful when the particles have a '
                           'directionality, but it is not very well defined. Values 0  :  no flip (default) 1  :  flips'
                           ' the coarsest level  in the multilevel grid   2  :  flips the full set (all levels).')
        form.addParam('inplane_check_peak', IntParam, label='Inplane check peak', default=0,
                      expertLevel=LEVEL_ADVANCED, help="Controls peak quality along the inplane rotation. Useful when "
                                                       "the particles are in close contact with some feature of high "
                                                       "intensity. Ensures that only angles yielding a real peak will "
                                                       "be considered. If the angle that yields the maximum of CC is in"
                                                       " the boundary, the particle will get aligned with its original "
                                                       "alignment parameters. Values: 0  : normal behaviour (no peak "
                                                       "quality control, default) any integer  : degrees that define "
                                                       "the peak. A maximum occurring within this distance to the "
                                                       "boundary of the defined range for inplane rotations will be "
                                                       "discarded.")
        form.addParam('rf', IntParam, label='Refine', default=5, expertLevel=LEVEL_ADVANCED,
                      help="How many refinement iterations are carried out on each single particle. This refinement "
                           "when comparing rotations of the reference against the data, takes the best orientation and "
                           "looks again with a finer sampling. The sampling in the refined search will be half of the "
                           "sampling used in the original one.  The range of the refined search encompasses all the "
                           "orientations that neighobur the best orientation found in the original search.")
        form.addParam('rff', IntParam, label='Refine factor', default=2, expertLevel=LEVEL_ADVANCED,
                      help="Controls the size of the angular neighborhood during the local refinement of the angular "
                           "grid.")

        form.addSection(label='Threshold')
        form.addParam('separation', IntParam, label='Separation in tomogram', default=0,
                      help='When tuned to  positive number, it will check the relative positions (positions in the '
                           'tomogram+shifts) of all the particles in each tomogram separately. Whenever two particles '
                           'are closer together than "separation_in_tomogram", only the particle with the higher '
                           'correlation will stay.')
        form.addParam('threshold', FloatParam, label='Threshold', default=0.2,
                      help='Different thresholding policies can be used in order to select which particles are averaged'
                           ' in view of their CC (cross correlation value) . The value of the thresholding parameter '
                           'defined here  will be interpreted differently depending on the "threshold_modus"')
        form.addParam('thresholdMode', IntParam, label='Threshold modus', default=1,
                      help='Possible values of the thresholding policy "threshold_modus": 0: no thesholding policy 1: '
                           'THRESHOLD is an absolute threshold (only particles with CC above this value are selected).'
                           ' 2: efective threshold = mean(CC) * THRESHOLD. 3: efective threshold = mean(CC) +std(CC)*'
                           'THRESHOLD. 4: THRESHOLD is the total number of particle (ordered by CC ). 5: THRESHOLD '
                           'ranges between 0 and 1  and sets the fraction of particles. * 11,21,31,34,41,51: select the'
                           ' same particles as 1,2,3,4 or 5, BUT non selected particles will be excluded: - from '
                           'averaging in the present iteration,and - ALSO from alignment in the next iteration (unlike '
                           '1,2,3,4,5).')
        form.addParam('threshold2', FloatParam, label='Second threshold', default=0.2, expertLevel=LEVEL_ADVANCED,
                      help="Thresholding II is operated against the average produced by the particles that survived the"
                           " first thresholding.")
        form.addParam('thresholdMode2', IntParam, label='Second threshold modus', default=0, expertLevel=LEVEL_ADVANCED,
                      help=" ")
        form.addParam('ccmatrix', BooleanParam, label='Compute  cross-correlation matrix', default=False,
                      help="Computation of a Cross-Correlation matrix among the aligned particles.")
        form.addParam('ccmatrixType', StringParam, label='Cross-correlation matrix type', default='align',
                      condition="ccmatrix", expertLevel=LEVEL_ADVANCED,
                      help="string with three characters, each position controling a different aspect: thresholding, "
                           "symmetrization, compensation")
        form.addParam('ccmatrixBatch', IntParam, label='Cross-correlation matrix batch', default=128,
                      expertLevel=LEVEL_ADVANCED,
                      help="Number of particles to be kept in memory simultaneously during the computation of the "
                           "ccmatrix. The larger this number, the more efficient the algorithm performance, as more "
                           "computations can be kept for reuse.However, trying to keep all the particles in memory can "
                           "lead to saturate it,blocking the CPU. Additionally, a small batch allows to divide the "
                           "matrix in more blocks. This might be useful in parallel computations.")
        form.addParam('low', IntParam, label='Low frequency', default=32,
                      expertLevel=LEVEL_ADVANCED, help='Cut off frequency for low pass filtering')
        form.addParam('high', IntParam, label='High frequency', default=2,
                      expertLevel=LEVEL_ADVANCED, help='Cut off frequency for high pass filtering')
        form.addParam('lim', StringParam, label='Area search', default='4 4 4', expertLevel=LEVEL_ADVANCED,
                      help='Restricts the search area to an ellipsoid centered and oriented in the last found position.'
                           ' The three parameters are the semiaxes of the ellipsoid. If a single parameter is '
                           'introduced, the ellipsoid collapses into a sphere. If no restriction should be imposed, put'
                           ' a zero on the "area search modus" parameter')
        form.addParam('limm', IntParam, label='Area search modus', default=0, expertLevel=LEVEL_ADVANCED,
                      help='States how exactly the shifts (area search) will be interpreted 0:  no limitations (can '
                           'easily produce artifacts if the initial reference is bad) 1:  limits are understood from '
                           'the center of the particle cube. 2:  limits are understood from the previous estimation on '
                           'the particle position (i.e., the shifts available in the table) With this option, the '
                           'origin of the shifts changes at every iteration. 3:  limits are understood from the '
                           'estimation provided for the first iteration of the round. The origin of the shifts will '
                           'change at each round. 4:  limits are understood from the estimation provided for the first '
                           'iteration')

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('alignStep')
        self._insertFunctionStep('createOutput')

    # --------------------------- STEPS functions -------------------------------

    def convertInputStep(self):
        fnDirData = self._getExtraPath("data")
        makePath(fnDirData)
        fnRoot = join(fnDirData, "particle_")
        inputVols = self.inputVolumes.get()
        writeSetOfVolumes(inputVols, fnRoot)
        pcaInt = int(self.pca.get())
        mraInt = int(self.mra.get())
        ccmatrixInt = int(self.ccmatrix.get())
        dim, _, _ = self.inputVolumes.get().getDimensions()
        self.projName = 'dynamoAlignmentProject'
        fnTable = self._getExtraPath("initial.tbl")
        fhTable = open(fnTable, 'w')
        writeDynTable(fhTable, inputVols)
        fhTable.close()
        fhCommands = open(self._getExtraPath("commands.doc"), 'w')
        content = "dcp.new('%s','data','data','table', 'initial.tbl','gui',0);" % self.projName + \
                  "dvput('%s', 'dim', '%s');" % (self.projName, dim) + \
                  "dvput('%s', 'sym', '%s');" % (self.projName, self.sym) + \
                  "dvput('%s', 'ite', '%s');" % (self.projName, self.numberOfIters) + \
                  "dvput('%s', 'mra', %s);" % (self.projName, mraInt) + \
                  "dvput('%s', 'pcas', %d);" % (self.projName, pcaInt) + \
                  "dvput('%s', 'cr', '%s');" % (self.projName, self.cr) + \
                  "dvput('%s', 'cs', '%s');" % (self.projName, self.cs) + \
                  "dvput('%s', 'cf', '%s');" % (self.projName, self.cf) + \
                  "dvput('%s', 'ccp', '%s');" % (self.projName, self.ccp) + \
                  "dvput('%s', 'rf', '%s');" % (self.projName, self.rf) + \
                  "dvput('%s', 'rff', '%s');" % (self.projName, self.rff) + \
                  "dvput('%s', 'ir', '%s');" % (self.projName, self.inplane_range) + \
                  "dvput('%s', 'is', '%s');" % (self.projName, self.inplane_sampling) + \
                  "dvput('%s', 'if', '%s');" % (self.projName, self.inplane_flip) + \
                  "dvput('%s', 'icp', '%s');" % (self.projName, self.inplane_check_peak) + \
                  "dvput('%s', 'thr', %s);" % (self.projName, self.threshold) + \
                  "dvput('%s', 'thrm', %s);" % (self.projName, self.thresholdMode) + \
                  "dvput('%s', 'thr2', %s);" % (self.projName, self.threshold2) + \
                  "dvput('%s', 'thr2m', %s);" % (self.projName, self.thresholdMode2) + \
                  "dvput('%s', 'ccms', %s);" % (self.projName, ccmatrixInt) + \
                  "dvput('%s', 'ccmt', '%s');" % (self.projName, self.ccmatrixType) + \
                  "dvput('%s', 'batch', '%s');" % (self.projName, self.ccmatrixBatch) +  \
                  "dvput('%s', 'stm', '%s');" % (self.projName, self.separation) + \
                  "dvput('%s', 'low', '%s');" % (self.projName, self.low) + \
                  "dvput('%s', 'high', '%s');" % (self.projName, self.high) + \
                  "dvput('%s', 'lim', '%s');" % (self.projName, self.lim) + \
                  "dvput('%s', 'limm', '%s');" % (self.projName, self.limm)

        # if self.mra.get():
            # fnDirMRA = self._getExtraPath("folder_multireference")
            # makePath(fnDirMRA)
            # for i in range(self.nref.get()):
            #     fnRootT = join(self._getExtraPath(), "table_initial_ref_00%s" % str(int(i)+1))
            #     fhTable = open(fnRootT, 'w')
            #     writeDynTable(fhTable, inputVols)
            #     fhTable.close()
            # content += "dvput('%s', 'nref', '%s');" % (projName, self.nref)

            # if self.templateSetRef.get() is not None and not self.generateTemplate.get():
            #     # fnDirTemps = self._getExtraPath("templates")
            #     makePath(self._getExtraPath())
            #     fnRoot = join(self._getExtraPath(), "template_initial_ref_")
            #     writeSetOfVolumes(self.templateSetRef.get(), fnRoot)
            #     content += "dvput('%s', 'template', 'templates');" % projName
            # else:
            #     pass
            #     # Generate templates

            # if self.setfmask.get() is not None:
            #     # fnDirFmasks = self._getExtraPath("fmasks")
            #     makePath(self._getExtraPath())
            #     fnRootf = join(self._getExtraPath(), "fmask_initial_ref_")
            #     writeSetOfVolumes(self.setfmask.get(), fnRootf)
            #     content += "dvput('%s', 'fmask', 'fmasks');" % projName

        template = self.templateRef.get()
        if template is not None:
            writeVolume(template, join(self._getExtraPath(), 'template'))
            content += "dvput('%s', 'template', 'template.mrc');" % self.projName
        else:
            # pass
            makePath(self._getExtraPath('templates'))
            content += "dynamo_write_multireference('table','initial.tbl','1:4','refs','data','subset','4/4');"

        if self.mask.get() is not None:
            writeVolume(self.mask.get(), join(self._getExtraPath(), 'mask'))
            content += "dvput('%s', 'mask', 'mask.mrc');" % self.projName
        if self.cmask.get() is not None:
            writeVolume(self.cmask.get(), join(self._getExtraPath(), 'cmask'))
            content += "dvput('%s', 'cmask', 'cmask.mrc');" % self.projName
        if self.fmask.get() is not None:
            writeVolume(self.fmask.get(), join(self._getExtraPath(), 'fmask'))
            content += "dvput('%s', 'fmask', 'fmask.mrc');" % self.projName
        if self.smask.get() is not None:
            writeVolume(self.smask.get(), join(self._getExtraPath(), 'smask'))
            content += "dvput('%s', 'smask', 'smask.mrc');" % self.projName

        # content += "dvcheck('%s');" % projName
        content += "dvcheck('%s');dvunfold('%s');dynamo_execute_project %s" % (self.projName, self.projName, self.projName)
        fhCommands.write(content)
        fhCommands.close()

    def alignStep(self):
        Plugin.runDynamo(self, 'commands.doc', cwd=self._getExtraPath())

    def createOutput(self):
        # self.fhTable = open(self._getExtraPath("initial.tbl"), 'r')  # Change to "real.tbl" or iter_x/...
        # self.subtomoSet = self._createSetOfSubTomograms()
        # inputSet = self.inputVolumes.get()
        # self.subtomoSet.copyInfo(inputSet)
        # self.subtomoSet.copyItems(inputSet, updateItemCallback=self._updateItem)
        # classesSubtomoSet = self._createSetOfClassesSubTomograms(self.subtomoSet)
        # classesSubtomoSet.classifyItems(updateClassCallback=self._updateClass)
        # self.fhTable.close()
        # self._defineOutputs(outputSubtomograms=self.subtomoSet)
        # self._defineSourceRelation(self.inputVolumes, self.subtomoSet)
        # self._defineOutputs(outputClassesSubtomo=classesSubtomoSet)
        # self._defineSourceRelation(self.inputVolumes, classesSubtomoSet)

        self.fhTable = open(self._getExtraPath('%s/results/ite_000%s/averages/refined_table_ref_001_ite_000%s.tbl') %
                            (self.projName, self.numberOfIters.get(), self.numberOfIters.get()), 'r')
        inputSet = self.inputVolumes.get()
        averageSubTomogram = AverageSubTomogram()
        readDynTable(self, averageSubTomogram)
        averageSubTomogram.setFileName(self._getExtraPath
                                       ('%s/results/ite_000%s/averages/average_ref_001_ite_000%s.em') %
                                       (self.projName, self.numberOfIters.get(), self.numberOfIters.get()))
        averageSubTomogram.setSamplingRate(inputSet.getSamplingRate())
        setOfAverageSubTomograms = self._createSet(SetOfAverageSubTomograms, 'subtomograms%s.sqlite', "")
        setOfAverageSubTomograms.copyInfo(inputSet)
        setOfAverageSubTomograms.setSamplingRate(inputSet.getSamplingRate())
        setOfAverageSubTomograms.append(averageSubTomogram)

        self._defineOutputs(averageSubTomogram=setOfAverageSubTomograms)
        self._defineSourceRelation(inputSet, setOfAverageSubTomograms)

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
