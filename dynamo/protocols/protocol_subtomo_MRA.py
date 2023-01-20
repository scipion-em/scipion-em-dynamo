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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os.path
from os import rename, remove
from os.path import join
from shutil import copy
from pwem import Domain
from pwem.objects.data import Volume, VolumeMask
from pyworkflow import BETA
from pyworkflow.object import Set
from pyworkflow.protocol import GPU_LIST, USE_GPU
from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam, StringParam, FloatParam, LEVEL_ADVANCED
from pyworkflow.utils.path import makePath, cleanPattern
from dynamo import Plugin
from dynamo.convert import writeVolume, writeSetOfVolumes, writeDynTable, readDynTable

ProtTomoSubtomogramAveraging = Domain.importFromPlugin("tomo.protocols.protocol_base", "ProtTomoSubtomogramAveraging")
AverageSubTomogram = Domain.importFromPlugin("tomo.objects", "AverageSubTomogram")
SetOfSubTomograms = Domain.importFromPlugin("tomo.objects", "SetOfSubTomograms")
SubTomogram = Domain.importFromPlugin("tomo.objects", "SubTomogram")


class DynamoSubTomoMRA(ProtTomoSubtomogramAveraging):
    """ This protocol will align subtomograms using Dynamo MRA Subtomogram Averaging"""

    _label = 'Subtomogram alignment'
    _devStatus = BETA

    def __init__(self, **args):
        ProtTomoSubtomogramAveraging.__init__(self, **args)

    # --------------------------- DEFINE param functions ------------------------

    def _defineParams(self, form):
        form.addSection(label='Input subtomograms')
        form.addHidden(USE_GPU, BooleanParam, default=True,
                       label="Use GPU for execution",
                       help="This protocol has both CPU and GPU implementation.\
                                   Select the one you want to use.")
        form.addHidden(GPU_LIST, StringParam, default='0',
                       expertLevel=LEVEL_ADVANCED,
                       label="Choose GPU IDs",
                       help="Add a list of GPU devices that can be used")
        form.addParam('inputVolumes', PointerParam, pointerClass="SetOfSubTomograms", label='Set of subtomograms',
                      help="Set of subtomograms to align with dynamo")
        form.addParam('sym', StringParam, default='c1', label='Symmetry group',
                      help="Specify the article's symmetry. Symmetrization is applied at the beginning of the round to "
                           "the input reference.")
        form.addParam('numberOfIters', IntParam, label='Iterations', default=5, help="Number of iterations per round")
        form.addParam('pca', BooleanParam, label='Perform PCA', default=False,
                      help="If selected, principal component analysis alignment is performed")

        form.addSection(label='Templates')
        form.addParam('mra', BooleanParam, label='Perform MRA', default=False,
                      help="If selected, multi-reference alignment (MRA) is performed")
        form.addParam('generateTemplate', BooleanParam, default=True, label='Generate reference template(s):',
                      help="Generate a reference template based on parameters")
        form.addParam('nref', IntParam, label='Number of references for MRA', default=2,
                      condition="mra and generateTemplate",
                      help="Number of references for multi-reference alignment (MRA)")
        form.addParam('templateRef', PointerParam, label="Template", condition="not generateTemplate",
                      pointerClass='Volume, SetOfVolumes',
                      help='The size of the template should be equal or smaller than the size of the particles. If you '
                           'pass a single file in multireference modus (MRA), Dynamo will just made copies of it.')
        form.addParam('useRandomTable', BooleanParam, label='Use a random table',
                      default=False, condition="generateTemplate", help="Rotate particles randomly")
        form.addParam('compensateMissingWedge', BooleanParam, label='Compensate for missing wedge', default=False,
                      condition="generateTemplate", help="")

        form.addSection(label='Masks')
        form.addParam('mask', PointerParam, label="Alignment mask", pointerClass='VolumeMask', allowsNull=True,
                      help='Needs to have the same dimmensionality as the template. It does NOT need to be smoothened. '
                           'A hard mask (only 0 and 1 values) is expected for the Roseman"s-like normalization')
        form.addParam('cmask', PointerParam, label="Classification mask", pointerClass='VolumeMask', allowsNull=True,
                      help='Needs to have the same dimmensionality as the template. It does NOT need to be smoothened. '
                           'A hard mask (only 0 and 1 values) is expected for the Roseman"s-like normalization')
        form.addParam('fmask', PointerParam, label="Fourier mask on template", pointerClass='Volume, SetOfVolumes',
                      allowsNull=True,
                      help='The fmask indicates which fourier coefficients are present at the starting reference volume'
                           '. The file should contain only ones or zeros.')
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
        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('alignStep')
        self._insertFunctionStep('createOutputStep')
        self._insertFunctionStep('closeSetsStep')

    # --------------------------- STEPS functions -------------------------------

    def convertInputStep(self):
        fnDirData = self._getExtraPath("data")
        makePath(fnDirData)
        fnRoot = join(fnDirData, "particle_")
        inputVols = self.inputVolumes.get()
        # writeSetOfEm(inputVols, fnRoot, self)
        writeSetOfVolumes(inputVols, fnRoot, 'id')
        pcaInt = int(self.pca.get())
        mraInt = int(self.mra.get())
        ccmatrixInt = int(self.ccmatrix.get())
        dim, _, _ = self.inputVolumes.get().getDimensions()
        self.projName = 'dynamoAlignmentProject'
        fnTable = self._getExtraPath("initial.tbl")
        fhTable = open(fnTable, 'w')
        writeDynTable(fhTable, inputVols)
        fhTable.close()
        fhCommands = open(self._getExtraPath("commands1.doc"), 'w')
        content = "dcp.new('%s','data','data','gui',0);" % self.projName + \
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
                  "dvput('%s', 'limm', '%s');" % (self.projName, self.limm) + \
                  "dvput('%s', 'destination', 'standalone');" % self.projName + \
                  "dvput('%s', 'cores', %d);" % (self.projName, self.numberOfThreads.get()) + \
                  "dvput('%s', 'mwa', 0);" % self.projName

        if self.useGpu.get():
            content += "dvput('%s', 'destination', 'standalone_gpu');" % (self.projName) + \
                       "dvput('%s', 'gpu_motor', 'spp');" % (self.projName) + \
                       "dvput('%s', 'gpu_identifier_set', %s);" % (self.projName, self.getGpuList())

        template = self.templateRef.get()
        fmask = self.fmask.get()
        if template is not None:
            if isinstance(template, Volume):
                writeVolume(template, join(self._getExtraPath(), 'template'))
                content += "dvput('%s', 'template', 'template.mrc');" % self.projName
                content += "dvput('%s', 'table', 'initial.tbl');" % self.projName
            else:
                makePath(self._getExtraPath('templates'))
                writeSetOfVolumes(template, join(self._getExtraPath(), 'templates/template_initial_ref_'), 'ix')
                for ix, template1 in enumerate(self.templateRef.get().iterItems()):
                    copy(join(self._getExtraPath(), 'initial.tbl'),
                         join(self._getExtraPath(), 'templates/table_initial_ref_%03d.tbl' % int(ix+1)))
                content += "dvput('%s', 'template', 'templates');" % self.projName
                content += "dvput('%s', 'table', 'templates');" % self.projName
                content += "dvput('%s', 'nref', '%d');" % (self.projName, self.nref)
                if fmask is None:
                        makePath(self._getExtraPath('fmasks'))
                        dimsf = inputVols.getFirstItem().getDim()
                        sizef = template.getSize()
                        content += "dynamo_write_multireference(ones(%d,%d,%d),'fmask','fmasks','refs',1:%d);" \
                                   % (dimsf[0], dimsf[1], dimsf[2], sizef)
                        content += "dvput('%s', 'fmask', 'fmasks');" % self.projName
        else:
            if not self.mra.get():
                if self.useRandomTable.get():
                    content += "dynamo_table_perturbation('initial.tbl','pshift',0,'paxis',0,'pnarot',360,'o'," \
                               "'initial.tbl');"
                if self.compensateMissingWedge.get():
                    content += "dynamo_table_randomize_azimuth('initial.tbl','o','initial.tbl');"
                content += "dynamo_average('data','table','initial.tbl','o','template.em');"
                content += "dvput('%s', 'template', 'template.em');" % self.projName
                content += "dvput('%s', 'table', 'initial.tbl');" % self.projName

            else:
                makePath(self._getExtraPath('templates'))
                for ix, template in enumerate(range(self.nref.get())):
                    if self.useRandomTable.get():
                        content += "dynamo_table_perturbation('initial.tbl','pshift',0,'paxis',0,'pnarot',360,'o'," \
                                   "'initial.tbl');"
                    if self.compensateMissingWedge.get():
                        content += "dynamo_table_randomize_azimuth('initial.tbl','o','initial.tbl');"
                    content += "dynamo_average('data','table','initial.tbl','o'," \
                               "'templates/template_initial_ref_%03d.em');" % int(ix+1)
                    copy(join(self._getExtraPath(), 'initial.tbl'),
                         join(self._getExtraPath(), 'templates/table_initial_ref_%03d.tbl' % int(ix+1)))
                content += "dvput('%s', 'template', 'templates');" % self.projName
                content += "dvput('%s', 'table', 'templates');" % self.projName
                content += "dvput('%s', 'nref', '%d');" % (self.projName, self.nref)
                if fmask is None:
                        makePath(self._getExtraPath('fmasks'))
                        dimsf = inputVols.getFirstItem().getDim()
                        sizef = self.nref.get()
                        content += "dynamo_write_multireference(ones(%d,%d,%d),'fmask','fmasks','refs',1:%d);" \
                                   % (dimsf[0], dimsf[1], dimsf[2], sizef)
                        content += "dvput('%s', 'fmask', 'fmasks');" % self.projName

        if self.mask.get() is not None:
            writeVolume(self.mask.get(), join(self._getExtraPath(), 'mask'))
            content += "dvput('%s', 'mask', 'mask.mrc');" % self.projName
        if self.cmask.get() is not None:
            writeVolume(self.cmask.get(), join(self._getExtraPath(), 'cmask'))
            content += "dvput('%s', 'cmask', 'cmask.mrc');" % self.projName
        if fmask is not None:
            if isinstance(fmask, Volume) or isinstance(fmask, SubTomogram) or isinstance(fmask, VolumeMask):
                writeVolume(fmask, join(self._getExtraPath(), 'fmask'))
                content += "dvput('%s', 'fmask', 'fmask.mrc');" % self.projName
            else:
                makePath(self._getExtraPath('fmasks'))
                writeSetOfVolumes(fmask, join(self._getExtraPath(), 'fmasks/fmask_initial_ref_'), 'ix')
                content += "dvput('%s', 'fmask', 'fmasks');" % self.projName
        if self.smask.get() is not None:
            writeVolume(self.smask.get(), join(self._getExtraPath(), 'smask'))
            content += "dvput('%s', 'smask', 'smask.mrc');" % self.projName

        content += "dynamo_data_format('templates/template_initial_ref_*.mrc'," \
                   "'templates','modus','convert','extension','.em');" \
                   "dynamo_data_format('data/particle_*.mrc'," \
                   "'data','modus','convert','extension','.em');"

        fhCommands.write(content)
        fhCommands.close()

        Plugin.runDynamo(self, 'commands1.doc', cwd=self._getExtraPath())

        cleanPattern(self._getExtraPath(join('data', 'particle_*.mrc')))

        if mraInt == 1:
            if self.generateTemplate.get() == 0:
                numberOfrefs = self.templateRef.get().getSize()
            else:
                numberOfrefs = self.nref.get()
            for iref in range(numberOfrefs):
                # TODO: Check if particle_%05d.em files are saved at some point in templates folder
                # TODO: or if we are trying to convert templates/template_initial_ref_%03d.mrc to em format
                if os.path.isfile(self._getExtraPath('templates/particle_%05d.em' % int(iref+1))):
                    rename(join(self._getExtraPath(), 'templates/particle_%05d.em' % int(iref+1)),
                           join(self._getExtraPath(), 'templates/template_initial_ref_%03d.em' % int(iref+1)))
                    remove(join(self._getExtraPath(), 'templates/template_initial_ref_%03d.mrc' % int(iref+1)))

    def alignStep(self):
        fhCommands2 = open(self._getExtraPath("commands2.doc"), 'w')
        content2 = "dvcheck('%s');" % self.projName
        content2 += "dvunfold('%s');dynamo_execute_project %s" % (self.projName, self.projName)
        fhCommands2.write(content2)
        fhCommands2.close()
        Plugin.runDynamo(self, 'commands2.doc', cwd=self._getExtraPath())

    def createOutputStep(self):
        self.projName = 'dynamoAlignmentProject'
        niters = int(self.numberOfIters.get())
        if self.mra.get() == 0:
            self.subtomoSet = self._createSetOfSubTomograms()
            inputSet = self.inputVolumes.get()
            self.subtomoSet.copyInfo(inputSet)
            self.fhTable = open(self._getExtraPath('%s/results/ite_%04d/averages/refined_table_ref_001_ite_%04d.tbl')
                                % (self.projName, niters, niters), 'r')
            self.subtomoSet.copyItems(inputSet, updateItemCallback=self._updateItem)
            self.fhTable.close()
            averageSubTomogram = AverageSubTomogram()
            averageSubTomogram.setFileName(
                self._getExtraPath('%s/results/ite_%04d/averages/average_ref_001_ite_%04d.em')
                % (self.projName, niters, niters))
            averageSubTomogram.setSamplingRate(inputSet.getSamplingRate())
            self._defineOutputs(outputSubtomograms=self.subtomoSet)
            self._defineSourceRelation(self.inputVolumes, self.subtomoSet)
            self._defineOutputs(averageSubTomogram=averageSubTomogram)
            self._defineSourceRelation(self.inputVolumes, averageSubTomogram)
        else:
            fhSurvivRefs = open(self._getExtraPath('%s/results/ite_%04d/currently_surviving_references_ite_%04d.txt')
                                % (self.projName, niters, niters), 'r')
            nline = next(fhSurvivRefs).rstrip()
            for nref in range(len(nline.split())):
                ref = int(nline.split()[nref])
                subtomoSet = self._createSetOfSubTomograms()
                inputSet = self.inputVolumes.get()
                subtomoSet.copyInfo(inputSet)
                self.fhTable = open(
                    self._getExtraPath('%s/results/ite_%04d/averages/refined_table_ref_%03d_ite_%04d.tbl') %
                    (self.projName, niters, ref, niters), 'r')
                subtomoSet.copyItems(inputSet, updateItemCallback=self._updateItem)
                self.fhTable.close()
                averageSubTomogram = AverageSubTomogram()
                averageSubTomogram.setFileName(
                    self._getExtraPath('%s/results/ite_%04d/averages/average_ref_%03d_ite_%04d.em')
                    % (self.projName, niters, ref, niters))
                averageSubTomogram.setSamplingRate(inputSet.getSamplingRate())

                name = 'outputSubtomogramsRef%s' % str(ref)
                args = {}
                args[name] = subtomoSet
                subtomoSet.setStreamState(Set.STREAM_OPEN)
                self._defineOutputs(**args)
                self._defineSourceRelation(inputSet, subtomoSet)

                name2 = 'AverageRef%s' % str(ref)
                args2 = {}
                args2[name2] = averageSubTomogram
                self._defineOutputs(**args2)
                self._defineSourceRelation(inputSet, averageSubTomogram)

    def closeSetsStep(self):
        if self.mra.get() == True:
            for outputset in self._iterOutputsNew():
                if isinstance(outputset, SetOfSubTomograms):
                    outputset[1].setStreamState(Set.STREAM_CLOSED)
            self._store()

    # --------------------------- UTILS functions --------------------------------
    def _updateItem(self, item, row):
        readDynTable(self, item)

    # --------------------------- INFO functions --------------------------------
    def _validate(self):
        validateMsgs = []
        if self.mra.get() == True:
            if self.generateTemplate.get():
                if self.nref.get() <= 1:
                    validateMsgs.append('If MRA is selected, number of references should be greater than one.')
            elif isinstance(self.templateRef.get(), Volume):
                validateMsgs.append('If MRA is selected, number of references should be greater than one.')
            elif self.templateRef.get().getSize() <= 1:
                validateMsgs.append('If MRA is selected, number of references should be greater than one.')
        return validateMsgs

    def _summary(self):
        summary = []
        summary.append("Input subtomograms: %d" % self.inputVolumes.get().getSize())
        if self.fmask.get() is not None:
            if isinstance(self.fmask.get(), Volume) or isinstance(self.fmask.get(), SubTomogram):
                summary.append("Input fmask: %s" % self.fmask.get())
            else:
                summary.append("Input fmasks: %d" % self.fmask.get().getSize())
        else:
            summary.append("Fmask(s) generated")
        if self.mra.get() == True:
            summary.append("Perform MRA with %s references" % self.nref.get())
        else:
            summary.append("No mra")
        if self.generateTemplate.get():
            summary.append("Template(s) generated")
        else:
            if isinstance(self.templateRef.get(), Volume) or isinstance(self.templateRef.get(), SubTomogram):
                summary.append("Provided template: %s" % self.templateRef.get())
            else:
                summary.append("Provided templates: %d" % self.templateRef.get().getSize())
        return summary

    def _methods(self):
        methods = []
        methods.append(
            'We aligned %d subtomograms from %s using Dynamo Subtomogram averaging.'
            % (self.inputVolumes.get().getSize(), self.getObjectTag('inputVolumes')))
        return methods

    def _citations(self):
        return ['CASTANODIEZ2012139']
