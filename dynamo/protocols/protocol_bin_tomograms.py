# **************************************************************************
# *
# * Authors:    David Herreros Calero (dherreros@cnb.csic.es)
# *             Scipion Team (scipion@cnb.csic.es)
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
from enum import Enum
from os.path import abspath
from dynamo.protocols.protocol_base_dynamo import DynamoProtocolBase, IN_TOMOS
from pwem.convert.headers import setMRCSamplingRate
from pwem.emlib.image import ImageHandler
from pyworkflow.object import Set
from pyworkflow.protocol import params, GT, STEPS_PARALLEL
from pyworkflow.utils import getExt, Message, createLink
from tomo.objects import Tomogram, SetOfTomograms
from dynamo import Plugin


class DynamoBinOuts(Enum):
    tomograms = SetOfTomograms


class DynamoBinTomograms(DynamoProtocolBase):
    """
    Reduces the size of tomographic datasets through controlled binning in
    order to decrease computational cost, simplify visualization, and prepare
    tomograms for downstream subtomogram averaging or particle picking tasks.

    AI Generated:

    Bin Tomograms (DynamoBinTomograms) - User Manual
        Overview

        The Bin Tomograms protocol performs spatial downsampling of a set of
        tomograms using Dynamo-based processing tools. In cryo-electron
        tomography workflows, binning is commonly used to reduce the size of
        volumetric datasets while preserving the overall structural content
        required for intermediate processing steps. This operation is
        especially important when working with large tomograms that would
        otherwise require excessive computational resources for visualization,
        alignment, segmentation, or particle extraction.

        From a biological perspective, tomogram binning is not intended to
        improve resolution or structural detail. Instead, it creates smaller
        and computationally lighter representations of the original data,
        enabling faster exploratory analyses and more efficient iterative
        processing. Researchers frequently apply binning during the early
        stages of tomographic workflows and later return to the original
        unbinned data for high-resolution refinement.

        Inputs and Dataset Preparation

        The protocol requires a set of tomograms as input. These tomograms are
        typically reconstructed volumes originating from cryo-electron
        tomography experiments and may contain cellular environments,
        macromolecular assemblies, membrane systems, or viral particles.

        Before processing, the protocol ensures that the tomograms are
        prepared in a format suitable for Dynamo-based operations. This
        preparation stage guarantees compatibility across different datasets
        and helps standardize the processing workflow. The resulting binned
        tomograms preserve the geometric identity and metadata of the original
        acquisitions while generating new volumes with reduced dimensions.

        Biological Role of Binning

        In practical cryo-ET workflows, binning serves several important
        biological and computational purposes. Large tomograms often contain
        substantial noise and extremely large voxel grids, making direct
        analysis computationally expensive. By reducing the sampling density,
        users can rapidly inspect cellular organization, identify regions of
        interest, and test processing strategies before committing to
        computationally demanding refinements.

        Binning is also commonly used during template matching and particle
        localization. Lower-resolution tomograms can accelerate the detection
        of candidate particles, after which the corresponding coordinates can
        later be refined against higher-resolution datasets. Similarly,
        subtomogram averaging pipelines often begin with heavily binned data
        to establish initial orientations and coarse structural organization.

        Choosing the Binning Factor

        The most important parameter in the protocol is the binning factor,
        which determines the degree of downsampling applied to the tomograms.
        Smaller binning values preserve more structural information but
        require more memory and processing time. Larger binning values produce
        significantly smaller datasets that are easier to manipulate but may
        remove fine structural features.

        In biological practice, moderate binning is frequently sufficient for
        visualization, manual inspection, or low-resolution alignment tasks.
        Strong binning should be applied carefully because excessive
        downsampling may obscure membranes, flexible domains, filamentous
        assemblies, or small macromolecular complexes that are biologically
        important.

        Users should also consider the final biological objective when
        selecting the binning level. For exploratory analyses or coordinate
        generation, aggressive binning may be acceptable. For quantitative
        interpretation or structural comparison, preserving adequate sampling
        becomes much more important.

        Memory Management and Large Tomograms

        Cryo-electron tomography datasets can be extremely large, especially
        when imaging thick cellular specimens or high-resolution tilt series.
        The protocol therefore supports memory-aware processing strategies
        that divide tomograms into manageable regions during execution.

        This approach is particularly valuable for facilities or laboratories
        handling large numbers of tomograms on shared computational resources.
        Users working with limited RAM can process large datasets more safely
        by adjusting memory-related execution parameters. Proper balancing
        between memory usage and processing speed is often essential for
        stable execution on workstation-class hardware.

        Parallel Processing Considerations

        The protocol supports multi-threaded execution in order to accelerate
        tomogram processing. Parallel execution can significantly reduce total
        runtime when multiple tomograms are processed simultaneously or when
        large tomograms are subdivided internally during computation.

        However, increased parallelization also increases memory consumption.
        In practical workflows, users should ensure that the available system
        memory is sufficient to support the chosen level of parallel
        execution. Excessive parallelization on memory-limited systems may
        lead to unstable execution or reduced overall performance.

        Outputs and Their Interpretation

        The protocol produces a new set of tomograms with reduced voxel
        dimensions and updated sampling information. These output tomograms
        maintain the structural organization of the original datasets while
        occupying less disk space and requiring fewer computational
        resources.

        The updated sampling rate is biologically important because all
        downstream processing steps must interpret the voxel size correctly.
        Proper sampling metadata ensures consistency during segmentation,
        coordinate assignment, subtomogram extraction, and averaging.

        The resulting tomograms are typically used as intermediate datasets
        for visualization, particle picking, classification, or subtomogram
        averaging workflows. They can also facilitate rapid quality control
        and exploratory analyses before returning to higher-resolution data.

        Practical Recommendations

        In most biological workflows, moderate binning provides the best
        balance between computational efficiency and preservation of
        structural context. Early exploratory processing often benefits from
        stronger binning, while later refinement stages generally require
        less aggressive downsampling or direct use of unbinned tomograms.

        Users working with crowded cellular environments should evaluate the
        impact of binning carefully because excessive reduction may merge
        neighboring densities or blur thin membrane structures. Visual
        inspection of representative outputs is strongly recommended before
        committing to large-scale downstream analyses.

        When processing many tomograms simultaneously, conservative memory
        settings and moderate parallelization usually provide the most stable
        execution environment. This is especially relevant for institutional
        clusters and shared tomography facilities.

        Final Perspective

        Tomogram binning is a foundational preprocessing step in many
        cryo-electron tomography pipelines. Although computational in nature,
        it strongly influences the speed, accessibility, and practicality of
        downstream structural analysis. Selecting an appropriate binning
        strategy allows researchers to efficiently explore complex biological
        datasets while preserving the structural information necessary for
        meaningful interpretation.
    """

    _label = 'bin tomograms'
    _possibleOutputs = DynamoBinOuts
    stepsExecutionMode = STEPS_PARALLEL

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.finalTomoNamesDict = {}
        self.ih = None
        self.sRate = None
        self.doConvertFiles = None

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam(IN_TOMOS, params.PointerParam,
                      pointerClass='SetOfTomograms',
                      label="Input Tomograms",
                      important=True)
        form.addParam('binning', params.IntParam,
                      default=2,
                      validators=[GT(0)],
                      label="Binning Factor",
                      help="A Binning Factor of 1 means that no binning will be carried out.")
        form.addParam('zChunk', params.IntParam,
                      default=300,
                      expertLevel=params.LEVEL_ADVANCED,
                      label="Number of slices kept in memory",
                      help="Maximum number of Z slices that are kept simultaneously in the memory during the "
                           "binning process. This parameter might be important for larger size tomograms, making "
                           "possible to process them in vertical slabs of thickness = value introduced in the  "
                           "current parameter. This procedure can be accelerated using the multiple threads to engage "
                           "several cores in parallel. However, this will only make sense if the total memory occupied "
                           "by all the slabs simultaneously in memory in a given time fits in the RAM of the machine.")
        self.insertBinThreads(form,
                              helpMsg='Number of threads used by Dynamo each time it is called in the protocol '
                                      'execution. For example, if 2 Scipion threads and 3 Dynamo threads are set, '
                                      'the tomograms will be processed in groups of 2 at the same time with a call '
                                      'of tomo3d with 3 threads each, so 6 threads will be used at the same time. '
                                      'Beware the memory of your machine has memory enough to load together the '
                                      'number of tomograms specified by Scipion threads.')
        form.addParallelSection(threads=1, mpi=0)

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._initialize()
        stepIds = []
        for tsId in self.tomoDict.keys():
            cInPid = self._insertFunctionStep(self.convertInputStep, tsId,
                                              prerequisites=[],
                                              needsGPU=False)
            binId = self._insertFunctionStep(self.binTomosStep, tsId,
                                             prerequisites=cInPid,
                                             needsGPU=False)
            cOutId = self._insertFunctionStep(self.createOutputStep, tsId,
                                              prerequisites=binId,
                                              needsGPU=False)
            stepIds.append(cOutId)
        self._insertFunctionStep(self._closeOutputSet,
                                 prerequisites=stepIds,
                                 needsGPU=False)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        inTomos = self.getInTomos()
        self.ih = ImageHandler()
        self.tomoDict = {tomo.getTsId(): tomo.clone() for tomo in inTomos}
        self.doConvertFiles = not self.isCompatibleFileFormat()
        self.sRate = inTomos.getSamplingRate() * self.getBinningFactor(fromDynamo=False)

    def convertInputStep(self, tsId: str):
        if self.doConvertFiles:
            tomo = self.tomoDict[tsId]
            origName = tomo.getFileName()
            finalName = self.getConvertedOrLinkedTsFn(tsId)
            self.ih.convert(origName, finalName)

    def binTomosStep(self, tsId: str):
        mFile = self.createMCodeFile(tsId)
        args = ' %s' % mFile
        self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    def createOutputStep(self, tsId: str):
        with self._lock:
            outTomos = self.getOutputSetOfTomograms()
            inTomo = self.tomoDict[tsId]
            outFn = self.getOutTsFn(tsId)
            setMRCSamplingRate(outFn, self.sRate)  # Update the apix value in file header
            tomo = Tomogram()
            tomo.copyInfo(inTomo)
            tomo.setSamplingRate(self.sRate)
            tomo.setFileName(outFn)
            outTomos.append(tomo)
            outTomos.update(tomo)
            outTomos.write()
            self._store(outTomos)

    # --------------------------- DEFINE utils functions ----------------------
    def isCompatibleFileFormat(self):
        """Compatible with MRC and em (MRC with that extension)"""
        compatibleExts = ['.em', '.mrc']
        return True if (getExt(self.getInTomos().getFirstItem().getFileName())
                        in compatibleExts) else False

    def getConvertedOrLinkedTsFn(self, tsId: str):
        return self._getExtraPath(f'in_{tsId}.mrc')

    def getOutTsFn(self, tsId: str):
        return self._getExtraPath(f'{tsId}.mrc')

    def createMCodeFile(self, tsId: str):
        # FROM DYNAMO:
        # ______________________________________________________________________________________________
        # bin(fileIn,fileOut,binFactor,varargin)
        # Create a Volume template object
        # p.addParamValue('slabSize',[],'short','ss');
        # p.addParamValue('matlabWorkers',0,'short','mw');
        # p.addParamValue('maximumMegaBytes',[]);
        # p.addParamValue('showStatistics',false,'short','sst');
        # ______________________________________________________________________________________________
        tomo = self.tomoDict[tsId]
        origName = abspath(tomo.getFileName())
        finalName = abspath(self.getOutTsFn(tsId))
        mFile = self._getExtraPath('binTomograms.m')
        with open(mFile, 'w') as codeFile:
            content = ("dpktomo.tools.bin('%s', '%s', %i, 'slabSize', %i, 'matlabWorkers', %i, "
                       "'showStatistics', true)\n") % (origName, finalName, super().getBinningFactor(),
                                                       self.zChunk.get(), self.binThreads.get())
            codeFile.write(content)
        return mFile

    def getOutputSetOfTomograms(self) -> SetOfTomograms:
        outTomograms = getattr(self, self._possibleOutputs.tomograms.name, None)
        if outTomograms:
            outTomograms.enableAppend()
            tomograms = outTomograms
        else:
            tomograms = SetOfTomograms.create(self._getPath(), template='tomograms%s.sqlite')
            tomograms.copyInfo(self.getInTomos())
            tomograms.setSamplingRate(self.sRate)
            tomograms.setStreamState(Set.STREAM_OPEN)
            setattr(self, self._possibleOutputs.tomograms.name, tomograms)
            self._defineOutputs(**{self._possibleOutputs.tomograms.name: tomograms})
            self._defineSourceRelation(self.inputTomos, tomograms)

        return tomograms

    # --------------------------- DEFINE info functions ----------------------
    def _methods(self):
        methodsMsgs = ["*Binning Factor*: %s" % self.binning.get()]
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
