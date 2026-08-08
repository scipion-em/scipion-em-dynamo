# **************************************************************************
# *
# * Authors:    Scipion Team (scipion@cnb.csic.es)
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
from pwem.protocols import EMProtocol
from pyworkflow.protocol import IntParam
from tomo.protocols import ProtTomoBase

IN_TOMOS = 'inputTomos'
IN_COORDS = 'inputCoords'
BIN_THREADS_MSG = 'Threads used by Dynamo each time it is called by Scipion'


class DynamoProtocolBase(EMProtocol, ProtTomoBase):
    """
    Provides a shared foundation for Dynamo-based cryo-electron tomography
    protocols. The class standardizes common tomography utilities, data access
    operations, and computational settings required for subtomogram processing
    workflows within Dynamo environments.

    AI Generated:

    Dynamo Protocol Base (DynamoProtocolBase) - User Manual
        Overview

        The Dynamo Protocol Base class serves as the common infrastructure layer
        for tomography protocols executed within Dynamo-based cryo-electron
        tomography workflows. Its purpose is to provide a consistent framework
        for managing tomographic datasets, coordinate information, binning
        conventions, and computational resources across multiple subtomogram
        analysis procedures.

        In practical cryo-electron tomography workflows, many processing stages
        share common requirements such as access to tomograms, particle
        coordinates, and standardized execution settings. By centralizing these
        operations, the protocol base ensures that downstream subtomogram
        alignment, averaging, classification, and reconstruction procedures can
        operate consistently and reproducibly within the same computational
        environment.

        Tomographic Data Management

        The protocol framework is designed to simplify the handling of tomograms
        and associated coordinate information throughout subtomogram analysis
        pipelines. Tomograms represent the reconstructed three-dimensional
        cellular or molecular environment, while coordinates define the spatial
        locations of particles or macromolecular assemblies extracted from those
        tomograms.

        From a biological perspective, maintaining consistent access to these
        datasets is essential because all subsequent subtomogram operations rely
        on accurate spatial relationships between particles and their original
        tomographic context. Errors in tomogram referencing or coordinate
        interpretation can propagate through the workflow and negatively affect
        averaging quality, alignment accuracy, and structural interpretation.

        Computational Resource Configuration

        Cryo-electron tomography processing is computationally demanding,
        particularly during subtomogram averaging and iterative alignment
        procedures. The protocol base therefore includes mechanisms for defining
        computational threading behavior in Dynamo executions.

        In routine biological workflows, selecting an appropriate number of
        computational threads can substantially improve processing speed while
        maintaining stable execution. Larger subtomogram datasets or high-
        resolution analyses generally benefit from increased parallelization,
        whereas smaller exploratory runs may require fewer computational
        resources.

        Binning and Scale Interpretation

        One of the important technical aspects addressed by the protocol base is
        the interpretation of tomogram binning conventions. Different tomography
        software environments may define binning factors differently, potentially
        leading to inconsistencies in scale interpretation and voxel size
        handling.

        The framework provides a standardized interpretation layer so that
        subtomogram workflows remain internally consistent regardless of the
        convention used during data preparation. This is biologically important
        because voxel scaling directly influences the interpretation of molecular
        dimensions, structural resolution, and coordinate accuracy.

        In practical workflows, proper handling of binning factors ensures that
        tomograms, subtomograms, and coordinate systems remain geometrically
        compatible throughout all stages of analysis.

        Integration Within Tomography Pipelines

        The protocol base is intended to support a wide range of Dynamo-based
        subtomogram analysis procedures, including particle extraction,
        alignment, classification, averaging, and refinement. By defining common
        interfaces and shared behaviors, it allows higher-level tomography
        protocols to focus on biological analysis rather than infrastructure
        management.

        This modular organization is particularly valuable in large cryo-electron
        tomography projects where multiple processing stages are chained together
        into reproducible workflows. Consistent data handling improves workflow
        stability and reduces the likelihood of mismatches between tomographic
        metadata and subtomogram processing parameters.

        Practical Recommendations

        In biological practice, users should ensure that tomograms and particle
        coordinates originate from compatible preprocessing pipelines before
        starting Dynamo-based analyses. Consistent voxel size calibration and
        careful binning interpretation are especially important when combining
        datasets from multiple microscopes, acquisition sessions, or processing
        environments.

        Computational thread settings should be adjusted according to available
        hardware resources and dataset size. Excessive parallelization may not
        always improve performance if memory availability becomes limiting during
        large subtomogram operations.

        Final Perspective

        The Dynamo Protocol Base framework provides the structural foundation
        required for robust and reproducible subtomogram analysis workflows. By
        standardizing tomographic data access, computational configuration, and
        scaling conventions, it supports reliable execution of downstream cryo-
        electron tomography procedures and contributes to more consistent
        biological interpretation of in situ structural data.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @staticmethod
    def insertBinThreads(form, helpMsg=BIN_THREADS_MSG):
        form.addParam('binThreads', IntParam,
                      label='Dynamo threads',
                      default=3,
                      help=helpMsg)

    def getBinningFactor(self, fromDynamo: bool = True) -> int:
        """From Dynamo: a binning Factor of 1 will decrease the size of the Tomograms by 2,
        a Binning Factor of 2 by 4... So Dynamo interprets the binning factor as 2**binFactor, while IMOD
        interprets it literally. Thus, this method will convert the binning introduced by the user in the
        Dynamo convention"""
        return (2 ** (self.binning.get() - 1)) / 2 if fromDynamo else self.binning.get()

    def getInTomos(self, isPointer: bool = False):
        inTomos = getattr(self, IN_TOMOS, None)
        return inTomos if isPointer else inTomos.get()

    def getInCoords(self, isPointer: bool = False):
        inCoords = getattr(self, IN_COORDS, None)
        return inCoords if isPointer else inCoords.get()

