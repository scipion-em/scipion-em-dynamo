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
from enum import Enum
from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.utils import getExt
from tomo.objects import SetOfSubTomograms, SetOfTomograms, AverageSubTomogram, SetOfCoordinates3D, SetOfMeshes
from tomo.protocols import ProtTomoBase


class DynamoOutputs(Enum):
    tomograms = SetOfTomograms
    subtomograms = SetOfSubTomograms
    average = AverageSubTomogram
    coordinates = SetOfCoordinates3D
    meshes = SetOfMeshes


class DynamoProtocolBase(EMProtocol, ProtTomoBase):

    _devStatus = BETA
    _possibleOutputs = DynamoOutputs

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def getBinningFactor(self, forDynamo=True):
        """From Dynamo: a binning Factor of 1 will decrease the size of the Tomograms by 2,
        a Binning Factor of 2 by 4... So Dynamo interprets the binning factor as 2**binFactor, while IMOD
        interprets it literally. Thus, this method will convert the binning introduced by the user in the
        Dynamo convention"""
        return (2 ** (self.binning.get() - 1)) / 2 if forDynamo else self.binning.get()

    def isCompatibleFileFormat(self):
        """Compatible with MRC and em (MRC with that extension)"""
        compatibleExts = ['.em', '.mrc']
        return True if getExt(self.inputTomos.getFirstItem().getFileName()) not in compatibleExts else False
