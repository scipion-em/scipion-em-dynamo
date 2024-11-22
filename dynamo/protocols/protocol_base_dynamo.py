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

