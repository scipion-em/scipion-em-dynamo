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
import math
from collections import OrderedDict

from pwem.protocols import EMProtocol
from pyworkflow import BETA
from pyworkflow.protocol import EnumParam
from tomo.protocols import ProtTomoBase

dynBinningChoices = OrderedDict()
dynBinningChoices[0] = 'none'
dynBinningChoices[2] = 'bin 2'
dynBinningChoices[4] = 'bin 4'
dynBinningChoices[8] = 'bin 8'


class DynamoProtocolBase(EMProtocol, ProtTomoBase):
    _devStatus = BETA

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def getBinningFactor(self, forDynamo=True):
        """From Dynamo: a binning Factor of 1 will decrease the size of the Tomograms by 2,
        a Binning Factor of 2 by 4... So Dynamo interprets the binning factor as 2**binFactor, while IMOD
        interprets it literally. Thus, this method will convert the binning introduced by the user in the
        Dynamo convention. Values that do not correspond to an integer value are floor-rounded internally by Dynamo,
        interpreting, for example, a literal bin factor of 6 as log2(6) = 2.5849, and then it's rewritten as 2,
        which corresponds to a literal binning factor of 2 ** 2 = 4."""
        binFactor = list(dynBinningChoices.keys())[self.binning.get()]
        return math.log2(binFactor) if forDynamo else binFactor

    @staticmethod
    def _addBinningParam(form):
        form.addParam('binning', EnumParam,
                      choices=dynBinningChoices.values(),
                      display=EnumParam.DISPLAY_HLIST,
                      default=2,  # Bin 4
                      label='Binning',
                      help='The tomogram reconstruction implies the generation of the interpolated tilt series. The '
                           'binning factor introduced will be applied to both.\n\n'
                           '*Note for Dynamo users*: the binning factor introduced must be interpreted literally, '
                           'not as in Dynamo (power of 2). This parameter will be transformed internally, so Dynamo '
                           'behaves as epexted. For example: a bin factor of 4 here will be passed to Dynamo as 2, a '
                           'bin factor of 2 will be passed as 1, and so on.\n\n')
