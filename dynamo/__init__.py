# **************************************************************************
# *
# * Authors:     you (you@yourinstitution.email)
# *
# * your institution
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

from os.path import join
import pyworkflow.em as pwem
import pyworkflow.utils as pwutils
from .constants import *

_logo = "icon.png"
_references = ['CASTANODIEZ2012139']

class Plugin(pwem.Plugin):
    _homeVar = DYNAMO_HOME
    _pathVars = [DYNAMO_HOME]
    # _supportedVersions =

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(DYNAMO_HOME, 'dynamo-{}'.format(DYNAMO_VERSION))

    @classmethod
    def getEnviron(cls):
        """ Create the needed environment for Dynamo programs. """
        environ = pwutils.Environ(os.environ)
        environ.update({
            'PATH': DYNAMO_HOME + '/matlab/bin:' + DYNAMO_HOME + '/matlab/src:' + DYNAMO_HOME + '/cuda/bin:' + DYNAMO_HOME + '/mpi:${PATH}',
            'LD_LIBRARY_PATH': DYNAMO_HOME + '/MCRLinux/runtime/glnxa64:' + DYNAMO_HOME + '/MCRLinux/bin/glnxa64:' + DYNAMO_HOME + '/MCRLinux/sys/os/glnxa64:${LD_LIBRARY_PATH}',
            'DYNAMO_ROOT': DYNAMO_HOME
        }, position=pwutils.Environ.BEGIN)
        return environ

    @classmethod
    def getDynamoProgram(cls):
        return join(DYNAMO_HOME, 'matlab', 'bin', DYNAMO_PROGRAM)


    # @classmethod
    # def defineBinaries(cls, env):
    #     env.addPackage(JANNI_GENMOD,
    #                    version=JANNI_GENMOD_20190703,
    #                    tar='void.tgz',
    #                    commands=[("wget https://github.com/MPI-Dortmund/sphire-janni/raw/master/janni_general_models/"
    #                               + JANNI_GENMOD_20190703_FN, JANNI_GENMOD_20190703_FN)],
    #                    neededProgs=["wget"],
    #                    default=True)

pwem.Domain.registerPlugin(__name__)
