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
import pwem
import pyworkflow.utils as pwutils
from .constants import *

__version__ = '3.1.7'
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
        dyn_home = cls.getHome()
        # For GPU, we need to add to LD_LIBRARY_PATH the path to Cuda/lib
        environ.update({
            'PATH': dyn_home + '/matlab/bin:' + dyn_home + '/matlab/src:' + dyn_home + '/cuda/bin:' + dyn_home + '/mpi',
            'LD_LIBRARY_PATH': dyn_home + '/MCRLinux/runtime/glnxa64:' + dyn_home + '/MCRLinux/bin/glnxa64:' +
                               dyn_home + '/MCRLinux/sys/os/glnxa64:' + pwem.Config.CUDA_LIB,
            'DYNAMO_ROOT': dyn_home
        }, position=pwutils.Environ.BEGIN)
        return environ

    @classmethod
    def getDynamoProgram(cls):
        return join(cls.getHome(), 'matlab', 'bin', DYNAMO_PROGRAM)

    @classmethod
    def runDynamo(cls, protocol, args, cwd=None):
        """ Run Dynamo command from a given protocol. """
        # args will be the .doc file which contains the MATLAB code
        program = cls.getDynamoProgram()
        protocol.runJob(program, args, env=cls.getEnviron(), cwd=cwd)

    @classmethod
    def defineBinaries(cls, env):
        # For GPU, we need to compile the CUDA binaries to ensure full compatibility
        SW_EM = env.getEmFolder()
        dyn_folder = 'dynamo-%s' % DYNAMO_VERSION
        compile_cuda = "cd %s/%s/cuda && make clean && ./config.sh && make all && touch cuda_compiled" % (SW_EM, dyn_folder)
        commands = [(compile_cuda, '%s/%s/cuda/cuda_compiled' % (SW_EM, dyn_folder))]

        env.addPackage('dynamo',
                       version='1.146',
                       tar='dynamo-1.146.tar.gz',
                       commands=commands,
                       default=True)

