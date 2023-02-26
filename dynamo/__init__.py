# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from os.path import join
import subprocess
import pwem
import pyworkflow.utils as pwutils
from .constants import *

__version__ = '3.1.14'
_logo = "icon.png"
_references = ['CASTANODIEZ2012139']

class Plugin(pwem.Plugin):
    _homeVar = DYNAMO_HOME
    _pathVars = [DYNAMO_HOME]
    _url = "https://github.com/scipion-em/scipion-em-dynamo"
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
        paths = (dyn_home + '/MCRLinux/runtime/glnxa64',
                 dyn_home + '/MCRLinux/bin/glnxa64',
                 dyn_home + '/MCRLinux/sys/os/glnxa64',
                 dyn_home + '/MCRLinux/sys/opengl/lib/glnxa64',
                 pwem.Config.CUDA_LIB)

        environ.update({
            'PATH': dyn_home + '/matlab/bin:' + dyn_home + '/matlab/src:' + dyn_home + '/cuda/bin:' + dyn_home + '/mpi',
            'LD_LIBRARY_PATH': ":".join(paths),
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
        # First, we check if CUDA is installed in the system
        preMsgs = []
        cudaMsgs = []
        nvidiaDriverVer = None
        if os.environ.get('CUDA', 'True') == 'True':
            try:
                nvidiaDriverVer = subprocess.Popen(["nvidia-smi",
                                                    "--query-gpu=driver_version",
                                                    "--format=csv,noheader"],
                                                   env=cls.getEnviron(),
                                                   stdout=subprocess.PIPE
                                                   ).stdout.read().decode('utf-8').split(".")[0]
            except Exception as e:
                preMsgs.append(str(e))

        if nvidiaDriverVer is not None:
            preMsgs.append("CUDA support find. Driver version: %s" % nvidiaDriverVer)
            msg = "Dynamo installed with CUDA SUPPORT."
            cudaMsgs.append(msg)
            useGpu = True
        else:
            preMsgs.append("CUDA will NOT be USED. (not found)")
            msg = ("Dynamo installed without GPU. Just CPU computations "
                   "enabled (slow computations). To enable CUDA, "
                   "set CUDA=True in 'scipion.conf' file")
            cudaMsgs.append(msg)
            useGpu = False

        SW_EM = env.getEmFolder()
        dyn_folder = 'dynamo-%s' % DYNAMO_VERSION

        compile_cuda = "echo ' > %s' && " % preMsgs
        if useGpu:
            compile_cuda += "cd %s/%s/cuda && make clean && ./config.sh && " \
                            "make all && " \
                            "touch cuda_compiled && " % (SW_EM, dyn_folder)
        compile_cuda += "echo ' > %s'" % cudaMsgs
        commands = [(compile_cuda, '%s/%s/cuda/cuda_compiled' % (SW_EM, dyn_folder))]

        env.addPackage('dynamo',
                       version='1.146',
                       tar='dynamo-1.146.tar.gz',
                       commands=commands,
                       default=False)

        # Dynamo 1.1.532
        commands = [("cd cuda && make clean && ./config.sh %s && make all && touch cuda_compiled" % os.path.dirname(pwem.Config.CUDA_LIB), 'cuda/cuda_compiled')]
        env.addPackage('dynamo', version='1.1.532',
                       tar='dynamo-v-1.1.532_MCR-9.9.0_GLNXA64_withMCR.tar',
                       createBuildDir = True,
                       commands=commands)

