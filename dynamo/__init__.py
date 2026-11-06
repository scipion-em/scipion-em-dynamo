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
from os.path import join, dirname
import subprocess
import pwem
import pyworkflow
import pyworkflow.utils as pwutils
from pyworkflow import VarTypes
from .constants import *

__version__ = '3.4.0'
_logo = "icon.png"
_references = ['CASTANODIEZ2012139']


class Plugin(pwem.Plugin):
    _homeVar = DYNAMO_HOME
    _pathVars = [DYNAMO_HOME]

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(DYNAMO_HOME, 'dynamo-{}'.format(DEFAULT_VERSION))
        cls._defineVar(DYNAMO_CUDA_LIB, pwem.Config.CUDA_LIB,
                       description="Path to the CUDA lib folder to use to build and run the "
                                   "Dynamo CUDA binaries. If unset, the global CUDA_LIB is used.",
                       var_type=VarTypes.FOLDER)
        cls._defineVar(DYNAMO_CUDA_BIN, pwem.Config.CUDA_BIN,
                       description="Path to the CUDA bin folder to use to build and run the "
                                   "Dynamo CUDA binaries. If unset, the global CUDA_BIN is used.",
                       var_type=VarTypes.FOLDER)

    @classmethod
    def getEnviron(cls, gpuId=0):
        """ Create the needed environment for Dynamo programs. """
        environ = pwutils.Environ(os.environ)
        dyn_home = cls.getHome()
        # For GPU, we need to add to LD_LIBRARY_PATH the path to Cuda/lib.
        # Use the plugin specific DYNAMO_CUDA_LIB variable (falls back to the
        # global CUDA_LIB) so build and runtime use the same CUDA.
        cudaLib = cls.getVar(DYNAMO_CUDA_LIB, pwem.Config.CUDA_LIB)
        paths = (dyn_home + '/MCRLinux/runtime/glnxa64',
                 dyn_home + '/MCRLinux/bin/glnxa64',
                 dyn_home + '/MCRLinux/sys/os/glnxa64',
                 dyn_home + '/MCRLinux/sys/opengl/lib/glnxa64',
                 cudaLib)

        environ.update({
            'PATH': dyn_home + '/matlab/bin:' + dyn_home + '/matlab/src:' + dyn_home + '/cuda/bin:' + dyn_home + '/mpi',
            'LD_LIBRARY_PATH': ":".join(paths),
            'DYNAMO_ROOT': dyn_home,
            'CUDA_VISIBLE_DEVICES': str(gpuId)
        }, position=pwutils.Environ.BEGIN)
        return environ

    @classmethod
    def getDynamoProgram(cls):
        return join(cls.getHome(), 'matlab', 'bin', DYNAMO_PROGRAM)

    @classmethod
    def runDynamo(cls, protocol, args, cwd=None, gpuId=0):
        """ Run Dynamo command from a given protocol. """
        # args will be the .doc file which contains the MATLAB code
        program = cls.getDynamoProgram()
        protocol.runJob(program, args, env=cls.getEnviron(gpuId=gpuId), cwd=cwd)

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

        # Dynamo 1.1.532
        commands = "./dynamo_setup_linux.sh "  # OpenMP commands
        # Clear executable-stack flag from bundled MCR 9.9 libs (required on
        # modern glibc >=2.41 / Debian 13, else 'dynamo' exits 251 with
        # "cannot enable executable stack as shared object requires")
        mcrPatch = join(os.path.dirname(os.path.abspath(__file__)), 'mcr_execstack_patch.py')
        commands += "&& python3 %s ." % mcrPatch
        if useGpu:
            # CUDA paths used to build the GPU binaries, resolved from the plugin
            # specific DYNAMO_CUDA_LIB / DYNAMO_CUDA_BIN variables (defined under
            # [PLUGINS] in scipion.conf). They follow the global CUDA_LIB / CUDA_BIN
            # convention and point to <CUDA_ROOT>/lib64 and <CUDA_ROOT>/bin, so the
            # CUDA *root* needed by the makefile is recovered with dirname(). If
            # neither was customized they default to the global CUDA_LIB / CUDA_BIN.
            cudaLib = cls.getVar(DYNAMO_CUDA_LIB, pwem.Config.CUDA_LIB)
            cudaBin = cls.getVar(DYNAMO_CUDA_BIN, pwem.Config.CUDA_BIN)
            if cudaLib and cudaLib != pwem.Config.CUDA_LIB:
                cudaHome = dirname(cudaLib)
            elif cudaBin and cudaBin != pwem.Config.CUDA_BIN:
                cudaHome = dirname(cudaBin)
            else:
                cudaHome = dirname(cudaLib)
            # Cuda commands.
            # The iteration binary is compiled with nvcc but linked with the system
            # g++ using -static-libstdc++ -static-libgcc: the protocol runtime puts
            # the old MCR libstdc++ (GLIBCXX_3.4.25) on LD_LIBRARY_PATH and the
            # system libstdc++ is needed by the modern build (GLIBCXX_3.4.32), which
            # the dynamic linker would otherwise resolve to the MCR one ("version
            # `GLIBCXX_3.4.32' not found").
            commands += ("&& cd cuda && ./config.sh " + cudaHome + " && " +
                         "sed -i 's/cuCtxCreate(&ctx, 0, dev);/cuCtxCreate(\\&ctx, NULL, 0, dev);/' src/infoGPU.cu && " +
                         "make clean && " +
                         cudaBin + "/nvcc -c src/dynamo_compute_iteration_spp.cu -o .spp.o " +
                         "-I" + cudaHome + "/include -Isrc -Isrc/classes -w && " +
                         "g++ -c src/classes/gpuAngles.cpp -o .gpuAngles.o -I" + cudaHome + "/include -w && " +
                         "g++ -static-libstdc++ -static-libgcc -o bin/dynamo_compute_iteration_spp " +
                         ".spp.o .gpuAngles.o -L" + cudaLib + " -lcudart -lcufft -lcublas && " +
                         "rm -f .spp.o .gpuAngles.o && make extended && touch cuda_compiled")
        commands = [(commands, 'cuda/cuda_compiled')]
        env.addPackage(DYNAMO_PROGRAM, version=DYNAMO_VERSION_1_1_532,
                       tar='dynamo-v-1.1.532_MCR-9.9.0_GLNXA64_withMCR.tar',
                       createBuildDir=True,
                       commands=commands,
                       default=True)

    @classmethod
    def checkDynamoVersion(cls):
        """Admitted versions must be higher or equal than version 1.1.532. The version number is extracted
        from the home variable, splitting by '-' and removing a possible 'v' for version. Finally the points are
        removed and the final numeric string is cast to an integer."""
        msg = []
        # Check the installed binary version
        dynamoVer = cls.getHome().split('-')[-1].replace('v', '')
        if int(dynamoVer.replace('.', '')) < MINIMUM_VERSION_NUM:
            msg = ['The Dynamo version pointed by variable %s '
                   '(%s) is not supported --> %s.\n\n'
                   'Please, update the variable value or comment it in %s' %
                   (DYNAMO_HOME, dynamoVer, DYNAMO_VERSION_1_1_532,
                    pyworkflow.Config.SCIPION_CONFIG)]
        return msg

    @classmethod
    def validateInstallation(cls):
        msg = super().validateInstallation()
        msg += cls.checkDynamoVersion()
        return msg
