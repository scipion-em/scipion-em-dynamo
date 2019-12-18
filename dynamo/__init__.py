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
import os
import pyworkflow
import pyworkflow.utils as pwutils

_logo = "icon.png"
_references = ['CASTANODIEZ2012139']

class Plugin(pyworkflow.em.Plugin):
    @classmethod
    def getEnviron(cls):
        """ Create the needed environment for Dynamo programs. """
        environ = pwutils.Environ(os.environ)
        environ.update({
            'PATH': '/home/davidh/dynamo_v1.146/matlab/bin:/home/davidh/dynamo_v1.146/matlab/src:/home/davidh/dynamo_v1.146/cuda/bin:/home/davidh/dynamo_v1.146/mpi:/home/davidh/dynamo_v1.146/examples:/home/davidh/dynamo_v1.146/doc${PATH}',
            'LD_LIBRARY_PATH': '/home/davidh/dynamo_v1.146/MCRLinux/runtime/glnxa64:/home/davidh/dynamo_v1.146/MCRLinux/bin/glnxa64:/home/davidh/dynamo_v1.146/MCRLinux/sys/os/glnxa64:${LD_LIBRARY_PATH}',
            'DYNAMO_ROOT': '/home/davidh/dynamo_v1.146/'
        }, position=pwutils.Environ.BEGIN)

        return environ

pyworkflow.em.Domain.registerPlugin(__name__)
