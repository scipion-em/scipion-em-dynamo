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

import os

from pwem.viewers import ObjectView

import pyworkflow.viewer as pwviewer

import tomo.objects
from tomo.viewers.views_tkinter_tree import MeshesTreeProvider

from dynamo.viewers.views_tkinter_tree import DynamoDialog



class DynamoDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        tomo.objects.SetOfMeshes
    ]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return ObjectView(
            self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, **kwargs):
        views = []
        cls = type(obj)

        if issubclass(cls, tomo.objects.SetOfMeshes):
            outputMeshes = obj

            # meshList = [item.clone() for item in outputMeshes.iterItems()]

            meshList = []
            for item in obj.iterItems():
                mesh = item.clone()
                mesh.setVolume(item.getVolume())
                meshList.append(mesh)
            # meshList.reverse()

            path = self.protocol._getTmpPath()

            meshProvider = MeshesTreeProvider(meshList,)

            setView = DynamoDialog(self._tkRoot, path, True, provider=meshProvider,)

        return views
