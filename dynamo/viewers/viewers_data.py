# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
# *              Scipion Team (scipion@cnb.csic.es)
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

from os.path import join, abspath, exists
from dynamo.viewers.DynamoTomoProvider import DynamoTomogramProvider
from pwem.viewers import ObjectView
import pyworkflow.viewer as pwviewer
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.utils.properties import Message
from tomo.objects import SetOfCoordinates3D, SetOfMeshes
from dynamo.viewers.views_tkinter_tree import DynamoTomoDialog
from dynamo import VLL_FILE


class DynamoDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [SetOfMeshes, SetOfCoordinates3D]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return ObjectView(
            self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, **kwargs):
        views = []

        if isinstance(obj, SetOfCoordinates3D) or isinstance(obj, SetOfMeshes):
            outputCoords = obj
            precedents = outputCoords.getPrecedents()
            tomoIdDict = {tomo.getTsId(): tomo.clone() for tomo in precedents}
            tomoList = list(tomoIdDict.values())
            nCoordsDict = {}
            for tomo in tomoList:
                tomoId = tomo.getTsId()
                nCoordsDict[tomoId] = len(list(outputCoords.iterItems(where='_tomoId=="%s"' % tomoId)))

            path = self.protocol._getExtraPath()
            tomoProvider = DynamoTomogramProvider(tomoList, path, nParticlesDict=nCoordsDict)
            listTomosFile = join(path, VLL_FILE)

            # Create list of tomos file (VLL file), only required if using the Dynamo viewer in a protocol from any
            # other plugin
            if not(exists(listTomosFile)):
                with open(listTomosFile, 'w') as tomoFid:
                    for tomo in tomoList:
                        tomoPath = abspath(tomo.getFileName())
                        tomoFid.write(tomoPath + '\n')

            self.dlg = DynamoTomoDialog(self._tkRoot, path,
                                        provider=tomoProvider,
                                        calledFromViewer=True)

            import tkinter as tk
            frame = tk.Frame()
            # TODO: check if the user has made changes and only ask in that case
            # TODO: Refactor this and make a correct output creating according to the last changes
            if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, frame):
                # textFile2Coords(self.protocol, outputCoords.getPrecedents(), path, directions=False)
                self.protocol._createOutput()

        return views
