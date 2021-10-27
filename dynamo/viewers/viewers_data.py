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
import numpy as np

from pwem.viewers import ObjectView

import pyworkflow.viewer as pwviewer
import pyworkflow.utils as pwutils
from pyworkflow.gui.dialog import askYesNo
from pyworkflow.utils.properties import Message
from pyworkflow.utils.process import runJob

import tomo.objects
from tomo.viewers.views_tkinter_tree import MeshesTreeProvider, TomogramsTreeProvider
import tomo.constants as const

from dynamo.viewers.views_tkinter_tree import DynamoTomoDialog
from dynamo.convert import textFile2Coords, matrix2eulerAngles
from dynamo import Plugin



class DynamoDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        tomo.objects.SetOfMeshes,
        tomo.objects.SetOfCoordinates3D
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

        if issubclass(cls, tomo.objects.SetOfCoordinates3D) or issubclass(cls, tomo.objects.SetOfMeshes):
            outputCoords = obj
            tomos = outputCoords.getPrecedents()

            volIds = outputCoords.aggregate(["MAX", "COUNT"], "_volId", ["_volId"])
            volIds = [(d['_volId'], d["COUNT"]) for d in volIds]

            tomoList = []
            for objId in volIds:
                tomogram = tomos[objId[0]].clone()
                tomogram.count = objId[1]
                tomoList.append(tomogram)

            path = self.protocol._getExtraPath()

            tomoProvider = TomogramsTreeProvider(tomoList, path, 'txt', )

            listTomosFile = os.path.join(path, "tomos.vll")

            # Create list of tomos file
            tomoFid = open(listTomosFile, 'w')
            for tomoFile in tomos.getFiles():
                tomoPath = os.path.abspath(tomoFile)
                tomoFid.write(tomoPath + '\n')
            tomoFid.close()

            for tomogram in tomoList:
                outFileCoord = os.path.join(path, pwutils.removeBaseExt(tomogram.getFileName())) + ".txt"
                # outFileAngle = os.path.join(path, 'angles_' + pwutils.removeBaseExt(tomogram.getFileName())) + ".txt"
                coords_tomo = []
                # angles_tomo = []
                for coord in outputCoords.iterCoordinates(tomogram):
                    coords_tomo.append(list(coord.getPosition(const.BOTTOM_LEFT_CORNER)) + [coord.getGroupId()])
                    # angles_shifts = matrix2eulerAngles(coord.getMatrix())
                    # angles_tomo.append(angles_shifts[:3])
                if coords_tomo:
                    np.savetxt(outFileCoord, np.asarray(coords_tomo), delimiter=' ')
                    # np.savetxt(outFileAngle, np.asarray(angles_tomo), delimiter=' ')


            # Create a catalogue with the Coordinates to be visualized
            catalogue = os.path.abspath(os.path.join(path, "tomos"))
            if not os.path.isdir(catalogue):
                codeFile = os.path.abspath(self.protocol._getExtraPath('writectlg.m'))
                contents = "dcm -create %s -fromvll %s\n" \
                           "path='%s'\n" \
                           "catalogue=dread(['%s' '.ctlg'])\n" \
                           "nVolumes=length(catalogue.volumes)\n" \
                           "for idv=1:nVolumes\n" \
                           "tomoPath=catalogue.volumes{idv}.fullFileName()\n" \
                           "tomoIndex=catalogue.volumes{idv}.index\n" \
                           "[~,tomoName,~]=fileparts(tomoPath)\n" \
                           "coordFile=[path '/' tomoName '.txt']\n" \
                           "if ~isfile(coordFile)\n" \
                           "continue\n" \
                           "end\n" \
                           "coords_ids=readmatrix(coordFile,'Delimiter',' ')\n" \
                           "idm_vec=unique(coords_ids(:,4))'\n" \
                           "for idm=idm_vec\n" \
                           "model_name=['model_',num2str(idm)]\n" \
                           "coords=coords_ids(coords_ids(:,4)==idm,1:4)\n" \
                           "general=dmodels.general()\n" \
                           "general.name=model_name\n" \
                           "addPoint(general,coords(:,1:3),coords(:,4))\n" \
                           "general.linkCatalogue('%s','i',tomoIndex,'s',1)\n" \
                           "general.saveInCatalogue()\n" \
                           "end\n" \
                           "end\n" \
                           "exit\n" % (catalogue, os.path.abspath(listTomosFile),
                                       os.path.abspath(path), catalogue, catalogue)
                codeFid = open(codeFile, 'w')
                codeFid.write(contents)
                codeFid.close()
                args = ' %s' % codeFile
                runJob(None, Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

            self.dlg = DynamoTomoDialog(self._tkRoot, path, provider=tomoProvider)

            import tkinter as tk
            frame = tk.Frame()
            if askYesNo(Message.TITLE_SAVE_OUTPUT, Message.LABEL_SAVE_OUTPUT, frame):
                textFile2Coords(self.protocol, outputCoords.getPrecedents(), path, directions=False)

        return views
