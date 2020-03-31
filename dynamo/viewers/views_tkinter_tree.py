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

import os, threading
import numpy as np

import pyworkflow.config as conf
from pyworkflow import utils as pwutils
from pyworkflow.utils.process import runJob
from pyworkflow.gui.dialog import ToolbarListDialog

from dynamo import Plugin


class DynamoTomoDialog(ToolbarListDialog):
    """
    This class extend from ListDialog to allow calling
    an Eman subprocess from a list of Tomograms.
    """

    def __init__(self, parent, path, **kwargs):
        self.path = path
        self.provider = kwargs.get("provider", None)
        ToolbarListDialog.__init__(self, parent,
                                   "Tomogram List",
                                    allowsEmptySelection=False,
                                    itemDoubleClick=self.doubleClickOnTomogram,
                                    **kwargs)

    def refresh_gui(self):
        if self.proc.isAlive():
            self.after(1000, self.refresh_gui)
        else:
            self.tree.update()

    def doubleClickOnTomogram(self, e=None):
        self.tomo = e
        self.proc = threading.Thread(target=self.lanchDynamoForTomogram, args=(self.tomo,))
        self.proc.start()
        self.after(1000, self.refresh_gui)

    def lanchDynamoForTomogram(self, tomo):
        commandsFile = self.writeMatlabCode(tomo)
        args = ' %s' % commandsFile
        runJob(None, Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    def writeMatlabCode(self, tomo):
        # Initialization params
        codeFilePath = os.path.join(os.getcwd(), "DynamoPicker.m")
        listTomosFile = os.path.join(self.path, "tomos.vll")
        catalogue = os.path.abspath(os.path.join(self.path, "tomos"))

        # Create list of tomos file
        tomoFid = open(listTomosFile, 'w')
        tomoPath = os.path.abspath(tomo.getFileName())
        tomoFid.write(tomoPath + '\n')
        tomoFid.close()

        # Write code to Matlab code file
        codeFid = open(codeFilePath, 'w')
        content = "catalogue_name='%s'\n" \
                  "c=dread(strcat(catalogue_name,'.ctlg'))\n" \
                  "n=cellfun(@(c) c.fullFileName,c.volumes,'UniformOutput',false)\n" \
                  "idt=find(cell2mat(cellfun(@(c) strcmp(c,'%s'),n,'UniformOutput',false)))\n" \
                  "dtmslice %s -c %s \n" \
                  "uiwait(dpkslicer.getHandles().figure_fastslicer)\n" \
                  "models = dread(dcmodels(catalogue_name,'i', idt))\n" \
                  "outFile='%s'\n" \
                  "savePath='%s'\n" \
                  "outPoints=[outFile '.txt']\n" \
                  "outAngles=['angles_' outFile '.txt']\n" \
                  "crop_points=[]\n" \
                  "crop_angles=[]\n" \
                  "for model=models\n" \
                  "if iscell(model)\n" \
                  "crop_points=[crop_points; model{end}.crop_points]\n" \
                  "crop_angles=[crop_angles; model{end}.crop_angles]\n" \
                  "else\n" \
                  "crop_points=[crop_points; model.crop_points]\n" \
                  "crop_angles=[crop_angles; model.crop_angles]\n" \
                  "end\n" \
                  "end\n" \
                  "if ~isempty(crop_points)\n" \
                  "writematrix(crop_points,fullfile(savePath,outPoints),'Delimiter',' ')\n" \
                  "writematrix(crop_angles,fullfile(savePath,outAngles),'Delimiter',' ')\n" \
                  "end\n" \
                  'exit' % (catalogue, tomo.getFileName(), tomo.getFileName(), catalogue,
                            pwutils.removeBaseExt(tomo.getFileName()), self.path)
        codeFid.write(content)
        codeFid.close()

        return codeFilePath