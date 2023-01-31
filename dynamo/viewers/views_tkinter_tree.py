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
from os.path import abspath, join, isdir

import numpy as np
from pyworkflow import utils as pwutils
from pyworkflow.utils.process import runJob
from pyworkflow.gui.dialog import ToolbarListDialog

from dynamo import Plugin, VLL_FILE, CATALOG_BASENAME, CATALOG_FILENAME


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
                                   allowSelect=False,
                                   **kwargs)

    def refresh_gui(self):
        if self.proc.is_alive():
            self.after(1000, self.refresh_gui)
        else:
            outPath = join(self.path, pwutils.removeBaseExt(self.tomo.getFileName()) + '.txt')
            self.tomo.count = np.loadtxt(outPath, delimiter=' ').shape[0]
            self.tree.update()

    def doubleClickOnTomogram(self, e=None):
        self.tomo = e
        # Create a catalogue with the Coordinates to be visualized
        mainPath = abspath(join(self.path))
        catalogue = join(mainPath, CATALOG_BASENAME)
        listTomosFile = join(mainPath, VLL_FILE)

        if not isdir(catalogue):
            codeFile = abspath(join(self.path, 'writectlg.m'))
            contents = "dcm -create %s -fromvll %s\n" % (catalogue, listTomosFile)
            contents += "catalogue=dread('%s')\n" % join(mainPath, CATALOG_FILENAME)
            contents += "nVolumes=length(catalogue.volumes)\n"
            contents += "for idv=1:nVolumes\n"
            contents += "tomoPath=catalogue.volumes{idv}.fullFileName()\n"
            contents += "tomoIndex=catalogue.volumes{idv}.index\n"
            contents += "[~,tomoName,~]=fileparts(tomoPath)\n"
            contents += "coordFile=fullfile('%s', tomoName '.txt')\n" % mainPath
            contents += "if ~isfile(coordFile)\n"
            contents += "continue\n"
            contents += "end\n"
            contents += "coords_ids=readmatrix(coordFile,'Delimiter',' ')\n"
            contents += "idm_vec=unique(coords_ids(:,4))'\n"
            contents += "for idm=idm_vec\n"
            contents += "model_name=['model_',num2str(idm)]\n"
            contents += "coords=coords_ids(coords_ids(:,4)==idm,1:4)\n"
            contents += "general=dmodels.general()\n"
            contents += "general.name=model_name\n"
            contents += "addPoint(general,coords(:,1:3),coords(:,4))\n"
            contents += "general.linkCatalogue('%s','i',tomoIndex,'s',1)\n" % catalogue
            contents += "general.saveInCatalogue()\n"
            contents += "end\n"
            contents += "end\n"
            contents += "exit\n"
            with open(codeFile, 'w') as codeFid:
                codeFid.write(contents)
            args = ' %s' % codeFile
            runJob(None, Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

        self.proc = threading.Thread(target=self.lanchDynamoForTomogram, args=(self.tomo,))
        self.proc.start()
        self.after(1000, self.refresh_gui)

    def lanchDynamoForTomogram(self, tomo):
        commandsFile = self.writeMatlabCode(tomo)
        args = ' %s' % commandsFile
        runJob(None, Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    def writeMatlabCode(self, tomo):
        # Initialization params
        codeFilePath = join(os.getcwd(), "DynamoPicker.m")
        catalogue = join(self.path, CATALOG_BASENAME)

        # Write code to Matlab code file
        codeFid = open(codeFilePath, 'w')
        content = "catalogue_name='%s'\n" \
                  "c=dread(strcat(catalogue_name,'.ctlg'))\n" \
                  "n=cellfun(@(c) c.fullFileName,c.volumes,'UniformOutput',false)\n" \
                  "idt=find(cell2mat(cellfun(@(c) strcmp(c,'%s'),n,'UniformOutput',false)))\n" \
                  "eval(strcat('dtmslice @{%s}',num2str(idt)))\n" \
                  "modeltrack.loadFromCatalogue('handles',c,'full',true,'select',false)\n" \
                  "uiwait(dpkslicer.getHandles().figure_fastslicer)\n" \
                  "modeltrack.saveAllInCatalogue\n" \
                  "models = dread(dcmodels(catalogue_name,'i', idt))\n" \
                  "outFile='%s'\n" \
                  "savePath='%s'\n" \
                  "outPoints=[outFile '.txt']\n" \
                  "outAngles=['angles_' outFile '.txt']\n" \
                  "crop_points=[]\n" \
                  "crop_angles=[]\n" \
                  "model_id=0\n" \
                  "for model=models\n" \
                  "model_id=model_id+1\n" \
                  "if iscell(model)\n" \
                  "crop_points=[crop_points; [model{end}.points model_id*ones(length(model{end}.points),1)]]\n" \
                  "else\n" \
                  "crop_points=[crop_points; [model.points model_id*ones(length(model.points),1)]]\n" \
                  "end\n" \
                  "end\n" \
                  "if ~isempty(crop_points)\n" \
                  "writematrix(crop_points,fullfile(savePath,outPoints),'Delimiter',' ')\n" \
                  "end\n" \
                  "exit\n" % (abspath(catalogue), abspath(tomo.getFileName()), catalogue,
                              pwutils.removeBaseExt(tomo.getFileName()), self.path)
        # "writematrix(crop_angles,fullfile(savePath,outAngles),'Delimiter',' ')\n" \
        # "crop_angles=[crop_angles; [model{end}.crop_angles model_id*ones(length(model{end}.crop_angles),1)]]\n" \
        # "crop_angles=[crop_angles; [model.crop_angles model_id*ones(length(model.crop_angles),1)]]\n" \
        codeFid.write(content)
        codeFid.close()

        return codeFilePath