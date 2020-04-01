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


class DynamoDialog(ToolbarListDialog):
    """
    This class extend from ListDialog to allow calling
    an Eman subprocess from a list of Tomograms.
    """

    def __init__(self, parent, path, viewer, **kwargs):
        self.path = path
        self.provider = kwargs.get("provider", None)
        if viewer:
            ToolbarListDialog.__init__(self, parent,
                                       "Tomogram List",
                                       allowsEmptySelection=False,
                                       itemDoubleClick=self.doubleClickViewer,
                                       **kwargs)
        else:
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
            pwutils.cleanPath(os.path.join(self.path, 'tomos'))
            pwutils.cleanPath(os.path.join(self.path, 'tomos.ctlg'))

    def refresh_gui_viewer(self):
        if self.proc.isAlive():
            self.after(1000, self.refresh_gui_viewer)
        else:
            pwutils.cleanPath(os.path.join(self.path, 'mesh.txt'))


    def doubleClickOnTomogram(self, e=None):
        self.tomo = e
        self.proc = threading.Thread(target=self.lanchDynamoForTomogram, args=(self.tomo,))
        self.proc.start()
        self.after(1000, self.refresh_gui)

    def doubleClickViewer(self, e=None):
        self.tomo = e
        self.proc = threading.Thread(target=self.lanchDynamoForViewing, args=(self.tomo,))
        self.proc.start()
        self.after(1000, self.refresh_gui_viewer)

    def lanchDynamoForTomogram(self, tomo):
        inputFilePath = self._writeInputFile(tomo)

        args = ' %s' % inputFilePath
        runJob(None, Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    def lanchDynamoForViewing(self, mesh):
        self.mesh = mesh
        inputFilePath = self._writeViewerFile(mesh.getVolume())

        args = ' %s' % inputFilePath
        runJob(None, Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    def _writeInputFile(self, tomo):
        inputFilePath = os.path.join(conf.Config.SCIPION_TMP, "commands.doc")
        tomoFile = os.path.join(conf.Config.SCIPION_TMP, "tomos.vll")
        catalogue = os.path.join(self.path, 'tomos')
        tomoFid = open(tomoFile, 'w')
        tomoName = tomo.getFileName()
        tomoName = os.path.basename(tomoName)
        tomoFid.write(tomoName + '\n')
        tomoFid.close()
        tomoBase = pwutils.removeBaseExt(tomo.getFileName())
        inputFid = open(inputFilePath, 'w')
        if not os.path.isfile(os.path.join(self.path, tomoBase + '.txt')):
            content = 'dcm -create %s -fromvll %s \n' \
                      'newClass = \'Yes\' \n' \
                      'class = 0 \n' \
                      'while strcmp(newClass, \'Yes\') \n' \
                      'class = class + 1 \n' \
                      'm = dmodels.membraneByLevels()\n' \
                      'm.linkCatalogue(\'%s\',\'i\',1,\'s\',class)\n' \
                      'm.saveInCatalogue()\n' \
                      'dtmslice %s -c %s \n' \
                      'modeltrack.addOne(\'model\',m)\n' \
                      'modeltrack.setActiveModel(class)\n' \
                      'uiwait(dpkslicer.getHandles().figure_fastslicer)\n' \
                      'm = dread(dcmodels(\'%s\',\'i\', 1))\n' \
                      'if class > 1 \n' \
                      'aux = [m{class}.points class*ones(length(m{class}.points), 1)] \n' \
                      'mat = [mat ; aux] \n' \
                      'else \n' \
                      'mat = [m.points class*ones(length(m.points), 1)] \n' \
                      'end \n' \
                      'newClass = questdlg(\'Create a new class\', \'Dialog\', \'Yes\', \'No\', \'Yes\') \n' \
                      'end \n' \
                      'writematrix(mat, \'%s\')\n' \
                      'exit' % (catalogue, tomoFile, catalogue, tomo.getFileName(), catalogue,
                                catalogue, os.path.join(self.path, tomoBase + '.txt'))
            # 'writematrix(m.mesh.tr.ConnectivityList, \'%s\')\n' \
        else:
            content = 'dcm -create %s -fromvll %s\n' \
                      'idx=1\n' \
                      'modelData=readmatrix(\'%s\')\n' \
                      'class=unique(modelData(:,4))\n' \
                      'for item=class\'\n' \
                      'm=dmodels.membraneByLevels()\n' \
                      'logic=repmat(modelData(:,4)==item,1,4)\n' \
                      'data=modelData(logic)\n' \
                      'data=reshape(data,length(data)/4,4)\n' \
                      'addPoint(m,data(:,1:3),data(:,3))\n' \
                      'm.linkCatalogue(\'%s\',\'i\',1,\'s\',item)\n' \
                      'm.saveInCatalogue()\n' \
                      'modeltrack.addOne(\'model\',m)\n' \
                      'end\n' \
                      'newClass=\'Yes\'\n' \
                      'while strcmp(newClass,\'Yes\')\n' \
                      'class=[class;class(end)+1]\n' \
                      'm=dmodels.membraneByLevels()\n' \
                      'm.linkCatalogue(\'%s\',\'i\',1,\'s\',class(end))\n' \
                      'm.saveInCatalogue()\n' \
                      'dtmslice %s -c %s\n' \
                      'modeltrack.addOne(\'model\',m)\n' \
                      'modeltrack.setActiveModel(class(end))\n' \
                      'uiwait(dpkslicer.getHandles().figure_fastslicer)\n' \
                      'm=dread(dcmodels(\'%s\',\'i\', 1))\n' \
                      'if idx==1 \n' \
                      'aux=[m.points class(end)*ones(length(m.points),1)]\n' \
                      'modelData=[modelData;aux] \n' \
                      'idx=idx+1\n' \
                      'else \n' \
                      'aux=[m{idx}.points class(end)*ones(length(m{idx}.points),1)]\n' \
                      'modelData=[modelData;aux] \n' \
                      'idx=idx+1\n' \
                      'end \n' \
                      'newClass=questdlg(\'Create a new class\',\'Dialog\',\'Yes\',\'No\',\'Yes\')\n' \
                      'end\n' \
                      'writematrix(modelData,\'%s\')\n' \
                      'exit' % (catalogue, tomoFile, os.path.join(self.path, tomoBase + '.txt'),
                                catalogue, catalogue, tomo.getFileName(), catalogue,
                                catalogue, os.path.join(self.path, tomoBase + '.txt'))
        inputFid.write(content)
        inputFid.close()
        return inputFilePath

    def _writeViewerFile(self, tomo):
        inputFilePath = os.path.join(conf.Config.SCIPION_TMP, "commands.doc")
        tomoFile = os.path.join(conf.Config.SCIPION_TMP, "tomos.vll")
        catalogue = os.path.join(self.path, 'tomos')
        tomoFid = open(tomoFile, 'w')
        tomoName = tomo.getFileName()
        tomoName = os.path.basename(tomoName)
        tomoFid.write(tomoName + '\n')
        tomoFid.close()
        mesh_group = self.mesh.getMesh()
        group = np.ones((1, len(mesh_group))) * self.mesh.getGroup()
        meshFile = 'mesh.txt'
        np.savetxt(os.path.join(self.path, meshFile), np.append(mesh_group, group.T, axis=1), fmt='%d', delimiter=',')
        inputFid = open(inputFilePath, 'w')
        content = 'dcm -create %s -fromvll %s\n' \
                  'modelData=readmatrix(\'%s\')\n' \
                  'class=unique(modelData(:,4))\n' \
                  'for item=class\'\n' \
                  'm=dmodels.membraneByLevels()\n' \
                  'logic=repmat(modelData(:,4)==item,1,4)\n' \
                  'data=modelData(logic)\n' \
                  'data=reshape(data,length(data)/4,4)\n' \
                  'addPoint(m,data(:,1:3),data(:,3))\n' \
                  'm.linkCatalogue(\'%s\',\'i\',1,\'s\',item)\n' \
                  'm.saveInCatalogue()\n' \
                  'modeltrack.addOne(\'model\',m)\n' \
                  'end\n' \
                  'dtmslice %s -c %s \n' \
                  'uiwait(dpkslicer.getHandles().figure_fastslicer)\n' \
                  'exit' % (catalogue, tomoFile, os.path.join(self.path, meshFile),
                            catalogue, tomo.getFileName(), catalogue)
        inputFid.write(content)
        inputFid.close()
        return inputFilePath

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
        catalogue = os.path.join(self.path, "tomos")

        # Write code to Matlab code file
        codeFid = open(codeFilePath, 'w')
        content = "catalogue_name='%s'\n" \
                  "c=dread(strcat(catalogue_name,'.ctlg'))\n" \
                  "n=cellfun(@(c) c.fullFileName,c.volumes,'UniformOutput',false)\n" \
                  "idt=find(cell2mat(cellfun(@(c) strcmp(c,'%s'),n,'UniformOutput',false)))\n" \
                  "dtmslice %s -c %s \n" \
                  "modeltrack.loadFromCatalogue('handles',c,'full',true,'select',false)\n" \
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
                  'exit\n' % (catalogue, tomo.getFileName(), tomo.getFileName(), catalogue,
                            pwutils.removeBaseExt(tomo.getFileName()), self.path)
        codeFid.write(content)
        codeFid.close()

        return codeFilePath