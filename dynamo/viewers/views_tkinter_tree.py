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

import threading
from os.path import abspath, join, isdir
from dynamo import Plugin, VLL_FILE, CATALOG_BASENAME, CATALOG_FILENAME
from dynamo.utils import getCurrentTomoTxtFile
from pyworkflow.gui.dialog import ToolbarListDialog
from pyworkflow.utils import makePath
from pyworkflow.utils.process import runJob


class DynamoTomoDialog(ToolbarListDialog):
    """
    This class extend from ListDialog to allow calling
    an Eman subprocess from a list of Tomograms.
    """

    def __init__(self, parent, path, **kwargs):
        self.proc = None
        self.currentTomoTxtFile = None
        self.tomo = None
        self.path = path
        self.coordsFileDict = kwargs.get('writtenCoordsFilesDict', None)
        self.calledFromViewer = kwargs.get('calledFromViewer', False)
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
            with open(self.currentTomoTxtFile, 'r') as fn:
                self.tomo.count = int(fn.read())
                self.tree.update()

    def doubleClickOnTomogram(self, e=None):
        self.tomo = e
        self.currentTomoTxtFile = getCurrentTomoTxtFile(self.path, e)
        # Create a catalogue with the Coordinates to be visualized
        extraPath = join(self.path)
        catalogue = join(extraPath, CATALOG_BASENAME)
        catalogueWithExt = join(extraPath, CATALOG_FILENAME)
        listTomosFile = join(extraPath, VLL_FILE)
        if self.calledFromViewer:
            makePath(catalogue)  # Needed for a correct catalog creation
            codeFile = abspath(join(self.path, 'writectlg.m'))
            contents = "if exist('%s', 'file')==2\n" % catalogueWithExt
            contents += "delete('%s')\n" % catalogueWithExt
            contents += "end\n"
            contents += "dcm -create %s -vll %s\n" % (catalogue, listTomosFile)  # create the catalog
            contents += "catalogue=dread('%s')\n" % catalogueWithExt  # read it
            contents += "for idv=1:length(catalogue.volumes)\n"
            contents += "cvolume = catalogue.volumes{idv}\n"
            contents += "tomoPath=cvolume.file\n"
            contents += "tomoIndex=cvolume.index\n"
            contents += "currentTomoModelsDir = fullfile('%s', 'tomograms', ['volume_', num2str(tomoIndex)], " \
                        "'models')\n" % catalogue
            contents += "if not(exist(currentTomoModelsDir)==7)\n"
            contents += "mkdir(currentTomoModelsDir)\n"
            contents += "end\n"
            contents += "[~,tomoName,~]=fileparts(tomoPath)\n"
            contents += "coordFile=fullfile('%s', [tomoName, '.txt'])\n" % extraPath
            contents += "if ~isfile(coordFile)\n"
            contents += "continue\n"
            contents += "end\n"
            contents += "coords_ids=readmatrix(coordFile,'Delimiter',' ')\n"
            contents += "idm_vec=unique(coords_ids(:,4))'\n"
            contents += "for idm=1:length(idm_vec)\n"
            contents += "model_name=['model_',num2str(idm)]\n"
            contents += "coords=coords_ids(coords_ids(:,4)==idm,1:4)\n"
            contents += "modelFilePath = fullfile(currentTomoModelsDir, [model_name, '.omd'])\n"
            contents += "model=dmodels.general()\n"  # Create a general model for each tomogram
            contents += "model.file = modelFilePath\n"
            contents += "model.name = model_name\n"
            contents += "model.cvolume = cvolume\n"
            contents += "addPoint(model, coords(:,1:3), coords(:,4))\n"  # Add the points to the model
            contents += "save(modelFilePath, 'model')\n"  # Save the model to a file in the expected location
            contents += "catalogue.models = [catalogue.models, model]\n"  # Add the current model to the catalog
            # contents += "general.linkCatalogue('%s','i',tomoIndex)\n" % catalogue
            # TODO: think about a method to check if the user has made any changes in the points, apart from just visualizing
            contents += "end\n"
            contents += "end\n"
            contents += "save('%s', 'catalogue')\n" % catalogueWithExt
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
        codeFilePath = join(self.path, "DynamoPicker.m")
        catalogue = join(self.path, CATALOG_BASENAME)

        # Write code to Matlab code file
        with open(codeFilePath, 'w') as codeFid:
            content = "ctlgNoExt = '%s'\n" % catalogue
            content += "ctlgName = [ctlgNoExt, '.ctlg']\n"
            content += "ctlg = dread(ctlgName)\n"  # Load the catalogue
            content += "tomoFiles = cellfun(@(x) x.file, ctlg.volumes, 'UniformOutput', false)\n"  # Cell with the tomo names
            content += "currentTomoInd = find(ismember(tomoFiles, '%s'))\n" % abspath(tomo.getFileName())  # Index of the current tomo in the catalog
            content += "loadDynSyntax = sprintf('dtmslice @{%s}%i', ctlgNoExt, currentTomoInd)\n"
            content += "eval(loadDynSyntax)\n"  # Launch Dynamo's picker GUI
            content += "uiwait(dpkslicer.getHandles().figure_fastslicer)\n"  # Wait until it's closed
            content += "modeltrack.saveAllInCatalogue\n"  # Save in the catalog
            content += "models = dcmodels(ctlgNoExt, 'i', currentTomoInd)\n"
            content += "nParticles = 0\n"
            content += "for i=1:length(models)\n"
            content += "model = dread(models{i})\n"
            content += "newParts = size(model.points, 1)\n"
            content += "nParticles = nParticles + newParts\n"  # Sum the no. of particles from all the models generated for the current tomogram
            content += "end\n"
            content += "fid = fopen('%s', 'w')\n" % self.currentTomoTxtFile
            content += "fprintf(fid, '%i', nParticles)\n"  # Save the number of particles to a text file
            content += "fclose(fid)\n"
            content += "nParticles"
            codeFid.write(content)
        return codeFilePath
