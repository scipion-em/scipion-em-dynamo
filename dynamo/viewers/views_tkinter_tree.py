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
import logging
logger = logging.getLogger(__name__)
import threading
from os import remove
from os.path import abspath, join, exists
from shutil import rmtree

from dynamo import Plugin, VLL_FILE, CATALOG_BASENAME, PROJECT_DIR, PRJ_FROM_VIEWER, TOMOGRAMS_DIR, \
    DATA_MODIFIED_FROM_VIEWER
from dynamo.utils import getCurrentTomoCountFile, getDynamoModels, getCatalogFile
from pyworkflow.gui.dialog import ToolbarListDialog, FloatingMessage
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
        self.catalgueMngCode = ''
        self.path = path
        self.coordsFileDict = kwargs.get('writtenCoordsFilesDict', None)
        self.calledFromViewer = kwargs.get('calledFromViewer', False)
        self.provider = kwargs.get("provider", None)
        ToolbarListDialog.__init__(self, parent,
                                   "Tomogram List",
                                   allowsEmptySelection=False,
                                   itemDoubleClick=self.doubleClickOnTomogram,
                                   allowSelect=False,
                                   cancelButton=True,
                                   **kwargs)

    def refresh_gui(self):
        if self.proc.is_alive():
            self.after(1000, self.refresh_gui)
        else:
            self.tree.update()

    def doubleClickOnTomogram(self, e=None):
        self.tomo = e
        self.currentTomoTxtFile = getCurrentTomoCountFile(self.path, e)
        # Create a catalogue with the Coordinates to be visualized
        extraPath = self.path
        catalogue = getCatalogFile(extraPath, withExt=False)
        catalogueWithExt = getCatalogFile(extraPath)
        listTomosFile = join(extraPath, VLL_FILE)
        if self.calledFromViewer and not self._isADynamoProj(extraPath):
            makePath(join(extraPath, PROJECT_DIR, PRJ_FROM_VIEWER))
            makePath(catalogue)  # Needed for a correct catalog creation
            # Dynamo fails if trying to create a catalog that already exists, so the previous one is deleted
            contents = "if exist('%s', 'file') == 2\n" % catalogueWithExt
            contents += "delete('%s')\n" % catalogueWithExt
            contents += "end\n"
            # Create the catalog with the tomograms involves, read it and use that info for a coherent indexation
            # of elements within Dynamo
            contents += "dcm -create %s -vll %s\n" % (catalogue, listTomosFile)  # create the catalog
            contents += "catalogue=dread('%s')\n" % catalogueWithExt  # read it
            # Create a model for each particles of the same groupId for each tomogram
            contents += "for idv=1:length(catalogue.volumes)\n"
            contents += "cvolume = catalogue.volumes{idv}\n"
            contents += "tomoPath = cvolume.file\n"
            contents += "tomoIndex = cvolume.index\n"
            contents += "currentTomoModelsDir = fullfile('%s', 'tomograms', ['volume_', num2str(tomoIndex)], " \
                        "'models')\n" % catalogue
            contents += "[~,tomoName,~] = fileparts(tomoPath)\n"
            contents += "coordFile = fullfile('%s', [tomoName, '.txt'])\n" % extraPath  # Get current coordinates file
            contents += "if exist(coordFile, 'file') == 2\n"
            contents += "s = dir(coordFile)\n"  # Check if the coordinates files is empty
            contents += "if s.bytes == 0\n"
            contents += "continue\n"
            contents += "end\n"
            contents += "else\n"
            contents += "continue\n"
            contents += "end\n"
            contents += "data = cellfun(@(x) regexp(x,',','Split'), importdata(coordFile), 'un', 0)\n"
            contents += "data = vertcat(data{:})\n"
            contents += "coordsMatrix = cell2mat(cellfun(@(x) str2double(x), data(:, 1:4), 'un', 0))\n"
            contents += "modelTypeList = data(:, 5)\n"
            contents += "idm_vec = unique(coordsMatrix(:,4))'\n"  # Get the groupIds
            contents += "for idm=1:length(idm_vec)\n"
            contents += "model_type = modelTypeList{idm}\n"
            contents += "model_name = [model_type, '_', num2str(idm)]\n"
            contents += "coords = coordsMatrix(coordsMatrix(:,4)==idm_vec(idm),:)\n"  # Use them for logical indexing of the coords
            contents += "modelFilePath = fullfile(currentTomoModelsDir, [model_name, '.omd'])\n"
            contents += "model=eval(['dmodels.', model_type, '()'])\n"  # Create a model of the same type as registered for each groupId in each tomogram
            contents += "model.file = modelFilePath\n"
            contents += "model.name = ['m', model_name]\n"
            contents += "model.cvolume = cvolume\n"
            contents += "nParticles = size(coords, 1)\n"
            contents += "model.individual_labels = 1:nParticles\n"
            contents += "addPoint(model, coords(:,1:3), coords(:,4))\n"  # Add the points to the model
            contents += "model.linkCatalogue('%s','i',tomoIndex,'s',1)\n" % abspath(catalogue)
            contents += "model.saveInCatalogue()\n"
            contents += "end\n"
            contents += "end\n"
            contents += "dynamo_write(catalogue, '%s')\n" % catalogueWithExt
        else:
            prjModifiedFile = join(extraPath, PROJECT_DIR, DATA_MODIFIED_FROM_VIEWER)
            if exists(prjModifiedFile):
                remove(prjModifiedFile)
            contents = "dcm -create %s -vll %s\n" % (catalogue, listTomosFile)

        self.catalgueMngCode = contents
        self.proc = threading.Thread(target=self.lanchDynamoForTomogram, args=(self.tomo,))
        self.proc.start()
        self.after(1000, self.refresh_gui)

    def lanchDynamoForTomogram(self, tomo):

        # TODO: Remove xPos and yPos when auto-centering released by pyworkflow
        msg = FloatingMessage(self, "Launching dynamo. Please wait!.", xPos=100, yPos=100)
        msg.show()
        #self.showFloatingMessage("Launching dynamo. Please wait!.")

        commandsFile = self.writeMatlabCode(tomo)
        args = ' %s' % commandsFile
        runJob(logger, Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

        msg.close()
        #self.closeFloatingMessage()

    def writeMatlabCode(self, tomo):
        # Initialization params
        codeFilePath = join(self.path, "DynamoPicker.m")
        catalogue = join(self.path, CATALOG_BASENAME)

        # Write code to Matlab code file
        with open(codeFilePath, 'w') as codeFid:
            content = self.catalgueMngCode
            content += "ctlgNoExt = '%s'\n" % catalogue
            content += "ctlgName = [ctlgNoExt, '.ctlg']\n"
            content += "ctlg = dread(ctlgName)\n"  # Load the catalogue
            content += "tomoFiles = cellfun(@(x) x.file, ctlg.volumes, 'UniformOutput', false)\n"  # Cell with the tomo names
            content += "currentTomoInd = find(ismember(tomoFiles, '%s'))\n" % abspath(tomo.getFileName())  # Index of the current tomo in the catalog
            content += "loadDynSyntax = sprintf('dtmslice @{%s}%i', ctlgNoExt, currentTomoInd)\n"
            content += "currentTomoModelsDir = fullfile('%s', ['volume_', num2str(currentTomoInd)], 'models')\n" % join(self.path, PROJECT_DIR, TOMOGRAMS_DIR)
            # Read models points and crop points before launching dynamo GUI
            content += "dirSt = dir(fullfile(currentTomoModelsDir, '*.omd'))\n"
            content += "prevModelList = cellfun(@(x) fullfile(currentTomoModelsDir, x), {dirSt.name}, 'un', 0)\n"
            content += "nModelsPrev = length(prevModelList)\n"
            content += "prevPointsStList = cell(1, nModelsPrev)\n"
            content += "for iModel=1:nModelsPrev\n"
            content += "iPrevModel = dread(prevModelList{iModel})\n"
            content += "prevPointsStList{iModel} = struct('points', iPrevModel.points, 'crop_points', iPrevModel.crop_points)\n"
            content += "end\n"
            # Launch Dynamo's picker GUI
            content += "eval(loadDynSyntax)\n"
            content += "modeltrack.loadFromCatalogue('handles', ctlg, 'full', true, 'select', false)\n"  # Load the models contained in the catalogue
            content += "uiwait(dpkslicer.getHandles().figure_fastslicer)\n"  # Wait until it's closed
            content += "modeltrack.saveAllInCatalogue\n"  # Save in the catalog
            # Read models points and crop points after having closed dynamo GUI
            content += "dirSt = dir(fullfile(currentTomoModelsDir, '*.omd'))\n"
            content += "postModelList = cellfun(@(x) fullfile(currentTomoModelsDir, x), {dirSt.name}, 'un', 0)\n"
            content += "nModelsPost = length(postModelList)\n"
            content += "postPointsStList = cell(1, nModelsPost)\n"
            content += "for iModel=1:nModelsPost\n"
            content += "iPostModel = dread(postModelList{iModel})\n"
            content += "postPointsStList{iModel} = struct('points', iPostModel.points, 'crop_points', iPostModel.crop_points)\n"
            content += "end\n"
            # Compare the points and cropped points stored from before and after running the viewer to check if the
            # they have changed
            content += "prjModified = false\n"
            content += "if nModelsPrev == nModelsPost\n"
            content += "for iModel=1:length(prevPointsStList)\n"
            content += "iPrevSt = prevPointsStList{iModel}\n"
            content += "iPostSt = postPointsStList{iModel}\n"
            content += "pickedPointsChanged = size(iPrevSt.points, 1) ~= size(iPostSt.points, 1)\n"
            content += "croppedPointsChanged = size(iPrevSt.crop_points, 1) ~= size(iPostSt.crop_points, 1)\n"
            content += "if pickedPointsChanged || croppedPointsChanged\n"
            content += "prjModified = true\n"
            content += "break\n"
            content += "end\n"
            content += "end\n"
            content += "else\n"
            content += "prjModified = true\n"
            content += "end\n"
            # If there was any data modification carried out using the viewer, the tree is updated consequently and
            # some text files are generated to indicate that the user operated that way
            content += "if prjModified == true\n"
            content += "models = dcmodels(ctlgNoExt, 'i', currentTomoInd)\n"
            content += "nParticles = 0\n"
            content += "for i=1:length(models)\n"
            content += "model = dread(models{i})\n"
            content += "newParts = size(model.points, 1)\n"
            content += "nParticles = nParticles + newParts\n"  # Sum the no. of particles from all the models generated for the current tomogram
            content += "end\n"
            # Save the number of particles to a text file
            content += "fid = fopen('%s', 'w')\n" % self.currentTomoTxtFile
            content += "fprintf(fid, '%i', nParticles)\n"
            content += "fclose(fid)\n"
            # Save a txt file to indicate that the project was modified using the viewer
            content += "fid = fopen('%s', 'w')\n" % join(self.path, PROJECT_DIR, DATA_MODIFIED_FROM_VIEWER)
            content += "fclose(fid)\n"
            content += "end\n"
            codeFid.write(content)
        return codeFilePath

    @staticmethod
    def _isADynamoProj(fpath):
        """Check the if a given path contains a valid Dynamo project:
            - Case 1: if the project was created by a previous use of the viewer, it will contain a file named
            prjFromViewer.txt. In that case, the project will be removed, as it may contain models corresponding to
            the coordinates of the valid model workflows instead of the failed meshes to be fixed, for example.
            - Case 2: Search recursively for Dynamo model files (omd) in a given directory (usually extra), to
            distinguish the usage of the viewer to visualize coordinates obtained with other plugins than Dynamo."""
        prjDir = join(fpath, PROJECT_DIR)
        prjFromViewer = join(prjDir, PRJ_FROM_VIEWER)
        if exists(prjFromViewer):
            rmtree(prjDir)
            return False
        modelFiles = getDynamoModels(fpath)
        return True if modelFiles else False
