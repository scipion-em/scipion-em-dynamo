# **************************************************************************
# *
# * Authors:    David Herreros Calero (dherreros@cnb.csic.es)
# *             Scipion Team (scipion@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
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
from enum import Enum

from pyworkflow import BETA
from pyworkflow.protocol import params
import pyworkflow.utils as pwutils

from pwem.protocols import EMProtocol
from tomo.objects import SetOfMeshes

from tomo.protocols import ProtTomoBase

from ..convert.convert import textFile2Coords
from dynamo import Plugin

# Model types mapping
MODEL_CHOICES = ["Ellipsoidal Vesicle", "Surface", "General"]

# Model types encoding
M_ELLIPSOIDAL = 0
M_SURFACE = 1
M_GENERAL = 2


class OutputsModelWf(Enum):
    meshes = SetOfMeshes


class DynamoModelWorkflow(EMProtocol, ProtTomoBase):
    """Apply a model workflow to a SetOfMeshes generated by Dynamo Boxing protocol. This workflow will use the
    models created by the user to create the corresponding cropping meshes needed to extract the crop points"""

    _label = 'model workflow'
    _devStatus = BETA
    _possibleOutputs = OutputsModelWf
    OUTPUT_PREFIX = 'output3DCoordinates'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMeshes', params.PointerParam,
                      pointerClass='SetOfMeshes',
                      label="Input Meshes",
                      important=True,
                      help="Input Meshes that will be used to create the cropping geometry and "
                           "to extract the crop points")
        form.addParam('boxSize', params.IntParam,
                      default=0,
                      label='Box size',
                      important=True)
        form.addParam('modelType', params.EnumParam,
                      choices=MODEL_CHOICES,
                      default=M_ELLIPSOIDAL,
                      label='Model type',
                      help='Select the type of model defined in the Tomograms.')

        form.addSection('Model parameters')
        form.addLine('Available models:')
        form.addLine('%s [EV]' % MODEL_CHOICES[M_ELLIPSOIDAL])
        form.addLine('%s [S]' % MODEL_CHOICES[M_SURFACE])
        form.addLine('%s [G]' % MODEL_CHOICES[M_GENERAL])
        form.addLine('%s model specific cases:' % MODEL_CHOICES[M_GENERAL])
        form.addLine('If oriented with an ellipsoidal vesicle model mesh [GE]')
        form.addLine('If oriented with a surface model mesh [GS]')
        form.addLine('Parameters corresponding to models not present in your picking data can be ignored')

        group = form.addGroup('Specific for %s model' % MODEL_CHOICES[M_GENERAL])
        group.addParam('orientMesh', params.PointerParam,
                       pointerClass='SetOfMeshes',
                       label="Orientation Meshes",
                       allowsNull=True,
                       # condition="modelType == %i" % M_GENERAL,
                       help="Specify if you will use a different SetOfMeshes to impart orientation to "
                            "the *Input Meshes* provided before. If not provided, no directional information "
                            "will be given to the output coordinates.")
        group.addParam('orientType', params.EnumParam,
                       choices=MODEL_CHOICES[:-1],
                       default=M_ELLIPSOIDAL,
                       # condition="modelType == %i" % M_GENERAL,
                       label='Model type for *Orientation Meshes*',
                       help='Select the type of model defined in the Tomograms.')

        group = form.addGroup('Common to multiple models')
        group.addParam('meshParameter', params.IntParam,
                       default=5,
                       label='Mesh parameter [EV][S][GE][GS]',
                       # condition='(modelType==%i or modelType==%i) or (modelType==%i and orientType==%i or '
                       #           'orientType==%i)' % (M_ELLIPSOIDAL, M_SURFACE, M_GENERAL, M_ELLIPSOIDAL, M_SURFACE),
                       help='Intended mesh parameter for the "mesh" that supports the depiction of the model')
        group.addParam('maxTr [EV][S][GE][GS]', params.IntParam,
                       default=100000,
                       # condition='(modelType==%i or modelType==%i) or (modelType==%i and orientType==%i or '
                       #           'orientType==%i)' % (M_ELLIPSOIDAL, M_SURFACE, M_GENERAL, M_ELLIPSOIDAL, M_SURFACE),
                       label="Maximun number of triangles",
                       help='Maximum number of triangles allowed during generation of a depiction mesh')
        group.addParam('auto', params.BooleanParam,
                       default=True,
                       # condition="modelType == %i or (modelType==%i and orientType==%i)" %
                       #           (M_ELLIPSOIDAL, M_GENERAL, M_ELLIPSOIDAL),
                       label='Auto detect geometry [EV][GE]')
        group.addParam('center', params.NumericListParam,
                       default='0 0 0',
                       label='Center (pix.) [EV][GE] only if auto = No',
                       # condition='(modelType==%i and auto==False) or (modelType==%i and orientType==%i and auto==False)'
                       #           % (M_ELLIPSOIDAL, M_GENERAL, M_ELLIPSOIDAL),
                       help='Center of globular models, or a point marking the interior part of a membrane')
        group.addParam('radius', params.NumericListParam,
                       default='10 10 10',
                       label='radius XYZ (pix.) [EV][GE] only if auto = No',
                       # condition='(modelType==%i and auto==False) or (modelType==%i and orientType==%i and auto==False)' %
                       #           (M_ELLIPSOIDAL, M_GENERAL, M_ELLIPSOIDAL),
                       help='Semi-axes for ellipsoidal vesicles or general models oriented with ellipsoidal vesicles.')
        group.addParam('cropping', params.IntParam,
                       default=10,
                       # condition='(modelType==%i or modelType==%i)' % (M_ELLIPSOIDAL, M_SURFACE),
                       label="Cropping parameter [EV][S]",
                       help='Intended mesh parameter for the "crop_mesh" that defined a cropping '
                            'geometry on a surface')
        group.addParam('subDivision', params.IntParam,
                       default=2,
                       # condition='modelType==%i or (modelType==%i and orientType==%i)' %
                       #           (M_SURFACE, M_GENERAL, M_SURFACE),
                       label="Number of Subdivision steps [S][GS]",
                       help="Specifiy the number of times the Mesh geometry will be subdivided. This will increase the "
                            "number of triangles in the mesh, making it smoother. However, it will also increase the "
                            "number of cropping points")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.applyWorkflowStep)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions -----------------------------
    def applyWorkflowStep(self):
        inputMeshes = self.inputMeshes.get()
        catalogue_path = os.path.abspath(pwutils.removeExt(inputMeshes._dynCatalogue.get()))
        self.volIds = inputMeshes.aggregate(["MAX"], "_volId", ["_volId"])
        for d in self.volIds:
            volId = d['_volId']
            commandsFile = self.writeMatlabFile(catalogue_path, volId)
            args = ' %s' % commandsFile
            self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    def createOutputStep(self):
        textFile2Coords(self, self.inputMeshes.get().getPrecedents(), self._getExtraPath())

    # --------------------------- DEFINE utils functions ----------------------
    def writeMatlabFile(self, catalogue_path, volId):
        codeFilePath = self._getExtraPath('modelWorkflow_Tomo_%d.m' % volId)
        outPath = pwutils.removeBaseExt(self.inputMeshes.get().getPrecedents()[volId].getFileName())
        if self.modelType.get() == 0:
            center = pwutils.getFloatListFromValues(self.center.get())
            radius = pwutils.getFloatListFromValues(self.radius.get())
            content = "path='%s'\n" \
                      "auto=%i\n" \
                      "crop_points=[]\n" \
                      "crop_angles=[]\n" \
                      "modelId=0\n" \
                      "outFile='%s'\n" \
                      "savePath='%s'\n" \
                      "outPoints=[outFile '.txt']\n" \
                      "outAngles=['angles_' outFile '.txt']\n" \
                      "modelFile=dir(fullfile(path,'tomograms','volume_%d','models'))\n" \
                      "modelFile=modelFile(~ismember({modelFile.name},{'.','..'}))\n" \
                      "for k=1:length(modelFile)\n" \
                      "modelId=modelId+1\n" \
                      "m=dread(fullfile(modelFile(k).folder,modelFile(k).name))\n" \
                      "if auto==1\n" \
                      "m.approximateGeometryFromPoints()\n" \
                      "else\n" \
                      "m.center=[%f %f %f]\n" \
                      "m.radius_x=%f\n" \
                      "m.radius_y=%f\n" \
                      "m.radius_z=%f\n" \
                      "end\n" \
                      "m.mesh_parameter=%d\n" \
                      "m.crop_mesh_parameter=%d\n" \
                      "m.mesh_maximum_triangles=%d\n" \
                      "m.createMesh\n" \
                      "m.createCropMesh\n" \
                      "m.grepTable()\n" \
                      "crop_points=[crop_points; [m.crop_points modelId*ones(length(m.crop_points),1)]]\n" \
                      "crop_angles=[crop_angles; [m.crop_angles modelId*ones(length(m.crop_angles),1)]]\n" \
                      "end\n" \
                      "writematrix(crop_points,fullfile(savePath,outPoints),'Delimiter',' ')\n" \
                      "writematrix(crop_angles,fullfile(savePath,outAngles),'Delimiter',' ')\n" \
                      "exit\n" % (catalogue_path, self.auto.get(), outPath, os.path.abspath(self._getExtraPath()),
                                  volId, center[0], center[1], center[2], radius[0], radius[1],
                                  radius[2], self.meshParameter.get(), self.cropping.get(),
                                  self.maxTr.get())
        elif self.modelType.get() == 1:
            content = "path='%s'\n" \
                      "crop_points=[]\n" \
                      "crop_angles=[]\n" \
                      "modelId=0\n" \
                      "outFile='%s'\n" \
                      "savePath='%s'\n" \
                      "outPoints=[outFile '.txt']\n" \
                      "outAngles=['angles_' outFile '.txt']\n" \
                      "modelFile=dir(fullfile(path,'tomograms','volume_%d','models'))\n" \
                      "modelFile=modelFile(~ismember({modelFile.name},{'.','..'}))\n" \
                      "for k=1:length(modelFile)\n" \
                      "modelId=modelId+1\n" \
                      "m=dread(fullfile(modelFile(k).folder,modelFile(k).name))\n" \
                      "m.mesh_parameter=%d\n" \
                      "m.crop_mesh_parameter=%d\n" \
                      "m.mesh_maximum_triangles=%d\n" \
                      "m.subdivision_iterations=%d\n" \
                      "if isempty(m.center)\n" \
                      "m.center=mean(m.points)\n" \
                      "end\n" \
                      "m.createMesh()\n" \
                      "m.refineMesh()\n" \
                      "m.createCropMesh()\n" \
                      "m.updateCrop()\n" \
                      "m.grepTable()\n" \
                      "crop_points=[crop_points; [m.crop_points modelId*ones(length(m.crop_points),1)]]\n" \
                      "crop_angles=[crop_angles; [m.crop_angles modelId*ones(length(m.crop_angles),1)]]\n" \
                      "end\n" \
                      "writematrix(crop_points,fullfile(savePath,outPoints),'Delimiter',' ')\n" \
                      "writematrix(crop_angles,fullfile(savePath,outAngles),'Delimiter',' ')\n" \
                      "exit\n" % (catalogue_path, outPath, os.path.abspath(self._getExtraPath()),
                                  volId, self.meshParameter.get(), self.cropping.get(),
                                  self.maxTr.get(), self.subDivision.get())
        elif self.modelType.get() == 2:
            if self.orientType.get() == 0 and self.orientMesh.get() is not None:
                center = pwutils.getFloatListFromValues(self.center.get())
                radius = pwutils.getFloatListFromValues(self.radius.get())
                orientation_mesh = os.path.abspath(pwutils.removeExt(self.orientMesh.get()._dynCatalogue.get()))
                content = "path='%s'\n" \
                          "path_o='%s'\n" \
                          "auto=%i\n" \
                          "crop_points=[]\n" \
                          "crop_angles=[]\n" \
                          "modelId=0\n" \
                          "outFile='%s'\n" \
                          "savePath='%s'\n" \
                          "outPoints=[outFile '.txt']\n" \
                          "outAngles=['angles_' outFile '.txt']\n" \
                          "modelFile=dir(fullfile(path,'tomograms','volume_%d','models'))\n" \
                          "modelFile=modelFile_o(~ismember({modelFile.name},{'.','..'}))\n" \
                          "modelFile_o=dir(fullfile(path_o,'tomograms','volume_%d','models'))\n" \
                          "modelFile_o=modelFile(~ismember({modelFile_o.name},{'.','..'}))\n" \
                          "for k=1:length(modelFile)\n" \
                          "modelId=modelId+1\n" \
                          "m=dread(fullfile(modelFile(k).folder,modelFile(k).name))\n" \
                          "m_o=dread(fullfile(modelFile_o(k).folder,modelFile_o(k).name))\n" \
                          "if auto==1\n" \
                          "m_o.approximateGeometryFromPoints()\n" \
                          "else\n" \
                          "m_o.center=[%f %f %f]\n" \
                          "m_o.radius_x=%f\n" \
                          "m_o.radius_y=%f\n" \
                          "m_o.radius_z=%f\n" \
                          "end\n" \
                          "m_o.mesh_parameter=%d\n" \
                          "m_o.mesh_maximum_triangles=%d\n" \
                          "m_o.createMesh\n" \
                          "t=m.grepTable()\n" \
                          "t=dpktbl.triangulation.fillTable(m_o.mesh,t)\n" \
                          "crop_points=[crop_points; [t(:,24:26) modelId*ones(length(m.crop_points),1)]]\n" \
                          "crop_angles=[crop_angles; [t(:,7:9) modelId*ones(length(m.crop_angles),1)]]\n" \
                          "end\n" \
                          "writematrix(crop_points,fullfile(savePath,outPoints),'Delimiter',' ')\n" \
                          "writematrix(crop_angles,fullfile(savePath,outAngles),'Delimiter',' ')\n" \
                          "exit\n" % (catalogue_path, orientation_mesh, self.auto.get(), outPath,
                                      os.path.abspath(self._getExtraPath()),
                                      volId, volId, center[0], center[1], center[2], radius[0], radius[1],
                                      radius[2], self.meshParameter.get(),
                                      self.maxTr.get())
            elif self.orientType.get() == 1 and self.orientMesh.get() is not None:
                orientation_mesh = os.path.abspath(pwutils.removeExt(self.orientMesh.get()._dynCatalogue.get()))
                content = "path='%s'\n" \
                          "path_o='%s'\n" \
                          "crop_points=[]\n" \
                          "crop_angles=[]\n" \
                          "modelId=0\n" \
                          "outFile='%s'\n" \
                          "savePath='%s'\n" \
                          "outPoints=[outFile '.txt']\n" \
                          "outAngles=['angles_' outFile '.txt']\n" \
                          "modelFile=dir(fullfile(path,'tomograms','volume_%d','models'))\n" \
                          "modelFile=modelFile(~ismember({modelFile.name},{'.','..'}))\n" \
                          "modelFile_o=dir(fullfile(path_o,'tomograms','volume_%d','models'))\n" \
                          "modelFile_o=modelFile_o(~ismember({modelFile_o.name},{'.','..'}))\n" \
                          "for k=1:length(modelFile)\n" \
                          "modelId=modelId+1\n" \
                          "m=dread(fullfile(modelFile(k).folder,modelFile(k).name))\n" \
                          "m_o=dread(fullfile(modelFile_o(k).folder,modelFile_o(k).name))\n" \
                          "m_o.mesh_parameter=%d\n" \
                          "m_o.mesh_maximum_triangles=%d\n" \
                          "m_o.subdivision_iterations=%d\n" \
                          "if isempty(m_o.center)\n" \
                          "m_o.center=mean(m_o.points)\n" \
                          "end\n" \
                          "m_o.createMesh()\n" \
                          "m_o.refineMesh()\n" \
                          "t=m.grepTable()\n" \
                          "t=dpktbl.triangulation.fillTable(m_o.mesh,t)\n" \
                          "crop_points=[crop_points; [t(:,24:26) modelId*ones(length(m.crop_points),1)]]\n" \
                          "crop_angles=[crop_angles; [t(:,7:9) modelId*ones(length(m.crop_angles),1)]]\n" \
                          "end\n" \
                          "writematrix(crop_points,fullfile(savePath,outPoints),'Delimiter',' ')\n" \
                          "writematrix(crop_angles,fullfile(savePath,outAngles),'Delimiter',' ')\n" \
                          "exit\n" % (catalogue_path, orientation_mesh, outPath,
                                      os.path.abspath(self._getExtraPath()),
                                      volId, volId, self.meshParameter.get(),
                                      self.maxTr.get(), self.subDivision.get())
            else:
                content = "path='%s'\n" \
                          "crop_points=[]\n" \
                          "crop_angles=[]\n" \
                          "modelId=0\n" \
                          "outFile='%s'\n" \
                          "savePath='%s'\n" \
                          "outPoints=[outFile '.txt']\n" \
                          "outAngles=['angles_' outFile '.txt']\n" \
                          "modelFile=dir(fullfile(path,'tomograms','volume_%d','models'))\n" \
                          "modelFile=modelFile(~ismember({modelFile.name},{'.','..'}))\n" \
                          "for k=1:length(modelFile)\n" \
                          "modelId=modelId+1\n" \
                          "m=dread(fullfile(modelFile(k).folder,modelFile(k).name))\n" \
                          "t=m.grepTable()\n" \
                          "crop_points=[crop_points; [t(:,24:26) modelId*ones(length(m.crop_points),1)]]\n" \
                          "crop_angles=[crop_angles; [t(:,7:9) modelId*ones(length(m.crop_angles),1)]]\n" \
                          "end\n" \
                          "writematrix(crop_points,fullfile(savePath,outPoints),'Delimiter',' ')\n" \
                          "writematrix(crop_angles,fullfile(savePath,outAngles),'Delimiter',' ')\n" \
                          "exit\n" % (catalogue_path, outPath,
                                      os.path.abspath(self._getExtraPath()),
                                      volId)
        codeFid = open(codeFilePath, 'w')
        codeFid.write(content)
        codeFid.close()
        return codeFilePath

    # --------------------------- DEFINE info functions ----------------------
    def _methods(self):
        methodsMsgs = ["*Model Type*: %s" % MODEL_CHOICES[self.modelType.get()]]
        if self.modelType.get() == 0:
            if self.auto.get():
                methodsMsgs.append("*Geometry Detection*: auto")
            else:
                methodsMsgs.append("*Geometry Detection*: manual")
                methodsMsgs.append("    *Center*: [%s]" % self.center.get())
                methodsMsgs.append("    *Semiaxes Radius*: [%s]" % self.radius.get())
            methodsMsgs.append("*Mesh parameter*: %d" % self.meshParameter.get())
            methodsMsgs.append("*Maximum number triangles*: %d" % self.maxTr.get())
            methodsMsgs.append("*Cropping parameter*: %d" % self.cropping.get())
        elif self.modelType.get() == 1:
            methodsMsgs.append("*Model Type*: Surface")
            methodsMsgs.append("*Mesh parameter*: %d" % self.meshParameter.get())
            methodsMsgs.append("*Maximum number triangles*: %d" % self.maxTr.get())
            methodsMsgs.append("*Cropping parameter*: %d" % self.cropping.get())
            methodsMsgs.append("*Number of subdivision steps*: %d" % self.subDivision.get())
        elif self.modelType.get() == 2:
            if self.orientMesh.get() is None:
                methodsMsgs.append("Particles extracted without orientation")
            else:
                if self.orientType.get() == 0:
                    methodsMsgs.append("*Orientation Model Type*: Ellipsoidal Vesicle")
                    if self.auto.get():
                        methodsMsgs.append("*Geometry Detection*: auto")
                    else:
                        methodsMsgs.append("*Geometry Detection*: manual")
                        methodsMsgs.append("    *Center*: [%s]" % self.center.get())
                        methodsMsgs.append("    *Semiaxes Radius*: [%s]" % self.radius.get())
                    methodsMsgs.append("*Mesh parameter*: %d" % self.meshParameter.get())
                    methodsMsgs.append("*Maximum number triangles*: %d" % self.maxTr.get())
                elif self.orientType.get() == 1:
                    methodsMsgs.append("*Orientation Model Type*: Surface")
                    methodsMsgs.append("*Mesh parameter*: %d" % self.meshParameter.get())
                    methodsMsgs.append("*Maximum number triangles*: %d" % self.maxTr.get())
                    methodsMsgs.append("*Number of subdivision steps*: %d" % self.subDivision.get())
        return methodsMsgs

    def _summary(self):
        summary = []
        if self.getOutputsSize() >= 1:
            for _, outCoords in self.iterOutputAttributes():
                summary.append("Output *%s*:" % outCoords.getNameId().split('.')[1])
                summary.append("    * Particle box size: *%s*" % self.boxSize.get())
                summary.append("    * Coordinates defined by geometry: *%s*" %
                               outCoords.getSize())
        else:
            summary.append("Output coordinates not ready yet.")
        return summary

    def _validate(self):
        errorMsg = []
        # Only sets of meshes generated using the Dynamo picking protocol are accepted (they must contain
        # an attribute named '_dynCatalogue')
        if not getattr(self.inputMeshes.get(), '_dynCatalogue', None):
            errorMsg.append('Only sets of meshes generated using the Dynamo picking protocol are accepted')
        return errorMsg
