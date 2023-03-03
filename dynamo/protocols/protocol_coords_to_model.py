# **************************************************************************
# *
# * Authors:    David Herreros Calero (dherreros@cnb.csic.es)
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
import numpy as np

from pyworkflow import BETA
from pyworkflow.protocol import params
import pyworkflow.utils as pwutils
import pyworkflow.object as pwobj

from pwem.protocols import EMProtocol

from tomo.protocols import ProtTomoBase
import tomo.constants as const
from tomo.objects import SetOfMeshes, MeshPoint

from dynamo import Plugin, VLL_FILE, CATALOG_FILENAME


class DynamoCoordsToModel(EMProtocol, ProtTomoBase):
    """Convert a SetOfCoordinates3D to a SetOfMeshes formatted to be appropriate to work with Dynamo
    protocols"""

    _label = 'coords to model'
    _devStatus = BETA

    OUTPUT_PREFIX = 'outputMeshes'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

   # --------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputCoords', params.PointerParam, pointerClass='SetOfCoordinates3D',
                      label="Input Coordinates", important=True,
                      help="Coordinates to be converted.")

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('modelFromCoordsStep')
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions -----------------------------
    def modelFromCoordsStep(self):
        inputCoords = self.inputCoords.get()
        tomos = inputCoords.getPrecedents()
        matlab_file = self.writeMatlabFile(tomos, inputCoords)
        args = ' %s' % matlab_file
        self.runJob(Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    def createOutputStep(self):
        inputCoords = self.inputCoords.get()
        tomos = inputCoords.getPrecedents()
        suffix = self._getOutputSuffix(SetOfMeshes)
        out_meshes = self._createSetOfMeshes(tomos, suffix)
        out_meshes.setName("meshesFromCoords")
        out_meshes.setSamplingRate(tomos.getSamplingRate())
        out_meshes.setBoxSize(inputCoords.getBoxSize())
        out_meshes._dynCatalogue = pwobj.String(os.path.join(self._getExtraPath(), CATALOG_FILENAME))

        for tomo in tomos.iterItems():
            outPoints = pwutils.join(self._getExtraPath(), pwutils.removeBaseExt(tomo.getFileName()) + '.txt')

            # Populate SetOfMeshes with Meshes
            points = np.loadtxt(outPoints, delimiter=' ')
            for idx in range(len(points)):
                mesh = MeshPoint()
                mesh.setVolume(tomo)
                mesh.setPosition(points[idx, 0], points[idx, 1], points[idx, 2], const.BOTTOM_LEFT_CORNER)
                mesh.setGroupId(points[idx, 3])
                out_meshes.append(mesh)

        name = self.OUTPUT_PREFIX + suffix
        args = {name: out_meshes}
        self._defineOutputs(**args)
        self._defineSourceRelation(inputCoords, out_meshes)
        self._updateOutputSet(name, out_meshes, state=out_meshes.STREAM_CLOSED)

    # --------------------------- DEFINE utils functions ----------------------
    def writeMatlabFile(self, tomos, coords):
        # Catalogue initialization files
        listTomosFile = self._getExtraPath(VLL_FILE)
        catalogue = os.path.abspath(self._getExtraPath("tomos"))
        tomoFid = open(listTomosFile, 'w')
        for tomo in tomos.iterItems():
            tomoPath = os.path.abspath(tomo.getFileName())
            tomoFid.write(tomoPath + '\n')
        tomoFid.close()

        # Save coordinates into .txt file for each tomogram
        for tomo in tomos.iterItems(iterate=False):
            outFileCoord = self._getExtraPath(pwutils.removeBaseExt(tomo.getFileName())) + ".txt"
            coords_tomo = []
            for coord in coords.iterCoordinates(tomo.getObjId()):
                if coord.getGroupId() is not None:
                    coords_tomo.append(list(coord.getPosition(const.BOTTOM_LEFT_CORNER)) + [coord.getGroupId()])
                else:
                    coords_tomo.append(list(coord.getPosition(const.BOTTOM_LEFT_CORNER)) + [1])
            if coords_tomo:
                np.savetxt(outFileCoord, np.asarray(coords_tomo), delimiter=' ')

        # MatLab program for coords to model conversion
        codeFile = self._getExtraPath('coords2model.m')
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
                   "coords=coords_ids(coords_ids(:,4)==idm,1:3)\n" \
                   "general=dmodels.general()\n" \
                   "general.name=model_name\n" \
                   "addPoint(general,coords(:,1:3),coords(:,3))\n" \
                   "general.linkCatalogue('%s','i',tomoIndex,'s',1)\n" \
                   "general.saveInCatalogue()\n" \
                   "end\n" \
                   "end\n" \
                   "exit\n" % (catalogue, listTomosFile,
                               os.path.abspath(self._getExtraPath()), catalogue, catalogue)

        codeFid = open(codeFile, 'w')
        codeFid.write(contents)
        codeFid.close()

        return codeFile

    # --------------------------- DEFINE info functions ----------------------
    def _methods(self):
        methods = []
        methods.append("Conversion of coordinates to a Dynamo General Model")
        return methods

    def _summary(self):
        summary = []
        if self.getOutputsSize() >= 1:
            summary.append("A total of %d coordinates where converted to "
                           "Dynamo format." % self.inputCoords.get().getSize())
        else:
            summary.append("Conversion of coordinates in progress")
        return summary