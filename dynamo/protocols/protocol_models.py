# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

import os

from pyworkflow.em import ProtAnalysis3D
from pyworkflow.protocol.params import PointerParam
from pyworkflow.utils import importFromPlugin, cleanPath
from pyworkflow.mapper.sqlite_db import SqliteDb

Mesh = importFromPlugin("tomo.objects", "Mesh")
SetOfMeshes = importFromPlugin("tomo.objects", "SetOfMeshes")


"""
Protocols to create models in Dynamo
"""

class DynamoModels(ProtAnalysis3D):
    """ It will align subtomograms using Dynamo"""
    _label = 'model manager'

    def __init__(self, **args):
        ProtAnalysis3D.__init__(self, **args)

        # --------------------------- DEFINE param functions ------------------------

    def _defineParams(self, form):
        form.addSection(label='Input Parameters')
        form.addParam('inputTomograms', PointerParam,
                      pointerClass="SetOfTomograms",
                      label='Set of tomograms',
                      help="Set of tomograms to create a model")


        # form.addParallelSection(threads=0, mpi=8)

        # --------------------------- INSERT steps functions --------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('writeInputFile')
        self._insertFunctionStep('modelStep')
        self._insertFunctionStep('createOutput')

        # --------------------------- STEPS functions -------------------------------

    def writeInputFile(self):
        self.inputFilePath = os.path.join(os.environ.get("SCIPION_HOME"), "software", "tmp", "commands.doc")
        tomoFile = os.path.join(os.environ.get("SCIPION_HOME"), "software", "tmp", "tomos.vll")
        catalogue = self._getExtraPath('tomos')
        tomoFid = open(tomoFile, 'w')
        for tomo in self.inputTomograms.get():
            tomoFile = tomo.getFileName()
            tomoName = os.path.basename(tomoFile)
            tomoFid.write(tomoName + '\n')
        tomoFid.close()
        inputFid = open(self.inputFilePath, 'w')
        content = 'dcm -create %s -fromvll %s \n' \
                  'dtmslice %s -c %s \n' \
                  'uiwait(msgbox(\'Click Ok when done\')) \n' \
                  'm = dread(dcmodels(\'%s\',\'i\',1)); \n' \
                  'writematrix(m.points, \'%s\'); \n' \
                  'exit' % (catalogue, tomoFile, tomoFile, catalogue, catalogue, self._getExtraPath('extra.txt'))
        inputFid.write(content)
        inputFid.close()


    def modelStep(self):
        program = '/home/davidh/dynamo_v1.146/matlab/bin/dynamo'
        args = ' %s' % self.inputFilePath
        dynamo = self.runJob(program, args)

    def createOutput(self):
        # TODO
        # Check of file corresponds to Tomo
        outSet = self._createSetOfMeshes()
        outFile = self._getExtraPath('extra.txt')
        roi = Mesh(outFile)
        for tomo in self.inputTomograms.get().iterItems():
                roi.setVolume(tomo)
        outSet.append(roi)
        outSet.setVolumes(self.inputTomograms.get())
        self._defineOutputs(outputMeshes=outSet)
        self._defineSourceRelation(self.inputTomograms.get(), outSet)

    def _createSetOfMeshes(self, suffix=''):
        return self.__createSet(SetOfMeshes,
                                'meshes%s.sqlite', suffix)

    def __createSet(self, SetClass, template, suffix, **kwargs):
        """ Create a set and set the filename using the suffix.
        If the file exists, it will be delete. """
        setFn = self._getPath(template % suffix)
        # Close the connection to the database if
        # it is open before deleting the file
        cleanPath(setFn)

        SqliteDb.closeConnection(setFn)
        setObj = SetClass(filename=setFn, **kwargs)
        return setObj


        # --------------------------- INFO functions --------------------------------

    def _summary(self):
        pass

    def _methods(self):
        pass

    def _citations(self):
        pass
