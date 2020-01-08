# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from tomo.protocols import ProtImportSubTomograms
from tomo.tests import DataSet
from xmipp3.protocols import XmippProtCreateMask3D
from dynamo.protocols import DynamoSubTomoMRA


class TestSubTomogramsAlignment(BaseTest):
    """ This class check if the protocol to import sub tomograms works
    properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.setOfSubtomograms = cls.dataset.getFile('basename.hdf')

    def _runAlignment(self):
        protImport = self.newProtocol(ProtImportSubTomograms,
                                      filesPath=self.setOfSubtomograms,
                                      samplingRate=5)
        self.launchProtocol(protImport)

        particles = protImport.outputSubTomograms

        protMask = self.newProtocol(XmippProtCreateMask3D,
                                    inputVolume=protImport,
                                    source=0, volumeOperation=0,
                                    threshold=0.4)
        protMask.inputVolume.setExtended("outputSubTomograms.1")
        protMask.setObjLabel('threshold mask')
        self.launchProtocol(protMask)
        self.assertIsNotNone(protMask.outputMask,
                             "There was a problem with create mask from volume")
        import time
        alignment = self.newProtocol(DynamoSubTomoMRA,
                                     projName=time.time(),
                                     inputVolumes=particles,
                                     numberOfIters=1,
                                     templateRef=protImport,
                                     alignmentMask=protMask.outputMask)
        alignment.templateRef.setExtended("outputSubTomograms.1")
        self.launchProtocol(alignment)
        self.assertIsNotNone(alignment.outputSubtomograms,
                             "There was a problem with SetOfSubtomograms output")
        self.assertIsNotNone(alignment.outputClassesSubtomo,
                             "There was a problem with outputClassesSubtomo output")
        return alignment

    def test_outputMRA(self):
        dynamoAlignment = self._runAlignment()
        outputSubtomos = getattr(dynamoAlignment, 'outputSubtomograms')
        outputClasses = getattr(dynamoAlignment, 'outputClassesSubtomo')
        self.assertTrue(outputSubtomos)
        self.assertTrue(outputSubtomos.getFirstItem().hasTransform())
        self.assertTrue(outputClasses)
        # self.assertTrue(outputClasses.hasRepresentatives())
        return dynamoAlignment
