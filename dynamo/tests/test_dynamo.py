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

from pyworkflow.utils import importFromPlugin
from pyworkflow.tests import BaseTest, setupTestProject, DataSet

from tomo.protocols import ProtImportSubTomograms
from tomo.tests import DataSet

from eman2.constants import TOMO_NEEDED_MSG

from xmipp3.protocols import XmippProtCreateMask3D

from dynamo.protocols import DynamoSubTomoMRA, DynamoExtraction


class TestDynamoBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @classmethod
    def setData(cls, projectData='tomo-em'):
        DataSet = importFromPlugin("tomo.tests", "DataSet", errorMsg=TOMO_NEEDED_MSG)
        cls.dataset = DataSet.getDataSet(projectData)
        cls.tomogram = cls.dataset.getFile('tomo1')
        cls.coords3D = cls.dataset.getFile('overview_wbp.txt')
        cls.angles = cls.dataset.getFile('overview_wbp.ang')
        cls.inputSetOfSubTomogram = cls.dataset.getFile('subtomo')
        cls.smallTomogram = cls.dataset.getFile('coremask_normcorona.mrc')

class TestSubTomogramsAlignment(BaseTest):
    """ This class check if the protocol to import sub tomograms works
    properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.setOfSubtomograms = cls.dataset.getFile('basename.hdf')

    def _runAlignment(self, randomInitialization=True, numberOfReferences=2,
                      numberOfIters=3, angularSampling=15):
        protImport = self.newProtocol(ProtImportSubTomograms,
                                      filesPath=self.setOfSubtomograms,
                                      samplingRate=5)
        self.launchProtocol(protImport)
        subtomos = protImport.outputSubTomograms

        protMask = self.newProtocol(XmippProtCreateMask3D,
                                    inputVolume=protImport,
                                    source=0, volumeOperation=0,
                                    threshold=0.4)
        protMask.inputVolume.setExtended("outputSubTomograms.1")
        protMask.setObjLabel('threshold mask')
        self.launchProtocol(protMask)
        self.assertIsNotNone(protMask.outputMask,
                             "There was a problem with create mask from volume")

        alignment = self.newProtocol(DynamoSubTomoMRA,
                                      inputVolumes=protImport.outputSubTomograms,
                                      alignmentMask=protMask.outputMask)
        self.launchProtocol(alignment)
        # self.assertIsNotNone(alignment.outputSubtomograms,
        #                      "There was a problem with SetOfSubtomograms output")
        # self.assertIsNotNone(alignment.outputClassesSubtomo,
        #                      "There was a problem with outputClassesSubtomo output")
        return alignment

    def test_outputMRA(self):
        dynamoAlignment = self._runAlignment()
        # outputSubtomos = getattr(dynamoAlignment, 'outputSubtomograms')
        # outputClasses = getattr(dynamoAlignment, 'outputClassesSubtomo')
        # self.assertTrue(outputSubtomos)
        # self.assertTrue(outputSubtomos.getFirstItem().hasTransform())
        # self.assertTrue(outputClasses)
        # self.assertTrue(outputClasses.hasRepresentatives())
        return dynamoAlignment

class TestDynamoExtraction(TestDynamoBase):
    '''This class checks if the protocol to extract subtomograms
    works properly'''

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestDynamoBase.setData()

    def _runExtraction(self, tomoSource = 0, downFactor = 1):
        ProtImportTomograms = importFromPlugin("tomo.protocols", "ProtImportTomograms",
                                               errorMsg=TOMO_NEEDED_MSG)
        ProtImportCoordinates3D = importFromPlugin("tomo.protocols", "ProtImportCoordinates3D",
                                                   errorMsg=TOMO_NEEDED_MSG)

        protImportTomogram = self.newProtocol(ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)

        self.launchProtocol(protImportTomogram)

        self.assertIsNotNone(protImportTomogram.outputTomograms,
                             "There was a problem with tomogram output")

        protImportCoordinatesNoAngles = self.newProtocol(ProtImportCoordinates3D,
                                                         auto=ProtImportCoordinates3D.IMPORT_FROM_EMAN,
                                                         filesPath=self.coords3D,
                                                         importTomograms=protImportTomogram.outputTomograms,
                                                         filesPattern='', boxSize=32,
                                                         samplingRate=5)

        protImportCoordinatesAngles = self.newProtocol(ProtImportCoordinates3D,
                                                       auto=ProtImportCoordinates3D.IMPORT_FROM_EMAN,
                                                       filesPath=self.coords3D,
                                                       importTomograms=protImportTomogram.outputTomograms,
                                                       filesPattern='', boxSize=32, importAngles=True,
                                                       samplingRate=5)

        self.launchProtocol(protImportCoordinatesNoAngles)
        self.launchProtocol(protImportCoordinatesAngles)

        coordsNoAngle = protImportCoordinatesNoAngles.outputCoordinates
        coordsAngle = protImportCoordinatesAngles.outputCoordinates

        self.assertIsNotNone(coordsNoAngle,
                             "There was a problem with coordinates 3d output")
        self.assertIsNotNone(coordsAngle,
                             "There was a problem with coordinates 3d output")

        protDynamoExtractionNoAngles = self.newProtocol(DynamoExtraction,
                                                        inputTomograms=protImportTomogram.outputTomograms,
                                                        inputCoordinates=coordsNoAngle,
                                                        boxSize=coordsNoAngle.getBoxSize(),
                                                        tomoSource=tomoSource,
                                                        downFactor=downFactor)

        protDynamoExtractionAngles = self.newProtocol(DynamoExtraction,
                                                      inputTomograms=protImportTomogram.outputTomograms,                                                      inputCoordinates=coordsAngle,
                                                      boxSize=coordsAngle.getBoxSize(),
                                                      tomoSource=tomoSource,
                                                      downFactor=downFactor)

        self.launchProtocol(protDynamoExtractionNoAngles)
        self.launchProtocol(protDynamoExtractionAngles)
        return protDynamoExtractionNoAngles, protDynamoExtractionAngles

    def test_Extraction_SameAsPicking(self):
        protExtraction = self._runExtraction()

        outputNoAngles = getattr( protExtraction[0], 'outputSetOfSubtomogram', None)
        self.assertTrue(outputNoAngles)
        self.assertTrue(outputNoAngles.hasCoordinates3D())
        self.assertTrue(outputNoAngles.getCoordinates3D().getObjValue())

        outputAngles = getattr( protExtraction[1], 'outputSetOfSubtomogram', None)
        self.assertTrue(outputAngles)
        self.assertTrue(outputAngles.hasCoordinates3D())
        self.assertTrue(outputAngles.getCoordinates3D().getObjValue())

        return protExtraction

    def test_Extraction_Other(self):
        protExtraction = self._runExtraction(tomoSource = 1)

        outputNoAngles = getattr(protExtraction[0], 'outputSetOfSubtomogram', None)
        self.assertTrue(outputNoAngles)
        self.assertTrue(outputNoAngles.hasCoordinates3D())
        self.assertTrue(outputNoAngles.getCoordinates3D().getObjValue())

        outputAngles = getattr(protExtraction[1], 'outputSetOfSubtomogram', None)
        self.assertTrue(outputAngles)
        self.assertTrue(outputAngles.hasCoordinates3D())
        self.assertTrue(outputAngles.getCoordinates3D().getObjValue())

        return protExtraction

    def test_Extraction_DownSampling(self):
        protExtraction = self._runExtraction(downFactor = 2)

        outputNoAngles = getattr(protExtraction[0], 'outputSetOfSubtomogram', None)
        self.assertTrue(outputNoAngles)
        self.assertTrue(outputNoAngles.hasCoordinates3D())
        self.assertTrue(outputNoAngles.getCoordinates3D().getObjValue())

        outputAngles = getattr(protExtraction[1], 'outputSetOfSubtomogram', None)
        self.assertTrue(outputAngles)
        self.assertTrue(outputAngles.hasCoordinates3D())
        self.assertTrue(outputAngles.getCoordinates3D().getObjValue())

        return protExtraction

    def test_Extraction_All(self):
        protExtraction = self._runExtraction(tomoSource = 1, downFactor = 2)

        outputNoAngles = getattr(protExtraction[0], 'outputSetOfSubtomogram', None)
        self.assertTrue(outputNoAngles)
        self.assertTrue(outputNoAngles.hasCoordinates3D())
        self.assertTrue(outputNoAngles.getCoordinates3D().getObjValue())

        outputAngles = getattr(protExtraction[1], 'outputSetOfSubtomogram', None)
        self.assertTrue(outputAngles)
        self.assertTrue(outputAngles.hasCoordinates3D())
        self.assertTrue(outputAngles.getCoordinates3D().getObjValue())

        return protExtraction