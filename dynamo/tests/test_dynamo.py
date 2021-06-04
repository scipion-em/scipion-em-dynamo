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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.tests import BaseTest, setupTestProject
from tomo.protocols import ProtImportSubTomograms, ProtImportCoordinates3D, ProtImportTomograms
from tomo.tests import DataSet
import tomo.constants as const
from xmipp3.protocols import XmippProtCreateMask3D
from dynamo.protocols import DynamoSubTomoMRA, DynamoExtraction, DynamoImportSubtomos


class TestDynamoBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @classmethod
    def setData(cls, projectData='tomo-em'):
        cls.dataset = DataSet.getDataSet(projectData)
        cls.tomogram = cls.dataset.getFile('tomo1')
        cls.emanCoords = cls.dataset.getFile('overview_wbp.txt')
        cls.dynCoords = cls.dataset.getFile('overview_wbp.tbl')
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

    def _runPreviousProtocols(self):
        protImport = self.newProtocol(ProtImportSubTomograms,
                                      filesPath=self.setOfSubtomograms,
                                      samplingRate=5)
        self.launchProtocol(protImport)
        return protImport

    def _runAlignment(self):
        protImport = self._runPreviousProtocols()
        particles = protImport.outputSubTomograms
        alignment = self.newProtocol(DynamoSubTomoMRA,
                                     inputVolumes=particles,
                                     generateTemplate=True)
        alignment.templateRef.setExtended("outputSubTomograms.1")
        self.launchProtocol(alignment)
        self.assertIsNotNone(alignment.averageSubTomogram,
                             "There was a problem with SetOfSubtomograms output")
        return alignment

    def _runAlignmentWithTemplate(self):
        protImport = self._runPreviousProtocols()
        particles = protImport.outputSubTomograms
        alignment = self.newProtocol(DynamoSubTomoMRA,
                                     inputVolumes=particles,
                                     templateRef=protImport)
        alignment.templateRef.setExtended("outputSubTomograms.1")
        self.launchProtocol(alignment)
        self.assertIsNotNone(alignment.averageSubTomogram,
                             "There was a problem with SetOfSubtomograms output")
        return alignment

    def _runAlignmentWithTemplatesMRA(self):
        protImport = self._runPreviousProtocols()
        particles = protImport.outputSubTomograms
        alignment = self.newProtocol(DynamoSubTomoMRA,
                                     inputVolumes=particles,
                                     templateRef=particles,
                                     mra=True,
                                     fmask=particles)
        alignment.fmask.setExtended("outputSubTomograms.1")
        self.launchProtocol(alignment)
        self.assertIsNotNone(alignment.AverageRef1,
                             "There was a problem with SetOfSubtomograms output")
        return alignment

    def _runAlignmentWithoutTemplatesMRA(self):
        protImport = self._runPreviousProtocols()
        particles = protImport.outputSubTomograms
        alignment = self.newProtocol(DynamoSubTomoMRA,
                                     inputVolumes=particles,
                                     mra=True,
                                     generateTemplate=True,
                                     nref=4)
        alignment.fmask.setExtended("outputSubTomograms.1")
        self.launchProtocol(alignment)
        self.assertIsNotNone(alignment.AverageRef1,
                             "There was a problem with SetOfSubtomograms output")
        return alignment

    def _runAlignmentWithMask(self):
        protImport = self._runPreviousProtocols()
        particles = protImport.outputSubTomograms
        protMask = self.newProtocol(XmippProtCreateMask3D,
                                    inputVolume=protImport,
                                    source=0,
                                    volumeOperation=0,
                                    threshold=0.4)
        protMask.inputVolume.setExtended("outputSubTomograms.1")
        protMask.setObjLabel('threshold mask')
        self.launchProtocol(protMask)
        self.assertIsNotNone(protMask.outputMask,
                             "There was a problem with create mask from volume")
        alignment = self.newProtocol(DynamoSubTomoMRA,
                                     inputVolumes=particles,
                                     templateRef=protImport,
                                     mask=protMask.outputMask)
        alignment.templateRef.setExtended("outputSubTomograms.1")
        self.launchProtocol(alignment)
        self.assertIsNotNone(alignment.averageSubTomogram,
                             "There was a problem with SetOfSubtomograms output")
        return alignment

    def _runAlignmentWithMasks(self):
        protImport = self._runPreviousProtocols()
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

        alignment = self.newProtocol(DynamoSubTomoMRA,
                                     inputVolumes=particles,
                                     templateRef=protImport,
                                     mask=protMask.outputMask,
                                     cmask=protMask.outputMask,
                                     fmask=protMask.outputMask,
                                     smask=protMask.outputMask)
        alignment.templateRef.setExtended("outputSubTomograms.1")
        self.launchProtocol(alignment)
        self.assertIsNotNone(alignment.averageSubTomogram,
                             "There was a problem with SetOfSubtomograms output")
        return alignment

    def _runAlignmentWithFmask(self):
        protImport = self._runPreviousProtocols()
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

        alignment = self.newProtocol(DynamoSubTomoMRA,
                                     inputVolumes=particles,
                                     generateTemplate=True,
                                     fmask=protMask.outputMask)
        alignment.templateRef.setExtended("outputSubTomograms.1")
        self.launchProtocol(alignment)
        self.assertIsNotNone(alignment.averageSubTomogram,
                             "There was a problem with SetOfSubtomograms output")
        return alignment

    def test_basicAlignment(self):
        dynamoAlignment = self._runAlignment()
        averageSubTomogram = getattr(dynamoAlignment, 'averageSubTomogram')
        outputSubtomograms = getattr(dynamoAlignment, 'outputSubtomograms')
        self.assertTrue(averageSubTomogram)
        self.assertTrue(outputSubtomograms)
        self.assertTrue(outputSubtomograms.getFirstItem().hasTransform())
        return dynamoAlignment

    def test_alignmentWithTemplate(self):
        dynamoAlignment = self._runAlignmentWithTemplate()
        averageSubTomogram = getattr(dynamoAlignment, 'averageSubTomogram')
        outputSubtomograms = getattr(dynamoAlignment, 'outputSubtomograms')
        self.assertTrue(averageSubTomogram)
        self.assertTrue(outputSubtomograms)
        self.assertTrue(outputSubtomograms.getFirstItem().hasTransform())
        return dynamoAlignment

    def test_alignmentWithMask(self):
        dynamoAlignment = self._runAlignmentWithMask()
        averageSubTomogram = getattr(dynamoAlignment, 'averageSubTomogram')
        outputSubtomograms = getattr(dynamoAlignment, 'outputSubtomograms')
        self.assertTrue(averageSubTomogram)
        self.assertTrue(outputSubtomograms)
        self.assertTrue(outputSubtomograms.getFirstItem().hasTransform())
        return dynamoAlignment

    def test_alignmentWithFask(self):
        dynamoAlignment = self._runAlignmentWithFmask()
        averageSubTomogram = getattr(dynamoAlignment, 'averageSubTomogram')
        outputSubtomograms = getattr(dynamoAlignment, 'outputSubtomograms')
        self.assertTrue(averageSubTomogram)
        self.assertTrue(outputSubtomograms)
        self.assertTrue(outputSubtomograms.getFirstItem().hasTransform())
        return dynamoAlignment

    def test_alignmentWithMasks(self):
        dynamoAlignment = self._runAlignmentWithMasks()
        averageSubTomogram = getattr(dynamoAlignment, 'averageSubTomogram')
        outputSubtomograms = getattr(dynamoAlignment, 'outputSubtomograms')
        self.assertTrue(averageSubTomogram)
        self.assertTrue(outputSubtomograms)
        self.assertTrue(outputSubtomograms.getFirstItem().hasTransform())
        return dynamoAlignment

    def test_alignmentWithTemplatesMRA(self):
        dynamoAlignment = self._runAlignmentWithTemplatesMRA()
        averageSubTomogram = getattr(dynamoAlignment, 'AverageRef1')
        outputSubtomograms = getattr(dynamoAlignment, 'outputSubtomogramsRef1')
        self.assertTrue(averageSubTomogram)
        self.assertTrue(outputSubtomograms)
        self.assertTrue(outputSubtomograms.getFirstItem().hasTransform())
        return dynamoAlignment

    def test_alignmentWithoutTemplatesMRA(self):
        dynamoAlignment = self._runAlignmentWithoutTemplatesMRA()
        averageSubTomogram = getattr(dynamoAlignment, 'AverageRef1')
        outputSubtomograms = getattr(dynamoAlignment, 'outputSubtomogramsRef1')
        self.assertTrue(averageSubTomogram)
        self.assertTrue(outputSubtomograms)
        self.assertTrue(outputSubtomograms.getFirstItem().hasTransform())
        return dynamoAlignment


class TestDynamoExtraction(TestDynamoBase):
    '''This class checks if the protocol to extract subtomograms
    works properly'''

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestDynamoBase.setData()

    def _runExtraction(self, tomoSource = 0, downFactor = 1, text=''):

        protImportTomogram = self.newProtocol(ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)

        self.launchProtocol(protImportTomogram)

        self.assertIsNotNone(protImportTomogram.outputTomograms,
                             "There was a problem with tomogram output")

        protImportCoordinatesNoAngles = self.newProtocol(ProtImportCoordinates3D,
                                                         objLabel='Coordinates without Angles',
                                                         auto=ProtImportCoordinates3D.IMPORT_FROM_EMAN,
                                                         filesPath=self.emanCoords,
                                                         importTomograms=protImportTomogram.outputTomograms,
                                                         filesPattern='', boxSize=32,
                                                         samplingRate=5)

        protImportCoordinatesAngles = self.newProtocol(ProtImportCoordinates3D,
                                                       objLabel='Coordinates with Angles',
                                                       auto=ProtImportCoordinates3D.IMPORT_FROM_EMAN,
                                                       filesPath=self.dynCoords,
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
                                                        objLabel='Extraction without Angles - %s' % text,
                                                        inputTomograms=protImportTomogram.outputTomograms,
                                                        inputCoordinates=coordsNoAngle,
                                                        boxSize=coordsNoAngle.getBoxSize(),
                                                        tomoSource=tomoSource,
                                                        downFactor=downFactor)

        protDynamoExtractionAngles = self.newProtocol(DynamoExtraction,
                                                      objLabel='Extraction with Angles - %s' % text,
                                                      inputTomograms=protImportTomogram.outputTomograms,
                                                      inputCoordinates=coordsAngle,
                                                      boxSize=coordsAngle.getBoxSize(),
                                                      tomoSource=tomoSource,
                                                      downFactor=downFactor)

        self.launchProtocol(protDynamoExtractionNoAngles)
        self.launchProtocol(protDynamoExtractionAngles)
        return protDynamoExtractionNoAngles, protDynamoExtractionAngles

    def test_Extraction_SameAsPicking(self):
        protExtraction = self._runExtraction(text='SameAsPicking')

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
        protExtraction = self._runExtraction(tomoSource = 1, text='OtherSource')

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
        protExtraction = self._runExtraction(downFactor = 2, text='DownSampling')

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
        protExtraction = self._runExtraction(tomoSource = 1, downFactor = 2, text='AllOptions')

        outputNoAngles = getattr(protExtraction[0], 'outputSetOfSubtomogram', None)
        self.assertTrue(outputNoAngles)
        self.assertTrue(outputNoAngles.hasCoordinates3D())
        self.assertTrue(outputNoAngles.getCoordinates3D().getObjValue())

        outputAngles = getattr(protExtraction[1], 'outputSetOfSubtomogram', None)
        self.assertTrue(outputAngles)
        self.assertTrue(outputAngles.hasCoordinates3D())
        self.assertTrue(outputAngles.getCoordinates3D().getObjValue())

        return protExtraction


class TestDynImportSubTomograms(BaseTest):
    """ This class check if the protocol to import subtomograms from Dynamo works
     properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.table = cls.dataset.getFile('initial.tbl')
        cls.path = cls.dataset.getPath()
        cls.subtomos = cls.dataset.getFile('basename.hdf')
        cls.tomo = cls.dataset.getFile('tomo1')

    def _runImportDynSubTomograms(self, tomos):
        protImport = self.newProtocol(DynamoImportSubtomos,
                                      filesPath=self.subtomos,
                                      samplingRate=1.35,
                                      tablePath=self.table,
                                      tomoSet=tomos)
        self.launchProtocol(protImport)
        return protImport

    def _runImportTomograms(self):
        protImportTomogram = self.newProtocol(ProtImportTomograms,
                                              filesPath=self.tomo,
                                              samplingRate=1.35)
        self.launchProtocol(protImportTomogram)
        return protImportTomogram

    def test_import_dynamo_subtomograms(self):
        protImportTomogram = self._runImportTomograms()
        protImport = self._runImportDynSubTomograms(protImportTomogram.outputTomograms)
        output = getattr(protImport, 'outputSubTomograms', None)
        self.assertTrue(output.getSamplingRate() == 1.35)
        self.assertTrue(output.getFirstItem().getSamplingRate() == 1.35)
        self.assertTrue(output.getDim()[0] == 32)
        self.assertTrue(output.getDim()[1] == 32)
        self.assertTrue(output.getDim()[2] == 32)
        self.assertTrue(output.getFirstItem().getObjId() == 4)
        self.assertTrue(output.getFirstItem().getClassId() == 1)
        self.assertTrue(output.getFirstItem().getAcquisition().getAngleMin() == -60)
        self.assertTrue(output.getFirstItem().getAcquisition().getAngleMax() == 60)
        self.assertAlmostEqual(output.getFirstItem().getCoordinate3D().getX(const.BOTTOM_LEFT_CORNER), 175, delta=1)
        self.assertAlmostEqual(output.getFirstItem().getCoordinate3D().getY(const.BOTTOM_LEFT_CORNER), 134, delta=1)
        self.assertAlmostEqual(output.getFirstItem().getCoordinate3D().getZ(const.BOTTOM_LEFT_CORNER), 115, delta=1)
