# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *              Scipion Team (scipion@cnb.csic.es)
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
import glob
from os.path import isfile

import numpy as np

from dynamo.protocols.protocol_base_dynamo import DynamoOutputs
from dynamo.protocols.protocol_extraction import SAME_AS_PICKING, OTHER
from dynamo.tests.test_dynamo_base import TestDynamoStaBase
from pyworkflow.object import String
from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.utils import magentaStr
from tomo.objects import SetOfSubTomograms
from tomo.protocols import ProtImportSubTomograms, ProtImportCoordinates3D, ProtImportTomograms
from tomo.tests import DataSet, EMD_10439, DataSetEmd10439
import tomo.constants as const
from xmipp3.protocols import XmippProtCreateMask3D
from dynamo.protocols import DynamoSubTomoMRA, DynamoExtraction, DynamoImportSubtomos, DynamoBinTomograms, \
    DynamoCoordsToModel

DYNAMOPARAMNAME = "dyn"
MYPROJECT = "myproject"


class TestDynamoSubTomogramsAlignment(TestDynamoStaBase):
    """ This class check if the protocol to import sub tomograms works
    properly."""

    samplingRate = 5

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        dataset = DataSet.getDataSet('tomo-em')
        cls.setOfSubtomograms = dataset.getFile('basename.hdf')

    def test_rounds(self):
        staDynamo = DynamoSubTomoMRA()

        def testValueCase(value, expectedresult):
            scipionParam = String(value)

            result = staDynamo.getRoundParams(DYNAMOPARAMNAME, scipionParam, projectName=MYPROJECT)

            self.assertEquals(expectedresult,
                              result, "getRounds not working for this value: %s." % value)

        testValueCase("1", "dvput('%s', '%s', '1');\n" % (MYPROJECT, DYNAMOPARAMNAME))

        testValueCase("5 3", "dvput('%s', '%s', '5');\ndvput('%s', '%s_r2', '3');\n" % (
            MYPROJECT, DYNAMOPARAMNAME, MYPROJECT, DYNAMOPARAMNAME))

    def _runPreviousProtocols(self):
        protImport = self.newProtocol(ProtImportSubTomograms,
                                      filesPath=self.setOfSubtomograms,
                                      samplingRate=self.samplingRate)
        self.launchProtocol(protImport)
        return protImport

    def _runAlignment(self):
        protImport = self._runPreviousProtocols()
        particles = protImport.outputSubTomograms
        alignment = self.newProtocol(DynamoSubTomoMRA,
                                     inputVolumes=particles,
                                     generateTemplate=True,
                                     useGpu=False)
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
                                     templateRef=protImport,
                                     useGpu=False)
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
                                     fmask=particles,
                                     useGpu=False)
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
                                     nref=4,
                                     useGpu=False)
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
                                     mask=protMask.outputMask,
                                     useGpu=False)
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
                                     smask=protMask.outputMask,
                                     useGpu=False)
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
                                     fmask=protMask.outputMask,
                                     useGpu=False)
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


class TestDynamoBinTomograms(TestDynamoStaBase):
    tomoFiles = None
    binningFactor = 2
    origSize = np.array([1024, 1024, 512])
    sRate = 5
    bin2SRate = 10
    importedTomos = None

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        ds = DataSet.getDataSet('tomo-em')
        cls.tomoFiles = ds.getFile('tomo1')
        cls.importedTomos = super().runImportTomograms(tomoFiles=cls.tomoFiles, sRate=cls.sRate)

    def testBinTomograms(self):
        binnedTomograms = super().runBinTomograms(self.importedTomos, binning=self.binningFactor)
        # Check the results
        self.assertSetSize(binnedTomograms, 1)
        self.assertEqual(binnedTomograms.getSamplingRate(), self.bin2SRate)
        self.assertTrue(np.array_equal(np.array(binnedTomograms.getDimensions()), self.origSize / self.binningFactor))


class TestDynamoCoordsToModel(TestDynamoStaBase):
    """This class checks if the protocol to convert a SetOfCoordinates3D to a Dynamo
    model works properly"""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        dataset = DataSet.getDataSet('tomo-em')
        cls.tomogram = dataset.getFile('tomo1')
        cls.dynCoords = dataset.getFile('overview_wbp.tbl')

    def _runCoordsToModel(self):
        protImportTomogram = self.newProtocol(ProtImportTomograms,
                                              filesPath=self.tomogram,
                                              samplingRate=5)

        self.launchProtocol(protImportTomogram)

        self.assertIsNotNone(protImportTomogram.outputTomograms,
                             "There was a problem with tomogram output")

        protImportCoordinates = self.newProtocol(ProtImportCoordinates3D,
                                                 objLabel='Initial coordinates',
                                                 filesPath=self.dynCoords,
                                                 importTomograms=protImportTomogram.outputTomograms,
                                                 filesPattern='', boxSize=32, importAngles=True,
                                                 samplingRate=5)

        self.launchProtocol(protImportCoordinates)

        coords = protImportCoordinates.outputCoordinates

        self.assertIsNotNone(coords,
                             "There was a problem with coordinates 3d output")

        protDynamoCoordsToModel = self.newProtocol(DynamoCoordsToModel,
                                                   objLabel='Coords to Dynamo model conversion',
                                                   inputCoords=coords)

        self.launchProtocol(protDynamoCoordsToModel)
        return protDynamoCoordsToModel

    def test_coords_to_model(self):
        protCoordsToModel = self._runCoordsToModel()

        outputMeshes = getattr(protCoordsToModel, 'outputMeshes', None)
        self.assertTrue(outputMeshes)
        self.assertTrue(isfile(outputMeshes._dynCatalogue.get()))

        return protCoordsToModel
