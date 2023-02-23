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
import os
from dynamo.protocols.protocol_bin_tomograms import DynBinTomosOutputs
from dynamo.protocols.protocol_extraction import DynExtractionOutputs, SAME_AS_PICKING, OTHER
from pyworkflow.object import String
from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.utils import magentaStr
from tomo.protocols import ProtImportSubTomograms, ProtImportCoordinates3D, ProtImportTomograms
from tomo.tests import DataSet
import tomo.constants as const
from xmipp3.protocols import XmippProtCreateMask3D
from dynamo.protocols import DynamoSubTomoMRA, DynamoExtraction, DynamoImportSubtomos, DynamoCoordsToModel, \
    DynamoBinTomograms

DYNAMOPARAMNAME = "dyn"
MYPROJECT = "myproject"


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


class TestDynamoSubTomogramsAlignment(BaseTest):
    """ This class check if the protocol to import sub tomograms works
    properly."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('tomo-em')
        cls.setOfSubtomograms = cls.dataset.getFile('basename.hdf')

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
                                      samplingRate=5)
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


class TestDynamoExtraction(TestDynamoBase):
    """This class checks if the protocol to extract subtomograms
    works properly"""

    nParticles = 4
    sRate = 5
    boxSize = 32
    binned2Tomos = None
    orientedCoords = None
    nonOrientedCoords = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestDynamoBase.setData()
        cls.runPreviousProtocols()

    @classmethod
    def runPreviousProtocols(cls):
        importedTomos = cls._runImportTomograms()
        cls.binned2Tomos = cls._runBinTomograms(importedTomos, binning=1)  # In dynamo, the binning is 2**binning
        cls.nonOrientedCoords = cls._runImportCoords(super().emanCoords, importedTomos)
        cls.orientedCoords = cls._runImportCoords(super().dynCoords, importedTomos, orientedCoords=True)

    @classmethod
    def _runImportTomograms(cls):
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=super().tomogram,
                                             samplingRate=cls.sRate)
        cls.launchProtocol(protImportTomogram)
        tomograms = protImportTomogram.outputTomograms
        cls.assertIsNotNone(tomograms, "There was a problem with tomogram output")
        return tomograms

    @classmethod
    def _runImportCoords(cls, coordsFile, tomograms, orientedCoords=False):
        print(magentaStr("\n==> Importing the 3D coordinates:"))
        objLabel = 'Oriented coords' if orientedCoords else 'Non-oriented coords'
        protImportCoordinates = cls.newProtocol(ProtImportCoordinates3D,
                                                objLabel=objLabel,
                                                filesPath=coordsFile,
                                                importTomograms=tomograms,
                                                filesPattern='',
                                                boxSize=cls.boxSize,
                                                samplingRate=cls.sRate)
        cls.launchProtocol(protImportCoordinates)
        importedCoords = protImportCoordinates.outputCoordinates
        cls.assertIsNotNone(importedCoords, "There was a problem with coordinates 3d output")
        return importedCoords

    @classmethod
    def _runBinTomograms(cls, tomoImported, binning=2):
        # Bin the tomogram to make it smaller
        print(magentaStr("\n==> Tomogram binning:"))
        protTomoNormalization = cls.newProtocol(DynamoBinTomograms,
                                                inputTomos=tomoImported,
                                                binning=binning)

        cls.launchProtocol(protTomoNormalization)
        tomosBinned = getattr(protTomoNormalization, DynBinTomosOutputs.tomograms.name, None)
        cls.assertIsNotNone(tomosBinned, 'No tomograms were genetated in tomo normalization.')
        return tomosBinned

    @classmethod
    def _runExtraction(cls, inputCoordinates, tomoSource=SAME_AS_PICKING, tomograms=None, downFactor=1, text=''):
        print(magentaStr("\n==> Extracting the subtomograms:"))
        argsDict = {
            'objLabel': 'Extraction - %s' % text,
            'inputCoordinates': inputCoordinates,
            'tomoSource': tomoSource,
            'boxSize': cls.boxSize,
            'downFactor': downFactor
        }
        if tomoSource == OTHER:
            argsDict['inputTomograms'] = tomograms
        protDynamoExtraction = cls.newProtocol(DynamoExtraction, **argsDict)
        cls.launchProtocol(protDynamoExtraction)
        return protDynamoExtraction

    def checkExtraction(self, protExtraction, binFactor=1):
        boxSize = self.boxSize / binFactor
        subtomograms = getattr(protExtraction[0], DynExtractionOutputs.subtomograms.name, None)
        self.assertIsNotNone(subtomograms)
        self.assertSetSize(subtomograms, self.nParticles)
        self.assertEqual(subtomograms.getSamplingRate, self.sRate * binFactor)
        self.assertEqual(subtomograms.getDimensions(), (boxSize, boxSize, boxSize))
        self.assertTrue(subtomograms.hasCoordinates3D())
        self.assertTrue(subtomograms._coordsPointer.hasExtended())

    def test_Extraction_SameAsPicking(self):
        protExtraction = self._runExtraction(self.orientedCoords,
                                             tomoSource=OTHER,
                                             tomograms=self.binned2Tomos,
                                             text='Same as Picking')
        self.checkExtraction(protExtraction, binFactor=2)

    def test_Extraction_Other(self):
        protExtraction = self._runExtraction(self.orientedCoords,
                                             tomoSource=OTHER,
                                             tomograms=self.binned2Tomos,
                                             text='Other Source')
        self.checkExtraction(protExtraction, binFactor=2)

    def test_Extraction_DownSampling(self):
        protExtraction = self._runExtraction(self.nonOrientedCoords,
                                             downFactor=2,
                                             text='Downsampling')
        self.checkExtraction(protExtraction)

    def test_Extraction_All(self):
        protExtraction = self._runExtraction(self.nonOrientedCoords,
                                             tomoSource=OTHER,
                                             tomograms=self.binned2Tomos,
                                             downFactor=2,
                                             text='All Options')
        self.checkExtraction(protExtraction, binFactor=2)


class TestDynamoImportSubTomograms(BaseTest):
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
        self.assertAlmostEqual(output.getFirstItem().getCoordinate3D().getX(const.SCIPION), -336.89, delta=1)
        self.assertAlmostEqual(output.getFirstItem().getCoordinate3D().getY(const.SCIPION), -377.65, delta=1)
        self.assertAlmostEqual(output.getFirstItem().getCoordinate3D().getZ(const.SCIPION), -140.26, delta=1)


class TestDynamoCoordsToModel(TestDynamoBase):
    '''This class checks if the protocol to convert a SetOfCoordinates3D to a Dynamo
    model works properly'''

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        TestDynamoBase.setData()

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
        self.assertTrue(os.path.isfile(outputMeshes._dynCatalogue.get()))

        return protCoordsToModel
