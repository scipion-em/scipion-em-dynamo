# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Scipion Team  (scipion@cnb.csic.es) [1]
# *
# * [1] Centro Nacional de Biotecnología (CSIC), Madrid, Spain
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
# from os.path import join
# from dynamo.tests import DYNAMO_TEST_DATASET, DataSetDynamo
# from pyworkflow.tests import BaseTest, setupTestProject, DataSet
# from pyworkflow.utils import magentaStr, makePath
#
#
# class TestDynamoAverageSubtomograms(BaseTest):
#
#     dsDynamo = None
#     dsReTomo = None
#
#     @classmethod
#     def setUpClass(cls):
#         setupTestProject(cls)
#         cls.dsDynamo = DataSet.getDataSet(DYNAMO_TEST_DATASET)
#         cls.dsReTomo = DataSet.getDataSet('reliontomo')
#         # Create a Dynamo project in /tmp with the Dynamo model files as expected by the sqlite that will
#         # be used to import the coordinates from a boxing protocol
#
#         cls.runPreviousProtocols()
#
#     @classmethod
#     def createDynProjectinTmp(cls):
#         tmpDir = '/tmp/dynamoTestProject/tomograms'
#         dirs2Create = [tmpDir,
#                        join(tmpDir, cls.dsDynamo.getFile(DataSetDynamo.vol1ModelsDir.name)),
#                        join(tmpDir, cls.dsDynamo.getFile(DataSetDynamo.vol2ModelsDir.name))]
#         makePath()
#
#
#     # @classmethod
#     # def runPreviousProtocols(cls):
#     #     importedTomos = cls.runImportTomograms()
#     #     binned2Tomos = cls.runBinTomograms(importedTomos)
#     #     importedCoords = cls.runImport3dCoords(binned2Tomos)
#     #     cls.extractedSubtomos = cls.runExtractSubtomograms(importedCoords, boxSize=cls.boxSize)
#     #
#     # @classmethod
#     # def runImportTomograms(cls):
#     #     # Import tomograms
#     #     print(magentaStr("\n==> Importing the tomograms:"))
#     #     protImportTomogram = cls.newProtocol(ProtImportTomograms,
#     #                                          filesPath=cls.ds.getFile(DataSetEmd10439.tomoEmd10439.name),
#     #                                          samplingRate=cls.sRate)
#     #
#     #     cls.launchProtocol(protImportTomogram)
#     #     tomoImported = protImportTomogram.outputTomograms
#     #     cls.assertIsNotNone(tomoImported, "There was a problem with tomogram output")
#     #     return tomoImported
#     #
#     # @classmethod
#     # def runImport3dCoords(cls, tomoImported):
#     #     # Import coordinates
#     #     print(magentaStr("\n==> Importing the 3D coordinates:"))
#     #     protImportCoordinates3d = cls.newProtocol(ProtImportCoordinates3DFromScipion,
#     #                                               sqliteFile=cls.ds.getFile(DataSetEmd10439.coords39Sqlite.name),
#     #                                               importTomograms=tomoImported,
#     #                                               boxSize=cls.boxSize)
#     #
#     #     cls.launchProtocol(protImportCoordinates3d)
#     #     coordsImported = getattr(protImportCoordinates3d, outputObjs.coordinates.name, None)
#     #     cls.assertIsNotNone(coordsImported, "There was a problem with the 3D coordinates output")
#     #     return coordsImported
#     #
#     # @classmethod
#     # def runBinTomograms(cls, tomoImported, binning=2):
#     #     # Bin the tomogram to make it smaller
#     #     print(magentaStr("\n==> Tomogram binning:"))
#     #     protTomoNormalization = cls.newProtocol(ProtImodTomoNormalization,
#     #                                             inputSetOfTomograms=tomoImported,
#     #                                             binning=binning)
#     #
#     #     cls.launchProtocol(protTomoNormalization)
#     #     tomosBinned = getattr(protTomoNormalization, OUTPUT_TOMOGRAMS_NAME, None)
#     #     cls.assertIsNotNone(tomosBinned, 'No tomograms were genetated in tomo normalization.')
#     #     return tomosBinned
#     #
#     # @classmethod
#     # def runExtractSubtomograms(cls, coordsImported, boxSize=None):
#     #     # Extract subtomograms
#     #     print(magentaStr("\n==> Extracting the subtomograms:"))
#     #     protLabel = 'Extraction - same as picking'
#     #     argsDict = {'inputCoordinates': coordsImported,
#     #                 'tomoSource': SAME_AS_PICKING,
#     #                 'boxSize': boxSize,
#     #                 'doInvert': True,
#     #                 'numberOfThreads': 4}
#     #     protTomoExtraction = cls.newProtocol(DynamoExtraction, **argsDict)
#     #     protTomoExtraction.setObjLabel(protLabel)
#     #     cls.launchProtocol(protTomoExtraction)
#     #     subtomosExtracted = getattr(protTomoExtraction, DynExtractionOutputs.subtomograms.name, None)
#     #     cls.assertIsNotNone(subtomosExtracted, "There was a problem with the subtomograms extraction")
#     #     return subtomosExtracted
#     #
#     # @classmethod
#     # def runAverageSubtomograms(cls):
#     #     print(magentaStr("\n==> Averaging the subtomograms:"))
#     #     protAvgSubtomo = cls.newProtocol(DynamoProtAvgSubtomograms,
#     #                                      inSubtomos=cls.extractedSubtomos)
#     #     cls.launchProtocol(protAvgSubtomo)
#     #     return getattr(protAvgSubtomo, DynAvgOutputs.average.name, None)
#     #
#     # def test_average(self):
#     #     testBoxSize = (self.boxSize, self.boxSize, self.boxSize)
#     #     avg = self.runAverageSubtomograms()
#     #     self.assertIsNotNone(avg)
#     #     self.assertTrue(exists(avg.getFileName()), "Average %s does not exists" % avg.getFileName())
#     #     self.assertTrue(avg.getFileName().endswith(".mrc"))
#     #     # Dynamo average protocol doesn't generate halves
#     #     self.assertFalse(avg.hasHalfMaps())
#     #     # The imported coordinates correspond to a binned 2 tomogram
#     #     self.assertEqual(avg.getSamplingRate(), self.Bin2SRate)
#     #     self.assertEqual(avg.getDimensions(), testBoxSize)
#     #
