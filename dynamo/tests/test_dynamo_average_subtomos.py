# ding=utf-8
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
from xmipp3.constants import MASK3D_CYLINDER
from xmipp3.protocols import XmippProtCreateMask3D
from xmipp3.protocols.protocol_preprocess.protocol_create_mask3d import SOURCE_GEOMETRY

from dynamo.protocols import DynamoProtAvgSubtomograms
from dynamo.protocols.protocol_extraction import SAME_AS_PICKING, DynamoExtraction, DynExtractionOutputs
from imod.protocols import ProtImodTomoNormalization
from imod.protocols.protocol_base import OUTPUT_TOMOGRAMS_NAME
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.utils import magentaStr
from tomo.protocols import ProtImportTomograms, ProtImportCoordinates3DFromScipion
from tomo.protocols.protocol_import_coordinates_from_scipion import outputObjs
from tomo.tests import EMD_10439, DataSetEmd10439


class TestDynamoAverageSubtomograms(BaseTest):

    ds = None
    boxSize = None
    sRate = None
    mask = None
    extractedSubtomos = None

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.ds = DataSet.getDataSet(EMD_10439)
        cls.sRate = 13.68
        cls.Bin2SRate = 2 * cls.sRate
        cls.boxSize = 44
        cls.runPreviousProtocols()

    @classmethod
    def runPreviousProtocols(cls):
        cls.mask = cls.runCreate3dMask()
        importedTomos = cls.runImportTomograms()
        binned2Tomos = cls.runBinTomograms(importedTomos)
        importedCoords = cls.runImport3dCoords(binned2Tomos)
        cls.extractedSubtomos = cls.runExtractSubtomograms(importedCoords, boxSize=cls.boxSize)

    @classmethod
    def runCreate3dMask(cls):
        print(magentaStr("\n==> Generating the 3D mask:"))
        protCreateParticleMask = cls.newProtocol(XmippProtCreateMask3D,
                                                 source=SOURCE_GEOMETRY,
                                                 samplingRate=cls.sRate,
                                                 size=cls.boxSize,
                                                 geo=MASK3D_CYLINDER,
                                                 radius=6,
                                                 shiftCenter=True,
                                                 centerZ=3,
                                                 height=15,
                                                 doSmooth=True)
        cls.launchProtocol(protCreateParticleMask)
        genMask = getattr(protCreateParticleMask, 'outputMask', None)
        cls.assertIsNotNone(genMask, 'the 3D mask was not generated')
        return genMask

    @classmethod
    def runImportTomograms(cls):
        # Import tomograms
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=cls.ds.getFile(DataSetEmd10439.tomoEmd10439.name),
                                             samplingRate=cls.sRate)

        cls.launchProtocol(protImportTomogram)
        tomoImported = protImportTomogram.outputTomograms
        cls.assertIsNotNone(tomoImported, "There was a problem with tomogram output")
        return tomoImported

    @classmethod
    def runImport3dCoords(cls, tomoImported):
        # Import coordinates
        print(magentaStr("\n==> Importing the 3D coordinates:"))
        protImportCoordinates3d = cls.newProtocol(ProtImportCoordinates3DFromScipion,
                                                  sqliteFile=cls.ds.getFile(DataSetEmd10439.coords39Sqlite.name),
                                                  importTomograms=tomoImported,
                                                  boxSize=cls.boxSize)

        cls.launchProtocol(protImportCoordinates3d)
        coordsImported = getattr(protImportCoordinates3d, outputObjs.coordinates.name, None)
        cls.assertIsNotNone(coordsImported, "There was a problem with the 3D coordinates output")
        return coordsImported

    @classmethod
    def runBinTomograms(cls, tomoImported, binning=2):
        # Bin the tomogram to make it smaller
        print(magentaStr("\n==> Tomogram binning:"))
        protTomoNormalization = cls.newProtocol(ProtImodTomoNormalization,
                                                inputSetOfTomograms=tomoImported,
                                                binning=binning)

        cls.launchProtocol(protTomoNormalization)
        tomosBinned = getattr(protTomoNormalization, OUTPUT_TOMOGRAMS_NAME, None)
        cls.assertIsNotNone(tomosBinned, 'No tomograms were genetated in tomo normalization.')
        return tomosBinned

    @classmethod
    def runExtractSubtomograms(cls, coordsImported, boxSize=None):
        # Extract subtomograms
        print(magentaStr("\n==> Extracting the subtomograms:"))
        protLabel = 'Extraction - same as picking'
        argsDict = {'inputCoordinates': coordsImported,
                    'tomoSource': SAME_AS_PICKING,
                    'boxSize': boxSize,
                    'doInvert': True,
                    'numberOfThreads': 4}
        protTomoExtraction = cls.newProtocol(DynamoExtraction, **argsDict)
        protTomoExtraction.setObjLabel(protLabel)
        cls.launchProtocol(protTomoExtraction)
        subtomosExtracted = getattr(protTomoExtraction, DynExtractionOutputs.subtomograms.name, None)
        cls.assertIsNotNone(subtomosExtracted, "There was a problem with the subtomograms extraction")
        return subtomosExtracted

    @classmethod
    def runAverageSubtomograms(cls):
        print(magentaStr("\n==> Averaging the subtomograms:"))
        protAvgSubtomo = cls.newProtocol(DynamoProtAvgSubtomograms,
                                         inSubtomos=cls.extractedSubtomos)
        cls.launchProtocol(protAvgSubtomo)