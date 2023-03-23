# **************************************************************************
# *
# * Authors:    Scipion Team (scipion@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
from dynamo.protocols import DynamoBinTomograms, DynamoProtAvgSubtomograms
from dynamo.protocols.protocol_extraction import SAME_AS_PICKING, OTHER, DynamoExtraction
from pyworkflow.tests import setupTestProject
from pyworkflow.utils import magentaStr
from tomo.protocols import ProtImportTomograms, ProtImportCoordinates3DFromScipion
from tomo.protocols.protocol_import_coordinates_from_scipion import outputObjs
from tomo.tests.test_base_centralized_layer import TestBaseCentralizedLayer


class TestDynamoStaBase(TestBaseCentralizedLayer):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    @classmethod
    def runImportTomograms(cls, tomoFiles=None, sRate=None):
        # Import tomograms
        print(magentaStr("\n==> Importing the tomograms:"))
        protImportTomogram = cls.newProtocol(ProtImportTomograms,
                                             filesPath=tomoFiles,
                                             samplingRate=sRate)

        cls.launchProtocol(protImportTomogram)
        tomoImported = protImportTomogram.outputTomograms
        cls.assertIsNotNone(tomoImported, "There was a problem while importing the tomograms")
        return tomoImported

    @classmethod
    def runBinTomograms(cls, inTomos=None, binning=None):
        # Bin the tomogram to make it smaller
        print(magentaStr("\n==> Tomogram binning:"))
        protBinTomos = cls.newProtocol(DynamoBinTomograms,
                                       inputTomos=inTomos,
                                       binning=binning)
        cls.launchProtocol(protBinTomos)
        tomosBinned = getattr(protBinTomos, protBinTomos._possibleOutputs.tomograms.name, None)
        cls.assertIsNotNone(tomosBinned, 'No tomograms were binned.')
        return tomosBinned

    @classmethod
    def runImport3dCoordsSqlite(cls, sqliteFile=None, inTomos=None, boxSize=None):
        # Import coordinates
        print(magentaStr("\n==> Importing the 3D coordinates:"))
        protImportCoordinates3d = cls.newProtocol(ProtImportCoordinates3DFromScipion,
                                                  sqliteFile=sqliteFile,
                                                  importTomograms=inTomos,
                                                  boxSize=boxSize)
        cls.launchProtocol(protImportCoordinates3d)
        coordsImported = getattr(protImportCoordinates3d, outputObjs.coordinates.name, None)
        cls.assertIsNotNone(coordsImported, "There was a problem with the 3D coordinates output")
        return coordsImported

    @classmethod
    def runExtractSubtomograms(cls, inCoords=None, tomoSource=SAME_AS_PICKING, tomograms=None, boxSize=None,
                               returnProtocol=False):
        print(magentaStr("\n==> Extracting the subtomograms:"))
        protLabel = 'Extraction - same as picking'
        argsDict = {'inputCoordinates': inCoords,
                    'tomoSource': tomoSource,
                    'boxSize': boxSize,
                    'doInvert': True}
        if tomoSource != SAME_AS_PICKING:
            argsDict['tomoSource'] = OTHER
            argsDict['inputTomograms'] = tomograms
            protLabel = 'Extraction - another tomo source'

        protTomoExtraction = cls.newProtocol(DynamoExtraction, **argsDict)
        protTomoExtraction.setObjLabel(protLabel)
        cls.launchProtocol(protTomoExtraction)
        if returnProtocol:
            return protTomoExtraction
        else:
            subtomosExtracted = getattr(protTomoExtraction, protTomoExtraction._possibleOutputs.subtomograms.name, None)
            cls.assertIsNotNone(subtomosExtracted, "There was a problem with the subtomogram extraction")
            return subtomosExtracted

    @classmethod
    def runAverageSubtomograms(cls, inSubtomos):
        print(magentaStr("\n==> Averaging the subtomograms:"))
        protAvgSubtomo = cls.newProtocol(DynamoProtAvgSubtomograms, inSubtomos=inSubtomos)
        cls.launchProtocol(protAvgSubtomo)
        avg = getattr(protAvgSubtomo, protAvgSubtomo._possibleOutputs.average.name, None)
        cls.assertIsNotNone(avg, "There was a problem with the subtomogram averaging")
        return avg




