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
from enum import Enum
from os.path import basename, dirname, join, normpath
import xmipp3
from pyworkflow.object import Float
from pyworkflow.protocol.params import PointerParam, FloatParam
from pyworkflow.utils.path import createLink, makePath
from pwem.emlib.image import ImageHandler
from tomo.protocols.protocol_base import ProtTomoImportFiles
from tomo.objects import SubTomogram, SetOfSubTomograms, SetOfCoordinates3D
from .protocol_base_dynamo import DynamoProtocolBase
from ..convert import dynTableLine2Subtomo


class DynImportSubtomosOuts(Enum):
    coordinates = SetOfCoordinates3D
    subtomograms = SetOfSubTomograms


class DynamoImportSubtomos(ProtTomoImportFiles, DynamoProtocolBase):
    """ This protocol imports subtomograms with metadata generated from Dynamo tables.
    The subtomograms files are generated in the same directory as the .tbl file, one for each tomogram."""

    _label = 'import subtomograms from tbl files'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.matchingFiles = None
        self.tblTomoDict = None
        self.tblSubtomoFilesDict = {}
        self.sRate = Float()

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        super()._defineImportParams(form)
        form.addParam('tomoSet', PointerParam,
                      pointerClass='SetOfTomograms',
                      allowsNull=True,
                      label="Tomograms (opt.)",
                      help="If not provided, the subtomograms won't be referred to any tomogram. "
                           "If provided, the sampling rate value will be read from them.")
        form.addParam('samplingRate', FloatParam,
                      condition='not tomoSet',
                      allowsNull=True,
                      label='Sampling rate [Å/px] (opt.)')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._initialize()
        self._insertFunctionStep(self.convertInputStep)
        self._insertFunctionStep(self.importSubTomogramsStep)

    # --------------------------- STEPS functions -----------------------------
    def _initialize(self):
        tomograms = self.tomoSet.get()
        self.matchingTblFiles = self.getMatchFiles()
        # If tomograms provided, the sampling rate value will be read from them. If not, the user
        # will have to introduce it manually
        if tomograms:
            self.tblTomoDict = {tblFile: tomo.clone() for tblFile, tomo in zip(self.matchingTblFiles, self.tomoSet.get())}
            self.sRate.set(tomograms.getSamplingRate())
        else:
            self.tblTomoDict = {tblFile: None for tblFile in self.matchingTblFiles}
            self.sRate.set(self.samplingRate.get())

    def convertInputStep(self):
        # Check the subtomograms file extension in the first directory expected to contain some (the directory in which
        # the corresponding tbl file is contained)
        for tblFile in self.matchingTblFiles:
            filesDir = dirname(tblFile)
            dirBaseName = self._getExtraPath(basename(normpath(filesDir)))
            makePath(dirBaseName)
            # Get the subtomo from the current tbl corresponding directory excluding the tbl file
            subtomoFiles = list(set(glob.glob(join(filesDir, '*'))) - set(glob.glob(join(filesDir, '*.tbl'))))
            files2store = []
            if subtomoFiles[0].endswith('.mrc'):
                # Create links
                for subtomoFile in sorted(subtomoFiles):
                    newFileName = join(dirBaseName, basename(subtomoFile))
                    createLink(subtomoFile, newFileName)
                    files2store.append(newFileName)
            else:
                # Convert the files
                for subtomoFile in subtomoFiles:
                    mrcTomoFile = join(dirBaseName, basename(subtomoFile) + '.mrc')
                    args = '-i %s -o %s -t vol ' % (subtomoFile, mrcTomoFile)
                    xmipp3.Plugin.runXmippProgram('xmipp_image_convert', args)
                    files2store.append(mrcTomoFile)
            self.tblSubtomoFilesDict[tblFile] = files2store

    def importSubTomogramsStep(self):
        imgh = ImageHandler()
        samplingRate = self.sRate.get()
        coordSet = None
        subtomoSet = SetOfSubTomograms.create(self._getPath(), template='subtomograms%s.sqlite')
        subtomoSet.setSamplingRate(samplingRate)
        subtomo = SubTomogram()
        subtomo.setSamplingRate(samplingRate)
        tomograms = self.tomoSet.get()
        if tomograms:
            coordSet = SetOfCoordinates3D.create(self._getPath(), template='coordinates3d%s.sqlite')
            coordSet.setPrecedents(tomograms)
            coordSet.setSamplingRate(samplingRate)
            coordSet.setBoxSize(20)
        for tblFile, tomo in self.tblTomoDict.items():
            with open(tblFile, 'r') as fhTable:
                lines = fhTable.readlines()
                # The subtomograms files are generated in the same directory as the .tbl file, one for each tomogram
                for line, fileName in zip(lines, self.tblSubtomoFilesDict[tblFile]):
                    _, _, _, n = imgh.getDimensions(fileName)
                    self._fillSubtomogram(line, subtomo, subtomoSet, fileName, tomo=tomo, coordSet=coordSet)
        # Needs to be registered before assigning it to the set of subtomograms
        self._defineOutputs(**{self._possibleOutputs.coordinates.name: coordSet})
        if tomograms:
            subtomoSet.setCoordinates3D(coordSet)
        self._defineOutputs(**{self._possibleOutputs.subtomograms.name: subtomoSet})
        if tomograms:
            self._defineSourceRelation(tomograms, subtomoSet)

    @staticmethod
    def _fillSubtomogram(line, subtomo, subtomoSet, newFileName, tomo=None, coordSet=None):
        """ adds a subtomogram to a set """
        subtomo.cleanObjId()
        subtomo.setFileName(newFileName)
        dynTableLine2Subtomo(line, subtomo, subtomoSet, tomo=tomo, coordSet=coordSet)

    # --------------------------- INFO functions ------------------------------
    def _validate(self):
        errors = []
        matchingFiles = self.getMatchFiles()
        tblFilesNotFoundMsg = 'No Dynamo .tbl files were detected in the introduced path/s.'
        if matchingFiles:
            if not matchingFiles[0].endswith('.tbl'):
                errors.append(tblFilesNotFoundMsg)
        else:
            errors.append(tblFilesNotFoundMsg)
        if not self.tomoSet.get() and not self.samplingRate.get():
            errors.append('If tomograms are not provided, a sampling rate value is required.')
        return errors
