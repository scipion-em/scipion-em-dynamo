# **************************************************************************
# *
# * Authors:    Scipion Team (scipion@cnb.csic.es)
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
import datetime
import glob
import pathlib
from os.path import join, basename, abspath, exists
from dynamo import CATALOG_FILENAME, CATALOG_BASENAME, SUFFIX_COUNT, Plugin, \
    BASENAME_CROPPED, BASENAME_PICKED, GUI_MW_FILE
from dynamo.convert import eulerAngles2matrix
from pyworkflow.object import String
from pyworkflow.utils import removeBaseExt
from tomo.constants import BOTTOM_LEFT_CORNER
from tomo.objects import SetOfCoordinates3D, Coordinate3D, SetOfMeshes


def getPickedFile(fPath, ext='.txt'):
    return join(fPath, BASENAME_PICKED + ext)


def getCroppedFile(fPath, ext='.txt'):
    return join(fPath, BASENAME_CROPPED + ext)


def getTomoPathAndBasename(filePath, tomo):
    """Path and base name of the txt file that will be generated with the corresponding
    extension using the methods below"""
    return join(filePath, tomo) if isinstance(tomo, str) else join(filePath, tomo.getTsId())


def getCurrentTomoCountFile(filePath, tomo, ext='.txt'):
    return getTomoPathAndBasename(filePath, tomo) + SUFFIX_COUNT + ext


def getFileMwFromGUI(outPath):
    return join(outPath, GUI_MW_FILE)


def getCatalogFile(fpath, withExt=True):
    return join(fpath, basename(CATALOG_FILENAME)) if withExt else join(fpath, CATALOG_BASENAME)


def genMCode4ReadDynModel(modelFile):
    """MATLAB code to read a model file from Dynamo"""
    return "m = dread('%s')\n" % abspath(modelFile)  # Load the model created in the boxing protocol


def readModels(prot, outPath, tmpPath, modelList, savePicked=True, saveCropped=True):
    """Read the models generated for each tomograms and write the info to a the corresponding file,
    depending if there was only a picking or a picking and a mesh calculation"""
    codeFile = join(outPath, 'parseModels.m')
    contents = genMCode4ReadAndSaveData(tmpPath, modelList, savePicked=savePicked, saveCropped=saveCropped)
    with open(codeFile, 'w') as codeFid:
        codeFid.write(contents)
    args = ' %s' % codeFile
    try:
        Plugin.runDynamo(prot, args)
    except:
        from pyworkflow.utils.process import runJob
        runJob(None, Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())


def genMCode4CheckModelWfFromGUI(modelFileList, outPath):
    """Create MATLAB code to read the generated models and write an empty file in case
    that the user carried out the model workflow from the boxing GUI"""
    if isinstance(modelFileList, str):
        modelFileList = [modelFileList]
    content = "modelList = {'%s'}\n" % "', '".join(modelFileList)
    content += "for i=1:length(modelList)\n"
    content += "modelFile = modelList{i}\n"
    content += "m = dread(modelFile)\n"  # Load the model
    content += "if not(isempty(m.crop_points))\n"
    content += "fid = fopen('%s', 'w')\n" % getFileMwFromGUI(outPath)
    content += "fclose(fid)\n"
    content += "break\n"
    content += "end\n"
    content += "end\n"
    return content


def genMCode4ReadAndSaveData(outPath, modelFileList, savePicked=True, saveCropped=True, modelWfCode=''):
    """MATLAB code to format and write the output data in a text file that will be read in the step create output.
    The column headers of the generated file are:
        - When no meshes were generated (only the clicked points, then): coordX, coordY, coordZ, vesicleId, modelName, modelFile, volumeFile
        - When meshes were generated (cropped points and angles, interpolation): coordX, coordY, coordZ, rot, tilt, psi, vesicleId, modelName, modelFile, volumeFile
    Input modelWfCode can be used to introduce specific code for a model workflow processing
    :param outPath: path where the results files will be generated, normally the directory 'extra'.
    :param modelFileList: list of the Dynamo model files to be processed.
    :param savePicked: flag to indicate if the manually picked points should be saved.
    :param saveCropped: the same as the previous one, but for the interpolated (model workflow result) points and angles.
    :param modelWfCode: string containing the MATLAB code which corresponds to the steps that have to be carried out for
    a specific Dynamo model.
    """
    # Generate the MATLAB code
    if isinstance(modelFileList, str):
        modelFileList = [modelFileList]
    content = "modelList = {'%s'}\n" % "', '".join(modelFileList)
    content += "for i=1:length(modelList)\n"
    content += "rng('shuffle')\n"
    content += "vesicleId = randi(1000)\n"  # To ensure the different vesicles are annotated with different id, a random number in [1, 1000] is generated for each
    content += "modelFile = modelList{i}\n"
    content += "m = dread(modelFile)\n"  # Load the model
    content += "modelVolume = m.cvolume.file\n"
    if savePicked:  # If a points file name is introduced, it means that the clicked points must be saves
        pointsFile = getPickedFile(outPath)
        content += "pointsClickedMatrix = m.points\n"
        content += "for row=1:size(pointsClickedMatrix, 1)\n"
        content += "cRow = pointsClickedMatrix(row, :)\n"
        content += "writecell({cRow(1), cRow(2), cRow(3), vesicleId, m.name, modelFile, modelVolume}, '%s', " \
                   "'WriteMode','append', 'Delimiter', 'tab')\n" % pointsFile
        content += "end\n"
    if saveCropped:
        croppedFile = getCroppedFile(outPath)
        # In the meshes were generated (in the Dynamo GUI or with the model workflow protocol), then there will
        # be cropped points and angles
        if modelWfCode:
            content += modelWfCode
        content += "coordsMatrix = m.crop_points\n"
        content += "anglesMatrix = m.crop_angles\n"
        content += "for row=1:size(coordsMatrix, 1)\n"
        content += "cRow = coordsMatrix(row, :)\n"  # Leave only two decimals for the coordinates
        content += "aRow = anglesMatrix(row, :)\n"  # The same for the angles
        content += "writecell({cRow(1), cRow(2), cRow(3), aRow(1), aRow(2), aRow(3), vesicleId, m.name, modelFile, modelVolume}, " \
                   "'%s', 'WriteMode','append', 'Delimiter', 'tab')\n" % croppedFile
        # content += "end\n"
        content += "end\n"
    content += "end\n"
    return content


def createSetOfOutputCoords(protPath, outPath, precedentsPointer, boxSize=20):
    # Create the output set
    precedents = precedentsPointer.get()
    outCoords = SetOfCoordinates3D.create(protPath, template='coordinates3d%s.sqlite')
    outCoords.setPrecedents(precedentsPointer)
    outCoords.setSamplingRate(precedents.getSamplingRate())
    outCoords.setBoxSize(boxSize)
    outCoords._dynCatalogue = String(getCatalogFile(outPath))  # Extended attribute
    return outCoords


def dynamoCroppingResults2Scipion(outCoords, croppedFile, tomoFileDict):
    with open(croppedFile, 'r') as coordFile:
        for line in coordFile:
            coord = Coordinate3D()
            values = line.replace('\n', '').split('\t')
            tomoFile = values[9]
            tomo = tomoFileDict[tomoFile]
            coord.setVolume(tomo)
            coord.setTomoId(tomo.getTsId())
            coordinates = float(values[0]), float(values[1]), float(values[2])
            angles = float(values[3]), float(values[4]), float(values[5])
            coord.setPosition(*coordinates, BOTTOM_LEFT_CORNER)
            matrix = eulerAngles2matrix(*angles, 0, 0, 0)  # There are no shifts at this point
            coord.setMatrix(matrix)
            coord.setGroupId(int(values[6]))
            # Extended attributes
            coord._dynModelName = String(values[7])
            coord._dynModelFile = String(values[8])
            outCoords.append(coord)


def getDynamoModels(fpath):
    """Search recursively for Dynamo model files (omd) in a given directory (usually extra)"""
    return glob.glob(join(fpath, '**/*.omd'), recursive=True)


def createBoxingOutputObjects(prot, precedentsPointer, boxSize=20, savePicked=True):
    meshes = None
    outCoords = None
    precedents = precedentsPointer.get()
    outPath = prot._getExtraPath()
    tomoFileDict = {abspath(tomo.getFileName()): tomo.clone() for tomo in precedents}
    # Create the output set of coordinates (only for the models to which the user carried out the model workflow
    # from the Dynamo GUI
    tmpPath = prot._getTmpPath()
    pickedFile = getPickedFile(tmpPath)
    croppedFile = getCroppedFile(tmpPath)
    dynamoModels = getDynamoModels(outPath)
    if dynamoModels:
        saveCropped = didUserMwInGui(prot, dynamoModels)
        readModels(prot, outPath, tmpPath, dynamoModels, savePicked=savePicked, saveCropped=saveCropped)
        if savePicked:
            # Create the output set of meshes (always produced)
            meshes = SetOfMeshes.create(prot._getPath(), template='meshes%s.sqlite')
            meshes.setPrecedents(precedents)
            meshes.setSamplingRate(precedents.getSamplingRate())
            meshes.setBoxSize(boxSize)
            meshes._dynCatalogue = String(getCatalogFile(outPath))  # Extended attribute
            # Save picked points to Scipion
            with open(pickedFile, 'r') as coordFile:
                for line in coordFile:
                    coord = Coordinate3D()
                    values = line.replace('\n', '').split('\t')
                    tomoFile = values[6]
                    tomo = tomoFileDict[tomoFile]
                    coord.setVolume(tomo)
                    coord.setTomoId(tomo.getTsId())
                    coord.setPosition(float(values[0]), float(values[1]), float(values[2]), BOTTOM_LEFT_CORNER)
                    coord.setGroupId(int(values[3]))
                    # Extended attributes
                    coord._dynModelName = String(values[4])
                    coord._dynModelFile = String(values[5])
                    meshes.append(coord)
        if saveCropped:
            outCoords = createSetOfOutputCoords(prot._getPath(), outPath, precedentsPointer, boxSize=prot.boxSize.get())
            dynamoCroppingResults2Scipion(outCoords, croppedFile, tomoFileDict)

    return meshes, outCoords


def didUserMwInGui(prot, dynamoModels):
    """Reads the generated models and creates an empty txt file in case the user carried out at least one
    model workflow from the boxing GUI"""
    tmpPath = prot._getTmpPath()
    mCode = genMCode4CheckModelWfFromGUI(dynamoModels, tmpPath)
    mCodeFile = prot._getExtraPath('checkModelWfFromGUI.m')
    with open(mCodeFile, 'w') as codeFid:
        codeFid.write(mCode)
    args = ' %s' % mCodeFile
    try:
        Plugin.runDynamo(prot, args)
    except:
        from pyworkflow.utils.process import runJob
        runJob(None, Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())
    return exists(getFileMwFromGUI(tmpPath))


def getNewestModelModDate(modelList):
    """Get the last modification datetime of the newest model file"""
    tSt = sorted([pathlib.Path(fname).stat().st_mtime for fname in modelList])[-1]
    return datetime.datetime.fromtimestamp(tSt)
