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
from os.path import join, basename, abspath
from dynamo import CATALOG_FILENAME, CATALOG_BASENAME, SUFFIX_COUNT, SUFFIX_ONLY_PICKED, SUFFIX_CROPPED


def getTomoPathAndBasename(filePath, tomo):
    """Path and base name of the txt file that will be generated with the corresponding
    extension using the methods below"""
    return join(filePath, tomo.getTsId())


def getCurrentTomoCountFile(filePath, tomo, ext='.txt'):
    return getTomoPathAndBasename(filePath, tomo) + SUFFIX_COUNT + ext


def getCurrentTomoPointsFile(filePath, tomo, ext='.txt'):
    return getTomoPathAndBasename(filePath, tomo) + SUFFIX_ONLY_PICKED + ext


def getCurrentTomoCroppedFile(filePath, tomo, ext='.txt'):
    return getTomoPathAndBasename(filePath, tomo) + SUFFIX_CROPPED + ext


def getCatalogFile(fpath, withExt=True):
    return join(basename(CATALOG_FILENAME)) if withExt else join(fpath, CATALOG_BASENAME)


def genMCode4ReadDynModel(modelFile):
    """MATLAB code to read a model file from Dynamo"""
    return "m = dread('%s')\n" % abspath(modelFile)  # Load the model created in the boxing protocol


def genMCode4ReadAndSaveData(vesicleId, modelFile, pointsFile=None, croppedFile=None):
    """MATLAB code to format and write the output data in a text file that will be read in the step create output.
    The column headers of the generated file are:
        - When no meshes were generated (pointsFile provided --> only the clicked points, then): coordX, coordY, coordZ, vesicleId, modelName, modelFile
        - When meshes were generated (croppedFile provided --> cropped points and angles, interpolation): coordX, coordY, coordZ, rot, tilt, psi, vesicleId, modelName, modelFile"""
    content = ''
    if pointsFile:
        content += "pointsClickedMatrix = m.points\n"
        content += "for row=1:size(pointsClickedMatrix, 1)\n"
        content += "cRow = pointsClickedMatrix(row, :)\n"
        content += "writecell({cRow(1), cRow(2), cRow(3), %i, m.name, '%s'}, '%s', " \
                   "'WriteMode','append', 'Delimiter', 'tab')\n" % (vesicleId, modelFile, pointsFile)
        content += "end\n"
    if croppedFile:
        # In the meshes were generated (in the Dynamo GUI or with the model workflow protocol), then there will
        # be cropped points and angles
        content += "coordsMatrix = m.crop_points\n"
        content += "anglesMatrix = m.crop_angles\n"
        content += "if not(isempty(coordsMatrix))"
        content += "for row=1:size(coordsMatrix, 1)\n"
        content += "cRow = coordsMatrix(row, :)\n"  # Leave only two decimals for the coordinates
        content += "aRow = anglesMatrix(row, :)\n"  # The same for the angles
        content += "writecell({cRow(1), cRow(2), cRow(3), aRow(1), aRow(2), aRow(3), %i, m.name, '%s'}, '%s', " \
                   "'WriteMode','append', 'Delimiter', 'tab')\n" % (vesicleId, modelFile, croppedFile)
        content += "end\n"
        content += "end\n"
    return content
