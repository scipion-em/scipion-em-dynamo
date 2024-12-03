# **************************************************************************
# *
# * Authors:  Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import logging
from dynamo import Plugin
from pwem.convert import transformations
from pwem.convert.transformations import euler_from_matrix, translation_from_matrix
import math, os
import numpy as np
from scipy.io import loadmat
import pyworkflow.utils as pwutils
from pyworkflow.utils.process import runJob
from pwem.convert.headers import getFileFormat, MRC
from pwem.emlib.image.image_handler import ImageHandler
from pwem.objects.data import Transform, Volume
from tomo.objects import Coordinate3D, TomoAcquisition
import tomo.constants as const

logger = logging.getLogger(__file__)


def convertOrLinkVolume(inVolume: Volume, outVolume: str):
    """Converts the inVolume into a compatible (mrc) dynamo volume named outVolume
    or links it if already compatible"""

    inFn = inVolume.getFileName()

    # If compatible with dynamo. Attention!! Assuming is not a stack of mrc volumes!!
    if getFileFormat(inFn) == MRC:
        pwutils.createLink(os.path.abspath(inFn), outVolume)
    else:
        ih = ImageHandler()
        ih.convert(inVolume, outVolume)


def writeSetOfVolumes(setOfVolumes, outputFnRoot, name):
    if name == 'id':  # write the ID of the object in the name
        for volume in setOfVolumes.iterSubtomos():
            convertOrLinkVolume(volume, "%s%03d.mrc" % (outputFnRoot, volume.getObjId()))
    if name == 'ix':  # write the INDEX of the object in the name
        for ix, volume in enumerate(setOfVolumes.iterSubtomos()):
            convertOrLinkVolume(volume, "%s%03d.mrc" % (outputFnRoot, int(ix + 1)))


def writeDynTable(fhTable, setOfSubtomograms):
    for subtomo in setOfSubtomograms.iterSubtomos():
        # Get 3d coordinates or 0, 0, 0
        if subtomo.hasCoordinate3D():
            x = subtomo.getCoordinate3D().getX(const.BOTTOM_LEFT_CORNER)
            y = subtomo.getCoordinate3D().getY(const.BOTTOM_LEFT_CORNER)
            z = subtomo.getCoordinate3D().getZ(const.BOTTOM_LEFT_CORNER)
        else:
            x = 0
            y = 0
            z = 0

        # Get alignment information
        if subtomo.hasTransform():
            tdrot, tilt, narot, shiftx, shifty, shiftz = matrix2eulerAngles(subtomo.getTransform().getMatrix())
        else:
            tilt = 0
            narot = 0
            tdrot = 0
            shiftx = 0
            shifty = 0
            shiftz = 0
        if subtomo.hasAcquisition():
            anglemin = subtomo.getAcquisition().getAngleMin()
            anglemax = subtomo.getAcquisition().getAngleMax()
        else:
            anglemin = 0
            anglemax = 0
        fhTable.write('%d 1 1 %d %d %d %d %d %d 0 0 0 1 %d %d 0 0 0 0 0 0 1 0 %d %d %d 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n'
                      % (subtomo.getObjId(), shiftx, shifty, shiftz, tdrot, tilt, narot, anglemin, anglemax, x, y, z))


def dynTableLine2Subtomo(inLine, subtomo, subtomoSet=None, tomo=None, coordSet=None):
    if type(inLine) != list:
        inLine = inLine.rstrip().split()
    # inLine = inLine.rstrip()
    subtomo.setObjId(int(inLine[0]))
    shiftx = inLine[3]
    shifty = inLine[4]
    shiftz = inLine[5]
    tdrot = inLine[6]
    tilt = inLine[7]
    narot = inLine[8]
    A = eulerAngles2matrix(tdrot, tilt, narot, shiftx, shifty, shiftz)
    transform = Transform()
    transform.setMatrix(A)
    subtomo.setTransform(transform)
    angleMin = inLine[13]
    angleMax = inLine[14]
    acquisition = TomoAcquisition()
    acquisition.setAngleMin(angleMin)
    acquisition.setAngleMax(angleMax)
    subtomo.setAcquisition(acquisition)
    volId = int(inLine[19])
    subtomo.setVolId(volId)
    classId = inLine[21]
    subtomo.setClassId(classId)
    if tomo:
        tomoOrigin = tomo.getOrigin()
        subtomo.setVolName(tomo.getFileName())
        subtomo.setOrigin(tomoOrigin)
        coordinate3d = Coordinate3D()
        coordinate3d.setVolId(tomo.getObjId())
        coordinate3d.setVolume(tomo)
        x = inLine[23]
        y = inLine[24]
        z = inLine[25]
        coordinate3d.setX(float(x), const.BOTTOM_LEFT_CORNER)
        coordinate3d.setY(float(y), const.BOTTOM_LEFT_CORNER)
        coordinate3d.setZ(float(z), const.BOTTOM_LEFT_CORNER)
        subtomo.setCoordinate3D(coordinate3d)
        coordSet.append(coordinate3d)
    if subtomoSet is not None:
        subtomoSet.append(subtomo)


def readDynCoord(tableFile, coord3DSet, tomo, scaleFactor=1):
    with open(tableFile) as fhTable:
        for nline in fhTable:
            coordinate3d = Coordinate3D()
            nline = nline.rstrip().split()
            shiftx = nline[3]
            shifty = nline[4]
            shiftz = nline[5]
            tdrot = nline[6]
            tilt = nline[7]
            narot = nline[8]
            A = eulerAngles2matrix(tdrot, tilt, narot, shiftx, shifty, shiftz)
            x = nline[23]
            y = nline[24]
            z = nline[25]
            groupId = nline[21]
            coordinate3d.setVolume(tomo)
            coordinate3d.setX(float(x) * scaleFactor, const.BOTTOM_LEFT_CORNER)
            coordinate3d.setY(float(y) * scaleFactor, const.BOTTOM_LEFT_CORNER)
            coordinate3d.setZ(float(z) * scaleFactor, const.BOTTOM_LEFT_CORNER)
            coordinate3d.setGroupId(int(groupId))
            coordinate3d.setMatrix(A)
            coord3DSet.append(coordinate3d)


# matrix2euler dynamo
def matrix2eulerAngles(matrix):
    # Relevant info:
    #   * Dynamo's transformation system is ZXZ
    #   * Sscipion = R * (-Sdynamo) ==> Sdynamo = Rinv * (-Sscipion)
    rotMatrix = matrix[:3, :3]
    rotMatrixInv = np.linalg.inv(rotMatrix)
    tdrot, tilt, narot = np.rad2deg(euler_from_matrix(rotMatrix, axes='szxz'))
    shiftsScipion = - translation_from_matrix(matrix)
    shiftsDynamo = np.dot(rotMatrixInv, shiftsScipion)
    shiftx, shifty, shiftz = shiftsDynamo
    return tdrot, tilt, narot, shiftx, shifty, shiftz


# euler2matrix dynamo
def eulerAngles2matrix(tdrot, tilt, narot, shiftx, shifty, shiftz):
    # Relevant info:
    #   * Dynamo's transformation system is ZXZ
    #   * Sscipion = R * (-Sdynamo) ==> Sdynamo = Rinv * (-Sscipion)
    M = np.eye(4)
    sx = float(shiftx)
    sy = float(shifty)
    sz = float(shiftz)
    Sdynamo = np.array([sx, sy, sz])
    tdrot = np.deg2rad(float(tdrot))
    narot = np.deg2rad(float(narot))
    tilt = np.deg2rad(float(tilt))
    R = transformations.euler_matrix(tdrot, tilt, narot, axes='szxz')
    R = R[:3, :3]
    Sscipion = - np.dot(R, Sdynamo)
    M[:3, :3] = R
    M[:3, 3] = Sscipion
    return M


def readDynCatalogue(ctlg_path, save_path):
    # MatLab script to convert an object into a structure
    matPath = os.path.join(save_path, 'structure.mat')
    codeFilePath = os.path.join(save_path, 'convert.m')
    codeFid = open(codeFilePath, 'w')
    content = "c=dread('%s')\n" \
              "s=struct(c)\n" \
              "volumes=c.volumes\n" \
              "for idv=1:length(volumes)\n" \
              "s.volumes{idv}=struct(volumes{idv})\n" \
              "end\n" \
              "save('%s','s','-v7');" % (os.path.abspath(ctlg_path),
                                         os.path.abspath(matPath))
    codeFid.write(content)
    codeFid.close()
    args = ' %s' % codeFilePath
    runJob(None, Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    # Read MatLab binary into Python
    return loadmat(matPath, struct_as_record=False, squeeze_me=True)['s']
