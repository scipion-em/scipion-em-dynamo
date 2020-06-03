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
import math
import numpy as np
from pwem import Domain
from pwem.emlib.image.image_handler import ImageHandler
from pwem.objects.data import Transform
Coordinate3D = Domain.importFromPlugin("tomo.objects", "Coordinate3D")
TomoAcquisition = Domain.importFromPlugin("tomo.objects", "TomoAcquisition")


def writeVolume(volume, outputFn):
    ih = ImageHandler()
    ih.convert(volume, "%s.mrc" % outputFn)


def writeSetOfVolumes(setOfVolumes, outputFnRoot, name):
    ih = ImageHandler()
    if name == 'id':  # write the ID of the object in the name
        for volume in setOfVolumes:
            ih.convert(volume, "%s%03d.mrc" % (outputFnRoot, volume.getObjId()))
    if name == 'ix':  # write the INDEX of the object in the name
        for ix, volume in enumerate(setOfVolumes):
            ih.convert(volume, "%s%03d.mrc" % (outputFnRoot, int(ix+1)))


def writeDynTable(fhTable, setOfSubtomograms):
    for subtomo in setOfSubtomograms.iterItems():
        if subtomo.hasCoordinate3D():
            x = subtomo.getCoordinate3D().getX()
            y = subtomo.getCoordinate3D().getY()
            z = subtomo.getCoordinate3D().getZ()
        else:
            x = 0
            y = 0
            z = 0
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


def readDynTable(self, item):
    nline = next(self.fhTable)
    nline = nline.rstrip()
    id = int(nline.split()[0])
    item.setObjId(id)
    shiftx = nline.split()[3]
    shifty = nline.split()[4]
    shiftz = nline.split()[5]
    tdrot = nline.split()[6]
    tilt = nline.split()[7]
    narot = nline.split()[8]
    A = eulerAngles2matrix(tdrot, tilt, narot, shiftx, shifty, shiftz)
    transform = Transform()
    transform.setMatrix(A)
    item.setTransform(transform)
    angleMin = nline.split()[13]
    angleMax = nline.split()[14]
    acquisition = TomoAcquisition()
    acquisition.setAngleMin(angleMin)
    acquisition.setAngleMax(angleMax)
    item.setAcquisition(acquisition)
    volId = nline.split()[19]
    item.setVolId(volId)
    x = nline.split()[23]
    y = nline.split()[24]
    z = nline.split()[25]
    coordinate3d = Coordinate3D()
    coordinate3d.setVolId(volId)
    coordinate3d.setX(float(x))
    coordinate3d.setY(float(y))
    coordinate3d.setZ(float(z))
    item.setCoordinate3D(coordinate3d)
    classId = nline.split()[21]
    item.setClassId(classId)


def readDynCoord(tableFile, coord3DSet, tomo):
    with open(tableFile) as fhTable:
        for nline in fhTable:
            coordinate3d = Coordinate3D()
            nline = nline.rstrip()
            shiftx = nline.split()[3]
            shifty = nline.split()[4]
            shiftz = nline.split()[5]
            tdrot = nline.split()[6]
            tilt = nline.split()[7]
            narot = nline.split()[8]
            A = eulerAngles2matrix(tdrot, tilt, narot, shiftx, shifty, shiftz)
            x = nline.split()[23]
            y = nline.split()[24]
            z = nline.split()[25]
            coordinate3d.setX(float(x))
            coordinate3d.setY(float(y))
            coordinate3d.setZ(float(z))
            coordinate3d.setMatrix(A)
            coordinate3d.setVolume(tomo)
            coord3DSet.append(coordinate3d)


# matrix2euler dynamo
def matrix2eulerAngles(A):
    tol = 1e-4
    if abs(A[2, 2] - 1) < tol:
        tilt = 0
        narot = math.atan2(A[1, 0], A[0, 0]) * 180 / math.pi
        tdrot = 0
    elif abs(A[2, 2] + 1) < tol:
        tdrot = 0
        tilt = 180
        narot = math.atan2(A[1, 0], A[0, 0]) * 180 / math.pi
    else:
        tdrot = math.atan2(A[2, 0], A[2, 1])
        tilt = math.acos(A[2, 2])
        narot = math.atan2(A[0, 2], -A[1, 2])
    tilt = tilt * 180 / math.pi
    narot = narot * 180 / math.pi
    tdrot = tdrot * 180 / math.pi
    return tdrot, tilt, narot, A[0, 3], A[1, 3], A[2, 3]


# euler2matrix dynamo
def eulerAngles2matrix(tdrot, tilt, narot, shiftx, shifty, shiftz):
    tdrot = float(tdrot)
    tilt = float(tilt)
    narot = float(narot)
    tdrot = tdrot * math.pi / 180
    narot = narot * math.pi / 180
    tilt = tilt * math.pi / 180
    costdrot = math.cos(tdrot)
    cosnarot = math.cos(narot)
    costilt = math.cos(tilt)
    sintdrot = math.sin(tdrot)
    sinnarot = math.sin(narot)
    sintilt = math.sin(tilt)
    A = np.empty([4, 4])
    A[3, 3] = 1
    A[3, 0:3] = 0
    A[0, 3] = float(shiftx)
    A[1, 3] = float(shifty)
    A[2, 3] = float(shiftz)
    A[0, 0] = costdrot*cosnarot - sintdrot*costilt*sinnarot
    A[0, 1] = -cosnarot*sintdrot - costdrot*costilt*sinnarot
    A[0, 2] = sinnarot*sintilt
    A[1, 0] = costdrot*sinnarot + cosnarot*sintdrot*costilt
    A[1, 1] = costdrot*cosnarot*costilt - sintdrot*sinnarot
    A[1, 2] = -cosnarot*sintilt
    A[2, 0] = sintdrot*sintilt
    A[2, 1] = costdrot*sintilt
    A[2, 2] = costilt
    return A
