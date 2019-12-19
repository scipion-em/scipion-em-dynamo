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
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data import Transform
from tomo.objects import TomoAcquisition, Coordinate3D


def writeSetOfVolumes(setOfVolumes, outputFnRoot):
    ih = ImageHandler()
    for volume in setOfVolumes:
        i = volume.getObjId()
        ih.convert(volume, "%s%06d.mrc" % (outputFnRoot, i))


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
        fhTable.write('%d 1 0 0 0 0 0 0 0 0 0 0 1 %d %d 0 0 0 0 0 0 1 0 %d %d %d 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n'
                      % (subtomo.getObjId(), subtomo.getAcquisition().getAngleMin(),
                         subtomo.getAcquisition().getAngleMax(), x, y, z))


def readDynTable(self, item):
    nline = self.fhTable.next()
    nline = nline.rstrip()
    id = int(nline.split()[0])
    item.setObjId(id)
    shiftx = nline.split()[3]
    shifty = nline.split()[4]
    shiftz = nline.split()[5]
    rot = nline.split()[6]      # Check if they are really the same angels (zxz??)
    tilt = nline.split()[7]
    psi = nline.split()[8]
    A = eulerAngles2matrix(rot, tilt, psi, shiftx, shifty, shiftz)
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
    x = nline.split()[23]
    y = nline.split()[24]
    z = nline.split()[25]
    coordinate3d = Coordinate3D()
    coordinate3d.setVolId(volId)
    coordinate3d.setX(x)
    coordinate3d.setY(y)
    coordinate3d.setZ(z)
    item.setCoordinate3D(coordinate3d)
    classId = nline.split()[21]
    item.setClassId(classId)


def eulerAngles2matrix(alpha, beta, gamma, shiftx, shifty, shiftz):  # Check func in dynamo
    A = np.empty([4, 4])
    A.fill(2)
    A[3, 3] = 1
    A[3, 0:3] = 0
    A[0, 3] = float(shiftx)
    A[1, 3] = float(shifty)
    A[2, 3] = float(shiftz)
    alpha = float(alpha)
    beta = float(beta)
    gamma = float(gamma)
    sa = np.sin(np.deg2rad(alpha))
    ca = np.cos(np.deg2rad(alpha))
    sb = np.sin(np.deg2rad(beta))
    cb = np.cos(np.deg2rad(beta))
    sg = np.sin(np.deg2rad(gamma))
    cg = np.cos(np.deg2rad(gamma))
    cc = cb * ca
    cs = cb * sa
    sc = sb * ca
    ss = sb * sa
    A[0, 0] = cg * cc - sg * sa
    A[0, 1] = cg * cs + sg * ca
    A[0, 2] = -cg * sb
    A[1, 0] = -sg * cc - cg * sa
    A[1, 1] = -sg * cs + cg * ca
    A[1, 2] = sg * sb
    A[2, 0] = sc
    A[2, 1] = ss
    A[2, 2] = cb
    return A


def matrix2eulerAngles(A):  # Check func in dynamo
    abs_sb = np.sqrt(A[0, 2] * A[0, 2] + A[1, 2] * A[1, 2])
    if abs_sb > 16 * np.exp(-5):
        gamma = math.atan2(A[1, 2], -A[0, 2])
        alpha = math.atan2(A[2, 1], A[2, 0])
        if abs(np.sin(gamma)) < np.exp(-5):
            sign_sb = np.sign(-A[0, 2] / np.cos(gamma))
        else:
            if np.sin(gamma) > 0:
                sign_sb = np.sign(A[1, 2])
            else:
                sign_sb = -np.sign(A[1, 2])
        beta = math.atan2(sign_sb * abs_sb, A[2, 2])
    else:
        if np.sign(A[2, 2]) > 0:
            alpha = 0
            beta = 0
            gamma = math.atan2(-A[1, 0], A[0, 0])
        else:
            alpha = 0
            beta = np.pi
            gamma = math.atan2(A[1, 0], -A[0, 0])
    gamma = np.rad2deg(gamma)
    beta = np.rad2deg(beta)
    alpha = np.rad2deg(alpha)
    return alpha, beta, gamma, A[0, 3], A[1, 3], A[2, 3]



