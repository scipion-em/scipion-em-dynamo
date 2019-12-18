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

from pyworkflow.em.convert import ImageHandler


def writeSetOfVolumes(setOfVolumes, outputFnRoot):
    ih = ImageHandler()
    for volume in setOfVolumes:
        i = volume.getObjId()
        ih.convert(volume, "%s%06d.mrc" % (outputFnRoot, i))

def writeTable(fhTable, setOfSubtomograms):
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
