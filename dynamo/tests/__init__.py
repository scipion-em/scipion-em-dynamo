# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Estrella Fernandez Gimenez (me.fernandez@cnb.csic.es)
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
from enum import Enum

from pyworkflow.tests import DataSet

DYNAMO_TEST_DATASET = 'Dynamo'


class DataSetDynamo(Enum):
    meshesFromBoxingProtSqlite = 'meshes.sqlite'
    projectCatalog = 'project.ctlg'
    listOfTomosFile = 'tomograms.vll'
    vol1ModelsDir = 'volume_1/models'
    vol2ModelsDir = 'volume_2/models'


DataSet(name=DYNAMO_TEST_DATASET, folder=DYNAMO_TEST_DATASET, files={el.name: el.value for el in DataSetDynamo})