# **************************************************************************
# *
# * Authors:    Scipion Team (scipion@cnb.csic.es)
# *
# *  BCU, Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
from os.path import exists
from dynamo.utils import getCurrentTomoCountFile
from tomo.viewers.views_tkinter_tree import TomogramsTreeProvider


class DynamoTomogramProvider(TomogramsTreeProvider):

    def __init__(self, tomoList, path, mode=None):
        super().__init__(tomoList, path, mode)

    def getColumns(self):
        return [('TomoId', 200), ("No. coords", 100), ('status', 150)]

    def getObjectInfo(self, tomogram):
        count = 0
        resDict = {'key': tomogram.getTsId(),
                   'parent': None,
                   'values': (count, 'TO DO'),
                   'tags': 'pending'}
        coordsInTomoCountFile = getCurrentTomoCountFile(self._path, tomogram)
        if exists(coordsInTomoCountFile):
            with open(coordsInTomoCountFile, 'r') as fn:
                count = int(fn.read())

        if count > 0:
            resDict['values'] = (count, 'DONE')
            resDict['tags'] = 'done'

        return resDict
