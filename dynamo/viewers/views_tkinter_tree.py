# **************************************************************************
# *
# * Authors:     David Herreros
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

import os, threading

from pyworkflow import utils as pwutils
from pyworkflow.utils.process import runJob
from pyworkflow.gui.dialog import ToolbarListDialog

from dynamo import Plugin


class DynamoDialog(ToolbarListDialog):
    """
    This class extend from ListDialog to allow calling
    an Eman subprocess from a list of Tomograms.
    """

    def __init__(self, parent, path, **kwargs):
        self.path = path
        self.provider = kwargs.get("provider", None)
        ToolbarListDialog.__init__(self, parent,
                                   "Tomogram List",
                                    allowsEmptySelection=False,
                                    itemDoubleClick=self.doubleClickOnTomogram,
                                    **kwargs)

    def refresh_gui(self):
        if self.proc.isAlive():
            self.after(1000, self.refresh_gui)
        else:
            self.tree.update()
            pwutils.cleanPath(os.path.join(self.path, 'tomos'))
            pwutils.cleanPath(os.path.join(self.path, 'tomos.ctlg'))


    def doubleClickOnTomogram(self, e=None):
        self.tomo = e
        self.proc = threading.Thread(target=self.lanchDynamoForTomogram, args=(self.tomo,))
        self.proc.start()
        self.after(1000, self.refresh_gui)

    def lanchDynamoForTomogram(self, tomo):
        inputFilePath = self._writeInputFile(tomo)

        args = ' %s' % inputFilePath
        runJob(None, Plugin.getDynamoProgram(), args, env=Plugin.getEnviron())

    def _writeInputFile(self, tomo):
        inputFilePath = os.path.join(os.environ.get("SCIPION_HOME"), "software", "tmp", "commands.doc")
        tomoFile = os.path.join(os.environ.get("SCIPION_HOME"), "software", "tmp", "tomos.vll")
        catalogue = os.path.join(self.path, 'tomos')
        tomoFid = open(tomoFile, 'w')
        tomoName = tomo.getFileName()
        tomoName = os.path.basename(tomoName)
        tomoFid.write(tomoName + '\n')
        tomoFid.close()
        tomoBase = pwutils.removeBaseExt(tomo.getFileName())
        inputFid = open(inputFilePath, 'w')
        if not os.path.isfile(os.path.join(self.path, tomoBase + '.txt')):
            content = 'dcm -create %s -fromvll %s \n' \
                      'm = dmodels.membraneByLevels()\n' \
                      'm.linkCatalogue(\'%s\',\'i\',1,\'s\',1)\n' \
                      'm.saveInCatalogue()\n' \
                      'dtmslice %s -c %s \n' \
                      'modeltrack.addOne(\'model\',m)\n' \
                      'modeltrack.setActiveModel(1)\n' \
                      'uiwait(dpkslicer.getHandles().figure_fastslicer)\n' \
                      'm = dread(dcmodels(\'%s\',\'i\',1))\n' \
                      'writematrix([m.points m.group_labels\'], \'%s\')\n' \
                      'writematrix(m.mesh.tr.ConnectivityList, \'%s\')\n' \
                      'exit' % (catalogue, tomoFile, catalogue, tomo.getFileName(), catalogue,
                                catalogue, os.path.join(self.path, tomoBase + '.txt'),
                                os.path.join(self.path, tomoBase + '_connectivity.txt'))
        else:
            content = 'dcm -create %s -fromvll %s \n' \
                      'm = dmodels.membraneByLevels()\n' \
                      'modelData = readmatrix(\'%s\')\n' \
                      'addPoint(m, modelData(:,1:3), modelData(:,4))\n' \
                      'm.linkCatalogue(\'%s\',\'i\',1,\'s\',1)\n' \
                      'm.saveInCatalogue()\n' \
                      'dtmslice %s -c %s \n' \
                      'modeltrack.addOne(\'model\',m)\n' \
                      'modeltrack.setActiveModel(1)\n' \
                      'uiwait(dpkslicer.getHandles().figure_fastslicer)\n' \
                      'm = dread(dcmodels(\'%s\',\'i\',1))\n' \
                      'writematrix([m.points m.group_labels\'], \'%s\')\n' \
                      'writematrix(m.mesh.tr.ConnectivityList, \'%s\')\n' \
                      'exit' % (catalogue, tomoFile, os.path.join(self.path, tomoBase + '.txt'),
                                catalogue, tomo.getFileName(), catalogue,
                                catalogue, os.path.join(self.path, tomoBase + '.txt'),
                                os.path.join(self.path, tomoBase + '_connectivity.txt'))
        inputFid.write(content)
        inputFid.close()
        return inputFilePath
