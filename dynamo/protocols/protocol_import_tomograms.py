# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
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
import numpy as np
import pyworkflow.protocol.params as params
from pyworkflow.utils import Message
from pwem.objects import Transform
from tomo.protocols.protocol_base import ProtTomoImportFiles
from tomo.objects import Tomogram
from tomo.protocols.protocol_base import ProtTomoImportAcquisition
from .. import Plugin
from ..convert import readDynCatalogue


class DynamoImportTomograms(ProtTomoImportFiles, ProtTomoImportAcquisition):
    """This protocols imports a series of Tomogram stored in a Dynamo catalogue into Scipion.
    In order to avoid handling a MatLab binary, the script relies on MatLab itself to turn a
    binary MatLab object into an Structure which can be afterwards read by Python."""

    _label = 'import tomograms from Dynamo'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection('Import')
        form.addParam('ctgPath', params.PathParam,
                      label="Dynamo catalogue file",
                      help='Path to Dynamo catalogues.')
        # form.addParam('filesPattern', params.StringParam,
        #               label='Pattern',
        #               help="Pattern of the files to be imported.\n\n"
        #                    "The pattern can contain standard wildcards such as\n"
        #                    "*, ?, etc, or special ones like ### to mark some\n"
        #                    "digits in the filename as ID.\n\n"
        #                    "NOTE: wildcards and special characters "
        #                    "('*', '?', '#', ':', '%') cannot appear in the "
        #                    "actual path.")
        form.addParam('samplingRate', params.FloatParam,
                      label=Message.LABEL_SAMP_RATE)
        form.addSection('Origin Info')
        form.addParam('setOrigCoord', params.BooleanParam,
                      label="Set origin of coordinates",
                      help="Option YES:\nA new volume will be created with "
                           "the "
                           "given ORIGIN of coordinates. This ORIGIN will be "
                           "set in the map file header.\nThe ORIGIN of "
                           "coordinates will be placed at the center of the "
                           "whole volume if you select n(x)/2, n(y)/2, "
                           "n(z)/2 as "
                           "x, y, z coordinates (n(x), n(y), n(z) are the "
                           "dimensions of the whole volume). However, "
                           "selecting "
                           "0, 0, 0 as x, y, z coordinates, the volume will be "
                           "placed at the upper right-hand corner.\n\n"
                           "Option NO:\nThe ORIGIN of coordinates will be "
                           "placed at the center of the whole volume ("
                           "coordinates n(x)/2, n(y)/2, n(z)/2 by default). "
                           "This "
                           "ORIGIN will NOT be set in the map file header.\n\n"
                           "WARNING: In case you want to process "
                           "the volume with programs requiring a specific "
                           "symmetry regarding the origin of coordinates, "
                           "for example the protocol extract unit "
                           "cell, check carefully that the coordinates of the "
                           "origin preserve the symmetry of the whole volume. "
                           "This is particularly relevant for loading "
                           "fragments/subunits of the whole volume.\n",
                      default=False)
        form.addLine('Offset',
                     help="A wizard will suggest you possible "
                          "coordinates for the ORIGIN. In MRC volume "
                          "files, the ORIGIN coordinates will be "
                          "obtained from the file header.\n "
                          "In case you prefer set your own ORIGIN "
                          "coordinates, write them here. You have to "
                          "provide the map center coordinates in "
                          "Angstroms (pixels x sampling).\n",
                     condition='setOrigCoord')
        # line.addParam would produce a nicer looking form
        # but them the wizard icon is drawn outside the visible
        # window. Until this bug is fixed form is a better option
        form.addParam('x', params.FloatParam, condition='setOrigCoord',
                      label="x", help="offset along x axis (Angstroms)")
        form.addParam('y', params.FloatParam, condition='setOrigCoord',
                      label="y", help="offset along y axis (Angstroms)")
        form.addParam('z', params.FloatParam, condition='setOrigCoord',
                      label="z", help="offset along z axis (Angstroms)")
        ProtTomoImportAcquisition._defineParams(self, form)

    # ----------------------- INSERT STEPS functions --------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importDynamoTomogramsStep', self.samplingRate.get())

    # --------------------------- STEPS functions -----------------------------
    def importDynamoTomogramsStep(self, samplingRate):
        """ Copy images matching the filename pattern. Register other parameters."""
        # self.info("Using pattern: '%s'" % pattern)
        self._parseAcquisitionData()

        # Create a Volume template object
        tomo = Tomogram()
        tomo.setSamplingRate(samplingRate)
        tomoSet = self._createSetOfTomograms()
        tomoSet.setSamplingRate(samplingRate)

        # for fileName, fileId in self.iterFiles():

        # Convert the catalogue bynary into a Python readable file
        tdb = readDynCatalogue(self.ctgPath.get(), self._getExtraPath())

        # Append Catalogue Tomograms to output set
        volumes_structs = tdb.volumes if isinstance(tdb.volumes, np.ndarray) else [tdb.volumes]
        for idv, volume in enumerate(volumes_structs):
            tomo.cleanObjId()
            tomo.setLocation(volume.file)
            x, y, z = tomo.getDim()
            origin = Transform()
            origin.setShifts(x / -2. * samplingRate,
                             y / -2. * samplingRate,
                             z / -2. * samplingRate)
            tomo.setOrigin(origin)
            tomo.setAcquisition(self._extractAcquisitionParameters(tomo.getFileName()))
            tomoSet.append(tomo)

        self._defineOutputs(outputTomograms=tomoSet)

    # --------------------------- INFO functions ------------------------------
    def _hasOutput(self):
        return self.hasAttribute('outputTomograms')

    def _getTomMessage(self):
        return "Tomograms %s" % self.getObjectTag('outputTomograms')

    def _summary(self):

        try:
            summary = []
            if self._hasOutput():
                summary.append("%s imported from:\n%s"
                               % (self._getTomMessage(), self.ctgPath.get()))

                if self.samplingRate.get():
                    summary.append(u"Sampling rate: *%0.2f* (Å/px)" % self.samplingRate.get())

                x, y, z = self.outputTomograms.getFirstItem().getShiftsFromOrigin()
                summary.append(u"Tomograms Origin (x,y,z):\n"
                               u"    x: *%0.2f* (Å/px)\n"
                               u"    y: *%0.2f* (Å/px)\n"
                               u"    z: *%0.2f* (Å/px)" % (x, y, z))

        except Exception as e:
            print(e)

        return summary

    def _methods(self):
        methods = []
        if self._hasOutput():
            methods.append(" %s imported with a sampling rate *%0.2f*" %
                           (self._getTomMessage(), self.samplingRate.get()),)
        return methods
