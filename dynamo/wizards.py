# **************************************************************************
# *
# * Authors:     David Herreros Calero (dherreros@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

from pwem.wizards import EmWizard

from dynamo.protocols import DynamoExtraction, DynamoModelWorkflow


# =============================================================================
# EXTRACTION
# =============================================================================


class DynamoExtractionWizard(EmWizard):
    _targets = [(DynamoExtraction, ['boxSize'])]

    def show(self, form):
        DynamoExtractProt = form.protocol
        inputCoordinates = DynamoExtractProt.inputCoordinates.get()
        if not inputCoordinates:
            print('You must specify input coordinates')
            return

        aux = inputCoordinates.getBoxSize()
        if not aux % 2 == 0:
            aux += 1
        boxSize = aux
        if not boxSize:
            print('These coordinates do not have box size. Please, enter box size manually.')
            return

        if DynamoExtractProt.downFactor.get() != 1:
            boxSize = float(boxSize/DynamoExtractProt.downFactor.get())

        form.setVar('boxSize', boxSize)


# =============================================================================
# MODEL WORKFLOW
# =============================================================================

class DynamoModelWorkflowWizard(EmWizard):
    _targets = [(DynamoModelWorkflow, ['boxSize'])]

    def show(self, form):
        DynamoModelWorkflowProt = form.protocol
        inputMeshes = DynamoModelWorkflowProt.inputMeshes.get()
        if not inputMeshes:
            print('You must specify input meshes')
            return

        aux = inputMeshes.getBoxSize()
        if not aux % 2 == 0:
            aux += 1
        boxSize = aux
        if not boxSize:
            print('These coordinates do not have box size. Please, enter box size manually.')
            return

        form.setVar('boxSize', boxSize)