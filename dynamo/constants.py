# **************************************************************************
# *
# * Authors:    David Herreros Calero (dherreros@cnb.csic.es)
# *             Scipion Team (scipion@cnb.csic.es)
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
"""
This file contains constants related to scipion-em-dynamo protocols
"""

DYNAMO_PROGRAM = 'dynamo'
DYNAMO_HOME = 'DYNAMO_HOME'
DYNAMO_VERSION_1_1_532 = '1.1.532'
DEFAULT_VERSION = DYNAMO_VERSION_1_1_532
MINIMUM_VERSION_NUM = int(DYNAMO_VERSION_1_1_532.replace('.', ''))
DYNAMO_SHIPPED_MCR = 'dynamo_activate_linux_shipped_MCR.sh'

# Dynamo files and dirs
PROJECT_DIR = 'project'
TOMOGRAMS_DIR = 'tomograms'
CATALOG_BASENAME = 'project'
CATALOG_FILENAME = '%s.ctlg' % CATALOG_BASENAME
VLL_FILE = 'tomograms.vll'
PRJ_FROM_VIEWER = 'prjFromViewer.txt'
DATA_MODIFIED_FROM_VIEWER = 'modified.txt'

# Test files created to pass data between Scipion and Dynamo
SUFFIX_COUNT = '_count'
BASENAME_PICKED = 'picked'
BASENAME_CROPPED = 'cropped'
GUI_MW_FILE = 'guiModelWf.txt'

# Tags of Dynamo models
MB_GENERAL = 'mgeneral'
MB_GENERAL_BOXES = 'mboxes'
DIPOLE_SET = 'mdipoleSet'
MB_BY_LEVELS = 'mmembraneByLevels'
MB_VESICLE = 'mvesicle'
MB_ELLIPSOIDAL = 'mellipsoidalVesicle'
MB_ELLIPSOIDAL_MARKED = 'mmarkedEllipsoidalVesicle'
FILAMENT = 'mfilament'
FILAMENT_WITH_TORSION = 'mfilamentWithTorsion'
FILAMENT_HELIX = 'mfilamentSubunitsInHelix'
FILAMENT_RINGS = 'mfilamentRings'
CUBIC_CRYSTAL = 'mcubicCrystal'
MODELS_NOT_PROCESSED_IN_MW = [DIPOLE_SET, FILAMENT, FILAMENT_WITH_TORSION,
                              FILAMENT_HELIX, FILAMENT_RINGS, CUBIC_CRYSTAL]

# Names of Dynamo models
M_GENERAL_NAME = "General"
M_GENERAL_WITH_BOXES_NAME = "General (shown in boxes)"
M_DIPOLE_SET_NAME = "Oriented particles (dipole set)"
M_SURFACE_NAME = "Surface"
M_VESICLE_NAME = "Vesicle"
M_SPH_VESICLE_NAME = "Vesicle (spherical)"
M_ELLIPSOIDAL_VESICLE_NAME = "Ellipsoidal vesicle"
M_MARKED_ELLIP_VESICLE_NAME = "Marked ellipsoidal vesicle"
M_FILAMENT_NAME = "Filament (crop along axis)"
M_FIL_WITH_TORSION_NAME = "Filament (crop on walls)"
M_FIL_HELIX_NAME = "Filament (crop on helical path)"
M_FIL_RINGS_NAME = "Filament (crop on rings along path)"
MODELS_ALLOWED_IN_MW_NAMES = [M_GENERAL_NAME, M_GENERAL_WITH_BOXES_NAME, M_SURFACE_NAME,
                              M_SPH_VESICLE_NAME, M_ELLIPSOIDAL_VESICLE_NAME, M_MARKED_ELLIP_VESICLE_NAME]

# Description of Dynamo models
M_GENERAL_DES = "Coordinates (x, y, z) without any specific property."
M_GENERAL_WITH_BOXES_DES = "Same as general, but allowing the creation of boxes when wisualized with " \
                           "the Dynamo dtmslice GUI.\n\n" \
                           "\tSpecific property: box size."
M_DIPOLE_SET_DES = "Isolated particles with orientation imparted individually.\n" \
                   "This is a simple model: it describes the picking geometry in which the user defines positions " \
                   "AND orientations of isolated particles.\n\nMore details here --> " \
                   "https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Dipole_set_models"
M_SURFACE_DES = "Model class for modelling of membranes patches by picking points on membranes boundaries. On " \
                "different x, y or z levels.\n\nMore details here --> " \
                "https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Membrane_models"
M_VESICLE_DES = "Basic container for vesicle-like models.\n\n" \
                "\tSpecific properties:\n" \
                "\t  - radius.\n" \
                "\t  - separation between points on the surface.\n\n" \
                "\tMore details here --> https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Vesicle_models"
M_ELLIPSOIDAL_VESICLE_DES = "Basic container for vesicle-like models of ellipsoidal shape.\n\n" \
                            "\tSpecific properties:\n" \
                            "\t  - radial directions (x, y, z)\n" \
                            "\t  - radius in x, y and z directions\n\n" \
                            "\tMore details here --> https://wiki.dynamo.biozentrum.unibas.ch/w/index.php/Vesicle_models"
M_MARKED_ELLIP_VESICLE_DES = "A vesicle marker with a single 3d point.\n\n" \
                             "\tSpecific property: exclusion radius around marker, used to exclude whatever is close " \
                             "to the marker."
