# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Alberto Manuel Parra Pérez (amparraperez@gmail.com)
# *
# * Biocomputing Unit, CNB-CSIC
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
This wizard will show the structure of the pdb using a matplotlib viewer 
to select the radius of the sphere that contains the protein or a desired zone.
"""

# Imports
from pwem.wizards.wizard import EmWizard
from pwem.wizards.wizards_3d.mask_structure_wizard import MaskStructureWizard
from ..protocols.protocol_generate_grid import Autodock_GridGeneration
from ..protocols.protocol_autodock import ProtChemAutodock


class GetDistance2Center(EmWizard):
    _targets = [(Autodock_GridGeneration, ['radius']),
                (ProtChemAutodock, ['radius'])]

    def show(self, form):
        protocol = form.protocol
        structure = protocol.inputAtomStruct.get()
        if not structure:
            print('You must specify input structure')
            return
        plt = MaskStructureWizard(structure.getFileName())
        plt.initializePlot()
        form.setVar('radius', plt.radius)
