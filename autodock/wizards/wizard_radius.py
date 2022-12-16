# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Alberto Manuel Parra Pérez (amparraperez@gmail.com)
# *           Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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
from pwchem.wizards import GetRadiusProtein
from ..protocols import *

for prot in [Autodock_GridGeneration, ProtChemAutodock, ProtChemVinaDocking, ProtChemAutoLigand, ProtChemAutodockGPU]:
    GetRadiusProtein().addTarget(protocol=prot,
                                 targets=['radius'],
                                 inputs=['inputAtomStruct'],
                                 outputs=['radius'])

GetRadiusProtein().addTarget(protocol=ProtChemAutodockScore,
                             targets=['radius'],
                             inputs=['inputSmallMolecules'],
                             outputs=['radius'])



