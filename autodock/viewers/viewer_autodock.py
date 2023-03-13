# **************************************************************************
# *
# * Authors:  Carlos Oscar Sorzano (coss@cnb.csic.es)
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

import os, re

from pwem.viewers.viewer_chimera import Chimera, ChimeraView

from pwchem.viewers import SmallMoleculesViewer
from autodock.protocols import ProtChemAutodock, ProtChemVinaDocking, ProtChemAutodockGPU

CHIMERA_ERROR = 'Chimera program is not found were it was expected: \n\n{}\n\n' \
                'Either install ChimeraX in this path or install our ' \
                'scipion-em-chimera plugin'.format(Chimera.getProgram())

SINGLE, MOLECULE, POCKET = 'single', 'molecule', 'pocket'

def chimeraInstalled():
  return Chimera.getHome() and os.path.exists(Chimera.getProgram())

class ProtAutodockDockingViewer(SmallMoleculesViewer):
    """ Visualize the output of protocol autodock """
    _label = 'Viewer autodock docking'
    _targets = [ProtChemAutodock, ProtChemAutodockGPU]

    def __init__(self, **args):
        super().__init__(**args)

class ProtVinaDockingViewer(SmallMoleculesViewer):
    """ Visualize the output of protocol autodock """
    _label = 'Viewer vina docking'
    _targets = [ProtChemVinaDocking]
    _viewerOptions = ['PyMol', 'ChimeraX', 'ViewDockX']

    def viewDockChimeraXMols(self, mols, ligandLabel, disable=True, e=None):
        if not chimeraInstalled():
            print(CHIMERA_ERROR)
            return [self.warnMessage(CHIMERA_ERROR, 'Chimera not found')]
        else:
            pmlsDir = self.getPmlsDir()
            chimScript = os.path.join(pmlsDir, '{}_chimeraX.py'.format(ligandLabel))

            molADFiles = set([])
            for mol in mols:
                poseFile = mol.getPoseFile()
                molADFiles.add(re.sub("_[1-9]*.pdbqt", ".pdbqt", poseFile))

            with open(chimScript, "w") as f:
                f.write("from chimerax.core.commands import run\n")

                f.write("run(session, 'cd %s')\n" % os.getcwd())
                f.write("run(session, 'cofr 0,0,0')\n")  # set center of coordinates

                _inputStruct = os.path.abspath(self.protocol.getOriginalReceptorFile())
                f.write("run(session, 'open %s')\n" % _inputStruct)

                i=2
                for mol in molADFiles:
                    f.write("run(session, 'open %s')\n" % mol)
                    if disable:
                        f.write("run(session, 'hide #%s models')\n" % i)
                    i += 1
                f.write("run(session, 'viewdockx')\n")

            view = ChimeraView(chimScript)
            return [view]


