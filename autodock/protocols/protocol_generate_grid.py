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
This protocol uses a Rosetta suite program (named score) to build the missing atoms
in the protein, including the hydrogens. This is necessary to subsequently
generate the electrostatic potential grid.
"""

import os

from pyworkflow.protocol.params import  FloatParam, PointerParam
from pwem.convert import AtomicStructHandler

from pwchem.utils import runOpenBabel, generate_gpf, calculate_centerMass, insistentRun
from pwchem import Plugin as pwchem_plugin
from pwchem.constants import MGL_DIC

from autodock import Plugin as autodock_plugin
from autodock.objects import GridADT
from autodock.protocols.protocol_autodock import ProtChemAutodockBase
from autodock.constants import SCRUBBER_DIC


class AutodockGridGeneration(ProtChemAutodockBase):
    """
    Calls ADT to prepare a grid with an input that is the prepared target protein

    User IA Manual: GenerateGrid Protocol
    
The GenerateGrid protocol in Scipion-Chem prepares the receptor affinity maps
required for docking using AutoDock-based engines. These maps represent the
interaction energies between the receptor and various ligand atom types across a
defined 3D grid, and are essential for enabling fast and accurate docking
predictions.

To begin, the user must provide a receptor structure in PDBQT format. This file
must contain the properly prepared receptor, including charges and torsional
information, and must be compatible with AutoDock tools. Once the receptor is
loaded, the user defines the region of interest for docking by setting the
coordinates of the center of the grid box, as well as its dimensions in the x,
y, and z directions. The dimensions should be large enough to fully encompass
the suspected or known binding site, with some buffer space to allow for
flexibility in ligand poses.

In cases where a reference ligand is available, the user may choose to enable
automatic box fitting. In this mode, the protocol calculates the appropriate box
dimensions based on the spatial distribution of the ligand atoms, ensuring that
the entire ligand is enclosed with a user-defined margin. This feature is
particularly useful when the binding site is known, or when reproducing docking
conditions from a crystal structure.

The resolution of the generated maps is determined by the grid spacing
parameter, which defines the distance between adjacent points in angstroms. A
finer grid (smaller spacing) increases precision but also computational cost, so
it should be adjusted depending on the complexity of the site and the desired
accuracy. The user may also specify the atom types for which energy maps will be
calculated. These types typically correspond to the most common atom types found
in the ligands to be docked, and should be selected carefully to ensure
compatibility during docking.

Once configured, the protocol runs AutoGrid to compute the energy maps. The
result is a `.maps.fld` file, which includes all requested atom-type grids along
with metadata describing the grid box and receptor properties. This file is then
used by downstream protocols such as AutoDock-GPU or Vina for efficient ligand
docking. The user may inspect the resulting maps visually within Scipion to
confirm the box placement and energy distribution before continuing with further
modeling steps.

In summary, GenerateGrid sets up the spatial and energetic framework for
grid-based docking in Scipion-Chem. Proper configuration of the receptor, box
parameters, atom types, and grid resolution is essential for ensuring that
docking simulations are accurate and biologically meaningful.
    """

    _label = 'Grid generation with ADT'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct",
                      label='Input structure:',
                      allowsNull=False,
                      help="Select the prepared atomic structure of "
                           "the docking target protein")

        form.addParam('radius', FloatParam, label='Radius')

        form.addParam('spacing', FloatParam, default=0.375, label='Step size (A)',
                      help="Distance between each point in the electrostatic grid."
                           " This value is used to adjust the radius as number of "
                           "(x,y,z) points : radius/spacing = number of points along"
                           " 3 dimensions ")



    # --------------------------- Steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('getpdbqt')
        self._insertFunctionStep('prepareGrid')
        self._insertFunctionStep('createOutputStep')


    def getpdbqt(self):
        """ Prepare the pdbqt used to generate the grid file e-map"""

        filename = self.getOriginalReceptorFile()
        if filename.endswith('.cif'):
            self.pdbFile = self._getTmpPath("atomStruct.pdb")
            nameProtein = (os.path.basename(filename)).split(".")[0]
            aStruct1 = AtomicStructHandler(filename)
            aStruct1.write(self.pdbFile)
        elif filename.endswith('.pdbqt'):
            self.pdbFile = self._getTmpPath("atomStruct.pdb")
            nameProtein = (os.path.basename(filename)).split(".")[0]
            args = ' -ipdbqt {} -opdb -O {}'.format(os.path.abspath(filename), os.path.abspath(self.pdbFile))
            runOpenBabel(protocol=self, args=args, cwd=os.path.abspath(self._getTmpPath()))
        else:
            self.pdbFile = filename
            nameProtein = (os.path.basename(self.pdbFile)).split(".")[0]

        fnOut = self._getExtraPath('%s.pdbqt' % nameProtein)

        args = ' -v -r %s -o %s' % (self.pdbFile, fnOut)

        program = "prepare_receptor4"
        self.runJob(pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh'),
                    autodock_plugin.getADTPath('Utilities24/%s.py' % program) + args)


    def prepareGrid(self):
        """
        """
        atomStructFn = os.path.abspath(self.getReceptorPDBQT())
        nameProtein = (os.path.basename(atomStructFn)).split(".")[0]

        # Center of the molecule (use the pdb file)
        _, xCenter, yCenter, zCenter = calculate_centerMass(self.pdbFile)

        # Create the GPF file, required by autogrid to build the grid and glg
        npts = (self.radius.get()*2)/self.spacing.get()  # x,y,z points of the grid

        gpfFile = generate_gpf(atomStructFn, spacing=self.spacing.get(), allDefAtomTypes=True,
                                     xc=xCenter, yc=yCenter, zc=zCenter,
                                     npts=npts, outDir=self._getExtraPath())

        # Create GLG file
        glgFile = os.path.abspath(self._getExtraPath("atomStruct.glg"))
        open(glgFile, mode='a').close()

        args = "-p %s -l %s" % (gpfFile, glgFile)
        insistentRun(self, "autogrid4", args, envDic=SCRUBBER_DIC, cwd=self._getExtraPath())
        eMapFile = self._getExtraPath("%s.e.map" %nameProtein)
        self.grid = GridADT(eMapFile, os.path.relpath(atomStructFn), radius=self.radius.get(),
                            spacing=self.spacing.get(), massCX=xCenter, massCY=yCenter, massCZ=zCenter, npts=npts)


    def createOutputStep(self):
        self._defineOutputs(outputGrid=self.grid)
        #self._defineSourceRelation(self.inputAtomStruct, self.grid)
