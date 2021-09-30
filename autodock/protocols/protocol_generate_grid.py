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
import pandas as pd

from pyworkflow.protocol.params import  FloatParam, PointerParam
from pyworkflow.protocol.constants import LEVEL_ADVANCED


from pwem.protocols import EMProtocol
from pwem.convert import AtomicStructHandler

from autodock import Plugin as autodock_plugin
from autodock.objects import GridADT
from autodock.utils import generate_gpf, calculate_centerMass


class Autodock_GridGeneration(EMProtocol):
    """
    Calls ADT to prepare a grid with an input that is the prepared target protein
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
                      expertLevel=LEVEL_ADVANCED,
                      help="Distance between each point in the electrostatic grid."
                           " This value is used to adjust the radius as number of "
                           "(x,y,z) points : radius/spacing = number of points along"
                           " 3 dimensions ")



    # --------------------------- Steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('getpdbqt')
        self._insertFunctionStep('prepareGrid')
        self._insertFunctionStep('createOutput')


    def getpdbqt(self):
        """ Prepare the pdbqt used to generate the grid file e-map"""

        filename = self.inputAtomStruct.get().getFileName()
        if filename.endswith('.cif'):
            fnIn = self._getTmpPath("atomStruct.pdb")
            name_protein = (os.path.basename(filename)).split(".")[0]
            aStruct1 = AtomicStructHandler(filename)
            aStruct1.write(fnIn)
        else:
            fnIn = filename
            name_protein = (os.path.basename(fnIn)).split(".")[0]

        fnOut = self._getExtraPath('%s.pdbqt' %name_protein)

        args = ' -v -r %s -o %s' % (fnIn, fnOut)

        program = "prepare_receptor4"
        self.runJob(autodock_plugin.getMGLPath('bin/pythonsh'),
                    autodock_plugin.getADTPath('Utilities24/%s.py' % program) + args)



    def prepareGrid(self):
        """
        """
        # Name of the grid
        inputAtomStruct = self.inputAtomStruct.get().getFileName()
        filename = os.path.abspath(inputAtomStruct)
        name_protein = (os.path.basename(inputAtomStruct)).split(".")[0]

        # pdbqt gasteiger path
        pdbqt = self._getExtraPath('%s.pdbqt' % name_protein)

        # Center of the molecule
        structure, x_center, y_center, z_center = calculate_centerMass(filename)

        # Create the GPF file, required by autogrid to build the grid and glg
        npts = (self.radius.get()*2)/self.spacing.get()  # x,y,z points of the grid

        gpf_file = generate_gpf(name_protein, spacing=self.spacing.get(),
                                     xc=x_center, yc=y_center, zc=z_center,
                                     npts=npts, outDir=self._getExtraPath())

        # Create GLG file
        glg_file = os.path.abspath(self._getExtraPath("library.glg"))
        open(glg_file, mode='a').close()

        args = "-p %s -l %s" % (gpf_file, glg_file)
        self.runJob(autodock_plugin.getAutodockPath("autogrid4"), args, cwd=self._getExtraPath())
        e_map_file = self._getExtraPath("%s.e.map" %name_protein)
        self.grid = GridADT(e_map_file, inputAtomStruct, radius=self.radius.get(), spacing=self.spacing.get(),
                            massCX=x_center, massCY=y_center, massCZ=z_center, npts=npts)


    def createOutput(self):
        self._defineOutputs(outputGrid=self.grid)
        #self._defineSourceRelation(self.inputAtomStruct, self.grid)
