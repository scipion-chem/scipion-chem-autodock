# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb
from ..protocols import Autodock_GridGeneration
import os


class TestGridADT(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')

        setupTestProject(cls)
        cls._runImportPDB()

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1,
            pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1.pdb'))
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    def _runCreateGrid(self, spacing, radius):
        protGrid = self.newProtocol(
            Autodock_GridGeneration,
            inputAtomStruct=self.protImportPDB.outputPdb,
            radius=radius,
            spacing=spacing)

        self.launchProtocol(protGrid)
        gridOut = getattr(protGrid, 'outputGrid', None)
        self.assertIsNotNone(gridOut)
        return protGrid

    def testGridGeneration(self):
        spacing, radius = 1.0, 37.0

        gridProt = self._runCreateGrid(spacing, radius)
        grid = gridProt.outputGrid
        grid_path = os.path.abspath(grid.getFileName())

        self.assertIsNotNone(grid,
                             "There was a problem with grid generation - Please check it")

        self.assertTrue(grid_path.endswith(".map"),
                        "Format grid is incorrect")

        print('Grid  stored radius: ', grid.getRadius())
        self.assertEqual(spacing, grid.getSpacing(), 'Parameters of the grid incorrectly saved')


