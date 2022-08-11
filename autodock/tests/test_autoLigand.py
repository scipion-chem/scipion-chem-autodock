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
from ..protocols import ProtChemAutoLigand, Autodock_GridGeneration


class TestAutoLigand(BaseTest):
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
            pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1_noHETATM.pdb'))
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    def _runCreateGrid(self):
        protGrid = self.newProtocol(
            Autodock_GridGeneration,
            inputAtomStruct=self.protImportPDB.outputPdb,
            radius=37.0,
            spacing=1.0)

        self.launchProtocol(protGrid)
        gridOut = getattr(protGrid, 'outputGrid', None)
        self.assertIsNotNone(gridOut)
        return protGrid

    def _runAutoLigandFind(self, protGrid=None):
        if protGrid == None:
            protAutoLigand = self.newProtocol(
                ProtChemAutoLigand,
                prevGrid=False,
                inputAtomStruct=self.protImportPDB.outputPdb,
                radius=35, spacing=1.0,
                nFillPoints=10)
        else:
            protAutoLigand = self.newProtocol(
                ProtChemAutoLigand,
                prevGrid=True,
                inputGrid=protGrid.outputGrid,
                nFillPoints=10)

        self.launchProtocol(protAutoLigand)
        pdbOut = getattr(protAutoLigand, 'outputStructROIs', None)
        self.assertIsNotNone(pdbOut)

    def testAutoLigand(self):
        self._runAutoLigandFind()

    def testAutoLigandGrid(self):
        protGrid = self._runCreateGrid()
        self._runAutoLigandFind(protGrid)



