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
from pwem.protocols import ProtImportPdb, ProtSetFilter
from ..protocols import ProtChemAutoLigand, ProtChemAutodock, \
    ProtChemADTPrepareLigands, ProtChemADTPrepareReceptor
from pwchem.protocols import ProtChemImportSmallMolecules, ProtChemOBabelPrepareLigands


class TestAutoDock(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls.dsLig = DataSet.getDataSet("smallMolecules")

        setupTestProject(cls)
        cls._runImportPDB()
        cls._runImportSmallMols()

        cls._runPrepareLigandsOBabel()
        cls._runPrepareLigandsADT()
        cls._runPrepareReceptorADT()

    @classmethod
    def _runImportSmallMols(cls):
        cls.protImportSmallMols = cls.newProtocol(
            ProtChemImportSmallMolecules,
            filesPath=cls.dsLig.getFile('mol2'))
        cls.launchProtocol(cls.protImportSmallMols)

    @classmethod
    def _runImportPDB(cls):
        cls.protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1,
            pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1_noHETATM.pdb'))
        cls.launchProtocol(cls.protImportPDB)

    @classmethod
    def _runPrepareLigandsOBabel(cls):
        cls.protOBabel = cls.newProtocol(
            ProtChemOBabelPrepareLigands,
            inputType=0, method_charges=0,
            inputSmallMols=cls.protImportSmallMols.outputSmallMolecules,
            doConformers=True, method_conf=0, number_conf=2,
            rmsd_cutoff=0.375)

        cls.launchProtocol(cls.protOBabel)

    @classmethod
    def _runPrepareLigandsADT(cls):
        cls.protPrepareLigand = cls.newProtocol(
            ProtChemADTPrepareLigands,
            inputSmallMols=cls.protImportSmallMols.outputSmallMolecules,
            doConformers=True, method_conf=0, number_conf=2,
            rmsd_cutoff=0.375)

        cls.launchProtocol(cls.protPrepareLigand)

    @classmethod
    def _runPrepareReceptorADT(cls):
        cls.protPrepareReceptor = cls.newProtocol(
            ProtChemADTPrepareReceptor,
            inputStructure=cls.protImportPDB.outputPdb,
            repair=3)

        cls.launchProtocol(cls.protPrepareReceptor)
      
    def _runSetFilter(self, inProt, number, property):
        protFilter = self.newProtocol(
            ProtSetFilter,
            operation=ProtSetFilter.CHOICE_RANKED,
            threshold=number, rankingField=property)
        protFilter.inputSet.set(inProt)
        protFilter.inputSet.setExtended('outputPockets')

        self.launchProtocol(protFilter)
        return protFilter

    def _runAutoLigandFind(self):
        protAutoLigand = self.newProtocol(
            ProtChemAutoLigand,
            inputAtomStruct=self.protPrepareReceptor.outputStructure,
            radius=36.0,
            nFillPoints=10)

        self.launchProtocol(protAutoLigand)
        pdbOut = getattr(protAutoLigand, 'outputPockets', None)
        self.assertIsNotNone(pdbOut)
        return protAutoLigand

    def _runAutoDock(self, pocketsProt=None):
        if pocketsProt == None:
            protAutoDock = self.newProtocol(
                ProtChemAutodock,
                wholeProt=True,
                inputAtomStruct=self.protPrepareReceptor.outputStructure,
                inputLibrary=self.protOBabel.outputSmallMolecules,
                radius=37, gaRun=2,
                numberOfThreads=4)
            self.launchProtocol(protAutoDock)
            smOut = getattr(protAutoDock, 'outputSmallMolecules', None)
            self.assertIsNotNone(smOut)

        else:
            protAutoDock = self.newProtocol(
                ProtChemAutodock,
                wholeProt=False,
                inputPockets=pocketsProt.outputPockets,
                inputLibrary=self.protPrepareLigand.outputSmallMolecules,
                pocketRadiusN=2, gaRun=2,
                numberOfThreads=4)
            self.launchProtocol(protAutoDock)
            smOut = getattr(protAutoDock, 'outputSmallMolecules_1', None)
            self.assertIsNotNone(smOut)

        return protAutoDock

    def testAutoDockPockets(self):
        autoLigProt = self._runAutoLigandFind()
        filterProt = self._runSetFilter(autoLigProt, 2, '_score')
        self._runAutoDock(filterProt)

    def testAutoDockProtein(self):
        self._runAutoDock()






