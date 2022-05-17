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
from ..protocols import ProtChemAutoLigand, ProtChemVina, \
    ProtChemADTPrepareLigands, ProtChemADTPrepareReceptor
from pwchem.protocols import ProtChemImportSmallMolecules, ProtChemOBabelPrepareLigands


class TestVina(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls.dsLig = DataSet.getDataSet("smallMolecules")

        setupTestProject(cls)
        cls._runImportPDB()
        cls._runImportSmallMols()
        cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)
        cls._waitOutput(cls.protImportSmallMols, 'outputSmallMolecules', sleepTime=5)

        cls._runPrepareLigandsOBabel()
        cls._runPrepareLigandsADT()
        cls._runPrepareReceptorADT()
        cls._waitOutput(cls.protOBabel, 'outputSmallMolecules', sleepTime=5)
        cls._waitOutput(cls.protPrepareLigandADT, 'outputSmallMolecules', sleepTime=5)
        cls._waitOutput(cls.protPrepareReceptor, 'outputStructure', sleepTime=5)

    @classmethod
    def _runImportSmallMols(cls):
        cls.protImportSmallMols = cls.newProtocol(
            ProtChemImportSmallMolecules,
            filesPath=cls.dsLig.getFile('mol2'))
        cls.proj.launchProtocol(cls.protImportSmallMols, wait=False)

    @classmethod
    def _runImportPDB(cls):
        cls.protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0,
            pdbId='4erf')
        cls.proj.launchProtocol(cls.protImportPDB, wait=False)

    @classmethod
    def _runPrepareLigandsOBabel(cls):
        cls.protOBabel = cls.newProtocol(
            ProtChemOBabelPrepareLigands,
            inputType=0, method_charges=0,
            inputSmallMols=cls.protImportSmallMols.outputSmallMolecules,
            doConformers=True, method_conf=0, number_conf=2,
            rmsd_cutoff=0.375)

        cls.proj.launchProtocol(cls.protOBabel, wait=False)

    @classmethod
    def _runPrepareLigandsADT(cls):
        cls.protPrepareLigandADT = cls.newProtocol(
            ProtChemADTPrepareLigands,
            doConformers=True, method_conf=0, number_conf=2,
            rmsd_cutoff=0.375)
        cls.protPrepareLigandADT.inputSmallMols.set(cls.protImportSmallMols)
        cls.protPrepareLigandADT.inputSmallMols.setExtended('outputSmallMolecules')

        cls.proj.launchProtocol(cls.protPrepareLigandADT, wait=False)

    @classmethod
    def _runPrepareReceptorADT(cls):
        cls.protPrepareReceptor = cls.newProtocol(
            ProtChemADTPrepareReceptor,
            inputAtomStruct=cls.protImportPDB.outputPdb,
            HETATM=True, rchains=True,
            chain_name='{"model": 0, "chain": "C", "residues": 93}',
            repair=3)

        cls.launchProtocol(cls.protPrepareReceptor)

    @classmethod
    def _runAutoLigandFind(cls):
        protPocketFinder = cls.newProtocol(
            ProtChemAutoLigand,
            inputAtomStruct=cls.protPrepareReceptor.outputStructure,
            radius=24,
            nFillPoints=10)

        cls.proj.launchProtocol(protPocketFinder, wait=False)
        return protPocketFinder

    @classmethod
    def _runSetFilter(cls, inProt, number, property):
        cls.protFilter = cls.newProtocol(
            ProtSetFilter,
            operation=ProtSetFilter.CHOICE_RANKED,
            threshold=number, rankingField=property)
        cls.protFilter.inputSet.set(inProt)
        cls.protFilter.inputSet.setExtended('outputPockets')

        cls.proj.launchProtocol(cls.protFilter, wait=False)
        return cls.protFilter


    def _runVina(self, pocketsProt=None):
        if pocketsProt == None:
            protAutoDock = self.newProtocol(
                ProtChemVina,
                wholeProt=True,
                inputAtomStruct=self.protPrepareReceptor.outputStructure,
                inputLibrary=self.protOBabel.outputSmallMolecules,
                radius=24, nPos=2,
                numberOfThreads=8)
            self.proj.launchProtocol(protAutoDock, wait=False)

        else:
            protAutoDock = self.newProtocol(
                ProtChemVina,
                wholeProt=False,
                inputPockets=pocketsProt.outputPockets,
                inputLibrary=self.protPrepareLigandADT.outputSmallMolecules,
                pocketRadiusN=5, nPos=2,
                mergeOutput=True,
                numberOfThreads=4)
            self.proj.launchProtocol(protAutoDock, wait=False)

        return protAutoDock

    def testVina(self):
        print('Docking with autodock vina in the whole protein')
        protVina1 = self._runVina()

        print('Docking with autodock vina in predicted pockets')
        protAutoLig = self._runAutoLigandFind()
        self._waitOutput(protAutoLig, 'outputPockets', sleepTime=5)
        self._runSetFilter(inProt=protAutoLig, number=2, property='_score')
        self._waitOutput(self.protFilter, 'outputPockets', sleepTime=5)

        protVina2 = self._runVina(self.protFilter)

        self._waitOutput(protVina1, 'outputSmallMolecules', sleepTime=10)
        self.assertIsNotNone(getattr(protVina1, 'outputSmallMolecules', None))
        self._waitOutput(protVina2, 'outputSmallMolecules', sleepTime=10)
        self.assertIsNotNone(getattr(protVina2, 'outputSmallMolecules', None))


