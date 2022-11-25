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

from pwchem.protocols import ProtChemImportSmallMolecules, ProtChemOBabelPrepareLigands

from ..protocols import *


class TestADPrepareReceptor(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        setupTestProject(cls)
        cls._runImportPDB()
        cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)

    @classmethod
    def _runImportPDB(cls):
        cls.protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0, pdbId='4erf')
        cls.proj.launchProtocol(cls.protImportPDB, wait=False)


    @classmethod
    def _runPrepareReceptorADT(cls):
        cls.protPrepareReceptor = cls.newProtocol(
            ProtChemADTPrepareReceptor,
            inputAtomStruct=cls.protImportPDB.outputPdb,
            HETATM=True, rchains=True,
            chain_name='{"model": 0, "chain": "C", "residues": 93}',
            repair=3)

        cls.launchProtocol(cls.protPrepareReceptor)

    def test(self):
        self._runPrepareReceptorADT()

        self._waitOutput(self.protPrepareReceptor, 'outputStructure', sleepTime=10)
        self.assertIsNotNone(getattr(self.protPrepareReceptor, 'outputStructure', None))


class TestAutoLigand(TestADPrepareReceptor):
    @classmethod
    def _runAutoLigandFind(cls):
        protPocketFinder = cls.newProtocol(
            ProtChemAutoLigand,
            inputAtomStruct=cls.protPrepareReceptor.outputStructure,
            radius=24,
            nFillPoints=10)

        cls.proj.launchProtocol(protPocketFinder, wait=False)
        return protPocketFinder

    def test(self):
        self._runPrepareReceptorADT()
        protAutoLig = self._runAutoLigandFind()
        self._waitOutput(protAutoLig, 'outputStructROIs', sleepTime=5)
        self.assertIsNotNone(getattr(protAutoLig, 'outputStructROIs', None))


class TestAutoSite(TestADPrepareReceptor):
    @classmethod
    def _runAutoSiteFind(cls):
        protPocketFinder = cls.newProtocol(
            ProtChemAutoSite,
            inputAtomStruct=cls.protPrepareReceptor.outputStructure)

        cls.proj.launchProtocol(protPocketFinder, wait=False)
        return protPocketFinder

    def test(self):
        self._runPrepareReceptorADT()
        protAutoSite = self._runAutoSiteFind()
        self._waitOutput(protAutoSite, 'outputStructROIs', sleepTime=5)
        self.assertIsNotNone(getattr(protAutoSite, 'outputStructROIs', None))


class TestADPrepareLigands(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.dsLig = DataSet.getDataSet("smallMolecules")
        setupTestProject(cls)

        cls._runImportSmallMols()
        cls._waitOutput(cls.protImportSmallMols, 'outputSmallMolecules', sleepTime=5)

    @classmethod
    def _runImportSmallMols(cls):
        cls.protImportSmallMols = cls.newProtocol(
            ProtChemImportSmallMolecules,
            filesPath=cls.dsLig.getFile('mol2'))
        cls.proj.launchProtocol(cls.protImportSmallMols, wait=False)

    @classmethod
    def _runPrepareLigandsADT(cls):
        cls.protPrepareLigandADT = cls.newProtocol(
            ProtChemADTPrepareLigands,
            doConformers=True, method_conf=0, number_conf=2,
            rmsd_cutoff=0.375)
        cls.protPrepareLigandADT.inputSmallMols.set(cls.protImportSmallMols)
        cls.protPrepareLigandADT.inputSmallMols.setExtended('outputSmallMolecules')

        cls.proj.launchProtocol(cls.protPrepareLigandADT, wait=False)


    def test(self):
        self._runPrepareLigandsADT()

        self._waitOutput(self.protPrepareLigandADT, 'outputSmallMolecules', sleepTime=10)
        self.assertIsNotNone(getattr(self.protPrepareLigandADT, 'outputSmallMolecules', None))


class TestAutoDock(TestAutoSite, TestADPrepareLigands):
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
    def _runPrepareLigandsOBabel(cls):
        cls.protOBabel = cls.newProtocol(
            ProtChemOBabelPrepareLigands,
            inputType=0, method_charges=0,
            inputSmallMols=cls.protImportSmallMols.outputSmallMolecules,
            doConformers=True, method_conf=0, number_conf=2,
            rmsd_cutoff=0.375)

        cls.proj.launchProtocol(cls.protOBabel, wait=False)

    @classmethod
    def _runSetFilter(cls, inProt, number, property):
        cls.protFilter = cls.newProtocol(
            ProtSetFilter,
            operation=ProtSetFilter.CHOICE_RANKED,
            threshold=number, rankingField=property)
        cls.protFilter.inputSet.set(inProt)
        cls.protFilter.inputSet.setExtended('outputStructROIs')

        cls.proj.launchProtocol(cls.protFilter, wait=False)
        return cls.protFilter


    def _runAutoDock(self, pocketsProt=None):
        if pocketsProt == None:
            protAutoDock = self.newProtocol(
                ProtChemAutodock,
                wholeProt=True,
                inputAtomStruct=self.protPrepareReceptor.outputStructure,
                inputLibrary=self.protOBabel.outputSmallMolecules,
                radius=24, gaRun=2,
                numberOfThreads=8)
            self.proj.launchProtocol(protAutoDock, wait=False)

        else:
            protAutoDock = self.newProtocol(
                ProtChemAutodock,
                wholeProt=False,
                inputStructROIs=pocketsProt.outputStructROIs,
                inputLibrary=self.protPrepareLigandADT.outputSmallMolecules,
                pocketRadiusN=5, gaRun=2,
                mergeOutput=True,
                numberOfThreads=4)
            self.proj.launchProtocol(protAutoDock, wait=False)

        return protAutoDock

    def test(self):
        print('Docking with autodock in the whole protein')
        protAutoDock1 = self._runAutoDock()

        print('Docking with autodock in predicted pockets')
        protAutoLig = self._runAutoSiteFind()
        self._waitOutput(protAutoLig, 'outputStructROIs', sleepTime=5)
        self._runSetFilter(inProt=protAutoLig, number=2, property='_score')
        self._waitOutput(self.protFilter, 'outputStructROIs', sleepTime=5)

        protAutoDock2 = self._runAutoDock(self.protFilter)

        self._waitOutput(protAutoDock1, 'outputSmallMolecules', sleepTime=10)
        self.assertIsNotNone(getattr(protAutoDock1, 'outputSmallMolecules', None))
        self._waitOutput(protAutoDock2, 'outputSmallMolecules', sleepTime=10)
        self.assertIsNotNone(getattr(protAutoDock2, 'outputSmallMolecules', None))


class TestVina(TestAutoDock):

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
                inputStructROIs=pocketsProt.outputStructROIs,
                inputLibrary=self.protPrepareLigandADT.outputSmallMolecules,
                pocketRadiusN=5, nPos=2,
                mergeOutput=True,
                numberOfThreads=4)
            self.proj.launchProtocol(protAutoDock, wait=False)

        return protAutoDock

    def test(self):
        print('Docking with autodock vina in the whole protein')
        protVina1 = self._runVina()

        print('Docking with autodock vina in predicted pockets')
        protAutoLig = self._runAutoSiteFind()
        self._waitOutput(protAutoLig, 'outputStructROIs', sleepTime=5)
        self._runSetFilter(inProt=protAutoLig, number=2, property='_score')
        self._waitOutput(self.protFilter, 'outputStructROIs', sleepTime=5)

        protVina2 = self._runVina(self.protFilter)

        self._waitOutput(protVina1, 'outputSmallMolecules', sleepTime=10)
        self.assertIsNotNone(getattr(protVina1, 'outputSmallMolecules', None))
        self._waitOutput(protVina2, 'outputSmallMolecules', sleepTime=10)
        self.assertIsNotNone(getattr(protVina2, 'outputSmallMolecules', None))


class TestAutoSitePharmacophore(TestAutoSite):
    @classmethod
    def _runAutoSitePharm(cls, protAutoSite):
        protPocketFinder = cls.newProtocol(
            ProtChemAutoSiteGenPharmacophore,
            inputStructROIs=protAutoSite.outputStructROIs,
            inputStructROISelect='Structural ROI 1, AutoSite class')

        cls.proj.launchProtocol(protPocketFinder, wait=False)
        return protPocketFinder

    def test(self):
        self._runPrepareReceptorADT()
        protAutoSite = self._runAutoSiteFind()
        self._waitOutput(protAutoSite, 'outputStructROIs', sleepTime=5)
        protPharm = self._runAutoSitePharm(protAutoSite)
        self._waitOutput(protPharm, 'outputPharmacophore', sleepTime=5)
        self.assertIsNotNone(getattr(protPharm, 'outputPharmacophore', None))




class TestGridADT(TestADPrepareReceptor):

    def _runCreateGrid(self, spacing, radius):
        protGrid = self.newProtocol(
            Autodock_GridGeneration,
            inputAtomStruct=self.protImportPDB.outputPdb,
            radius=radius,
            spacing=spacing)

        self.launchProtocol(protGrid)
        return protGrid

    def test(self):
        spacing, radius = 1.0, 37.0
        protGrid = self._runCreateGrid(spacing, radius)
        self.assertIsNotNone(getattr(protGrid, 'outputGrid', None))
