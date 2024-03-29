# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

# Scipion em imports
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb, ProtSetFilter

# Scipion chem imports
from pwchem.protocols import ProtChemImportSmallMolecules, ProtChemOBabelPrepareLigands
from pwchem.utils import assertHandle

# Plugin imports
from ..protocols import ProtChemADTPrepareReceptor, ProtChemADTPrepareLigands, ProtChemMeekoLigands
from ..protocols import ProtChemAutoLigand, ProtChemAutoSite, ProtChemAutodock
from ..protocols import ProtChemAutodockGPU, ProtChemVinaDocking, ProtChemAutoSiteGenPharmacophore
from ..protocols import AutodockGridGeneration, ProtChemAutodockScore

# Receptor and ligands preparations
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
			HETATM=True, rchains=True, repair=3,
			chain_name='{"model": 0, "chain": "C", "residues": 93}')

		cls.launchProtocol(cls.protPrepareReceptor)

	def test(self):
		self._runPrepareReceptorADT()

		self._waitOutput(self.protPrepareReceptor, 'outputStructure', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(self.protPrepareReceptor, 'outputStructure', None), cwd=self.protPrepareReceptor.getWorkingDir())

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
			doConformers=True, method_conf=0, number_conf=2, rmsd_cutoff=0.375)
		cls.protPrepareLigandADT.inputSmallMolecules.set(cls.protImportSmallMols)
		cls.protPrepareLigandADT.inputSmallMolecules.setExtended('outputSmallMolecules')

		cls.proj.launchProtocol(cls.protPrepareLigandADT, wait=False)


	def test(self):
		self._runPrepareLigandsADT()

		self._waitOutput(self.protPrepareLigandADT, 'outputSmallMolecules', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(self.protPrepareLigandADT, 'outputSmallMolecules', None), cwd=self.protPrepareLigandADT.getWorkingDir())

class TestADMeekoLigands(TestADPrepareLigands):
	@classmethod
	def _runPrepareLigandsMeeko(cls):
		cls.protMeeko = cls.newProtocol(
			ProtChemMeekoLigands)
		cls.protMeeko.inputSmallMolecules.set(cls.protImportSmallMols)
		cls.protMeeko.inputSmallMolecules.setExtended('outputSmallMolecules')

		cls.proj.launchProtocol(cls.protMeeko, wait=False)


	def test(self):
		self._runPrepareLigandsMeeko()

		self._waitOutput(self.protMeeko, 'outputSmallMolecules', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(self.protMeeko, 'outputSmallMolecules', None), cwd=self.protMeeko.getWorkingDir())

# Binding site predictions
class TestAutoLigand(TestADPrepareReceptor):
	@classmethod
	def _runAutoLigandFind(cls):
		protPocketFinder = cls.newProtocol(
			ProtChemAutoLigand,
			inputAtomStruct=cls.protPrepareReceptor.outputStructure,
			radius=24, nFillPoints=10)

		cls.proj.launchProtocol(protPocketFinder, wait=False)
		return protPocketFinder

	def test(self):
		self._runPrepareReceptorADT()
		protAutoLig = self._runAutoLigandFind()
		self._waitOutput(protAutoLig, 'outputStructROIs', sleepTime=5)
		assertHandle(self.assertIsNotNone, getattr(protAutoLig, 'outputStructROIs', None), cwd=protAutoLig.getWorkingDir())

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
		assertHandle(self.assertIsNotNone, getattr(protAutoSite, 'outputStructROIs', None), cwd=protAutoSite.getWorkingDir())

# Docking
class TestAutoDock(TestAutoSite, TestADMeekoLigands):
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
		cls._runPrepareLigandsMeeko()
		cls._waitOutput(cls.protOBabel, 'outputSmallMolecules', sleepTime=5)
		cls._waitOutput(cls.protPrepareLigandADT, 'outputSmallMolecules', sleepTime=5)
		cls._waitOutput(cls.protMeeko, 'outputSmallMolecules', sleepTime=5)
		cls._waitOutput(cls.protPrepareReceptor, 'outputStructure', sleepTime=5)

	@classmethod
	def _runPrepareLigandsOBabel(cls):
		cls.protOBabel = cls.newProtocol(
			ProtChemOBabelPrepareLigands,
			inputType=0, method_charges=0,
			inputSmallMolecules=cls.protImportSmallMols.outputSmallMolecules,
			doConformers=True, method_conf=0, number_conf=2, rmsd_cutoff=0.375)

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
				inputAtomStruct=self.protPrepareReceptor.outputStructure,
				inputSmallMolecules=self.protOBabel.outputSmallMolecules,
				fromReceptor=0, radius=24, nRuns=2,
				numberOfThreads=4)
			self.proj.launchProtocol(protAutoDock, wait=False)

		else:
			protAutoDock = self.newProtocol(
				ProtChemAutodock,
				inputStructROIs=pocketsProt.outputStructROIs,
				inputSmallMolecules=self.protPrepareLigandADT.outputSmallMolecules,
				fromReceptor=1, pocketRadiusN=1.2, nRuns=2,
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
		assertHandle(self.assertIsNotNone, getattr(protAutoDock1, 'outputSmallMolecules', None), cwd=protAutoDock1.getWorkingDir())
		self._waitOutput(protAutoDock2, 'outputSmallMolecules', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(protAutoDock2, 'outputSmallMolecules', None), cwd=protAutoDock2.getWorkingDir())

class TestAutoDockGPU(TestAutoDock):
	def _runAutoDock(self, pocketsProt=None):
		if pocketsProt == None:
			protAutoDock = self.newProtocol(
				ProtChemAutodockGPU,
				inputAtomStruct=self.protPrepareReceptor.outputStructure,
				inputSmallMolecules=self.protMeeko.outputSmallMolecules,
				fromReceptor=0, radius=24, nRuns=2,
				numberOfThreads=4, useGpu=True)
			self.proj.launchProtocol(protAutoDock, wait=False)
		else:
			protAutoDock = self.newProtocol(
				ProtChemAutodockGPU,
				inputStructROIs=pocketsProt.outputStructROIs,
				inputSmallMolecules=self.protMeeko.outputSmallMolecules,
				fromReceptor=1, pocketRadiusN=1.2, nRuns=2,
				numberOfThreads=4, useGpu=True)
			self.proj.launchProtocol(protAutoDock, wait=False)

		return protAutoDock

	def test(self):
		print('Docking with autodock-GPU in the whole protein')
		protAutoDock1 = self._runAutoDock()

		print('Docking with autodock-GPU in predicted pockets')
		protAutoLig = self._runAutoSiteFind()
		self._waitOutput(protAutoLig, 'outputStructROIs', sleepTime=5)
		self._runSetFilter(inProt=protAutoLig, number=2, property='_score')
		self._waitOutput(self.protFilter, 'outputStructROIs', sleepTime=5)

		protAutoDock2 = self._runAutoDock(self.protFilter)

		self._waitOutput(protAutoDock1, 'outputSmallMolecules', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(protAutoDock1, 'outputSmallMolecules', None), cwd=protAutoDock1.getWorkingDir())
		self._waitOutput(protAutoDock2, 'outputSmallMolecules', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(protAutoDock2, 'outputSmallMolecules', None), cwd=protAutoDock2.getWorkingDir())

class TestVina(TestAutoDock):
	def _runVina(self, pocketsProt=None):
		if pocketsProt == None:
			protAutoDock = self.newProtocol(
				ProtChemVinaDocking,
				fromReceptor=0,
				inputAtomStruct=self.protPrepareReceptor.outputStructure,
				inputSmallMolecules=self.protOBabel.outputSmallMolecules,
				radius=24, nRuns=2, exhaust=3,
				numberOfThreads=4)
			self.proj.launchProtocol(protAutoDock, wait=False)
		else:
			protAutoDock = self.newProtocol(
				ProtChemVinaDocking,
				fromReceptor=1,
				inputStructROIs=pocketsProt.outputStructROIs,
				inputSmallMolecules=self.protPrepareLigandADT.outputSmallMolecules,
				pocketRadiusN=1.2, nRuns=2, exhaust=3,
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
		assertHandle(self.assertIsNotNone, getattr(protVina1, 'outputSmallMolecules', None), cwd=protVina1.getWorkingDir())
		self._waitOutput(protVina2, 'outputSmallMolecules', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(protVina2, 'outputSmallMolecules', None), cwd=protVina2.getWorkingDir())

# Pharmacophore generation
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
		assertHandle(self.assertIsNotNone, getattr(protPharm, 'outputPharmacophore', None), cwd=protPharm.getWorkingDir())

# Grid generation
class TestGridADT(TestADPrepareReceptor):
	def _runCreateGrid(self, spacing, radius):
		protGrid = self.newProtocol(
			AutodockGridGeneration,
			inputAtomStruct=self.protImportPDB.outputPdb,
			radius=radius,
			spacing=spacing)

		self.launchProtocol(protGrid)
		return protGrid

	def test(self):
		spacing, radius = 1.0, 37.0
		protGrid = self._runCreateGrid(spacing, radius)
		assertHandle(self.assertIsNotNone, getattr(protGrid, 'outputGrid', None), cwd=protGrid.getWorkingDir())

class TestAutoDockScoring(TestAutoDockGPU):
	def _runScoring(self, inputProt):
		protScore = self.newProtocol(ProtChemAutodockScore, radius=20)
		protScore.inputSmallMolecules.append(inputProt)
		protScore.inputSmallMolecules[-1].setExtended('outputSmallMolecules')

		self.proj.launchProtocol(protScore, wait=True)

		return protScore

	def test(self):
		protAutoDock1 = self._runAutoDock()
		self._waitOutput(protAutoDock1, 'outputSmallMolecules', sleepTime=10)
		protScore = self._runScoring(protAutoDock1)
		assertHandle(self.assertIsNotNone, getattr(protScore, 'outputSmallMolecules', None), cwd=protScore.getWorkingDir())
