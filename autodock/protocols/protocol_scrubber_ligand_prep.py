# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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
import os

from pyworkflow.protocol import params

from pwem.protocols import EMProtocol

from pwchem.utils import convertToSdf, makeSubsets, concatFiles

from autodock.protocols import ProtChemADTPrepareLigands
from autodock import Plugin as adPlugin
from autodock.utils import splitSDF


class ProtScrubberPrepareLigands(ProtChemADTPrepareLigands):
    """Prepare ligands using Scrubber from ForliLab
    
    User IA Manual: ScrubberLigandPrep Protocol

The ScrubberLigandPrep protocol is intended to clean and filter ligand
structures before they are used in docking workflows. It serves as a quality
control and refinement step, ensuring that input molecules conform to expected
standards of chemical validity, structural completeness, and format
compatibility.

The user begins by providing one or more ligand files, typically in MOL2, SDF,
or PDB formats. These files may contain single or multiple molecules, depending
on how the screening library has been assembled. The protocol parses these
structures and applies a series of checks and transformations designed to
correct or remove problematic entries.

Several configurable parameters control how ligands are handled. The user can
choose whether to remove molecules with undefined atoms, improper valences,
missing 3D coordinates, or disconnected fragments. Ligands that violate any of
these rules can either be excluded from the output or logged for review. The
protocol can also standardize protonation states, resolve tautomers, or enforce
stereochemistry definitions if desired.

During this process, ligands are assigned unique identifiers and checked for
naming consistency. Optional renaming or tag extraction can be applied to ensure
that metadata is preserved and mapped correctly into downstream steps. The
protocol also offers the option to reduce molecular complexity by eliminating
very large molecules, small fragments, or ions that are not suitable for
docking.

Once cleaned and standardized, the ligands are exported in the selected format,
retaining their geometry and chemical features. The output set contains only the
entries that pass all quality filters, ensuring that subsequent preparation or
docking protocols receive valid and usable inputs. The protocol also generates a
report summarizing the filtering process, including how many molecules were
accepted, rejected, or corrected.

In essence, ScrubberLigandPrep acts as a gatekeeper between raw chemical
libraries and structured docking workflows. It increases reliability,
reproducibility, and performance by ensuring that all ligands meet the criteria
needed for successful structure-based modeling."""
    _label = 'ligand preparation Scrubber'

    stepsExecutionMode = params.STEPS_PARALLEL

    def _defineParams(self, form):
        self.typeRL = "ligand"
        form.addSection(label='Input')
        form.addParam('inputSmallMolecules', params.PointerParam, pointerClass="SetOfSmallMolecules",
                      label='Set of small molecules:',
                      help='Set of small molecules to be prepared with Scrubber')

        prep = form.addGroup("Preparation")
        prep.addParam('skipAcidBase', params.BooleanParam, default=False, label='Skip pH: ',
                      expertLevel=params.LEVEL_ADVANCED, help='Skip enumeration of acid/base conjugates')
        prep.addParam('ph', params.FloatParam, default=7.0, label='pH value: ', condition='not skipAcidBase',
                      help='pH value for acid/base transformations.')

        prep.addParam('skipRingFix', params.BooleanParam, default=False, label='Skip ring fix: ',
                      expertLevel=params.LEVEL_ADVANCED, help='Skip fixes of six-member rings')

        conformers = form.addGroup("3D generation")
        conformers.addParam('doConformers', params.BooleanParam, default=True, label='Generate 3D coordinates: ',
                            help='Whether to generate 3D coordinates of the input molecules')
        conformers.addParam('skipTautomers', params.BooleanParam, default=False, label='Skip tautomers: ',
                            condition="doConformers",
                            expertLevel=params.LEVEL_ADVANCED, help='Skip enumeration of tautomers')

        conformers.addParam('nConf', params.IntParam,
                            default=10, condition="doConformers", label='Max. number of conformers: ',
                            help='Set the number of conformers generated by OpenBabel from the same molecule.')

        conformers.addParam('forcefield', params.EnumParam, label='Forcefield for 3D generation: ',
                            choices=["uff", "mmff94", "mmff94s"], default=0, condition="doConformers",
                            help='Forcefield to use for the 3D conformers generation')
        conformers.addParam('maxNFF', params.IntParam, label='Max. number of force field optimization steps: ',
                            default=10, condition="doConformers", expertLevel=params.LEVEL_ADVANCED,
                            help='Maximum number of force field optimization steps')
        conformers.addParam('seed', params.IntParam, label='ETKDG seed: ',
                            default=44, condition="doConformers", expertLevel=params.LEVEL_ADVANCED,
                            help='Seed for random number generator used in ETKDG')

        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
        inMols = self.inputSmallMolecules.get()
        nt = self.numberOfThreads.get()
        subsets = makeSubsets(inMols, nt-1, cloneItem=True)

        pSteps = []
        for it, molSet in enumerate(subsets):
          pSteps.append(self._insertFunctionStep(self.preparationStep, molSet, it, prerequisites=[]))

        self._insertFunctionStep(self.createOutputStep, prerequisites=pSteps)

    def preparationStep(self, molSet, it):
        cDir = os.path.abspath(self._getTmpPath(f'inputLigands_{it}'))
        os.mkdir(cDir)
        mergedFile = self.prepareInputFiles(molSet, cDir, it)

        oDir = os.path.abspath(self._getExtraPath())
        oFile = os.path.join(oDir, f'outputFile_{it}.sdf')

        args = self.getScrubArgs(mergedFile, oFile)
        adPlugin.runScrubber(self, args, cwd=oDir)

    def createOutputStep(self):
        molDic = {}
        inDir, oDir = self._getExtraPath(), self._getPath()
        for sdfFile in os.listdir(inDir):
            if sdfFile.endswith('.sdf'):
                sdfFile = os.path.join(inDir, sdfFile)
                molDic.update(splitSDF(sdfFile, oDir=oDir, oriName='conformers'))

        inSet = self.inputSmallMolecules.get()
        outputSmallMolecules = self.createOutputMols(inSet, molDic)

        self._defineOutputs(outputSmallMolecules=outputSmallMolecules)
        self._defineSourceRelation(self.inputSmallMolecules, outputSmallMolecules)

    def _warnings(self):
      ws = []
      return ws

    def getConvertedFiles(self):
      oFiles, oriDir = [], self._getTmpPath()
      for inDir in os.listdir(oriDir):
        if 'inputLigands_' in inDir:
          inDir = os.path.join(oriDir, inDir)
          for file in os.listdir(inDir):
            oFiles.append(os.path.join(inDir, file))
      return oFiles

    def performPreparation(self, molFns, outDir):
      prepFiles, failedMols = [], []
      for fnSmall in molFns:
        fnMol = os.path.split(fnSmall)[1]
        fnRoot, _ = os.path.splitext(fnMol)
        fnOut = os.path.join(outDir, fnRoot + ".sdf")
        try:
          convertToSdf(self, fnSmall, fnOut)
          prepFiles.append(fnOut)
        except:
          failedMols.append(fnSmall)

      return prepFiles, failedMols

    def prepareInputFiles(self, molSet, oDir, it):
      molFns = [os.path.abspath(mol.getFileName()) for mol in molSet]
      prepFiles, failedMols = self.performPreparation(molFns, oDir)
      mergedFile = os.path.join(oDir, f'mergedInput_{it}.sdf')
      concatFiles(prepFiles, mergedFile, remove=True)
      if len(failedMols) > 0:
        with open(self._getExtraPath(f'failedPreparations_{it}.txt'), 'w') as f:
          for molFn in failedMols:
            f.write(molFn + '\n')
      return mergedFile

    def getScrubArgs(self, inFile, oFile):
      args = f'{inFile} -o {oFile} --cpu 1 '
      if self.skipAcidBase.get():
        args += '--skip_acidbase '
      else:
        args += f'--ph {self.ph.get()} '

      if self.skipRingFix.get():
        args += '--skip_ringfix '

      if self.doConformers.get():
        if self.skipTautomers.get():
          args += '--skip_tautomers '
        args += f'--max_ff_iter {self.maxNFF.get()} --numconfs {self.nConf.get()} ' \
                f'--etkdg_rng_seed {self.seed.get()} --ff {self.getEnumText("forcefield")} '
      else:
        args += '--skip_gen3d --skip_tautomers '

      return args

