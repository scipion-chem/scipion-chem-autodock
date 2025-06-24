# **************************************************************************
# *
# * Authors:    Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

from pyworkflow.protocol.params import PointerParam, BooleanParam, EnumParam, IntParam, FloatParam, LEVEL_ADVANCED

from pwchem.utils import makeSubsets
from pwchem.constants import RDKIT_DIC

from autodock import Plugin
from autodock.protocols import ProtChemADTPrepareLigands

scriptName = 'meeko_preparation.py'

class ProtChemMeekoLigands(ProtChemADTPrepareLigands):
    """Prepare ligands using Meeko from Autodock 
    
    User IA Manual: MeekoLigands Protocol

The MeekoLigands protocol is used to prepare ligand molecules for docking by
converting them from common chemical structure formats, such as SDF or MOL2,
into the PDBQT format required by AutoDock-based engines. This preparation
includes steps such as 3D coordinate generation, torsion detection, charge
assignment, and atom typing, ensuring the ligand is properly configured for
structure-based virtual screening.

To begin, the user must provide a file containing one or more ligands. These may
include two-dimensional or three-dimensional structures. If only 2D coordinates
are present, the protocol attempts to generate a 3D conformation automatically.
Proper 3D geometry is essential for successful docking, so ligand files must be
checked to ensure their suitability prior to execution.

The protocol handles key steps such as the assignment of Gasteiger partial
charges, the addition of hydrogens at physiological pH, and the detection of
rotatable bonds. It offers the option to enable or disable automatic torsion
assignment, which can be useful for molecules that require fixed conformations
or have specific torsional constraints. If needed, users may instruct the
protocol to preserve original atom names or chemical features that are sensitive
to stereochemistry and tautomeric forms.

Each processed ligand is written as an individual PDBQT file, with unique
identifiers to maintain traceability throughout the workflow. These output files
are fully compatible with downstream protocols such as AutoDock-GPU, AutoDock
Vina, or EncoderDockScoring. All relevant conformational and chemical
information is embedded into the output, allowing immediate use for docking or
rescoring.

Errors encountered during the conversion process, such as missing 3D
coordinates, undefined stereochemistry, or failures in charge assignment, are
reported clearly in the log. This ensures that problematic molecules can be
quickly identified and corrected.

In summary, this protocol provides an automated and reliable interface to Meeko
for the preparation of ligands in docking-ready format. It integrates smoothly
into Scipion-Chem workflows and serves as the standard entry point for molecular
screening campaigns involving structure-based methods."""
    _label = 'meeko ligand preparation'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                      label='Set of small molecules:', allowsNull=False,
                      help='Input small molecules to be prepared with Meeko')

        form.addParam('hydrate', BooleanParam, default=False,
                      label='Hydrate molecule:', expertLevel=LEVEL_ADVANCED,
                      help='Hydrate molecule for hydrated-docking')

        conformers = form.addGroup("Conformers generation", condition="not hydrate")
        conformers.addParam('doConformers', BooleanParam, default=False, condition="not hydrate",
                            label='Do you want to generate conformers? ', allowsNull=False,
                            help='You can produce conformers of the ligand in order to do a better rigid docking')
        conformers.addParam('method_conf', EnumParam, condition="not hydrate and doConformers",
                            choices=["OpenBabel Genetic Algorithm", "OpenBabel Confab"], default=0,
                            label='Method of conformers generation',
                            help='Method of conformers generation. If Confab fails due to the impossibility '
                                 'of assigning a force fields (there is a possibility that it may occur), you should'
                                 'use Genetic Algorithm generator ')

        conformers.addParam('number_conf', IntParam,
                            default=10, condition="not hydrate and doConformers",
                            label='Max. number of conformers:',
                            help='Set the number of conformers generated by OpenBabel from the same molecule.')

        conformers.addParam('rmsd_cutoff', FloatParam, condition="not hydrate and method_conf != 0 and doConformers",
                            default=0.5, label='RMSD cutoff:',
                            help='Set the number of conformers generated by OpenBabel from the same molecule.',
                            expertLevel=LEVEL_ADVANCED)

        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
      inMols = self.inputSmallMolecules.get()
      nt = self.numberOfThreads.get()
      subsets = makeSubsets(inMols, nt-1, cloneItem=True)

      pSteps, cSteps = [], []
      for it, molSet in enumerate(subsets):
        pSteps.append(self._insertFunctionStep('preparationStep', molSet, it, prerequisites=[]))
        if self.doConformers.get() and not self.hydrate.get():
          cSteps.append(self._insertFunctionStep('conformerGenerationStep', it, prerequisites=pSteps[-1]))

      self._insertFunctionStep('createOutputStep', prerequisites=pSteps+cSteps)

    def preparationStep(self, molSet, it):
        molFiles = [mol.getFileName() for mol in molSet]
        self.performMeekoPrep(molFiles, it)

    def performMeekoPrep(self, molFiles, it):
      oDir = self.getPreparedDirPath(it)
      if not os.path.exists(oDir):
        os.makedirs(oDir)
      paramsFile = self.writeParamsFile(molFiles, oDir, it)
      Plugin.runScript(self, scriptName, paramsFile, RDKIT_DIC)

    def writeParamsFile(self, molFiles, oDir, it=None):
        paramsFile = os.path.abspath(self._getExtraPath('inputParams.txt'))
        if it is not None:
          paramsFile = paramsFile.replace('.txt', f'_{it}.txt')

        with open(paramsFile, 'w') as f:
          f.write(f'ligandFiles:: {" ".join(molFiles)}\n')
          f.write(f'hydrate:: {self.hydrate.get()}\n')

          f.write(f'outDir:: {oDir}\n')

        return paramsFile

    def _warnings(self):
        return []
