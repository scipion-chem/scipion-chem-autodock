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
import shutil
import os

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, BooleanParam, EnumParam, IntParam, FloatParam, LEVEL_ADVANCED
import pyworkflow.object as pwobj

from pwchem.utils import runOpenBabel, splitConformerFile, appendToConformersFile
from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pwchem.constants import RDKIT_DIC

from autodock import Plugin

scriptName = 'meeko_preparation.py'

class ProtChemMeekoLigands(EMProtocol):
    """Prepare ligands using Meeko from Autodock """
    _label = 'meeko ligand preparation'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                      label='Set of small molecules:', allowsNull=False,
                      help='Input small molecules to be prepared with Meeko')

        form.addParam('keepHs', BooleanParam, default=False,
                      label='Keep non-polar hydrogens:',
                      help='Keep non-polar hydrogens:')

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

    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('preparationStep')
        if self.doConformers.get() and not self.hydrate.get():
            self._insertFunctionStep('conformer_generation')
        self._insertFunctionStep('createOutput')

    def preparationStep(self):
        mols = self.inputSmallMolecules.get()
        paramsFile = self.writeParamsFile(mols)

        Plugin.runScript(self, scriptName, paramsFile, envDict=RDKIT_DIC, cwd=self._getExtraPath())

    def conformer_generation(self):
      """ Generate a number of conformers of the same small molecule in pdbqt format with
          openbabel using two different algorithm
      """
      with open(self._getPath('meeko_files.txt')) as fIn:
          prepFiles = fIn.read().split()

      for file in prepFiles:
        fnRoot = os.path.split(file)[1].split(".pdbqt")[0]  # ID or filename without _meeko.pdbqt

        if self.method_conf.get() == 0:  # Genetic algorithm
          args = " %s --conformer --nconf %s --score rmsd --writeconformers -O %s_conformers.pdbqt" % \
                 (os.path.abspath(file), str(self.number_conf.get() - 1), fnRoot)
        else:  # confab
          args = " %s --confab --original --verbose --conf %s --rcutoff %s -O %s_conformers.pdbqt" % \
                 (os.path.abspath(file), str(self.number_conf.get() - 1), str(self.rmsd_cutoff.get()), fnRoot)

        runOpenBabel(protocol=self, args=args, cwd=os.path.abspath(self._getExtraPath()))

    def createOutput(self):
        with open(self._getPath('meeko_files.txt')) as fIn:
          prepFiles = fIn.read().split()

        outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='')
        for file in prepFiles:
            fnRoot = os.path.split(file)[1].split('.pdbqt')[0]
            if self.doConformers.get():
                outDir = self._getExtraPath(fnRoot)
                os.mkdir(outDir)
                firstConfFile = self._getTmpPath('{}-{}.pdbqt'.format(fnRoot, 1))
                shutil.copy(file, firstConfFile)
                confFile = self._getExtraPath("{}_conformers.pdbqt".format(fnRoot))
                confFile = appendToConformersFile(confFile, firstConfFile, beginning=True)
                confDir = splitConformerFile(confFile, outDir=outDir)
                for molFile in os.listdir(confDir):
                    molFile = os.path.join(confDir, molFile)
                    confId = molFile.split('-')[-1].split('.')[0]

                    newSmallMol = SmallMolecule(smallMolFilename=molFile, type='AutoDock')
                    newSmallMol.setMolName(fnRoot)
                    newSmallMol._ConformersFile = pwobj.String(confFile)
                    newSmallMol.setConfId(confId)
                    outputSmallMolecules.append(newSmallMol)
            else:
                newSmallMol = SmallMolecule(smallMolFilename=file, type='AutoDock')
                newSmallMol.setMolName(fnRoot)
                outputSmallMolecules.append(newSmallMol)

        self._defineOutputs(outputSmallMolecules=outputSmallMolecules)
        self._defineSourceRelation(self.inputSmallMolecules, outputSmallMolecules)


    def writeParamsFile(self, molsScipion):
        paramsFile = os.path.abspath(self._getExtraPath('inputParams.txt'))

        molFiles = []
        f = open(paramsFile, 'w')
        for mol in molsScipion:
            molFiles.append(os.path.abspath(mol.getFileName()))

        f.write('ligandFiles:: {}\n'.format(' '.join(molFiles)))
        f.write('keepNonPolar:: {}\n'.format(self.keepHs.get()))
        f.write('hydrate:: {}\n'.format(self.hydrate.get()))

        f.write('outDir:: {}\n'.format(os.path.abspath(self._getPath())))

        return paramsFile