# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
import shutil, os, glob

from pyworkflow.protocol.params import PointerParam, BooleanParam, EnumParam, IntParam, FloatParam, LEVEL_ADVANCED
from pyworkflow.utils.path import createLink
import pyworkflow.object as pwobj

from pwchem.utils import runOpenBabel, splitConformerFile, appendToConformersFile, makeSubsets, getBaseName

from autodock.protocols.protocol_preparation_receptor import ProtChemADTPrepare


class ProtChemADTPrepareLigands(ProtChemADTPrepare):
    """Prepare ligands using Autodocking Tools from MGL
    
    User IA Manual: PreparationLigands Protocol

The PreparationLigands protocol is used to convert ligand structures into a
docking-compatible format, typically PDBQT, ensuring they are properly prepared
for use with AutoDock-based protocols. It is designed for ligands that have
already been processed externally and are stored in formats such as MOL2, PDB,
or SDF, and need to undergo conversion or refinement within the Scipion-Chem
workflow.

To start, the user must provide one or more ligand files that contain the
chemical structures to be prepared. These can be single-molecule files or
multi-ligand files, depending on the source and use case. The protocol will
apply a series of steps to each input ligand to ensure that it meets the
structural and formatting requirements of downstream docking engines.

Among the key operations performed are the generation or optimization of 3D
coordinates, detection of rotatable bonds, addition of hydrogens, and the
calculation of partial charges. The protocol uses Meeko internally to handle
this conversion process, leveraging RDKit and OpenBabel-based functions to
assign bond orders and perceive chemical features.

The user can choose whether to apply automatic protonation at a standard pH,
and whether to include all tautomers or stereoisomers when ambiguous structures
are detected. It is also possible to restrict torsional flexibility by
preventing the protocol from defining certain bonds as rotatable. These options
are useful for controlling the conformational space that will be explored during
docking.

Each ligand is then exported in PDBQT format, ready to be used by protocols such
as AutoDock-GPU or Vina. The resulting files preserve chemical integrity,
torsional flexibility, and partial charge information. Any issues encountered
during the preparation process, such as invalid valences or missing coordinates,
are reported in the output log to facilitate manual correction.

In essence, this protocol serves as a bridge between raw ligand structures and
the docking engine, offering a reproducible and customizable pipeline for
ligand standardization. It ensures that all ligands entering virtual screening
campaigns are chemically sound, geometrically valid, and formatted for
compatibility with the AutoDock family of tools.
    """
    _label = 'ligand preparation ADT'
    _program = ""

    def _defineParams(self, form):
        self.typeRL = "ligand"
        form.addSection(label='Input')
        form.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                      label='Set of small molecules:', allowsNull=False,
                      help='It must be in pdb or mol2 format, you may use Schrodinger convert to change it')
        ProtChemADTPrepare._defineParamsBasic(self, form)

        conformers = form.addGroup("Conformers generation")
        conformers.addParam('doConformers', BooleanParam, default=False,
                            label='Do you want to generate conformers? ', allowsNull=False,
                            help='You can produce conformers of the ligand in order to do a better rigid docking')
        conformers.addParam('method_conf', EnumParam, condition="doConformers",
                            choices=["OpenBabel Genetic Algorithm", "OpenBabel Confab"], default=0,
                            label='Method of conformers generation',
                            help='Method of conformers generation. If Confab fails due to the impossibility '
                                 'of assigning a force fields (there is a possibility that it may occur), you should'
                                 'use Genetic Algorithm generator ')

        conformers.addParam('number_conf', IntParam,
                            default=10, condition="doConformers", label='Max. number of conformers:',
                            help='Set the number of conformers generated by OpenBabel from the same molecule.')

        conformers.addParam('rmsd_cutoff', FloatParam, condition="method_conf != 0 and doConformers",
                            default=0.5, label='RMSD cutoff:',
                            help='Set the number of conformers generated by OpenBabel from the same molecule.',
                            expertLevel=LEVEL_ADVANCED)

        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
        inMols = self.inputSmallMolecules.get()
        nt = self.numberOfThreads.get()
        subsets = makeSubsets(inMols, nt - 1, cloneItem=True)

        pSteps, cSteps = [], []
        for it, molSet in enumerate(subsets):
          pSteps.append(self._insertFunctionStep('preparationStep', molSet, it, prerequisites=[]))
          if self.doConformers.get():
            cSteps.append(self._insertFunctionStep('conformerGenerationStep', it, prerequisites=pSteps[-1]))

        self._insertFunctionStep('createOutputStep', prerequisites=pSteps + cSteps)

    def preparationStep(self, molSet, it):
        molFns = [os.path.abspath(mol.getFileName()) for mol in molSet]
        failedMols = self.performPreparation(molFns, it)

        if len(failedMols) > 0:
          with open(os.path.abspath(self._getExtraPath(f'failedPreparations_{it}.txt')), 'w') as f:
            for molFn in failedMols:
              f.write(molFn + '\n')

    def conformerGenerationStep(self, it):
      """ Generate a number of conformers of the same small molecule in pdbqt format with
          openbabel using two different algorithm
      """
      inDir = self.getPreparedDirPath(it)
      molFns = [molFile for molFile in glob.glob(os.path.join(inDir, '*.pdbqt'))]
      failedMols = self.performConfGeneration(molFns, it)

      if len(failedMols) > 0:
        with open(os.path.abspath(self._getExtraPath(f'failedConfomerGeneration_{it}.txt')), 'w') as f:
          for molFn in failedMols:
            f.write(molFn + '\n')

    def createOutputStep(self):
      outMolDic = {}
      outDir, outConfDir = self._getPath('outputLigands'), self._getExtraPath()
      for inDir in self.getPreparedDirs():
        for file in os.listdir(inDir):
          if file.endswith('.pdbqt') and 'conformers.pdbqt' not in file:
            file = os.path.join(inDir, file)
            outMolDic.update(self.indOutputCreation(file, outDir, outConfDir=outConfDir))

      outputSmallMolecules = self.createOutputMols(self.inputSmallMolecules.get(), outMolDic)
      self._defineOutputs(outputSmallMolecules=outputSmallMolecules)
      self._defineSourceRelation(self.inputSmallMolecules, outputSmallMolecules)

    #################### MAIN FUNCTIONS ################

    def performPreparation(self, molFns, it):
      oDir = self.getPreparedDirPath(it)
      if not os.path.exists(oDir):
        os.makedirs(oDir)

      failedMols = []
      for fnSmall in molFns:
        fnMol = os.path.split(fnSmall)[1]
        fnRoot, ext = os.path.splitext(fnMol)
        fnOut = os.path.join(oDir, fnRoot + ".pdbqt")
        try:
          if ext == '.sdf' and (self.repair.get() == 3 or self.repair.get() == 1):
            # AUTODOCK: Cannot handle add hydrogens to files coming from 2D sdf files
            args = f'-isdf {os.path.abspath(fnSmall)} -h -opdbqt -O {os.path.abspath(fnOut)} '
            if self.preserveCharges.get() == 0:
              args += '--partialcharge gasteiger '
            runOpenBabel(protocol=self, args=args, cwd=oDir)
          else:
            if ext == '.sdf':
              auxPDB = os.path.abspath(self._getTmpPath(os.path.basename(fnOut).replace('.pdbqt', '.pdb')))
              args = f'-isdf {os.path.abspath(fnSmall)} -h -opdb -O {auxPDB} '
              runOpenBabel(protocol=self, args=args, cwd=oDir)
              fnSmall = auxPDB

            # Neccessary to have a local copy of ligandFile from mgltools 1.5.7
            molLink = os.path.join(oDir, os.path.basename(fnSmall))
            createLink(fnSmall, molLink)
            args = f'-l {fnSmall} -o {fnOut} '
            ProtChemADTPrepare.callPrepare(self, "prepare_ligand4", args, outDir=oDir)
            os.remove(molLink)
        except:
          failedMols.append(fnRoot)

      return failedMols

    def performConfGeneration(self, molFns, it):
      failedMols, inDir = [], self.getPreparedDirPath(it)
      for file in molFns:
        if file.endswith('.pdbqt'):
          fnRoot = getBaseName(file)
          if self.method_conf.get() == 0:  # Genetic algorithm
            args = f" {os.path.abspath(file)} --conformer --nconf {self.number_conf.get() - 1} --score rmsd " \
                   f"--writeconformers -O {fnRoot}_conformers.pdbqt"
          else:  # confab
            args = f" {os.path.abspath(file)} --confab --original --conf {self.number_conf.get() - 1} " \
                   f"--rcutoff {self.rmsd_cutoff.get()} -O {fnRoot}_conformers.pdbqt"

          try:
            runOpenBabel(protocol=self, args=args, cwd=os.path.abspath(inDir))
          except:
            failedMols.append(fnRoot)

      return failedMols

    def indOutputCreation(self, file, outDir, molName=None, outConfDir=None):
      '''Returns a dict as {molName: [(molFile, confFile), ....]} from the prepared molecules
      '''
      if not os.path.exists(outDir):
        os.mkdir(outDir)
      if outConfDir is None:
        outConfDir = outDir

      inDir = os.path.dirname(file)
      molName = getBaseName(file) if not molName else molName
      outMolDic = {molName: []}
      if self.doConformers.get():
        firstConfFile = self._getTmpPath('{}-{}.pdbqt'.format(molName, 1))
        shutil.copy(file, firstConfFile)
        confFile = os.path.join(inDir, "{}_conformers.pdbqt".format(molName))
        outConfFile = os.path.join(outConfDir, "{}_conformers.pdbqt".format(molName))

        confFile = appendToConformersFile(confFile, firstConfFile, beginning=True, outConfFile=outConfFile)
        molFiles = splitConformerFile(confFile, outDir=outDir)
        for molFile in molFiles:
          outMolDic[molName].append((molFile, confFile))
      else:
        oFile = os.path.join(outDir, os.path.split(file)[-1])
        os.rename(file, oFile)
        outMolDic[molName].append((oFile, None))
      return outMolDic

    def createOutputMols(self, inMols, outMolDic):
      '''Creates the output SetOFSmallMolecules from a copy of the input.
      Original set and mols attributes are kept in the new prepared molecules
      '''
      objId = 1
      outputSmallMolecules = inMols.createCopy(self._getPath(), copyInfo=True)
      for inMol in inMols:
        nMol = inMol.clone()
        molName = nMol.getMolName()
        if molName in outMolDic:
          for molFile, confFile in outMolDic[molName]:
            nMol.setFileName(molFile)
            nMol.setMolClass('AutoDock')
            nMol._ConformersFile = pwobj.String(confFile)
            if confFile:
              nMol.setConfId(molFile.split('-')[-1].split('.')[0])
            nMol.setObjId(objId)
            objId += 1
            outputSmallMolecules.append(nMol)

      outputSmallMolecules.updateMolClass()
      return outputSmallMolecules

    #################### UTILS FUNCTIONS ################

    def getPreparedDirs(self):
      pDir = os.path.dirname(self.getPreparedDirPath(1))
      return [os.path.join(pDir, d) for d in os.listdir(pDir) if 'thread_' in d]

    def getPreparedDirPath(self, it):
      return os.path.abspath(self._getTmpPath(f'thread_{it}'))


#################### VALIDATION FUNCTIONS ################

    def _warnings(self):
      ws = []
      for mol in self.inputSmallMolecules.get():
          if mol.getFileName().endswith('.sdf') and (self.repair.get()==3 or self.repair.get()==1):
              ws.append('The Autodock4 script prepare_ligand4.py cannot handle to add hydrogens to '
                        'molecules coming from 2D sdf files. Do you want to continue performing a '
                        'similar operation with openBabel? (Note that the resulting molecules might be different)')
              break
      return ws

