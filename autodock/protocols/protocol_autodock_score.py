# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

import os, glob, shutil

from pyworkflow.protocol.params import FloatParam, MultiPointerParam
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath, createLink

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import generate_gpf, calculate_centerMass, getBaseFileName, makeSubsets, insistentRun
from pwchem import Plugin as pwchem_plugin
from pwchem.constants import MGL_DIC

from autodock.protocols.protocol_autodock import ProtChemAutodockBase
from autodock import Plugin as autodock_plugin


LGA, GA, LS, SA = 0, 1, 2, 3
searchDic = {LGA: 'Lamarckian Genetic Algorithm', GA: 'Genetic Algorithm', LS: 'Local Search',
             SA: 'Simulated annealing'}

class ProtChemAutodockScore(ProtChemAutodockBase):
  """Perform a docking experiment with autodock. Grid must be generated in this protocol in order to take into
       account ligands atom types. See the help at
       http://autodock.scripps.edu/faqs-help/manual/autodock-4-2-user-guide/AutoDock4.2_UserGuide.pdf"""
  _label = 'Autodock scoring'
  _program = ""


  def _defineParams(self, form):
    form.addSection(label='Input')
    group = form.addGroup('Input')
    group.addParam('inputSmallMolecules', MultiPointerParam, pointerClass="SetOfSmallMolecules",
                   label='Input small molecules: ',
                   help="Input small molecules to be scored with AutoDock4")
    group.addParam('radius', FloatParam, label='Grid radius for whole protein: ', allowsNull=False,
                   help='Radius of the Autodock grid for the whole protein')
    group.addParam('spacing', FloatParam, label='Spacing of the grids: ', default=0.5, allowsNull=False,
                   help='Spacing of the generated Autodock grids')

    form.addParallelSection(threads=4, mpi=1)

  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
    self.ligandFileNames = []

    convSteps, scoreSteps = [], []
    subsets = makeSubsets(self.getInputMols(), self.numberOfThreads.get(), cloneItem=True)

    for mols in subsets:
      cId = self._insertFunctionStep('convertStep', mols, prerequisites=[])
      convSteps.append(cId)

    gridId = self._insertFunctionStep('generateGridsStep', prerequisites=convSteps)

    for mols in subsets:
        dockId = self._insertFunctionStep('scoreStep', mols, prerequisites=[gridId])
        scoreSteps.append(dockId)

    self._insertFunctionStep('createOutputStep', prerequisites=scoreSteps)

  def getInputMols(self):
      mols = []
      for inputPointer in self.inputSmallMolecules:
          for mol in inputPointer.get():
              mols.append(mol.clone())
      return mols

  def convertStep(self, mols):
    for mol in mols:
        fnSmall, smallDir = self.convertLigand2PDBQT(mol, self._getExtraPath('conformers'), pose=True)
        self.ligandFileNames.append(fnSmall)

    receptorFile = self.getOriginalReceptorFile()
    if receptorFile.endswith('.pdb'):
        self.convertReceptor2PDBQT(receptorFile)
        shutil.copy(receptorFile, self.getReceptorPDB())
    elif receptorFile.endswith('.pdbqt'):
        self.convertReceptor2PDB(receptorFile)
        shutil.copy(receptorFile, self.getReceptorPDBQT())

  def generateGridsStep(self):
    fnReceptor = self.getReceptorPDBQT()
    outDir = self.getScoringDir()

    radius = self.radius.get()
    # Use the original pdb for mass center
    pdbFile = self.getOriginalReceptorFile()
    if not os.path.splitext(pdbFile)[1] == '.pdb':
      pdbFile = self.convertReceptor2PDB(pdbFile)
    structure, x_center, y_center, z_center = calculate_centerMass(pdbFile)

    makePath(outDir)

    npts = (radius * 2) / self.spacing.get()
    gpf_file = generate_gpf(fnReceptor, spacing=self.spacing.get(),
                            xc=x_center, yc=y_center, zc=z_center,
                            npts=npts, outDir=outDir, ligandFns=self.ligandFileNames)

    args = "-p {} -l {}.glg".format(gpf_file, self.getReceptorName())
    insistentRun(self, autodock_plugin.getPackagePath(package='AUTODOCK', path="autogrid4"), args, cwd=outDir)

  def scoreStep(self, mols):
    flexReceptorFn, receptorFn = None, self.getReceptorPDBQT()
    outDir = self.getScoringDir()

    for mol in mols:
        molFn = self.getMolLigandName(mol)
        dpfFile = self.writeDPF(outDir, molFn, receptorFn, flexReceptorFn)

        fnDLG = dpfFile.replace('.dpf', '.dlg')
        args = "-p %s -l %s" % (dpfFile, fnDLG)
        self.runJob(autodock_plugin.getPackagePath(package='AUTODOCK', path="autodock4"), args, cwd=outDir)


  def createOutputStep(self):
    outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())

    scoresDic = {}
    for dlgFile in self.getDockedLigandsFiles():
        molName = getBaseFileName(dlgFile)
        scoresDic[molName] = self.parseDockedMolsDLG(dlgFile)


    for smallMol in self.getInputMols():
      molName = smallMol.getUniqueName()
      if molName in scoresDic:
        molDic = scoresDic[molName]

        newSmallMol = SmallMolecule()
        newSmallMol.copy(smallMol, copyId=False)
        newSmallMol._ADEnergy = pwobj.Float(molDic['energy'])
        if 'ki' in molDic:
          newSmallMol._ligandEfficiency = pwobj.String(molDic['ki'])
        else:
          newSmallMol._ligandEfficiency = pwobj.String(None)

        outputSet.append(newSmallMol)

    # todo: manage several receptor conformations when flexible docking
    outputSet.proteinFile.set(self.getOriginalReceptorFile())
    outputSet.setDocked(True)
    self._defineOutputs(outputSmallMolecules=outputSet)
    self._defineSourceRelation(self.inputSmallMolecules, outputSet)

  ########################### Utils functions ############################

  def reorderIds(self, inSet):
    '''Return the set with the reordered ids and a mapper dictionary {newId: oldId}'''
    idsDic = {}
    for i, item in enumerate(inSet):
      idsDic[i + 1] = item.getObjId()
      item.setObjId(i + 1)
    return inSet, idsDic

  def getScoringDir(self):
      return self._getExtraPath('scoring')

  def getDockedLigandsFiles(self):
    outDir = self.getScoringDir()
    return list(glob.glob(os.path.join(outDir, '*.dlg')))

  def getOriginalReceptorFile(self):
      return self.inputSmallMolecules[0].get().getProteinFile()

  def getMolLigandName(self, mol):
    molName = mol.getUniqueName()
    for ligName in self.ligandFileNames:
      if molName + '.pdbqt' in ligName:
        return ligName

  def writeDPF(self, outDir, molFn, receptorFn, flexFn=None):
      molName, _ = os.path.splitext(os.path.basename(molFn))
      createLink(molFn, os.path.join(outDir, molName + '.pdbqt'))
      makePath(outDir)

      fnDPF = os.path.abspath(os.path.join(outDir, molName + ".dpf"))
      baseDPF = os.path.abspath(os.path.join(outDir, molName + "_base.dpf"))
      # General parameters
      args = " -l %s -r %s -o %s" % (molFn, receptorFn, baseDPF)
      if flexFn:
          args += ' -x ' + flexFn

      self.runJob(pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh'),
                  autodock_plugin.getADTPath('Utilities24/prepare_dpf42.py') + args, cwd=outDir)


      myDPFstr, cont = '', True
      with open(baseDPF) as f:
        for line in f:
          if not cont:
              break
          elif line.startswith('move'):
              cont = False
              line = line.replace('#', ' #')

          myDPFstr += line

      myDPFstr += '\n\n# EVALUATE LIGAND ENERGY in RECEPTOR SECTION\n' \
                  'epdb\n'

      with open(fnDPF, 'w') as f:
          f.write(myDPFstr)
      return fnDPF

  def parseDockedMolsDLG(self, fnDlg):
    molDic = {}
    with open(fnDlg) as fRes:
      for line in fRes:
        if line.startswith('epdb: USER    Estimated Free Energy'):
          molDic['energy'] = line.split('=')[1].split('kcal/mol')[0]
        elif line.startswith('epdb: USER    Estimated Inhibition'):
          molDic['ki'] = ' '.join(line.split()[7:9])
    return molDic



