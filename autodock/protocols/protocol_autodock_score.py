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

import os, glob

from pyworkflow.protocol.params import FloatParam, MultiPointerParam
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath, createLink

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import runOpenBabel, generate_gpf, calculate_centerMass, getBaseFileName, performBatchThreading
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
    self.receptorFile = self.getOriginalReceptorFile()

    convSteps, scoreSteps = [], []

    mols = self.getInputMols()
    for mol in mols:
        cId = self._insertFunctionStep('convertStep', mol.clone(), prerequisites=[])
        convSteps.append(cId)

    gridId = self._insertFunctionStep('generateGridsStep', prerequisites=convSteps)

    for mol in mols:
        dockId = self._insertFunctionStep('scoreStep', mol.clone(), prerequisites=[gridId])
        scoreSteps.append(dockId)

    self._insertFunctionStep('createOutputStep', mols, prerequisites=scoreSteps)

  def getInputMols(self):
      mols = []
      for inputPointer in self.inputSmallMolecules:
          for mol in inputPointer.get():
              mols.append(mol.clone())
      return mols

  def convertStep(self, mol):
    fnSmall, smallDir = self.convert2PDBQT(mol, self._getExtraPath('conformers'), pose=True)
    self.ligandFileNames.append(fnSmall)

    if not '.pdbqt' in self.receptorFile:
      self.receptorFile = self.convertReceptor2PDBQT(self.receptorFile)

  def generateGridsStep(self):
    fnReceptor = self.receptorFile
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
    self.runJob(autodock_plugin.getAutodockPath("autogrid4"), args, cwd=outDir)

  def scoreStep(self, mol, pocket=None):
    # Prepare grid
    flexReceptorFn, receptorFn = None, self.receptorFile
    outDir = self.getScoringDir()

    molFn = self.getMolLigandName(mol)
    dpfFile = self.writeDPF(outDir, molFn, receptorFn, flexReceptorFn)

    fnDLG = dpfFile.replace('.dpf', '.dlg')
    args = "-p %s -l %s" % (dpfFile, fnDLG)
    self.runJob(autodock_plugin.getAutodockPath("autodock4"), args, cwd=outDir)


  def createOutputStep(self, mols):
    outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())

    scoresDic = {}
    for dlgFile in self.getDockedLigandsFiles():
        molName = getBaseFileName(dlgFile)
        scoresDic[molName] = self.parseDockedMolsDLG(dlgFile)


    for smallMol in mols:
      molName = smallMol.getUniqueName()
      if molName in scoresDic:
        molDic = scoresDic[molName]

        newSmallMol = SmallMolecule()
        newSmallMol.copy(smallMol, copyId=True)
        newSmallMol._energy = pwobj.Float(molDic['energy'])
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

  def convertReceptor2PDB(self, proteinFile):
    inName, inExt = os.path.splitext(os.path.basename(proteinFile))
    oFile = os.path.abspath(os.path.join(self._getTmpPath(inName + '.pdb')))

    args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(proteinFile), oFile)
    runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())

    return oFile

  def convertReceptor2PDBQT(self, proteinFile):
    inName, inExt = os.path.splitext(os.path.basename(proteinFile))
    oFile = os.path.abspath(os.path.join(self._getExtraPath(inName + '.pdbqt')))

    args = ' -v -r %s -o %s' % (proteinFile, oFile)
    self.runJob(pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh'),
                autodock_plugin.getADTPath('Utilities24/prepare_receptor4.py') + args)

    return oFile

  def getScoringDir(self):
      return self._getExtraPath('scoring')

  def getDockedLigandsFiles(self):
    outDir = self.getScoringDir()
    return list(glob.glob(os.path.join(outDir, '*.dlg')))

  def getOriginalReceptorFile(self):
      return self.inputSmallMolecules[0].get().getProteinFile()

  def getReceptorName(self):
      atomStructFn = self.getOriginalReceptorFile()
      return atomStructFn.split('/')[-1].split('.')[0]

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



