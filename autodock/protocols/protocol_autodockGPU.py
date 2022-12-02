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

from pyworkflow.protocol.params import IntParam, FloatParam, BooleanParam, \
  LEVEL_ADVANCED, USE_GPU, GPU_LIST, StringParam, EnumParam
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import runOpenBabel, generate_gpf, calculate_centerMass, getBaseFileName
from pwchem import Plugin as pwchem_plugin
from pwchem.constants import MGL_DIC

from autodock import Plugin as autodock_plugin
from autodock.protocols.protocol_autodock import ProtChemAutodockBase

SW, SD, FIRE, AD, ADAM = 0, 1, 2, 3, 4
searchDic = {SW: 'Solis-Wets', SD: 'Steepest-Descent', FIRE: 'FIRE',
             AD: 'ADADELTA', ADAM: 'ADAM'}
searchKeys = {SW: 'sw', SD: 'sd', FIRE: 'fire', AD: 'ad', ADAM: 'adam'}


class ProtChemAutodockGPU(ProtChemAutodockBase):
  """Perform a docking experiment with AutoDock-GPU https://github.com/ccsb-scripps/AutoDock-GPU"""
  _label = 'AutoDock-GPU docking'
  _program = ""

  def _defineParams(self, form):
    form.addHidden(USE_GPU, BooleanParam, default=True,
                   label="Use GPU for execution: ",
                   help="This protocol has both CPU and GPU implementation.\
                                         Select the one you want to use.")

    form.addHidden(GPU_LIST, StringParam, default='0', label="Choose GPU IDs",
                   help="Add a list of GPU devices that can be used")

    super()._defineParams(form)

    form.addSection(label="Search")
    group = form.addGroup('Heuristics')
    group.addParam('heuristics', BooleanParam, label='Use heuristics to guess parameters: ', default=True,
                   help='Use heuristics implemented in the AutoDock-GPU software to guess the optimal number of '
                        'evaluations and the local search method')
    group.addParam('heurmax', IntParam, label='Limit number of evaluations (heurMax): ',
                   default=12000000, condition='heuristics', expertLevel=LEVEL_ADVANCED,
                   help='Asymptotic heuristics # evals limit (smooth limit)')
    group.addParam('autostop', BooleanParam, label='Perform autoStop: ', default=True,
                   help='Automatic stopping criterion based on convergence')
    line = group.addLine('AutoStop parameters: ', expertLevel=LEVEL_ADVANCED, condition='autostop',
                         help='AutoStop testing frequency (in # of generations) / '
                              'energy standard deviation tolerance (kcal/mol)')
    line.addParam('asfreq', IntParam, label='Frequency: ', default=5)
    line.addParam('stopstd', FloatParam, label='Energy tolerance: ', default=0.15)

    group = form.addGroup('Global search', condition='not heuristics')
    group.addParam('gaPop', IntParam, label='Population size', default=150,
                   help='This is the number of individuals in the population. Each individual is a coupling of a '
                        'genotype and its associated phenotype')
    group.addParam('gaNumEvals', IntParam, label='Number of evaluations: ', default=2500000,
                   help='Set the maximum number of energy evaluations performed during each GA, LGA, or LS run')
    group.addParam('gaNumGens', IntParam, label='Number of generations: ', default=42000,
                   help='This is the maximum number of generations simulated during each GA or LGA run')

    line = group.addLine('Genetic algorithm (LGA) rates: ', expertLevel=LEVEL_ADVANCED,
                         help='Probability that a particular gene is mutated / pairs in the population will exchange '
                              'genetic material / selection rates')
    line.addParam('mrat', FloatParam, label='Mutation: ', default=0.02)
    line.addParam('crat', FloatParam, label='Crossover: ', default=0.8)
    line.addParam('trat', FloatParam, label='Selection: ', default=0.6)

    line = group.addLine('Maximum LGA deltas: ', expertLevel=LEVEL_ADVANCED,
                         help='Maximum LGA movement (�) / angle (�) delta')
    line.addParam('dmov', FloatParam, label='Movement: ', default=2.0)
    line.addParam('dang', FloatParam, label='Angle: ', default=90)

    group = form.addGroup("Local search", condition='not heuristics')
    group.addParam('lsFreq', FloatParam, label='Local search frequency: ', default=1.00,
                   help='This is the probability of any particular phenotype being subjected to local search')
    group.addParam('lsType', EnumParam, label='Local search method: ', choices=list(searchDic.values()), default=3,
                   help='Method for local search to use.')
    group.addParam('lsMaxIts', IntParam, label='Number of iterations: ', default=300,
                   help='This is the maximum number of iterations that the local search procedure applies to the '
                        'phenotype of any given individual, per generation')

    line = group.addLine('Solis-Wets variance: ', expertLevel=LEVEL_ADVANCED, condition='lsType=={}'.format(SW),
                         help='Maximum Solis-Wets movement (�) / angle (�) variance for making changes to genes '
                              '(i.e. translations, orientation and torsions)')
    line.addParam('swLbRho', FloatParam, label='Lower bound: ', default=0.01)
    line.addParam('swMovD', FloatParam, label='Movement: ', default=2)
    line.addParam('swAngD', FloatParam, label='Angle: ', default=75)

    group.addParam('swMaxSucc', IntParam, label='Successes/failures in a raw: ', default=4, expertLevel=LEVEL_ADVANCED,
                   condition='lsType=={}'.format(SW), help='Solis-Wets consecutive success/failure limit to adjust rho')

    form.addParallelSection(threads=4, mpi=1)

  def getADGPUArgs(self):
    args = ''
    if self.heuristics:
      args += '-H 1 -E {} '.format(self.heurmax.get())
    else:
      args += '-e {} -l {} '.format(self.gaNumEvals.get(), searchKeys[self.lsType.get()])

      args += '-g {} -p {} '. \
        format(self.gaNumGens.get(), self.lsMaxIts.get(), self.gaPop.get())
      args += '--mrat {} --crat {} --trat {} '. \
        format(self.mrat.get(), self.crat.get(), self.trat.get())

      args += '--dmov {} --dang {} '.format(self.dang.get(), self.dang.get())
      args += '--lsit {} --lsrat {} '.format(self.lsMaxIts.get(), self.lsFreq.get())

      if self.lsType.get() == SW:
        args += '--rholb {} --lsmov {} --lsang {} --cslim {} '. \
          format(self.swLbRho.get(), self.swMovD.get(), self.swAngD.get(), self.swMaxSucc.get())

    if self.autostop:
      args += '-A 1 -a {} --stopstd {} '.format(self.asfreq.get(), self.stopstd.get())
    return args

  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
    cId = self._insertFunctionStep('convertStep', prerequisites=[])

    mols = self.inputSmallMolecules.get()
    gpuList = self.getGPU_Ids()

    dockSteps = []
    if self.fromReceptor.get() == 0:
      gridId = self._insertFunctionStep('generateGridsStep', prerequisites=[cId])
      dockId = self._insertFunctionStep('dockStep', mols, gpuList, prerequisites=[gridId])
      dockSteps.append(dockId)
    else:
      for pocket in self.inputStructROIs.get():
        gridId = self._insertFunctionStep('generateGridsStep', pocket.clone(), prerequisites=[cId])
        dockId = self._insertFunctionStep('dockStep', mols, gpuList, pocket.clone(), prerequisites=[gridId])
        dockSteps.append(dockId)

    self._insertFunctionStep('createOutputStep', prerequisites=dockSteps)

  def convertStep(self):
    self.ligandFileNames = []
    for mol in self.inputSmallMolecules.get():
      fnSmall, smallDir = self.convert2PDBQT(mol.clone(), self._getExtraPath('ligands'))
      self.ligandFileNames.append(fnSmall)

    self.receptorFile = self.getOriginalReceptorFile()
    if not '.pdbqt' in self.receptorFile:
      self.receptorFile = self.convertReceptor2PDBQT(self.receptorFile)

  def generateGridsStep(self, pocket=None):
    fnReceptor = self.receptorFile
    outDir = self.getOutputPocketDir(pocket)
    if self.fromReceptor.get() == 0:
      radius = self.radius.get()
      # Use the original pdb for mass center
      pdbFile = self.getOriginalReceptorFile()
      if not os.path.splitext(pdbFile)[1] == '.pdb':
        pdbFile = self.convertReceptor2PDB(pdbFile)
      structure, x_center, y_center, z_center = calculate_centerMass(pdbFile)
    else:
      radius = (pocket.getDiameter() / 2) * self.pocketRadiusN.get()
      x_center, y_center, z_center = pocket.calculateMassCenter()

    makePath(outDir)

    if self.doFlexRes:
      flexFn, fnReceptor = self.buildFlexReceptor(fnReceptor)

    npts = (radius * 2) / self.spacing.get()
    gpf_file = generate_gpf(fnReceptor, spacing=self.spacing.get(),
                            xc=x_center, yc=y_center, zc=z_center,
                            npts=npts, outDir=outDir, ligandFns=self.ligandFileNames)

    args = "-p {} -l {}.glg".format(gpf_file, self.getReceptorName())
    self.runJob(autodock_plugin.getAutodockPath("autogrid4"), args, cwd=outDir)

  def dockStep(self, mols, gpuIdxs, pocket=None):
    # Prepare grid
    if self.doFlexRes:
      flexReceptorFn, receptorFn = self.getFlexFiles()
    else:
      flexReceptorFn, receptorFn = None, self.receptorFile
    outDir = self.getOutputPocketDir(pocket)

    if self.doFlexRes:
      fldFile = 'receptor_rig.maps.fld'
    else:
      fldFile = '{}.maps.fld'.format(self.getReceptorName())

    batchFile = os.path.abspath(os.path.join(outDir, 'batchFile.txt'))
    with open(batchFile, 'w') as f:
      f.write('{}\n'.format(fldFile))
      for mol in mols:
        molFn = self.getMolLigandFileName(mol)
        molBase = molFn.split('/')[-1]
        molLink = os.path.join(outDir, molBase)
        os.link(molFn, molLink)

        f.write('{}\n{}\n'.format(molBase, getBaseFileName(molBase)))

    args = '-B {} -D {} -n {} --rmstol {} -C 1 --output-cluster-poses auto '. \
      format(batchFile, ','.join(gpuIdxs), self.nRuns.get(), self.rmsTol.get())
    if self.doFlexRes:
      args += '-F {} '.format(flexReceptorFn)
    args += self.getADGPUArgs()
    autodock_plugin.runAutodockGPU(self, args, outDir)

  def createOutputStep(self):
    outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())
    for pocketDir in self.getPocketDirs():
      pocketDic = {}
      gridId = self.getGridId(pocketDir)
      for dlgFile in self.getDockedLigandsFiles(pocketDir):
        molName = getBaseFileName(dlgFile)
        pocketDic[molName] = self.parseDockedMolsDLG(dlgFile)

        for modelId in pocketDic[molName]:
          pdbqtFile = self._getPath('g{}_{}_{}.pdbqt'.format(gridId, molName, modelId))
          with open(pdbqtFile, 'w') as f:
            f.write(pocketDic[molName][modelId]['pdb'])
          pocketDic[molName][modelId]['file'] = pdbqtFile

      for smallMol in self.inputSmallMolecules.get():
        molFile = smallMol.getFileName()
        molName = getBaseFileName(molFile)
        if molName in pocketDic:
          molDic = pocketDic[molName]

          for posId in molDic:
            newSmallMol = SmallMolecule()
            newSmallMol.copy(smallMol, copyId=False)
            newSmallMol._energy = pwobj.Float(molDic[posId]['energy'])
            if 'ki' in molDic[posId]:
              newSmallMol._ligandEfficiency = pwobj.Float(molDic[posId]['ki'])
            else:
              newSmallMol._ligandEfficiency = pwobj.Float(None)
            if os.path.getsize(molDic[posId]['file']) > 0:
              newSmallMol.poseFile.set(molDic[posId]['file'])
              newSmallMol.setPoseId(posId)
              newSmallMol.gridId.set(gridId)
              newSmallMol.setMolClass('Autodock4')
              newSmallMol.setDockId(self.getObjId())

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

  def getReceptorDir(self):
    atomStructFn = self.getOriginalReceptorFile()
    return '/'.join(atomStructFn.split('/')[:-1])

  def getGPU_Ids(self):
    gpus = []
    for gp in getattr(self, GPU_LIST).get().split(','):
      gpus.append(str(int(gp) + 1))
    return gpus

  def getLigandsFileNames(self):
    ligFns = []
    for mol in self.inputSmallMolecules.get():
      ligFns.append(mol.getFileName())
    return ligFns

  def getMolLigandFileName(self, mol):
    molName = mol.getUniqueName()
    for ligFileName in self.ligandFileNames:
      if molName + '.pdbqt' in ligFileName:
        return ligFileName

  def parseDockedMolsDLG(self, fnDlg):
    molDic = {}
    i = 1
    with open(fnDlg) as fRes:
      for line in fRes:
        if line.startswith('DOCKED: MODEL'):
          posId = line.split()[-1]
          molDic[posId] = {'pdb': ''}
        elif line.startswith('DOCKED: USER    Estimated Free Energy'):
          molDic[posId]['energy'] = line.split('=')[1].split('kcal/mol')[0]
        elif line.startswith('DOCKED: USER    Estimated Inhibition'):
          molDic[posId]['ki'] = line.split()[7]
        elif line.startswith('DOCKED: REMARK') or line.startswith('TER'):
          molDic[posId]['pdb'] += line[8:]
        elif line.startswith('DOCKED: ATOM'):
          molDic[posId]['pdb'] += line[8:]
          i += 1
    return molDic

  def _validate(self):
    vals = []
    return vals

  def commentFirstLine(self, fn):
    with open(fn) as f:
      all = '# ' + f.read()
    with open(fn, 'w') as f:
      f.write(all)
    return fn

  def getDockedLigandsFiles(self, outDir):
    return list(glob.glob(os.path.join(outDir, '*.dlg')))
