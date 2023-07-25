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
from pwchem.utils import *

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
    group.addParam('gaPop', IntParam, label='Population size: ', default=150,
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

      args += '-g {} -i {} -p {} '. \
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

    args += '-s 44 '
    return args

  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
      cId = self._insertFunctionStep('convertStep', prerequisites=[])
      gpuList = self.getGPU_Ids()
      self.receptorName = self.getReceptorName()

      dockSteps = []
      if self.fromReceptor.get() == 0:
        gridId = self._insertFunctionStep('generateGridsStep', None, False, prerequisites=[cId])
        dockId = self._insertFunctionStep('dockStep', gpuList, prerequisites=[gridId])
        dockSteps.append(dockId)
      else:
        for pocket in self.inputStructROIs.get():
          gridId = self._insertFunctionStep('generateGridsStep', pocket.clone(), False, prerequisites=[cId])
          dockId = self._insertFunctionStep('dockStep', gpuList, pocket.clone(), prerequisites=[gridId])
          dockSteps.append(dockId)

      self._insertFunctionStep('createOutputStep', prerequisites=dockSteps)

  def dockStep(self, gpuIdxs, pocket=None):
      molFns = self.getInputPDBQTFiles()
      if self.doFlexRes:
        flexReceptorFn, receptorFn = self.getFlexFiles()
      else:
        flexReceptorFn, receptorFn = None, self.getReceptorPDBQT()
      outDir = self.getOutputPocketDir(pocket)

      fldFile = '{}.maps.fld'.format(self.getReceptorName())
      self.fixFldFile(os.path.join(outDir, fldFile))

      batchFile = os.path.abspath(os.path.join(outDir, 'batchFile.txt'))
      with open(batchFile, 'w') as f:
        f.write('{}\n'.format(fldFile))
        for molFn in molFns:
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
      outDir = self._getPath('outputLigands')
      makePath(outDir)
      outputSet = SetOfSmallMolecules().create(outputPath=outDir)

      for pocketDir in self.getPocketDirs():
        pocketDic = {}
        gridId = self.getGridId(pocketDir)
        for dlgFile in self.getDockedLigandsFiles(pocketDir):
          molName = getBaseFileName(dlgFile)
          pocketDic[molName] = self.parseDockedMolsDLG(dlgFile)

          for modelId in pocketDic[molName]:
            pdbqtFile = os.path.join(outDir, 'g{}_{}_{}.pdbqt'.format(gridId, molName, modelId))
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

  def fixFldFile(self, fldFile):
    s = ''
    with open(fldFile) as f:
      for line in f:
        if line.startswith('#MACROMOLECULE'):
          sline = line.split()
          line = '{} ../{}\n'.format(sline[0], getBaseFileName(sline[1]).strip()+'.pdbqt')
        s += line
    with open(fldFile, 'w') as f:
      f.write(s)


  def getGPU_Ids(self):
    gpus = []
    for gp in getattr(self, GPU_LIST).get().split(','):
      gpus.append(str(int(gp) + 1))
    return gpus

  def parseDockedMolsDLG(self, fnDlg):
    molDic = {}
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
    return molDic

  def commentFirstLine(self, fn):
    with open(fn) as f:
      all = '# ' + f.read()
    with open(fn, 'w') as f:
      f.write(all)
    return fn

  def getDockedLigandsFiles(self, outDir):
    return list(glob.glob(os.path.join(outDir, '*.dlg')))
