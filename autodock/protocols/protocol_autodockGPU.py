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
from pwchem.utils import getBaseName, performBatchThreading, replaceInFiles, makeSubsets

from autodock import Plugin as autodockPlugin
from autodock.protocols.protocol_autodock import ProtChemAutodockBase
from autodock.objects import RingtailDatabase

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
    form.addParam('ringtailOutput', BooleanParam, label='Create ringtail output: ', default=False,
                  help='Create a ringtail database as output of the docking execution')
    form.addParam('remTmp', BooleanParam, label='Remove intermediate files: ', default=True,
                  expertLevel=LEVEL_ADVANCED, condition='not ringtailOutput',
                  help='Whether to remove the intermediate files generate by AutoDock to reduce the memory usage')

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

  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
      inMols = self.inputSmallMolecules.get()
      nt = self.numberOfThreads.get()
      gpuList = self.getGPU_Ids()
      subsets = makeSubsets(inMols, nt - 1, cloneItem=True)

      cRStep = self._insertFunctionStep(self.convertReceptorStep, prerequisites=[], needsGPU=False)

      cSteps = []
      for it, molSet in enumerate(subsets):
        cSteps.append(self._insertFunctionStep(self.convertLigandsStep, molSet, it, prerequisites=[], needsGPU=False))

      dockSteps = []
      gridReqs = [cRStep] + cSteps
      if self.fromReceptor.get() == 0:
        gridId = self._insertFunctionStep(self.generateGridsStep, prerequisites=gridReqs, needsGPU=False)
        for it, _ in enumerate(subsets):
          dockId = self._insertFunctionStep(self.dockStep, it, gpuList, prerequisites=[gridId])
          dockSteps.append(dockId)
      else:
        for pocket in self.inputStructROIs.get():
          gridId = self._insertFunctionStep(self.generateGridsStep, pocket.clone(), prerequisites=gridReqs, needsGPU=False)
          for it, _ in enumerate(subsets):
            dockId = self._insertFunctionStep(self.dockStep, it, gpuList, pocket.clone(), prerequisites=[gridId])
            dockSteps.append(dockId)

      self._insertFunctionStep(self.createOutputStep, prerequisites=dockSteps, needsGPU=False)

  def dockStep(self, it, gpuIdxs, pocket=None):
      molFns = self.getConvertedLigandsFiles(it)
      flexReceptorFn = self.getFlexFiles()[0] if self.doFlexRes else None
      outDir = self.getOutputPocketDir(pocket)

      fldFile = f'{self.getReceptorName()}.maps.fld'
      self.fixFldFile(os.path.join(outDir, fldFile))

      batchFile = self.writeBatchFile(fldFile, molFns, outDir, it)
      args = f"-B {batchFile} -D {','.join(gpuIdxs)} -n {self.nRuns.get()} --rmstol {self.rmsTol.get()} -C 1 " \
             f"--output-cluster-poses auto "
      if self.doFlexRes:
        args += f'-F {flexReceptorFn} '
      args += self.getADGPUArgs()
      autodockPlugin.runAutodockGPU(self, args, outDir)

  def getBatchFile(self, outDir, it):
    return os.path.abspath(os.path.join(outDir, f'batchFile_{it}.txt'))

  def createOutputStep(self):
      nt = self.numberOfThreads.get()
      recFile = self.getReceptorPDBQT()
      if self.ringtailOutput.get():

        self.fixDLGReceptor()
        outDir = os.path.abspath(self._getExtraPath())
        args = f'write --file_path {outDir} --recursive -o ringtail.db -mpr {nt} --overwrite -sr -rf {recFile}'
        autodockPlugin.runRingtail(self, args, cwd=self._getPath())

        outputDB = RingtailDatabase(filename=self._getPath('ringtail.db'))
        outputDB.setReceptorFile(recFile)
        outputDB.createSumFile(self.getSumPath())
        outputDB.performBaseFilter()
        self._defineOutputs(outputRingtail=outputDB)
      else:
        outDir = self._getPath('outputLigands')
        makePath(outDir)

        inputMols = self.inputSmallMolecules.get()
        outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())
        for pocketDir in self.getPocketDirs():
          dlgFiles = self.getDockedLigandsFiles(pocketDir)
          gridId = self.getGridId(pocketDir)
          pocketDics = performBatchThreading(self.performOutputParsing, dlgFiles, nt, cloneItem=False,
                                             gridId=gridId, outDir=outDir)
          pocketDic = {k: v for pDic in pocketDics for (k, v) in pDic.items()}
  
          outputMols = performBatchThreading(self.performOutputCreation, inputMols, nt,
                                             gridId=gridId, pocketDic=pocketDic, recFile=recFile)

          print(f'------output Mols: {outputMols}')

          if self.remTmp.get():
            self.removeTmpFiles(pocketDir)

          for smallMol in outputMols:
            smallMol.guessMolName()
            outputSet.append(smallMol)
  
        outputSet.setProteinFile(recFile)
        outputSet.setDocked(True)
        outputSet.saveGroupIndexes()
        self._defineOutputs(outputSmallMolecules=outputSet)
        self._defineSourceRelation(self.inputSmallMolecules, outputSet)

      self.cleanTmpFiles()


  def performOutputParsing(self, dlgFiles, molLists, it, gridId, outDir):
    pocketDic = {}
    for dlgFile in dlgFiles:
      molName = getBaseName(dlgFile)
      pocketDic[molName] = self.parseDockedMolsDLG(dlgFile)

      for modelId in pocketDic[molName]:
        pdbqtFile = os.path.join(outDir, 'g{}_{}_{}.pdbqt'.format(gridId, molName, modelId))
        with open(pdbqtFile, 'w') as f:
          f.write(pocketDic[molName][modelId]['pdb'])
        pocketDic[molName][modelId]['file'] = pdbqtFile

    molLists[it] = [pocketDic]

  def performOutputCreation(self, mols, molLists, it, pocketDic, gridId, recFile):
    outMols = []
    for smallMol in mols:
      molFile = smallMol.getFileName()
      molName = getBaseName(molFile)
      print(f'----mol name:{molName}')
      print(f'----pocket dic: {pocketDic}')
      if molName in pocketDic:
        molDic = pocketDic[molName]

        for posId in molDic:
          newSmallMol = SmallMolecule()
          newSmallMol.copy(smallMol, copyId=False)
          newSmallMol._energy = pwobj.Float(molDic[posId]['energy'])
          newSmallMol._ligandEfficiency = pwobj.Float(molDic[posId]['ligEfficiency'])
          ki = molDic[posId]['ki'] if 'ki' in molDic[posId] else None
          newSmallMol._ki = pwobj.String(ki)

          poseFile = molDic[posId]['file']
          if os.path.getsize(poseFile) > 0:
            if self.doFlexRes:
              poseFile, curRecFile = self.makeFlexPoseFiles(poseFile, recFile)
              newSmallMol.setProteinFile(os.path.relpath(curRecFile))

            newSmallMol.poseFile.set(os.path.relpath(poseFile))
            newSmallMol.setPoseId(posId)
            newSmallMol.gridId.set(gridId)
            newSmallMol.setMolClass('Autodock4')
            newSmallMol.setDockId(self.getObjId())

            outMols.append(newSmallMol)
    molLists[it] = outMols

  ########################### Utils functions ############################

  def fixFldFile(self, fldFile):
    s = ''
    with open(fldFile) as f:
      for line in f:
        if line.startswith('#MACROMOLECULE'):
          sline = line.split()
          line = '{} ../{}\n'.format(sline[0], getBaseName(sline[1]).strip()+'.pdbqt')
        s += line
    with open(fldFile, 'w') as f:
      f.write(s)

  def writeBatchFile(self, fldFile, molFns, outDir, it):
    batchFile = self.getBatchFile(outDir, it)
    with open(batchFile, 'w') as f:
      f.write(f'{fldFile}\n')
      for molFn in molFns:
        molBase = molFn.split('/')[-1]
        molLink = os.path.join(outDir, molBase)
        if not os.path.exists(molLink):
          os.link(molFn, molLink)

        f.write(f'{molBase}\n{getBaseName(molBase)}\n')
    return batchFile

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

  def getGPU_Ids(self):
    gpus = []
    for gp in getattr(self, GPU_LIST).get().split(','):
      gpus.append(str(int(gp) + 1))
    return gpus

  def commentFirstLine(self, fn):
    with open(fn) as f:
      allStr = '# ' + f.read()
    with open(fn, 'w') as f:
      f.write(allStr)
    return fn

  def getDockedLigandsFiles(self, outDir):
    return list(glob.glob(os.path.join(outDir, '*.dlg')))
  
  def getSumPath(self):
    return os.path.abspath(self._getExtraPath('ringSum.txt'))

  def fixDLGReceptor(self):
    recName = self.getReceptorName()
    for pocketDir in self.getPocketDirs():
      replaceInFiles(os.path.abspath(pocketDir), f'..\/{recName}', recName, file_extension='.dlg')

  def removeTmpFiles(self, pDir):
    for file in os.listdir(pDir):
      if file.split('.')[-1] in ['dlg', 'xml', 'pdbqt', 'map']:
        os.remove(os.path.join(pDir, file))

  def _summary(self):
    s = []
    if os.path.exists(self.getSumPath()):
      with open(self.getSumPath()) as f:
        s.append(f.read())
    return s
