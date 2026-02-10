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

# todo: try whether it works on big GPU memories (on laptop it crashed unknown error)

import os, glob

from pyworkflow.protocol.params import IntParam, FloatParam, BooleanParam, \
  LEVEL_ADVANCED, USE_GPU, GPU_LIST, StringParam, EnumParam
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import getBaseName, performBatchThreading, replaceInFiles, makeSubsets, \
  calculate_centerMass, calculateCoordLimits

from autodock import Plugin as autodockPlugin
from autodock.protocols.protocol_autodock import ProtChemAutodockBase
from autodock.utils import parseDockedPDBQT

MGL = 'MGLTools'

ADT, QUICK2, QUICKW = 0, 1, 2
DOCKSOFT = ['AutoDock-Vina', 'QuickVina2', 'QuickVina-W']


class ProtChemVinaGPU(ProtChemAutodockBase):
  """Perform a docking experiment with Vina-GPU https://github.com/DeltaGroupNJUPT/Vina-GPU-2.1"""
  _label = 'Vina-GPU docking'
  _program = ""

  def _defineParams(self, form):
    form.addSection(label='Input')
    form.addHidden(USE_GPU, BooleanParam, default=True,
                   label="Use GPU for execution: ",
                   help="This protocol has both CPU and GPU implementation.\
                                         Select the one you want to use.")

    form.addHidden(GPU_LIST, StringParam, default='0', label="Choose GPU IDs",
                   help="Add a list of GPU devices that can be used")

    super()._defineInput(form)
    flexGroup = self._defineFlexParams(form)
    form.addParam('remTmp', BooleanParam, label='Remove intermediate files: ', default=True,
                  expertLevel=LEVEL_ADVANCED,
                  help='Whether to remove the intermediate files generate by AutoDock to reduce the memory usage')

    dGroup = form.addGroup('Docking parameters')
    dGroup.addParam('dockSoft', EnumParam, label='Docking software: ', default=ADT, choices=DOCKSOFT,
                    help='Vina-GPU software to perform the docking')
    dGroup.addParam('nRuns', IntParam, label='Number of docking runs: ', default=10,
                    help='Number of independent runs using the selected strategy. \n'
                         'Different docking positions will be found for each of them.')
    dGroup.addParam('searchDepth', IntParam, label='Search depth: ', default=-1, expertLevel=LEVEL_ADVANCED,
                    help='The number of search depth in monte carlo. If < 0, heuristically determined')
    dGroup.addParam('energyRange', IntParam, label='Energy range: ', default=3, expertLevel=LEVEL_ADVANCED,
                    help='Maximum energy difference between the best binding mode and the worst one '
                         'displayed (kcal/mol)')
    dGroup.addParam('lbfgs', BooleanParam, label='Use RILC-BFGS: ', default=True, expertLevel=LEVEL_ADVANCED,
                    help='Use RILC-BFGS optimization algorithm for docking')

    dGroup.addParam('seed', IntParam, label='Random seed: ', default=44, expertLevel=LEVEL_ADVANCED,
                    help='Random seed for the initial random generation')

    form.addParallelSection(threads=4, mpi=1)

  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
      inMols = self.inputSmallMolecules.get()
      nt = self.numberOfThreads.get()
      subsets = makeSubsets(inMols, max(nt - 1, 1), cloneItem=True)

      cRStep = self._insertFunctionStep(self.convertReceptorStep, prerequisites=[], needsGPU=False)

      cSteps = []

      for it, molSet in enumerate(subsets):
        cSteps.append(self._insertFunctionStep(self.convertLigandsStep, molSet, it,
                                               self.getLigConvertedDir(), prerequisites=[], needsGPU=False))

      dockSteps = []
      gridReqs = [cRStep] + cSteps
      if self.fromReceptor.get() == 0:
          dockId = self._insertFunctionStep(self.dockStep, prerequisites=gridReqs)
          dockSteps.append(dockId)
      else:
        for pocket in self.inputStructROIs.get():
            dockId = self._insertFunctionStep(self.dockStep, pocket.clone(), prerequisites=gridReqs)
            dockSteps.append(dockId)

      self._insertFunctionStep(self.createOutputStep, prerequisites=dockSteps, needsGPU=False)

  def dockStep(self, pocket=None):
      confFile, valid = self.writeConfigFile(pocket)
      if valid:
          program = self.getEnumText('dockSoft')
          args = f'--config {confFile}'

          autodockPlugin.runVinaGPU(self, program, args)

  def createOutputStep(self):
    recFile = self.getReceptorPDBQT()
    outDir = self._getPath('outputLigands')
    makePath(outDir)
    outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())

    for pocketDir in self.getPocketDirs():
      pocketDic = {}
      gridId = self.getGridId(pocketDir)
      dockFiles = self.getDockedLigandsFiles(pocketDir)
      for dockFile in dockFiles:
        molName = getBaseName(dockFile)
        pocketDic[molName] = parseDockedPDBQT(dockFile)

      for smallMol in self.inputSmallMolecules.get():
        molName = smallMol.getUniqueName(conf=True)
        if molName in pocketDic:
          molDic = pocketDic[molName]

          for posId in molDic:
              newSmallMol = SmallMolecule()
              newSmallMol.copy(smallMol, copyId=False)
              newSmallMol._energy = pwobj.Float(molDic[posId]['energy'])

              poseFile = molDic[posId]['file']
              if os.path.getsize(poseFile) > 0:

                filename = f'g{gridId}_{os.path.split(poseFile)[-1]}'
                newPoseFile = os.path.join(outDir, filename)
                os.rename(poseFile, newPoseFile)

                newSmallMol.poseFile.set(newPoseFile)
                newSmallMol.setPoseId(posId)
                newSmallMol.gridId.set(gridId)
                newSmallMol.setMolClass('VinaGPU')
                newSmallMol.setDockId(self.getObjId())

                outputSet.append(newSmallMol)

        else:
          print(f'Molecule {molName} was not found in the docking results')

      if self.remTmp.get():
        self.removeTmpFiles(pocketDir)

    outputSet.updateMolClass()
    outputSet.setProteinFile(recFile)
    outputSet.setDocked(True)
    self._defineOutputs(outputSmallMolecules=outputSet)
    self._defineSourceRelation(self.inputSmallMolecules, outputSet)

    self.cleanTmpFiles()

  ########################### Utils functions ############################

  def getLigConvertedDir(self):
      return os.path.abspath(self._getTmpPath('convertedFiles'))

  def getConfigFile(self, pocket):
    it = 1 if pocket is None else pocket.getObjId()
    return os.path.abspath(self._getExtraPath(f'configFile_{it}.txt'))

  def getPocketArgs(self, pocket):
    if self.fromReceptor.get() == 0:
      pdbFile = self.getReceptorPDB()
      minMaxCoords = calculateCoordLimits(pdbFile)
      _, xCenter, yCenter, zCenter = calculate_centerMass(pdbFile)
    else:
      minMaxCoords = pocket.getLimits()
      xCenter, yCenter, zCenter = pocket.calculateMassCenter()
    diams = [(minMax[1] - minMax[0]) * self.pocketRadiusN.get() for minMax in minMaxCoords]
    return (xCenter, yCenter, zCenter), diams

  def writeConfigFile(self, pocket):
    fnReceptor = self.getReceptorPDBQT()
    flexFn = None
    if self.doFlexRes:
      flexFn, fnReceptor = self.buildFlexReceptor(fnReceptor)

    ligDir = self.getLigConvertedDir()
    outDir = self.getOutputPocketDir(pocket)

    center, diams = self.getPocketArgs(pocket)
    valid = self.checkBoxSize(diams)

    confFile = self.getConfigFile(pocket)
    with open(confFile, 'w') as f:
        f.write(f'receptor = {fnReceptor}\n')
        if flexFn is not None:
          f.write(f'flex = {flexFn}\n')
        f.write(f'ligand_directory = {ligDir}\n')
        f.write(f'output_directory = {outDir}\n')

        f.write(f'center_x = {center[0]}\n')
        f.write(f'center_y = {center[1]}\n')
        f.write(f'center_z = {center[2]}\n')
        f.write(f'size_x = {diams[0]}\n')
        f.write(f'size_y = {diams[1]}\n')
        f.write(f'size_z = {diams[2]}\n')

        f.write(f'num_modes = {self.nRuns.get()}\n')
        f.write(f'energy_range = {self.energyRange.get()}\n')

        f.write(f'rilc_bfgs = {1 if self.lbfgs.get() else 0}\n')
        f.write(f'seed = {self.seed.get()}\n')
        f.write('thread = 8000\n')

        if self.searchDepth.get() > 1:
          f.write(f'search_depth = {self.searchDepth.get()}\n')

        # if self.doFlexRes:
        #   args += f'-F {flexReceptorFn} '

    return confFile, valid

  def getVinaGPUArgs(self):
    args = f'--num_modes {self.nRuns.get()} --energy_range {self.energyRange.get()} '
    lbfgs = 1 if self.lbfgs.get() else 0
    args += f'--search_depth {self.searchDepth.get()} --rilc_bfgs {lbfgs} '

    args += f'--seed {self.seed.get()} '
    return args

  def getDockedLigandsFiles(self, outDir):
    oriFiles = list(glob.glob(os.path.join(outDir, '*.pdbqt')))
    newFiles = []
    for file in oriFiles:
      if '_out' in file:
        nFile = file.replace('_out.pdbqt', '.pdbqt')
        os.rename(file, nFile)
      else:
        nFile = file
      newFiles.append(nFile)
    return newFiles
  
  def getSumPath(self):
    return os.path.abspath(self._getExtraPath('ringSum.txt'))

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

  def checkBoxSize(self, diams):
      check = True
      vol = diams[0] * diams[1] * diams[2]
      if vol > 27000:
        check = False
      return check

  def _warnings(self):
      warns = []
      if self.fromReceptor.get() == 0:
          warns.append('The size of the box for the full receptor must be below 27000 A³ for VinaGPU to manage.')
      else:
          for pocket in self.inputStructROIs.get():
              center, diams = self.getPocketArgs(pocket)
              if not self.checkBoxSize(diams):
                  warns.append(f'The size of the pocket {pocket.getObjId()} is too big for VinaGPU to manage.')

      return warns

