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

from pyworkflow.protocol import params
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath

from pwchem import Plugin as pwchem_plugin
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import calculate_centerMass, generate_gpf, insistentRun, getBaseName, \
  makeSubsets, calculateCoordLimits

from autodock import Plugin as autodockPlugin
from autodock.protocols.protocol_autodock import ProtChemAutodockBase
from autodock.constants import VINA_DIC, VINA, SCRUBBER_DIC
from autodock.objects import RingtailDatabase
from autodock.utils import parseDockedPDBQT

meekoScript = 'meeko_preparation.py'
scriptName = 'vina_docking.py'

PDBext, PDBQText = '.pdb', '.pdbqt'

class ProtChemVinaDocking(ProtChemAutodockBase):
    """Dock ligands using Vina"""
    _label = 'Vina docking'
    _program = ""

    _scoreNames = ['Vina', 'AD4']

    def _defineParams(self, form):
        dockGroup = super()._defineParams(form)[1]
        dockGroup.addParam('scoreName', params.EnumParam, label='Score function: ',
                       choices=self._scoreNames, default=0,
                       help='Score function to use for docking')
        dockGroup.addParam('exhaust', params.IntParam, label='Exhaustiveness number: ', default=8,
                       help='Number that controls the amount of runs or trials to use in the vina algorithm.\n'
                            'The higher, the more the search pace is explored, but also it will be slower')
        dockGroup.addParam('maxEvals', params.IntParam, label='Maximum number of evaluation: ',
                       expertLevel=params.LEVEL_ADVANCED, default=0, help='Maximum number of evaluation')

        form.addParam('ringtailOutput', params.BooleanParam, label='Create ringtail output: ', default=False,
                      help='Create a ringtail database as output of the docking execution')

        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
      inMols = self.inputSmallMolecules.get()
      nt = self.numberOfThreads.get()
      subsets = makeSubsets(inMols, nt - 1, cloneItem=True)

      cRStep = self._insertFunctionStep(self.convertReceptorStep, prerequisites=[], needsGPU=False)

      cSteps, dockSteps = [], []
      for it, molSet in enumerate(subsets):
        cSteps.append(self._insertFunctionStep(self.convertLigandsStep, molSet, it, prerequisites=[], needsGPU=False))

        convReqs = [cRStep, cSteps[-1]]
        if self.fromReceptor.get() == 0:
            dockId = self._insertFunctionStep(self.dockStep, it, prerequisites=convReqs, needsGPU=False)
            dockSteps.append(dockId)
        else:
          for pocket in self.inputStructROIs.get():
              dockId = self._insertFunctionStep(self.dockStep, it, pocket.clone(), prerequisites=convReqs, needsGPU=False)
              dockSteps.append(dockId)
      self._insertFunctionStep(self.createOutputStep, prerequisites=dockSteps, needsGPU=False)


    def dockStep(self, it, pocket=None):
      fnReceptor = self.getReceptorPDBQT()
      pdbqtFiles = self.getConvertedLigandsFiles(it)

      outDir = self.getOutputPocketDir(pocket)
      if not os.path.exists(outDir):
          makePath(outDir)

      if self.fromReceptor.get() == 0:
          pdbFile = self.getReceptorPDB()
          minMaxCoords = calculateCoordLimits(pdbFile)
          _, xCenter, yCenter, zCenter = calculate_centerMass(self.getReceptorPDB())
      else:
          minMaxCoords = pocket.getLimits()
          xCenter, yCenter, zCenter = pocket.calculateMassCenter()

      radius = [((minMax[1] - minMax[0]) / 2) * self.pocketRadiusN.get() for minMax in minMaxCoords]
      npts = [(r*2) / self.spacing.get() for r in radius]

      znFFfile = autodockPlugin.getPackagePath(package='VINA', path='AutoDock-Vina/data/AD4Zn.dat') \
        if self.doZnDock.get() else None
      gpfFile = generate_gpf(fnReceptor, spacing=self.spacing.get(), allDefAtomTypes=True,
                              xc=xCenter, yc=yCenter, zc=zCenter,
                              npts=npts, outDir=outDir, ligandFns=pdbqtFiles, znFFfile=znFFfile)

      if self.doFlexRes:
          flexFn, fnReceptor = self.buildFlexReceptor(fnReceptor, cleanZn=self.doZnDock.get())
      else:
          flexFn = None

      nThreads = self.getnThreads()
      scoreFunc = self.getEnumText('scoreName').lower()
      if not self.doZnDock.get() and scoreFunc == 'vina':
        self.performScriptDocking(pdbqtFiles, it, outDir=outDir, fnReceptor=fnReceptor, radius=radius, 
                                  xCenter=xCenter, yCenter=yCenter, zCenter=zCenter,
                                  gpfFile=gpfFile, flexFn=flexFn)
        
      else:
          scoreFunc = scoreFunc if not self.doZnDock.get() and not flexFn else 'ad4'
          args = "-p {} -l {}.glg".format(gpfFile, self.getReceptorName())
          insistentRun(self, "autogrid4", args, envDic=SCRUBBER_DIC, cwd=outDir)

          batchDirs = self.getBatchDirs(pdbqtFiles)
          for molDir in batchDirs:
              paramsFile = self.writeConfigFile(radius, [xCenter, yCenter, zCenter],
                                                outDir, nThreads, scoreFunc, flexFn)
              args = "--batch {}/*.pdbqt --maps {} --config {}".format(molDir, self.getReceptorName(), paramsFile)  # batch cannot be read from config
              self.runJob(pwchem_plugin.getEnvPath(VINA_DIC, 'bin/vina'), args, cwd=outDir)

    def createOutputStep(self):
      recFile = self.getReceptorPDBQT()
      if self.ringtailOutput.get():
        nt = self.numberOfThreads.get()
        outDir = os.path.abspath(self._getExtraPath())
        args = f'write --file_path {outDir} --recursive -o ringtail.db -m vina -sr -rf {os.path.abspath(recFile)} ' \
               f'-mpr {nt} --overwrite'
        autodockPlugin.runRingtail(self, args, cwd=self._getPath())

        outputDB = RingtailDatabase(filename=self._getPath('ringtail.db'))
        outputDB.setReceptorFile(recFile)
        outputDB.setType(VINA)
        outputDB.createSumFile(self.getSumPath())
        outputDB.performBaseFilter()
        self._defineOutputs(outputRingtail=outputDB)
      else:
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
                    if self.doFlexRes:
                      poseFile, curRecFile = self.makeFlexPoseFiles(poseFile, recFile)
                      newSmallMol.setProteinFile(os.path.relpath(curRecFile))

                    filename = f'g{gridId}_{os.path.split(poseFile)[-1]}'
                    newPoseFile = os.path.join(outDir, filename)
                    os.rename(poseFile, newPoseFile)

                    newSmallMol.poseFile.set(newPoseFile)
                    newSmallMol.setPoseId(posId)
                    newSmallMol.gridId.set(gridId)
                    newSmallMol.setMolClass('AutodockVina')
                    newSmallMol.setDockId(self.getObjId())

                    outputSet.append(newSmallMol)
              else:
                print(f'Molecule {molName} was not found in the docking results')

        outputSet.updateMolClass()
        outputSet.setProteinFile(recFile)
        outputSet.setDocked(True)
        self._defineOutputs(outputSmallMolecules=outputSet)
        self._defineSourceRelation(self.inputSmallMolecules, outputSet)

      self.cleanTmpFiles()

    ########################### Parameters functions ############################

    def getnThreads(self):
        '''Get the number of threads available for each pocket execution'''
        nThreads = self.numberOfThreads.get()
        if self.fromReceptor.get() != 0:
            nPockets = len(self.inputStructROIs.get())
            nThreads = nThreads // nPockets
        nThreads = 1 if nThreads == 0 else nThreads
        return nThreads

    def writeParamsFile(self, fnReceptor, molFiles, radius, center, gpfFile, outDir, nCPUs, it, flexFn=None):
        paramsFile = os.path.join(outDir, f'inputParams_{it}.txt')

        f = open(paramsFile, 'w')
        f.write('ligandFiles:: {}\n'.format(' '.join(molFiles)))
        f.write('receptorFile:: {}\n'.format(fnReceptor))
        if flexFn:
            f.write('flexRecFile:: {}\n'.format(flexFn))
        f.write('mapsName:: {}\n'.format(self.getReceptorName()))
        f.write('gpfFile:: {}\n'.format(gpfFile))

        f.write('boxSize:: {}\n'.format([r*2 for r in radius]))
        f.write('boxCenter:: {}\n'.format(center))

        f.write('scoreName:: {}\n'.format(self.getEnumText('scoreName')))
        f.write('exhaust:: {}\n'.format(self.exhaust.get()))
        f.write('nPoses:: {}\n'.format(self.nRuns.get()))
        f.write('minRMSD:: {}\n'.format(self.rmsTol.get()))
        f.write('maxEvals:: {}\n'.format(self.maxEvals.get()))

        f.write('nCPUs:: {}\n'.format(nCPUs))
        f.write('outDir:: {}\n'.format(outDir))
        f.write('it:: {}\n'.format(it))

        return paramsFile

    def writeConfigFile(self, radius, center, outDir, nCPUs, scoring, flexFn=None):
        paramsFile = os.path.join(outDir, 'inputParams.txt')

        f = open(paramsFile, 'w')
        if flexFn:
            f.write('flex = {}\n'.format(flexFn))
        f.write('scoring = {}\n'.format(scoring))

        f.write('center_x = {}\ncenter_y = {}\ncenter_z = {}\n'.format(*center))
        f.write('size_x = {}\nsize_y = {}\nsize_z = {}\n'.format(*3*[radius*2]))

        f.write('dir = {}\n'.format(outDir))

        f.write('cpu = {}\n'.format(nCPUs))
        f.write('exhaustiveness = {}\n'.format(self.exhaust.get()))
        f.write('max_evals = {}\n'.format(self.maxEvals.get()))
        f.write('num_modes = {}\n'.format(self.nRuns.get()))
        f.write('min_rmsd = {}\n'.format(self.rmsTol.get()))
        f.write('verbosity = 2\n')

        return paramsFile

########################### Utils functions ############################

    def performScriptDocking(self, pdbqtFiles, it, outDir, **kwargs):
      k = kwargs
      paramsFile = self.writeParamsFile(k['fnReceptor'], pdbqtFiles, k['radius'],
                                        [k['xCenter'], k['yCenter'], k['zCenter']], k['gpfFile'],
                                        outDir, 1, it, k['flexFn'])
      autodockPlugin.runScript(self, scriptName, paramsFile, envDict=VINA_DIC, cwd=outDir)


    def getBatchDirs(self, molFiles):
        ds = []
        for mf in molFiles:
            ds.append(os.path.dirname(mf))
        return list(set(ds))

    def getDockedLigandsFiles(self, outDir):
        dockFiles = []
        for file in os.listdir(outDir):
            if 'docked_files' in file and not "failed" in file:
                dockedFilesFile = os.path.join(outDir, file)
                with open(dockedFilesFile) as fIn:
                    dockFiles += fIn.read().split()
        return dockFiles

    def getSumPath(self):
      return os.path.abspath(self._getExtraPath('ringSum.txt'))

    def _summary(self):
      s = []
      if os.path.exists(self.getSumPath()):
        with open(self.getSumPath()) as f:
          s.append(f.read())
      return s

    def _warnings(self):
      ws = []
      if self.doZnDock.get() and self.getEnumText('scoreName').lower() != 'ad4':
          ws.append('Zn docking can only be performed with autodock scoring, we will use it here.')
      return ws

