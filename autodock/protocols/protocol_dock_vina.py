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

import os, shutil

from pyworkflow.protocol import params
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath

from pwchem import Plugin as pwchem_plugin
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import calculate_centerMass, generate_gpf, insistentRun, getBaseFileName, performBatchThreading

from autodock import Plugin as autodock_plugin
from autodock.protocols.protocol_autodock import ProtChemAutodockBase
from autodock.constants import VINA_DIC

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

        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
      cId = self._insertFunctionStep('convertStep', prerequisites=[])
      self.receptorName = self.getReceptorName()

      dockSteps = []
      if self.fromReceptor.get() == 0:
          dockId = self._insertFunctionStep('dockStep', prerequisites=[cId])
          dockSteps.append(dockId)
      else:
        for pocket in self.inputStructROIs.get():
            dockId = self._insertFunctionStep('dockStep', pocket.clone(), prerequisites=[cId])
            dockSteps.append(dockId)

      self._insertFunctionStep('createOutputStep', prerequisites=dockSteps)

    def convertStep(self):
      inputMols, nt = self.inputSmallMolecules.get(), self.numberOfThreads.get()
      molLists = performBatchThreading(self.performLigConversion, inputMols, nt)
      with open(self.getConvertedLigandsFile(), 'w') as f:
        f.write('\n'.join(molLists))

      receptorFile = self.getOriginalReceptorFile()
      if receptorFile.endswith(PDBext):
        self.convertReceptor2PDBQT(receptorFile)
        shutil.copy(receptorFile, self.getReceptorPDB())
      elif receptorFile.endswith(PDBQText):
        self.convertReceptor2PDB(receptorFile)
        shutil.copy(receptorFile, self.getReceptorPDBQT())

    def dockStep(self, pocket=None):
      fnReceptor = self.getReceptorPDBQT()
      pdbqtFiles = self.getConvertedLigandsFiles()

      outDir = self.getOutputPocketDir(pocket)
      if not os.path.exists(outDir):
          makePath(outDir)

      if self.fromReceptor.get() == 0:
          radius = self.radius.get()
          # Use the original pdb for mass center
          _, xCenter, yCenter, zCenter = calculate_centerMass(self.getReceptorPDB())
      else:
          radius = (pocket.getDiameter() / 2) * self.pocketRadiusN.get()
          xCenter, yCenter, zCenter = pocket.calculateMassCenter()

      spacing = self.spacing.get()
      npts = (radius * 2) / spacing

      znFFfile = autodock_plugin.getPackagePath(package='VINA', path='AutoDock-Vina/data/AD4Zn.dat') \
        if self.doZnDock.get() else None
      gpfFile = generate_gpf(fnReceptor, spacing=spacing,
                              xc=xCenter, yc=yCenter, zc=zCenter,
                              npts=npts, outDir=outDir, ligandFns=pdbqtFiles, znFFfile=znFFfile)

      if self.doFlexRes:
          flexFn, fnReceptor = self.buildFlexReceptor(fnReceptor, cleanZn=self.doZnDock.get())
      else:
          flexFn = None

      nThreads = self.getnThreads()
      scoreFunc = self.getEnumText('scoreName').lower()
      if not self.doZnDock.get() and scoreFunc == 'vina':
        performBatchThreading(self.performScriptDocking, pdbqtFiles, nThreads, cloneItem=False, outDir=outDir,
                              fnReceptor=fnReceptor, radius=radius, xCenter=xCenter, yCenter=yCenter, zCenter=zCenter,
                              gpfFile=gpfFile, nThreads=nThreads, flexFn=flexFn)

      else:
          scoreFunc = scoreFunc if not self.doZnDock.get() and not flexFn else 'ad4'
          args = "-p {} -l {}.glg".format(gpfFile, self.getReceptorName())
          insistentRun(self, autodock_plugin.getPackagePath("AUTOSITE", path='bin/autogrid4'), args, cwd=outDir)

          batchDirs = self.getBatchDirs(pdbqtFiles)
          for molDir in batchDirs:
              paramsFile = self.writeConfigFile(radius, [xCenter, yCenter, zCenter],
                                                outDir, nThreads, scoreFunc, flexFn)
              args = "--batch {}/*.pdbqt --maps {} --config {}".format(molDir, self.getReceptorName(), paramsFile)  # batch cannot be read from config
              self.runJob(pwchem_plugin.getCondaEnvPath('vina', path='bin/vina'), args, cwd=outDir)

    def performScriptDocking(self, pdbqtFiles, molLists, it, **kwargs):
      k = kwargs
      paramsFile = self.writeParamsFile(k['fnReceptor'], pdbqtFiles, k['radius'],
                                        [k['xCenter'], k['yCenter'], k['zCenter']], k['gpfFile'],
                                        k['outDir'], k['nThreads'], k['flexFn'])
      autodock_plugin.runScript(self, scriptName, paramsFile, envDict=VINA_DIC, cwd=k['outDir'], popen=True)


    def getBatchDirs(self, molFiles):
        ds = []
        for mf in molFiles:
            ds.append(os.path.dirname(mf))
        return list(set(ds))

    def createOutputStep(self):
      outDir = self._getPath('outputLigands')
      makePath(outDir)
      outputSet = SetOfSmallMolecules().create(outputPath=outDir)

      for pocketDir in self.getPocketDirs():
        pocketDic = {}
        gridId = self.getGridId(pocketDir)
        for dockFile in self.getDockedLigandsFiles(pocketDir):
            molName = os.path.split(dockFile)[1].split(PDBQText)[0]
            pocketDic[molName] = self.parseDockedPDBQT(dockFile)

        for smallMol in self.inputSmallMolecules.get():
            molName = smallMol.getUniqueName(conf=True)
            molDic = pocketDic[molName]

            for posId in molDic:
              newSmallMol = SmallMolecule()
              newSmallMol.copy(smallMol, copyId=False)
              newSmallMol._energy = pwobj.Float(molDic[posId]['energy'])
              if os.path.getsize(molDic[posId]['file']) > 0:
                newSmallMol.poseFile.set(os.path.relpath(molDic[posId]['file']))
                newSmallMol.setPoseId(posId)
                newSmallMol.gridId.set(gridId)
                newSmallMol.setMolClass('AutodockVina')
                newSmallMol.setDockId(self.getObjId())

                outputSet.append(newSmallMol)

      outputSet.proteinFile.set(self.getOriginalReceptorFile())
      outputSet.setDocked(True)
      self._defineOutputs(outputSmallMolecules=outputSet)
      self._defineSourceRelation(self.inputSmallMolecules, outputSet)

    ########################### Parameters functions ############################

    def getnThreads(self):
        '''Get the number of threads available for each pocket execution'''
        nThreads = self.numberOfThreads.get()
        if self.fromReceptor.get() != 0:
            nPockets = len(self.inputStructROIs.get())
            nThreads = nThreads // nPockets
        nThreads = 1 if nThreads == 0 else nThreads
        return nThreads

    def performLigConversion(self, inMols, molLists, it):
      for mol in inMols:
        molFile = mol.getFileName()
        if not molFile.endswith(PDBQText):
          fnSmall = self.convertLigand2PDBQT(mol, self._getTmpPath(), popen=True)[0]
        else:
          fnSmall = os.path.abspath(self._getTmpPath(getBaseFileName(molFile + PDBQText)))
          shutil.copy(molFile, fnSmall)
        molLists[it].append(fnSmall)

    def writeParamsFile(self, fnReceptor, molFiles, radius, center, gpfFile, outDir, nCPUs, flexFn=None):
        paramsFile = os.path.join(outDir, 'inputParams.txt')

        f = open(paramsFile, 'w')
        f.write('ligandFiles:: {}\n'.format(' '.join(molFiles)))
        f.write('receptorFile:: {}\n'.format(fnReceptor))
        if flexFn:
            f.write('flexRecFile:: {}\n'.format(flexFn))
        f.write('mapsName:: {}\n'.format(self.getReceptorName()))
        f.write('gpfFile:: {}\n'.format(gpfFile))

        f.write('boxSize:: {}\n'.format(3*[radius*2]))
        f.write('boxCenter:: {}\n'.format(center))

        f.write('scoreName:: {}\n'.format(self.getEnumText('scoreName')))
        f.write('exhaust:: {}\n'.format(self.exhaust.get()))
        f.write('nPoses:: {}\n'.format(self.nRuns.get()))
        f.write('minRMSD:: {}\n'.format(self.rmsTol.get()))
        f.write('maxEvals:: {}\n'.format(self.maxEvals.get()))

        f.write('nCPUs:: {}\n'.format(nCPUs))
        f.write('outDir:: {}\n'.format(outDir))

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

    def parseDockedPDBQT(self, pdbqtFile):
        dockedDic = {}
        towrite = ''
        with open(pdbqtFile) as fIn:
            for line in fIn:
                if line.startswith('MODEL'):
                    if towrite != '':
                        newFile = pdbqtFile.replace(PDBQText, '_{}.pdbqt'.format(modelId))
                        dockedDic[modelId] = {'file': newFile, 'energy': energy}
                        with open(newFile, 'w') as f:
                            f.write(towrite)
                    towrite = ''
                    modelId = line.strip().split()[1]
                elif line.startswith('REMARK VINA RESULT:'):
                    energy = line.split()[3]
                else:
                    towrite += line
        if towrite:
            newFile = pdbqtFile.replace(PDBQText, '_{}.pdbqt'.format(modelId))
            dockedDic[modelId] = {'file': newFile, 'energy': energy}
            with open(newFile, 'w') as f:
                f.write(towrite)
        return dockedDic

    def getMeekoFiles(self, oDir):
        with open(os.path.join(oDir, 'meeko_files.txt')) as fIn:
          prepFiles = fIn.read().split()
        return prepFiles
    
    def getConvertedLigandsFile(self):
        return os.path.abspath(self._getExtraPath('inputLigands.txt'))

    def getConvertedLigandsFiles(self):
        with open(self.getConvertedLigandsFile()) as f:
            a = f.read()
        return a.split('\n')

    def getDockedLigandsFiles(self, outDir):
        dockedFilesFile = os.path.join(outDir, 'docked_files.txt')
        if os.path.exists(dockedFilesFile):
            with open(dockedFilesFile) as fIn:
                files = fIn.read().split()
        else:
            files = []
            for f in os.listdir(outDir):
                if '_out.pdbqt' in f:
                    inFile, outFile = os.path.join(outDir, f), os.path.join(outDir, f.replace('_out.pdbqt', PDBQText))
                    os.rename(inFile, outFile)
                    files.append(outFile)
        return files

    def _warnings(self):
      ws = []
      if self.doZnDock.get() and self.getEnumText('scoreName').lower() != 'ad4':
          ws.append('Zn docking can only be performed with autodock scoring, we will use it here.')
      return ws

