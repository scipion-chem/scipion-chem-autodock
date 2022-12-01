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

from pyworkflow.protocol.params import PointerParam, EnumParam, IntParam, FloatParam, \
  LEVEL_ADVANCED
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import runOpenBabel, calculate_centerMass, generate_gpf
from pwchem import Plugin as pwchem_plugin
from pwchem.constants import MGL_DIC

from autodock import Plugin as autodock_plugin
from autodock.protocols.protocol_autodock import ProtChemAutodockBase

meekoScript = 'meeko_preparation.py'
scriptName = 'vina_docking.py'

class ProtChemVinaDocking(ProtChemAutodockBase):
    """Dock ligands using Vina"""
    _label = 'Vina docking'
    _program = ""

    _scoreNames = ['Vina', 'Vinardo', 'AD4']

    def _defineParams(self, form):
        super()._defineParams(form)
        group = form.addGroup('Docking')
        group.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Input small molecules: ', allowsNull=False,
                       help="Input small molecules to be docked with AutoDock")
        group.addParam('nRuns', IntParam, label='Number of docking runs: ', default=50,
                       help='Number of independent runs using the selected strategy. \nDifferent docking positions will be '
                            'found for each of them.')
        group.addParam('rmsTol', FloatParam, label='Cluster tolerance (A): ', default=2.0, expertLevel=LEVEL_ADVANCED,
                       help='Maximum RMSD for 2 docked structures to be in the same cluster')

        group.addParam('scoreName', EnumParam, label='Score function: ',
                       choices=self._scoreNames, default=0,
                       help='Score function to use for docking')
        group.addParam('exhaust', IntParam, label='Exhaustiveness number: ', default=8,
                       help='Number that controls the amount of runs or trials to use in the vina algorithm.\n'
                            'The higher, the more the search pace is explored, but also it will be slower')
        group.addParam('maxEvals', IntParam, label='Maximum number of evaluation: ',
                       expertLevel=LEVEL_ADVANCED, default=0, help='Maximum number of evaluation')

        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
      cId = self._insertFunctionStep('convertStep', prerequisites=[])

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
      mols = self.inputSmallMolecules.get()
      for mol in mols:
        self.convert2PDBQT(mol, self._getTmpPath())

      receptorFile = self.getOriginalReceptorFile()
      if receptorFile.endswith('.pdb'):
        self.convertReceptor2PDBQT(receptorFile)
        shutil.copy(receptorFile, self.getReceptorPDB())
      elif receptorFile.endswith('.pdbqt'):
        self.convertReceptor2PDB(receptorFile)
        shutil.copy(receptorFile, self.getReceptorPDBQT())

    def dockStep(self, pocket=None):
      outDir = self.getOutputPocketDir(pocket)
      if not os.path.exists(outDir):
          makePath(outDir)

      if self.fromReceptor.get() == 0:
        radius = self.radius.get()
        # Use the original pdb for mass center
        structure, x_center, y_center, z_center = calculate_centerMass(self.getReceptorPDB())
      else:
        radius = (pocket.getDiameter() / 2) * self.pocketRadiusN.get()
        x_center, y_center, z_center = pocket.calculateMassCenter()

      mols = self.inputSmallMolecules.get()
      molFiles = []
      for mol in mols:
        molFiles.append(os.path.abspath(mol.getFileName()))

      spacing = self.spacing.get()
      npts = (radius * 2) / spacing
      gpf_file = generate_gpf(self.getReceptorPDBQT(), spacing=spacing,
                              xc=x_center, yc=y_center, zc=z_center,
                              npts=npts, outDir=outDir, ligandFns=molFiles)

      pdbqtFiles = self.getConvertedLigandsFiles()
      paramsFile = self.writeParamsFile(pdbqtFiles, radius, [x_center, y_center, z_center], gpf_file, outDir)
      autodock_plugin.runScript(self, scriptName, paramsFile, env='vina', cwd=outDir)

    def createOutputStep(self):
      outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())

      for pocketDir in self.getPocketDirs():
        pocketDic = {}
        gridId = self.getGridId(pocketDir)
        for dockFile in self.getDockedLigandsFiles(pocketDir):
            molName = os.path.split(dockFile)[1].split('.pdbqt')[0]
            pocketDic[molName] = self.parseDockedPDBQT(dockFile)

        for smallMol in self.inputSmallMolecules.get():
            molName = smallMol.getUniqueName(conf=True)
            molDic = pocketDic[molName]

            for posId in molDic:
              newSmallMol = SmallMolecule()
              newSmallMol.copy(smallMol, copyId=False)
              newSmallMol._energy = pwobj.Float(molDic[posId]['energy'])
              if os.path.getsize(molDic[posId]['file']) > 0:
                newSmallMol.poseFile.set(molDic[posId]['file'])
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

    def writeParamsFile(self, molFiles, radius, center, gpfFile, outDir):
        paramsFile = os.path.join(outDir, 'inputParams.txt')

        f = open(paramsFile, 'w')
        f.write('ligandFiles:: {}\n'.format(' '.join(molFiles)))

        fnReceptor = self.getReceptorPDBQT()
        if self.doFlexRes:
            flexFn, fnReceptor = self.buildFlexReceptor(fnReceptor)

        f.write('receptorFile:: {}\n'.format(fnReceptor))
        if self.doFlexRes:
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

        f.write('nCPUs:: {}\n'.format(self.numberOfThreads.get()))
        f.write('outDir:: {}\n'.format(outDir))

        return paramsFile

########################### Utils functions ############################

    def convert2PDBQT(self, smallMol, oDir):
        oFile, oDir = super().convert2PDBQT(smallMol, oDir)

        write = 'w'
        if os.path.exists(self.getConvertedLigandsFile()):
            write = 'a'
        with open(self.getConvertedLigandsFile(), write) as f:
            f.write(oFile + '\n')

        return oFile, oDir

    def convertReceptor2PDB(self, proteinFile):
        inName, inExt = os.path.splitext(os.path.basename(proteinFile))
        oFile = self.getReceptorPDB()
        if not os.path.exists(oFile):
            args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(proteinFile), oFile)
            runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())

        return oFile

    def convertReceptor2PDBQT(self, proteinFile):
        oFile = self.getReceptorPDBQT()

        args = ' -v -r %s -o %s' % (proteinFile, oFile)
        self.runJob(pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh'),
                    autodock_plugin.getADTPath('Utilities24/prepare_receptor4.py') + args)

        return oFile

    def parseDockedPDBQT(self, pdbqtFile):
        dockedDic = {}
        towrite = ''
        with open(pdbqtFile) as fIn:
            for line in fIn:
                if line.startswith('MODEL'):
                    if towrite != '':
                        newFile = pdbqtFile.replace('.pdbqt', '_{}.pdbqt'.format(modelId))
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
            newFile = pdbqtFile.replace('.pdbqt', '_{}.pdbqt'.format(modelId))
            dockedDic[modelId] = {'file': newFile, 'energy': energy}
            with open(newFile, 'w') as f:
                f.write(towrite)
        return dockedDic

    def getMeekoFiles(self, oDir):
        with open(os.path.join(oDir, 'meeko_files.txt')) as fIn:
          prepFiles = fIn.read().split()
        return prepFiles
    
    def getReceptorPDBQT(self):
        return os.path.abspath(os.path.join(self._getExtraPath('receptor.pdbqt')))
    
    def getReceptorPDB(self):
        return os.path.abspath(os.path.join(self._getExtraPath('receptor.pdb')))

    def getConvertedLigandsFile(self):
        return os.path.abspath(self._getExtraPath('inputLigands.txt'))

    def getConvertedLigandsFiles(self):
        with open(self.getConvertedLigandsFile()) as f:
            a = f.read()
        return a.split('\n')

    def getDockedLigandsFiles(self, outDir):
        with open(os.path.join(outDir, 'docked_files.txt')) as fIn:
            files = fIn.read().split()
        return files

