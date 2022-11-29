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

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, BooleanParam, EnumParam, IntParam, FloatParam, \
  LEVEL_ADVANCED, STEPS_PARALLEL
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath, createLink

from autodock import Plugin as autodock_plugin
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import runOpenBabel, calculate_centerMass, generate_gpf
from pwchem import Plugin as pwchem_plugin
from pwchem.constants import MGL_DIC

meekoScript = 'meeko_preparation.py'
scriptName = 'vina_docking.py'

class ProtChemVinaDocking(EMProtocol):
    """Dock ligands using Vina"""
    _label = 'Vina docking'
    _program = ""

    _scoreNames = ['Vina', 'Vinardo', 'AD4']

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addHidden('nCPUs', IntParam, label='Number of CPUs to use', default=1)
        group = form.addGroup('Receptor specification')
        group.addParam('wholeProt', BooleanParam, label='Dock on whole protein: ', default=True,
                       help='Whether to dock on a whole protein surface or on specific regions')

        # Docking on whole protein
        group.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct",
                       label='Input atomic structure:', condition='wholeProt',
                       help="The atom structure to use as receptor in the docking")
        group.addParam('radius', FloatParam, label='Grid radius for whole protein: ',
                       condition='wholeProt', allowsNull=False,
                       help='Radius of the Autodock grid for the whole protein')

        # Docking on pockets
        group.addParam('inputStructROIs', PointerParam, pointerClass="SetOfStructROIs",
                       label='Input StructROIs:', condition='not wholeProt',
                       help="The protein StructROIs to dock in")
        group.addParam('pocketRadiusN', FloatParam, label='Grid radius vs StructROI radius: ',
                       condition='not wholeProt', default=1.1, allowsNull=False,
                       help='The radius * n of each StructROI will be used as grid radius')
        group.addParam('spacing', FloatParam, label='Spacing of affinity maps: ',
                       default=0.375, help='Affinity maps spacing')

        group = form.addGroup('Docking')
        group.addParam('inputSmallMols', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Set of small molecules:', allowsNull=False,
                       help='Input small molecules to be docked with Vina')

        group.addParam('scoreName', EnumParam, label='Score function: ',
                       choices=self._scoreNames, default=0,
                       help='Score function to use for docking')
        group.addParam('exhaust', IntParam, label='Exhaustiveness number: ', default=8,
                       help='Number that controls the amount of runs or trials to use in the vina algorithm.\n'
                            'The higher, the more the search pace is explored, but also it will be slower')
        group.addParam('nPoses', IntParam, label='Number of positions per ligand: ', default=20,
                       help='Number of positions outputed for each of the ligands provided')
        group.addParam('minRMSD', FloatParam, label='Minimum RMSD difference between poses: ',
                       default=1.0, help='Minimum RMSD difference between poses')
        group.addParam('maxEvals', IntParam, label='Maximum number of evaluation: ',
                       expertLevel=LEVEL_ADVANCED, default=0, help='Maximum number of evaluation')

        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
      cId = self._insertFunctionStep('convertStep', prerequisites=[])

      dockSteps = []
      if self.wholeProt:
          dockId = self._insertFunctionStep('dockStep', prerequisites=[cId])
          dockSteps.append(dockId)
      else:
        for pocket in self.inputStructROIs.get():
            dockId = self._insertFunctionStep('dockStep', pocket.clone(), prerequisites=[cId])
            dockSteps.append(dockId)

      self._insertFunctionStep('createOutputStep', prerequisites=dockSteps)

    def convertStep(self):
      mols = self.inputSmallMols.get()
      self.convert2PDBQT(mols, self._getTmpPath())

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

      if self.wholeProt:
        radius = self.radius.get()
        # Use the original pdb for mass center
        structure, x_center, y_center, z_center = calculate_centerMass(self.getReceptorPDB())
      else:
        radius = (pocket.getDiameter() / 2) * self.pocketRadiusN.get()
        x_center, y_center, z_center = pocket.calculateMassCenter()

      mols = self.inputSmallMols.get()
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

        for smallMol in self.inputSmallMols.get():
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
      self._defineSourceRelation(self.inputSmallMols, outputSet)

    ########################### Parameters functions ############################

    def writeMeekoParamsFile(self, molsScipion, oDir):
        paramsFile = os.path.abspath(self._getExtraPath('inputParams.txt'))

        molFiles = []
        f = open(paramsFile, 'w')
        for mol in molsScipion:
            molFiles.append(os.path.abspath(mol.getFileName()))

        f.write('ligandFiles:: {}\n'.format(' '.join(molFiles)))
        f.write('keepNonPolar:: False\n')
        f.write('hydrate:: False\n')

        f.write('outDir:: {}\n'.format(os.path.abspath(oDir)))

        return paramsFile

    def writeParamsFile(self, molFiles, radius, center, gpfFile, outDir):
        paramsFile = os.path.join(outDir, 'inputParams.txt')

        f = open(paramsFile, 'w')
        f.write('ligandFiles:: {}\n'.format(' '.join(molFiles)))
        f.write('receptorFile:: {}\n'.format(self.getReceptorPDBQT()))
        f.write('mapsName:: {}\n'.format(self.getReceptorName()))
        f.write('gpfFile:: {}\n'.format(gpfFile))

        f.write('boxSize:: {}\n'.format(3*[radius*2]))
        f.write('boxCenter:: {}\n'.format(center))

        f.write('scoreName:: {}\n'.format(self.getEnumText('scoreName')))
        f.write('exhaust:: {}\n'.format(self.exhaust.get()))
        f.write('nPoses:: {}\n'.format(self.nPoses.get()))
        f.write('minRMSD:: {}\n'.format(self.minRMSD.get()))
        f.write('maxEvals:: {}\n'.format(self.maxEvals.get()))

        f.write('nCPUs:: {}\n'.format(self.numberOfThreads.get()))
        f.write('outDir:: {}\n'.format(outDir))

        return paramsFile

########################### Utils functions ############################

    def getGridId(self, outDir):
        return outDir.split('_')[-1]

    def convert2PDBQT(self, mols, oDir):
        '''Convert ligand to pdbqt using meeko'''
        if not os.path.exists(oDir):
            os.mkdir(oDir)

        write = 'w'
        noPDBQT_mols, oFiles = [], []
        for mol in mols:
            inFile = mol.getFileName()
            inName, inExt = os.path.splitext(os.path.basename(inFile))
            oFile = os.path.abspath(os.path.join(oDir, mol.getUniqueName() + '.pdbqt'))

            if inExt != '.pdbqt':
                noPDBQT_mols.append(mol)
            else:
                createLink(inFile, oFile)
                oFiles.append(oFile)

        if len(noPDBQT_mols) > 0:
            paramsFile = self.writeMeekoParamsFile(noPDBQT_mols, oDir)
            autodock_plugin.runScript(self, meekoScript, paramsFile, env='rdkit', cwd=oDir)
            meekoResults = os.path.join(oDir, 'meeko_files.txt')
            shutil.copy(meekoResults, self.getConvertedLigandsFile())
            write = 'a'

        with open(self.getConvertedLigandsFile(), write) as f:
            for oFile in oFiles:
                f.write('\n' + oFile)


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

    def getOriginalReceptorFile(self):
        if self.wholeProt:
            return self.inputAtomStruct.get().getFileName()
        else:
            return self.inputStructROIs.get().getProteinFile()

    def getOutputPocketDir(self, pocket=None):
        if pocket==None:
            outDir = self._getExtraPath('pocket_1')
        else:
            outDir = self._getExtraPath('pocket_{}'.format(pocket.getObjId()))
        return os.path.abspath(outDir)

    def getPocketDirs(self):
        dirs = []
        for file in os.listdir(self._getExtraPath()):
            d = self._getExtraPath(file)
            if os.path.isdir(d) and 'pocket' in file:
                dirs.append(d)
        dirs.sort()
        return dirs

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

    def getReceptorName(self):
        if self.wholeProt:
            atomStructFn = self.getOriginalReceptorFile()
            return atomStructFn.split('/')[-1].split('.')[0]
        else:
            return self.inputStructROIs.get().getProteinName()
