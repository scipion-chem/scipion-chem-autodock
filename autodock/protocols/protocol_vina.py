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

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, IntParam, FloatParam, STEPS_PARALLEL, EnumParam, LEVEL_ADVANCED
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath, createLink

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import runOpenBabel, calculate_centerMass, getBaseFileName
from pwchem import Plugin as pwchem_plugin
from pwchem.constants import MGL_DIC

from autodock import Plugin as autodock_plugin
from autodock.protocols.protocol_autodock import _defineFlexParams

class ProtChemVina(EMProtocol):
    """Perform a docking experiment with autodock vina."""
    _label = 'Autodock Vina docking'
    _program = ""

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Receptor specification')
        group.addParam('fromReceptor', EnumParam, label='Dock on : ', default=1,
                       choices=['Whole protein', 'SetOfStructROIs'],
                       help='Whether to dock on a whole protein surface or on specific regions')

        #Docking on whole protein
        group.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct",
                      label='Input atomic structure:', condition='fromReceptor == 0',
                      help="The atom structure to use as receptor in the docking")
        group.addParam('radius', FloatParam, label='Grid radius for whole protein: ',
                       condition='fromReceptor == 0', allowsNull=False,
                       help='Radius of the Autodock grid for the whole protein')

        #Docking on pockets
        group.addParam('inputStructROIs', PointerParam, pointerClass="SetOfStructROIs",
                      label='Input StructROIs:', condition='not fromReceptor == 0',
                      help="The protein StructROIs to dock in")
        group.addParam('pocketRadiusN', FloatParam, label='Grid radius vs StructROI radius: ',
                       condition='not fromReceptor == 0', default=1.1, allowsNull=False,
                       help='The radius * n of each StructROI will be used as grid radius')

        _defineFlexParams(form)

        group = form.addGroup('Docking')
        group.addParam('inputLibrary', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Ligand library:', allowsNull=False,
                       help="The library must be prepared for autodock")

        group.addParam('nPos', IntParam, label='Number of positions per ligand: ', default=9,
                       help='Number of positions outputed for each of the ligands provided')
        group.addParam('maxEDiff', FloatParam, label='Maximum energy difference: ',
                       default=3.0, expertLevel=LEVEL_ADVANCED,
                       help='Maximum energy difference between the best and the worst ligand positions')
        group.addParam('exhaust', IntParam, label='Exhaustiveness number: ', default=8,
                       help='Number that controls the amount of runs or trials to use in the vina algorithm.\n'
                            'The higher, the more the search pace is explored, but also it will be slower')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        cId = self._insertFunctionStep('convertStep', prerequisites=[])

        dockSteps = []
        if self.fromReceptor.get() == 0:
            for mol in self.inputLibrary.get():
                dockId = self._insertFunctionStep('dockStep', mol.clone(), prerequisites=[cId])
                dockSteps.append(dockId)
        else:
            for pocket in self.inputStructROIs.get():
                for mol in self.inputLibrary.get():
                    dockId = self._insertFunctionStep('dockStep', mol.clone(), pocket.clone(), prerequisites=[cId])
                    dockSteps.append(dockId)

        self._insertFunctionStep('createOutputStep', prerequisites=dockSteps)

    def convertStep(self):
        self.ligandFileNames = []
        for mol in self.inputLibrary.get():
            fnSmall, smallDir = self.convert2PDBQT(mol.clone(), self._getExtraPath('conformers'))
            self.ligandFileNames.append(fnSmall)

        receptorFile = self.getOriginalReceptorFile()
        if receptorFile.endswith('.pdb'):
            self.pdbqtFile = self.convertReceptor2PDBQT(receptorFile)
            self.pdbFile = receptorFile
        elif receptorFile.endswith('.pdbqt'):
            self.pdbqtFile = receptorFile
            self.pdbFile = self.convertReceptor2PDB(receptorFile)


    def dockStep(self, mol, pocket=None):
        fnReceptor = os.path.abspath(self.pdbqtFile)
        molFn = self.getMolLigandName(mol)
        molName, _ = os.path.splitext(os.path.basename(molFn))

        outDir = self.getOutputPocketDir(pocket)
        smallDir = os.path.join(outDir, molName)
        makePath(smallDir)

        if self.fromReceptor.get() == 0:
            radius = self.radius.get()
            # Use the original pdb for mass center
            structure, x_center, y_center, z_center = calculate_centerMass(self.pdbFile)
        else:
            radius = (pocket.getDiameter() / 2) * self.pocketRadiusN.get()
            x_center, y_center, z_center = pocket.calculateMassCenter()

        if self.doFlexRes:
            flexFn, fnReceptor = self.buildFlexReceptor(fnReceptor)

        args = '--receptor {} --ligand {} '.format(fnReceptor, molFn)
        if self.doFlexRes:
            args += '--flexRec {} '.format(flexFn)
        args += '--center_x {} --center_y {} --center_z {} '.format(x_center, y_center, z_center)
        args += '--size_x {} --size_y {} --size_z {} '.format(radius*2, radius*2, radius*2)
        args += '--exhaustiveness {} --energy_range {} --num_modes {} --cpu {} '.\
            format(self.exhaust.get(), self.maxEDiff.get(), self.nPos.get(), self.numberOfThreads.get())

        outName = os.path.abspath(os.path.join(smallDir, molName))
        args += '--out {}.pdbqt --log {}.log '.format(outName, outName)

        autodock_plugin.runVina(self, args=args, cwd=smallDir)

    def createOutputStep(self):
        outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())

        for pocketDir in self.getPocketDirs():
            gridId = self.getGridId(pocketDir)

            for smallMol in self.inputLibrary.get():
                smallFn = smallMol.getFileName()
                molName = os.path.splitext(os.path.basename(smallFn))[0]
                fnSmallDir = os.path.join(pocketDir, molName)

                molDic = self.parseDockedPDBQT(os.path.join(fnSmallDir, molName+'.pdbqt'))

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
        self._defineSourceRelation(self.inputLibrary, outputSet)
      
########################### Utils functions ############################

    def getGridId(self, outDir):
        return outDir.split('_')[-1]

    def convert2PDBQT(self, smallMol, oDir):
        '''Convert ligand to pdbqt using obabel'''
        inFile = smallMol.getFileName()
        if os.path.splitext(inFile)[1] not in ['.pdb', '.mol2', '.pdbq']:
            # Convert to formats recognized by ADT
            outName, outDir = os.path.splitext(os.path.basename(inFile))[0], os.path.abspath(self._getTmpPath())
            args = ' -i "{}" -of mol2 --outputDir "{}" --outputName {}'.format(os.path.abspath(inFile),
                                                                               os.path.abspath(outDir), outName)
            pwchem_plugin.runScript(self, 'obabel_IO.py', args, env='plip', cwd=outDir)
            inFile = self._getTmpPath(outName + '.mol2')

        if not os.path.exists(oDir):
            os.mkdir(oDir)

        inName, inExt = os.path.splitext(os.path.basename(inFile))
        oFile = os.path.abspath(os.path.join(oDir, smallMol.getUniqueName() +'.pdbqt'))

        if inExt != '.pdbqt':
            args = ' -l {} -o {}'.format(inFile, oFile)
            self.runJob(pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh'),
                        autodock_plugin.getADTPath('Utilities24/prepare_ligand4.py') + args)
        else:
            createLink(inFile, oFile)
        return oFile, oDir

    def convertReceptor2PDB(self, proteinFile):
        inName, inExt = os.path.splitext(os.path.basename(proteinFile))
        oFile = os.path.abspath(os.path.join(self._getTmpPath(inName + '.pdb')))
        if not os.path.exists(oFile):
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

    def getReceptorName(self):
        if self.fromReceptor.get() == 0:
            atomStructFn = self.getOriginalReceptorFile()
            return atomStructFn.split('/')[-1].split('.')[0]
        else:
            return self.inputStructROIs.get().getProteinName()

    def getOriginalReceptorFile(self):
        if self.fromReceptor.get() == 0:
            return self.inputAtomStruct.get().getFileName()
        else:
            return self.inputStructROIs.get().getProteinFile()

    def getReceptorDir(self):
        atomStructFn = self.getOriginalReceptorFile()
        return '/'.join(atomStructFn.split('/')[:-1])

    def getLigandsFileNames(self):
        ligFns = []
        for mol in self.inputLibrary.get():
            ligFns.append(mol.getFileName())
        return ligFns

    def getOutputPocketDir(self, pocket=None):
        if pocket==None:
            outDir = self._getExtraPath('pocket_1')
        else:
            outDir = self._getExtraPath('pocket_{}'.format(pocket.getObjId()))
        return outDir

    def buildFlexReceptor(self, receptorFn):
        molName = getBaseFileName(receptorFn)
        allFlexRes = self.parseFlexRes(molName)
        flexFn, rigFn = self.getFlexFiles()
        args = ' -r {} -s {} -g {} -x {}'.format(receptorFn, allFlexRes, rigFn, flexFn)
        self.runJob(pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh'),
                    autodock_plugin.getADTPath('Utilities24/prepare_flexreceptor4.py') + args, cwd=self._getExtraPath())
        return flexFn, rigFn

    def getFlexFiles(self):
        return os.path.abspath(self._getExtraPath('receptor_flex.pdbqt')), \
               os.path.abspath(self._getExtraPath('receptor_rig.pdbqt'))

    def parseFlexRes(self, molName):
        allFlexDic = {}
        for line in self.flexList.get().split('\n'):
            if line.strip():
                chain, res = line.strip().split(':')
                if chain in allFlexDic:
                    allFlexDic[chain] += res + '_'
                else:
                    allFlexDic[chain] = res + '_'

        allFlexStr = ''
        for chain in allFlexDic:
            allFlexStr += '{}:{}:{},'.format(molName, chain, allFlexDic[chain][:-1])
        return allFlexStr[:-1]

    def getPocketDirs(self):
        dirs = []
        for file in os.listdir(self._getExtraPath()):
            d = self._getExtraPath(file)
            if os.path.isdir(d) and 'pocket' in file:
                dirs.append(d)
        dirs.sort()
        return dirs

    def getMolLigandName(self, mol):
        molName = mol.getUniqueName()
        for ligName in self.ligandFileNames:
            if molName + '.pdbqt' in ligName:
                return ligName

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


    def _citations(self):
        return ['Morris2009']
