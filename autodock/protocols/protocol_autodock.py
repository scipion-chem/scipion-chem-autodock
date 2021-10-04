# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
from pyworkflow.protocol.params import PointerParam, IntParam, FloatParam, STEPS_PARALLEL, BooleanParam
import pyworkflow.object as pwobj
from autodock import Plugin as autodock_plugin
from pyworkflow.utils.path import makePath, createLink, cleanPattern
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import writePDBLine, splitPDBLine
from autodock.utils import generate_gpf, calculate_centerMass

class ProtChemAutodock(EMProtocol):
    """Perform a docking experiment with autodock. See the help at
       http://autodock.scripps.edu/faqs-help/manual/autodock-4-2-user-guide/AutoDock4.2_UserGuide.pdf"""
    _label = 'autodock'
    _program = ""

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Grids generation')
        group.addParam('wholeProt', BooleanParam, label='Dock on whole protein: ', default=True,
                      help='Whether to dock on a whole protein surface or on specific regions')
        group.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct",
                      label='Input atomic structure:', condition='wholeProt',
                      help="The atom structure to dock in")
        group.addParam('inputPockets', PointerParam, pointerClass="SetOfPockets",
                      label='Input pockets:', condition='not wholeProt',
                      help="The protein pockets to dock in")
        group.addParam('mergeOutput', BooleanParam, default=False,
                       label='Merge outputs from pockets:', condition='not wholeProt',
                       help="Merge the outputs from the different pockets")
        group.addParam('spacing', FloatParam, label='Spacing of the grids: ', default=0.5, allowsNull=False,
                       help='Spacing of the generated Autodock grids')
        group.addParam('radius', FloatParam, label='Grid radius for whole protein: ',
                       condition='wholeProt', allowsNull=False,
                       help='Radius of the Autodock grid for the whole protein')
        group.addParam('pocketRadiusN', FloatParam, label='Grid radius vs pocket radius: ',
                       condition='not wholeProt', default=1.1, allowsNull=False,
                       help='The radius * n of each pocket will be used as grid radius')

        group = form.addGroup('Docking')
        group.addParam('inputLibrary', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Ligand library:', allowsNull=False,
                       help="The library must be prepared for autodock")
        group.addParam('rmsTol', FloatParam, label='Cluster tolerance (A)', default=2.0)
        group.addParam('gaRun', IntParam, label='Number of positions per ligand', default=10)

        form.addSection(label="Genetic algorithm")
        form.addParam('gaPop', IntParam, label='Population size', default=150)
        form.addParam('gaNumEvals', IntParam, label='Number of evaluations', default=2500000)
        form.addParam('gaNumGens', IntParam, label='Number of generations', default=27000)
        form.addParam('gaElitism', IntParam, label='Elitism', default=1)
        form.addParam('gaMutationRate', FloatParam, label='Mutation rate', default=0.02)
        form.addParam('gaCrossOverRate', FloatParam, label='Crossover rate', default=0.8)
        form.addParam('gaWindowSize', IntParam, label='Window size', default=10)
        form.addParam('lsFreq', FloatParam, label='Local search frequency', default=0.06)

        form.addSection(label="Solis & Wets algorithm")
        form.addParam('swMaxIts', IntParam, label='Max. Number of iterations', default=300)
        form.addParam('swMaxSucc', IntParam, label='Successes in a raw', default=4,
                      help='Number of successes before changing rho')
        form.addParam('swMaxFail', IntParam, label='Failures in a raw', default=4,
                      help='Number of failures before changing rho')
        form.addParam('swRho', FloatParam, label='Initial variance', default=1.0,
                      help='It defines the size of the local search')
        form.addParam('swLbRho', FloatParam, label='Variance lower bound', default=0.01)
        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        cId = self._insertFunctionStep('convertStep', prerequisites=[])

        dockSteps = []
        if self.wholeProt:
            gridId = self._insertFunctionStep('generateGridsStep', prerequisites=[cId])
            for mol in self.inputLibrary.get():
                dockId = self._insertFunctionStep('dockStep', mol.clone(), prerequisites=[gridId])
                dockSteps.append(dockId)
        else:
            for pocket in self.inputPockets.get():
                gridId = self._insertFunctionStep('generateGridsStep', pocket.clone(), prerequisites=[cId])
                for mol in self.inputLibrary.get():
                    dockId = self._insertFunctionStep('dockStep', mol.clone(), pocket.clone(), prerequisites=[gridId])
                    dockSteps.append(dockId)

        self._insertFunctionStep('createOutputStep', prerequisites=dockSteps)

    def convertStep(self):
        self.ligandFileNames = []
        for mol in self.inputLibrary.get():
            fnSmall, smallDir = self.convert2PDBQT(mol.clone(), self._getExtraPath())
            self.ligandFileNames.append(fnSmall)

        self.receptorFile = self.convertReceptor2PDBQT(self.getOriginalReceptorFile())

    def generateGridsStep(self, pocket=None):
        self.gridDirs = []
        fnReceptor = self.receptorFile
        if self.wholeProt:
            radius = self.radius.get()
            #Use the original pdb for mass center
            structure, x_center, y_center, z_center = calculate_centerMass(self.getOriginalReceptorFile())
        else:
            radius = (pocket.getDiameter() / 2) * self.pocketRadiusN.get()
            x_center, y_center, z_center = pocket.calculateMassCenter()

        outDir = self.getOutputGridDir(pocket)
        makePath(outDir)

        npts = (radius * 2) / self.spacing.get()
        gpf_file = generate_gpf(fnReceptor, spacing=self.spacing.get(),
                                xc=x_center, yc=y_center, zc=z_center,
                                npts=npts, outDir=outDir, ligandFns=self.ligandFileNames)

        args = "-p {} -l {}.glg".format(gpf_file, self.getReceptorName())
        self.runJob(autodock_plugin.getAutodockPath("autogrid4"), args, cwd=outDir)


    def dockStep(self, mol, pocket=None):
        #Prepare grid
        fnReceptor = self.receptorFile
        nameReceptor, extReceptor = os.path.splitext(os.path.basename(fnReceptor))
        outDir = self.getOutputGridDir(pocket)

        molFn = self.getMolLigandName(mol)
        molName, _ = os.path.splitext(os.path.basename(molFn))
        createLink(molFn, os.path.join(outDir, molName+'.pdbqt'))
        smallDir = os.path.join(outDir, molName)
        makePath(smallDir)

        fnDPF = os.path.abspath(os.path.join(smallDir, nameReceptor+".dpf"))
        args = " -l %s -r %s -o %s"%(molFn, fnReceptor, fnDPF)
        args += " -p ga_pop_size=%d"%self.gaPop.get()
        args += " -p ga_num_evals=%d"%self.gaNumEvals.get()
        args += " -p ga_num_generations=%d"%self.gaNumGens.get()
        args += " -p ga_elitism=%d"%self.gaElitism.get()
        args += " -p ga_mutation_rate=%f"%self.gaMutationRate.get()
        args += " -p ga_crossover_rate=%f"%self.gaCrossOverRate.get()
        args += " -p ga_window_size=%d"%self.gaWindowSize.get()
        args += " -p sw_max_its=%d"%self.swMaxIts.get()
        args += " -p sw_max_succ=%d"%self.swMaxSucc.get()
        args += " -p sw_max_fail=%d"%self.swMaxFail.get()
        args += " -p sw_rho=%f"%self.swRho.get()
        args += " -p sw_lb_rho=%d"%self.swLbRho.get()
        args += " -p ls_search_freq=%f"%self.lsFreq.get()
        args += " -p ga_run=%d"%self.gaRun.get()
        args += " -p rmstol=%f"%self.rmsTol.get()

        self.runJob(autodock_plugin.getMGLPath('bin/pythonsh'),
                    autodock_plugin.getADTPath('Utilities24/prepare_dpf42.py')+args,
                    cwd=outDir)

        fnDLG = fnDPF.replace('.dpf', '.dlg')
        args = "-p %s -l %s" % (fnDPF, fnDLG)
        self.runJob(autodock_plugin.getAutodockPath("autodock4"), args, cwd=outDir)

    def createOutputStep(self):
        if self.checkSingleOutput():
            outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())

        for i, gridDir in enumerate(self.getGridDirs()):
            if not self.checkSingleOutput():
                outputSet = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix=i+1)

            for smallMol in self.inputLibrary.get():
                smallFn = smallMol.getFileName()
                smallName = os.path.splitext(os.path.basename(smallFn))[0]
                fnSmallDir = os.path.join(gridDir, smallName)

                fnDlg = os.path.join(fnSmallDir, self.getReceptorName()+".dlg")
                if os.path.exists(fnDlg):
                    molDic = self.parseDockedMolsDLG(fnDlg)
                    for posId in molDic:
                        pdbFile = '{}/{}_{}.pdb'.format(fnSmallDir, smallName, posId)
                        with open(pdbFile, 'w') as f:
                            f.write(molDic[posId]['pdb'])

                        newSmallMol = SmallMolecule()
                        newSmallMol.copy(smallMol)
                        newSmallMol.cleanObjId()
                        newSmallMol.dockingScoreLE = pwobj.Float(molDic[posId]['energy'])
                        newSmallMol.ligandEfficiency = pwobj.Float(molDic[posId]['ki'])
                        newSmallMol.poseFile.set(pdbFile)
                        newSmallMol.gridId.set(i+1)

                        outputSet.append(newSmallMol)

            if not self.checkSingleOutput():
                self._defineOutputs(**{'outputSmallMolecules_{}'.format(i+1): outputSet})
                self._defineSourceRelation(self.inputLibrary, outputSet)

        if self.checkSingleOutput():
            self._defineOutputs(outputSmallMolecules = outputSet)
            self._defineSourceRelation(self.inputLibrary, outputSet)
      
########################### Utils functions ############################

    def convert2PDBQT(self, smallMol, oDir):
        '''Convert ligand to pdbqt using obabel'''
        inFile = smallMol.getFileName()
        inName, inExt = os.path.splitext(os.path.basename(inFile))
        oFile = os.path.abspath(os.path.join(oDir, inName+'.pdbqt'))

        if inExt != '.pdbqt':
            args = ' -l {} -o {}'.format(inFile, oFile)
            self.runJob(autodock_plugin.getMGLPath('bin/pythonsh'),
                        autodock_plugin.getADTPath('Utilities24/prepare_ligand4.py') + args)
        else:
            createLink(inFile, oFile)
        return oFile, oDir

    def convertReceptor2PDBQT(self, proteinFile):
        inName, inExt = os.path.splitext(os.path.basename(proteinFile))
        oFile = os.path.abspath(os.path.join(self._getExtraPath(inName + '.pdbqt')))

        args = ' -v -r %s -o %s' % (proteinFile, oFile)
        self.runJob(autodock_plugin.getMGLPath('bin/pythonsh'),
                    autodock_plugin.getADTPath('Utilities24/prepare_receptor4.py') + args)

        return oFile

    def getReceptorName(self):
        if self.wholeProt:
            atomStructFn = self.inputAtomStruct.get().getFileName()
            return atomStructFn.split('/')[-1].split('.')[0]
        else:
            return self.inputPockets.get().getProteinName()

    def getOriginalReceptorFile(self):
        if self.wholeProt:
            return self.inputAtomStruct.get().getFileName()
        else:
            return self.inputPockets.get().getProteinFile()

    def getReceptorDir(self):
        if self.wholeProt:
            atomStructFn = self.inputAtomStruct.get().getFileName()
        else:
            atomStructFn = self.inputPockets.get().getProteinFile()
        return '/'.join(atomStructFn.split('/')[:-1])

    def getLigandsFileNames(self):
        ligFns = []
        for mol in self.inputLibrary.get():
            ligFns.append(mol.getFileName())
        return ligFns

    def getLigandsNames(self):
        ligNames = []
        for mol in self.inputLibrary.get():
            ligNames.append(mol.getMolName())
        return ligNames

    def getOutputGridDir(self, pocket=None):
        if pocket==None:
            outDir = self._getExtraPath('grid')
        else:
            outDir = self._getExtraPath('grid_{}'.format(pocket.getObjId()))
        return outDir

    def getGridDirs(self):
        dirs = []
        for file in os.listdir(self._getExtraPath()):
            d = self._getExtraPath(file)
            if os.path.isdir(d) and 'grid' in file:
                dirs.append(d)
        dirs.sort()
        return dirs

    def getMolLigandName(self, mol):
        molName = mol.getMolName()
        for ligName in self.ligandFileNames:
            if molName in ligName:
                return ligName

    def parseDockedMolsDLG(self, fnDlg):
        molDic = {}
        i=1
        with open(fnDlg) as fRes:
            for line in fRes:
                if line.startswith('DOCKED: MODEL'):
                    posId = line.split()[-1]
                    molDic[posId] = {'pdb': ''}
                elif line.startswith('DOCKED: USER    Estimated Free Energy'):
                    molDic[posId]['energy'] = line.split()[8]
                elif line.startswith('DOCKED: USER    Estimated Inhibition'):
                    molDic[posId]['ki'] = line.split()[7]
                elif line.startswith('DOCKED: REMARK') or line.startswith('TER'):
                    molDic[posId]['pdb'] += line[8:]
                elif line.startswith('DOCKED: ATOM'):
                    molDic[posId]['pdb'] += line[8:]
                    i+=1
        return molDic

    def checkSingleOutput(self):
        return self.mergeOutput.get() or len(self.getGridDirs()) == 1

    def _citations(self):
        return ['Morris2009']
