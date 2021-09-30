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
from pyworkflow.protocol.params import PointerParam, IntParam, FloatParam, STEPS_PARALLEL
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
        form.addParam('input', PointerParam, pointerClass="GridADT",
                       label='Input grid:', allowsNull=False,
                       help="The grid must be prepared for autodock")
        form.addParam('inputLibrary', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Ligand library:', allowsNull=False,
                       help="The library must be prepared for autodock")
        form.addParam('rmsTol', FloatParam, label='Cluster tolerance (A)', default=2.0)
        form.addParam('gaRun', IntParam, label='Number of positions per ligand', default=10)

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
        grid = self.inputGrid.get()
        dockSteps = []
        for smallMol in self.inputLibrary.get():
            stepId = self._insertFunctionStep('dockStep', grid, smallMol.clone(), prerequisites=[])
            dockSteps.append(stepId)
        self._insertFunctionStep('createOutputStep', prerequisites=dockSteps)

    def convert2PDBQT(self, smallMol):
        '''Convert ligand to pdbqt using obabel'''
        inFile = smallMol.getFileName()
        inName, inExt = os.path.splitext(os.path.basename(inFile))
        oDir = self._getExtraPath(inName)
        if not os.path.exists(oDir):
            makePath(oDir)
        oFile = os.path.abspath(os.path.join(oDir, inName+'.pdbqt'))

        if inExt != '.pdbqt':
            args = '-i{} {} -opdbqt -O {}'.format(inExt[1:], inFile, oFile)
            self.runJob('obabel', args)
        else:
            shutil.copy(inFile, oFile)
        return oFile, oDir


    def dockStep(self, grid, smallMol):
        #Prepare grid
        name_protein = (os.path.basename(grid.getProteinFile())).split(".")[0]
        dirPdbqt = '/'.join(grid.getFileName().split('/')[:-1])
        fnReceptor = os.path.abspath('{}/{}.pdbqt'.format(dirPdbqt, name_protein))

        fnSmall, smallDir = self.convert2PDBQT(smallMol)
        smallName, inExt = os.path.splitext(os.path.basename(fnSmall))

        fnDPF = os.path.join(smallDir, smallName+".dpf")
        args = " -l %s -r %s -o %s"%(fnSmall, fnReceptor, fnDPF)
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
                    autodock_plugin.getADTPath('Utilities24/prepare_dpf42.py')+args)

        createLink(fnReceptor, os.path.join(smallDir, "{}.pdbqt".format(name_protein)))

        args = " -r {}.pdbqt -l {} -o library.gpf".format(name_protein, fnSmall)
        self.runJob(autodock_plugin.getMGLPath('bin/pythonsh'),
                    autodock_plugin.getADTPath('Utilities24/prepare_gpf4.py') + args,
                    cwd=smallDir)

        args = "-p library.gpf -l library.glg"
        self.runJob(autodock_plugin.getAutodockPath("autogrid4"), args, cwd=smallDir)

        args = "-p %s.dpf -l %s.dlg"%(smallName, smallName)
        self.runJob(autodock_plugin.getAutodockPath("autodock4"), args, cwd=smallDir)

        # Clean a bit
        cleanPattern(os.path.join(smallDir,"atomStruct.*.map"))

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

    def createOutputStep(self):
        outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())
        for smallMol in self.inputLibrary.get():
            fnSmall = smallMol.getFileName()
            fnBase = os.path.splitext(os.path.split(fnSmall)[1])[0]
            fnSmallDir = self._getExtraPath(fnBase)
            fnDlg = os.path.join(fnSmallDir, fnBase+".dlg")
            if os.path.exists(fnDlg):
                molDic = self.parseDockedMolsDLG(fnDlg)
                for posId in molDic:
                    pdbFile = '{}/{}_{}.pdb'.format(fnSmallDir, fnBase, posId)
                    with open(pdbFile, 'w') as f:
                        f.write(molDic[posId]['pdb'])

                    newSmallMol = SmallMolecule()
                    newSmallMol.copy(smallMol)
                    newSmallMol.cleanObjId()
                    newSmallMol.dockingScoreLE = pwobj.Float(molDic[posId]['energy'])
                    newSmallMol.ligandEfficiency = pwobj.Float(molDic[posId]['ki'])
                    newSmallMol.poseFile.set(pdbFile)

                    outputSet.append(newSmallMol)

        self._defineOutputs(outputSmallMolecules=outputSet)
        self._defineSourceRelation(self.inputGrid, outputSet)
        self._defineSourceRelation(self.inputLibrary, outputSet)

    def _citations(self):
        return ['Morris2009']
