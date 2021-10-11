# **************************************************************************
# *
# * Authors:  Carlos Oscar Sorzano (coss@cnb.csic.es)
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
from pyworkflow.viewer import DESKTOP_TKINTER, ProtocolViewer
from pyworkflow.protocol.params import IntParam, EnumParam, BooleanParam
from pwchem.viewers import PyMolViewer

from autodock.protocols.protocol_autodock import ProtChemAutodock

class ProtAutodockDockingViewer(ProtocolViewer):
    """ Visualize the output of protocol autodock """
    _label = 'viewer autodock docking'
    _targets = [ProtChemAutodock]
    _environments = [DESKTOP_TKINTER]

    def __init__(self, **args):
        ProtocolViewer.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Docking visualization')
        group = form.addGroup('Pocket ligands')
        group.addParam('displayPymol', EnumParam,
                      choices=self.getChoices(pymol=True), default=0,
                      label='Display docking on pocket result: ',
                      help='Docking results are grouped by their pocket, choose the one to visualize')

        group = form.addGroup('Each ligand')
        group.addParam('displayPymolEach', EnumParam,
                       choices=self.getChoicesEach(), default=0,
                       label='Display one ligand type: ',
                       help='Display each ligand position with the target')

        group = form.addGroup('Single Ligand')
        group.addParam('displayPymolSingle', EnumParam,
                      choices=self.getChoicesSingle(), default=0,
                      label='Display single ligand: ',
                      help='Display this single ligand with the target')


    def getChoicesSingle(self):
        self.outputLigands = {}
        for oAttr in self.protocol.iterOutputAttributes():
            if 'outputSmallMolecules' in oAttr[0]:
                molSet = getattr(self.protocol, oAttr[0])
                for mol in molSet:
                    curMol = mol.clone()
                    pName = curMol.getUniqueName()
                    self.outputLigands[pName] = curMol

        outputLabels = list(self.outputLigands.keys())
        outputLabels.sort()
        return outputLabels

    def getChoices(self, pymol=True):
        outputLabels = []
        for oAttr in self.protocol.iterOutputAttributes():
            outputLabels.append(oAttr[0])

        outputLabels.sort()
        if pymol and len(outputLabels) > 1:
            outputLabels = ['All'] + outputLabels
        return outputLabels

    def getChoicesEach(self):
        self.outputLigandsC = {}
        for oAttr in self.protocol.iterOutputAttributes():
            if 'outputSmallMolecules' in oAttr[0]:
                molSet = getattr(self.protocol, oAttr[0])
                for mol in molSet:
                    curMol = mol.clone()
                    pName = curMol.getMolBase()
                    if not pName in self.outputLigandsC:
                        self.outputLigandsC[pName] = [curMol]
                    else:
                        self.outputLigandsC[pName] += [curMol]

        outputLabels = list(self.outputLigandsC.keys())
        outputLabels.sort()
        return outputLabels


    def _getVisualizeDict(self):
        return {'displayPymol': self._viewPymol,
                'displayPymolSingle': self._viewSinglePymol,
                'displayPymolEach': self._viewEachPymol,
                }

    def _viewSinglePymol(self, e=None):
        ligandLabel = self.getEnumText('displayPymolSingle')
        mol = self.outputLigands[ligandLabel].clone()
        pmlFile = self.protocol._getExtraPath(mol.getMolName())+'/{}.pml'.format(ligandLabel)
        self.writePmlFile(pmlFile, self.buildPMLDockingSingleStr(mol, ligandLabel, addTarget=True))

        pymolV = PyMolViewer(project=self.getProject())
        pymolV.visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

    def _viewEachPymol(self, e=None):
        pmlStr = ''
        ligandLabel = self.getEnumText('displayPymolEach')
        pmlFile = self.protocol._getExtraPath('{}.pml'.format(ligandLabel))

        mols = self.outputLigandsC[ligandLabel]
        mols = self.sortMolsByUnique(mols)
        for mol in mols:
            pmlStr += self.buildPMLDockingSingleStr(mol.clone(), mol.getUniqueName(), addTarget=True)

        self.writePmlFile(pmlFile, pmlStr)

        pymolV = PyMolViewer(project=self.getProject())
        pymolV.visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

    def _viewPymol(self, e=None):
        if self.getEnumText('displayPymol') != 'All':
            outName = self.getEnumText('displayPymol')
            pmlFile = self.protocol._getPath('{}.pml'.format(outName))
            pmlStr = self.buildPMLDockingsStr(outName)
        else:
            pmlFile = self.protocol._getPath('allOutputMols.pml')
            outName = self.getChoices()[1]
            pmlStr = self.buildPMLDockingsStr(outName)
            for outName in self.getChoices()[2:]:
                pmlStr += self.buildPMLDockingsStr(outName, addTarget=False)

        self.writePmlFile(pmlFile, pmlStr)
        pymolV = PyMolViewer(project=self.getProject())
        pymolV.visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))


    def buildPMLDockingSingleStr(self, mol, molName, addTarget=True):
        pmlStr = ''
        if addTarget:
            pmlStr = 'load {}\n'.format(os.path.abspath(self.protocol.getOriginalReceptorFile()))

        pdbFile = os.path.abspath(mol.getPoseFile())
        pmlStr += 'load {}, {}\ndisable {}\nhide spheres, {}\nshow sticks, {}\n'.\
            format(pdbFile, molName, molName, molName, molName)
        return pmlStr

    def buildPMLDockingsStr(self, outName, addTarget=True):
        molSet = getattr(self.protocol, outName)
        pmlStr = ''
        if addTarget:
            pmlStr = 'load {}\n'.format(os.path.abspath(self.protocol.getOriginalReceptorFile()))
        if not os.path.exists(self.getPDBsDir(outName)):
            molDic = {}
            for i, mol in enumerate(molSet):
                pFile = os.path.abspath(mol.getPoseFile())
                pName = mol.getUniqueName()
                molDic[pName] = pFile

            molKeys = list(molDic.keys())
            molKeys.sort()
            for molName in molKeys:
                pFile = molDic[molName]
                pmlStr += 'load {}, {}\ndisable {}\nhide spheres, {}\nshow sticks, {}\n'.format(
                    pFile, molName, molName, molName, molName).format(pFile, molName, molName)

        else:
            molFiles = list(os.listdir(self.getPDBsDir(outName)))
            molFiles.sort()
            for pdbFile in molFiles:
                pFile = os.path.abspath(os.path.join(self.getPDBsDir(outName), pdbFile))
                pName = pdbFile.split('.')[0]
                pmlStr += 'load {}, {}\ndisable {}\n'.format(pFile, pName, pName)
        return pmlStr

    def writePmlFile(self, pmlFile, pmlStr):
        with open(pmlFile, 'w') as f:
            f.write(pmlStr)
            f.write('zoom')
        return pmlFile, pmlStr

    def getPDBsDir(self, outName):
        return self.protocol._getExtraPath(outName)

    def sortMolsByUnique(self, mols):
        uniques = []
        for mol in mols:
            uniques.append(mol.getUniqueName())

        zipped_lists = sorted(zip(uniques, mols))
        uniques, mols = zip(*zipped_lists)
        return mols

