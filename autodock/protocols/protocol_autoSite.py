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

"""
This protocol is used to perform a pocket search on a autodock grid from
protein structure using the AutoSite software

"""

import os, shutil

from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam, FloatParam, LEVEL_ADVANCED
from pyworkflow.protocol import params

from pwchem.objects import SetOfStructROIs, StructROI
from pwchem.utils import insistentRun

from autodock import Plugin as autodock_plugin
from autodock.protocols.protocol_autodock import ProtChemAutodockBase



class ProtChemAutoSite(ProtChemAutodockBase):
    """Perform a pocket find experiment with AutoSite. See the help at
       http://autodock.scripps.edu/faqs-help/manual/autodock-4-2-user-guide/AutoDock4.2_UserGuide.pdf"""
    _label = 'AutoSite'

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Input')
        group.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct",
                      label='Input atomic structure:', allowsNull=False,
                      help="The atom structure to search pockets in")

        group = form.addGroup('Parameters')
        group.addParam('spacing', FloatParam, default=1, label='Step size (A): ',
                      help="Distance between each point in the electrostatic grid. This value is used to adjust the "
                           "radius as number of (x,y,z) points : radius/spacing = number of points along 3 dimensions ")

        group.addParam('nneighbors', IntParam, default=14, label='Number of neighbors: ',
                      help="Minimum number of neighbor grid points for a cluster to be considered a pocket")

        line = group.addLine('Energy cutoffs: ', expertLevel=LEVEL_ADVANCED,
                             help='Energy cutoffs for the C, O, and H atom types points (kcal/mol) to be considered of '
                                  'high affinity')
        line.addParam('cutC', params.FloatParam, label='Cutoff C: ', default=-0.3)
        line.addParam('cutO', params.FloatParam, label='Cutoff O: ', default=-0.66)
        line.addParam('cutH', params.FloatParam, label='Cutoff H: ', default=-0.5)

        group.addParam('pepScore', BooleanParam, default=False,
                       label='Score for peptides: ', expertLevel=LEVEL_ADVANCED,
                       help="Whether to tune the scoring function to be directed to peptides")


    # --------------------------- Steps functions --------------------
    def _insertAllSteps(self):
        iniDep = [self._insertFunctionStep('convertInputStep')]
        clustDeps = [self._insertFunctionStep('predictPocketStep', prerequisites=iniDep)]

        self._insertFunctionStep('createOutputStep', prerequisites=clustDeps)

    def convertInputStep(self):
        '''Moves necessary files to current extra path'''
        oriReceptorFile = self.getOriginalReceptorFile()
        if os.path.splitext(oriReceptorFile)[1] != '.pdbqt':
            self.convertReceptor2PDBQT(oriReceptorFile)
        else:
            shutil.copy(oriReceptorFile, self.getReceptorPDBQT())
        self.cleanPDBQT(self.getReceptorPDBQT())

    def predictPocketStep(self):
        fnReceptor = self.getReceptorPDBQT()
        argsSite = self.getAutoSiteArgs(fnReceptor)
        insistentRun(self, autodock_plugin.getPackagePath(package='AUTOSITE', path='bin/autosite'), argsSite,
                     nMax=5, cwd=self._getExtraPath(), sleepTime=5)

    def createOutputStep(self):
        outFiles = self.getOutFiles()
        scDic = self.getScoresDic()

        outPockets = SetOfStructROIs(filename=self._getPath('pockets.sqlite'))
        for oFile in outFiles:
            oId = self.getIdFromFile(oFile)
            featFile = oFile.replace('_cl_', '_fp_')
            pock = StructROI(oFile, self.getOriginalReceptorFile(), pClass='AutoSite', objId=oId, extraFile=featFile,
                             score=scDic[oId]['score'], energy=scDic[oId]['energy'], nPoints=scDic[oId]['nPoints'],
                             volume=scDic[oId]['volume'])
            outPockets.append(pock)

        if len(outPockets) > 0:
            outHETMFile = outPockets.buildPDBhetatmFile()
            self._defineOutputs(outputStructROIs=outPockets)


    # --------------------------- Utils functions --------------------
    def getOriginalReceptorFile(self):
        return self.inputAtomStruct.get().getFileName()


    def getAutoSiteArgs(self, fnReceptor):
        args = ' -r {} --spacing {} --nneighbors {} -bc user {} {} {}'.\
            format(fnReceptor, self.spacing.get(), self.nneighbors.get(), self.cutC.get(), self.cutO.get(),
                   self.cutH.get())
        if self.pepScore:
            args += ' --pep'
        return args


    def getIdFromFile(self, file):
        return int(os.path.basename(file).split('_')[-1].split('.')[0])

    def getOutFiles(self, key='_cl_'):
        outFiles, pdbName = [], self.getReceptorName()
        allFiles = os.listdir(self._getExtraPath(pdbName))
        for file in allFiles:
            if key in file:
                outFiles.append(os.path.abspath(self._getExtraPath(pdbName, file)))
        return outFiles

    def getScoresDic(self):
        scDic, resultsFile = {}, self._getExtraPath(self.getReceptorName() + '_AutoSiteSummary.log')
        inScores = False
        with open(resultsFile) as f:
            for line in f:
                if not inScores and line.startswith('----'):
                    inScores = True
                elif line.startswith('AutoSite identified'):
                    inScores = False

                elif inScores:
                    sline = line.split()
                    scDic[int(sline[0])] = {'energy': sline[1], 'nPoints': sline[2], 'volume': sline[2],
                                            'score': sline[-4]}
        return scDic

    def getTmpSizePath(self, pocketSize, file=''):
        return self._getTmpPath('Size_{}/{}'.format(pocketSize, file))

    # --------------------------- INFO functions -----------------------------------

    def _citations(self):
        return []

    def _warnings(self):
        warns = []
        if self.nneighbors.get() < 14:
            warns.append('If the number of neighbors is too low (<14), AutoSite might fail')
        return warns
