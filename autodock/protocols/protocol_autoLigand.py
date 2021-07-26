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

import os

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, IntParam, EnumParam, FloatParam
from autodock import Plugin as autodock_plugin
from autodock.objects import AutoLigandPocket
from pwchem.objects import SetOfPockets
from pyworkflow.utils.path import copyTree, cleanPath
from pyworkflow.protocol import params
from ..constants import PML_STR

NUMBER, RANGE = 0, 1

class ProtChemAutoLigand(EMProtocol):
    """Perform a pocket find experiment with autoLigand. See the help at
       http://autodock.scripps.edu/faqs-help/manual/autodock-4-2-user-guide/AutoDock4.2_UserGuide.pdf"""
    _label = 'autoLigand'
    fillChoices = ['Number', 'Range']

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.stepsExecutionMode = params.STEPS_PARALLEL
        self.finalPockets = []

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputGrid', PointerParam, pointerClass="GridADT",
                       label='Input grid:', allowsNull=False,
                       help="The grid must be prepared for autodock")

        form.addParam('fillType', EnumParam, default=NUMBER,
                       label='Parameters to write',
                       choices=self.fillChoices,
                       help='List of parameters which will be written from the selected protocol')

        form.addParam('nFillPoints', IntParam, default=100, label='Number of fill points',
                       condition="fillType==0",
                       help="Number of fill points to use. The resulting grids will have this size")

        form.addParam('iniFillPoints', IntParam, default=50, label='Initial fill points',
                       condition="fillType==1",
                       help="Number of fill points to use. The resulting grids will have this size")
        form.addParam('endFillPoints', IntParam, default=100, label='Final fill points',
                       condition="fillType==1",
                       help="Number of fill points to use. The resulting grids will have this size")
        form.addParam('stepFillPoints', IntParam, default=10, label='Step of fill points',
                       condition="fillType==1",
                       help="Number of fill points to use. The resulting grids will have this size")

        form.addParam('propShared', FloatParam, default=0.5, label='Proportion of points for overlapping',
                      condition="fillType==1",
                      help="Min proportion of points (from the smaller) of two pockets to be considered overlapping")

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- Steps functions --------------------
    def _insertAllSteps(self):
        iniDep = [self._insertFunctionStep('convertInputStep')]
        predDeps, clustDeps = [], iniDep.copy()
        if self.fillType.get() == NUMBER:
            predDeps += [self._insertFunctionStep('predictPocketStep', self.nFillPoints.get(), prerequisites=iniDep)]
        else:
            for pocketSize in range(self.iniFillPoints.get(), self.endFillPoints.get()+1, self.stepFillPoints.get()):
                predDeps += [self._insertFunctionStep('predictPocketStep', pocketSize, prerequisites=iniDep)]
                #Cluster pockets when each is predicted and when the previous clustering is done (iniDep for the first)
                clustDeps += [self._insertFunctionStep('clusterPocketsStep', pocketSize,
                                                       prerequisites=[predDeps[-1], clustDeps[-1]])]

        self._insertFunctionStep('createOutputStep', prerequisites=clustDeps)

    def convertInputStep(self):
        '''Moves necessary files to current extra path'''
        gridPath = self.inputGrid.get().getFileName()
        prevExtraPath = '/'.join(gridPath.split('/')[:-1])
        copyTree(prevExtraPath, self._getExtraPath())

    def predictPocketStep(self, pocketSize):
        pdbName = self.getPDBName()
        program = "AutoLigand"
        args = ' -r {} -p {}'.format(pdbName, pocketSize)

        copyTree(self._getExtraPath(), self.getTmpSizePath(pocketSize))
        self.runJob(autodock_plugin.getMGLPath('bin/pythonsh'),
                    autodock_plugin.getADTPath('%s.py' % program) + args,
                    cwd=self.getTmpSizePath(pocketSize))

    def clusterPocketsStep(self, pocketSize):
        outFiles, resFile = self.getOutFiles(pocketSize)
        for oFile in outFiles:
            overlaps = []
            curPockets = self.finalPockets.copy()
            newPock = AutoLigandPocket(oFile, resFile)
            for curPock in curPockets:
              propOverlap = self.calculateOverlap(newPock, curPock)
              if propOverlap > self.propShared.get():
                  overlaps.append(curPock)

            #New envelope is accepted if the average energy of the overlaps is higher (worse)
            if len(overlaps) > 0:
              if self.getMinEnergy(overlaps) > newPock.getEnergyPerVol():
                for overPock in overlaps:
                  self.finalPockets.remove(overPock)
                self.finalPockets += [newPock]
            else:
              self.finalPockets += [newPock]


    def createOutputStep(self):
        outFiles, resultsFile = self.organizeOutput()
        self.createPML(outFiles)

        outPockets = SetOfPockets(filename=self._getExtraPath('pockets.sqlite'))
        for oFile in outFiles:
            pock = AutoLigandPocket(oFile, resultsFile)
            outPockets.append(pock)
        self._defineOutputs(outputPockets=outPockets)

    # --------------------------- Utils functions --------------------
    def getIdFromFile(self, file):
        return int(os.path.basename(file).split('out')[1].split('.')[0])

    def getSizeFromFile(self, file, results = False):
        if results:
            return int(os.path.basename(file).split('_')[1].split('Result')[0])
        else:
            return int(os.path.basename(file).split('_')[1].split('out')[0])

    #Do in setofpockets defs
    def getAverageEnergy(self, pockets):
        return sum(self.getPocketEnergies(pockets)) / len(pockets)

    def getMinEnergy(self, pockets):
        return min(self.getPocketEnergies(pockets))

    def getPocketEnergies(self, pockets):
      es = []
      for p in pockets:
        es.append(p.getEnergyPerVol())
      return es

    def getPocketsFiles(self, pockets):
        fs = []
        for p in pockets:
          fs.append(p.getFileName())
        return fs

    def getPocketsSizes(self, pockets):
        sizes = set([])
        for p in pockets:
          sizes.add(int(p.getNumberOfPoints()))
        return sizes


    def calculateOverlap(self, pock1, pock2):
        p1, p2 = pock1.getPoints(), pock2.getPoints()
        minSize = min(pock1.getNumberOfPoints(), pock2.getNumberOfPoints())
        numSharedPoints = len(set(p1).intersection(set(p2)))
        return numSharedPoints / int(minSize)

    def getOutFiles(self, pocketSize):
        outFiles = []
        for file in os.listdir(self.getTmpSizePath(pocketSize)):
          if 'FILL' in file:
            outFiles.append(self.getTmpSizePath(pocketSize, file))
          elif 'Result' in file:
            resultsFile = self.getTmpSizePath(pocketSize, file)
        return outFiles, resultsFile

    def organizeOutput(self):
        '''Scan the Tmp for the output and cleans it'''
        if self.fillType.get() == NUMBER:
            pocketSize = self.nFillPoints.get()
            outFiles, resultsFile = self.getOutFiles(pocketSize)
        else:
            resultsFiles = {}
            outFiles = self.getPocketsFiles(self.finalPockets)
            for pocketSize in self.getPocketsSizes(self.finalPockets):
                _, resFile = self.getOutFiles(pocketSize)
                resultsFiles[pocketSize] = resFile
            outFiles, resultsFile = self.manageIds(outFiles, resultsFiles)

        return outFiles, resultsFile

    def manageIds(self, outFiles, resultsFiles):
        newResFile = self._getPath('{}_Results.txt'.format(self.getPDBName()))
        newOutFiles, i = [], 1
        with open(newResFile, 'w') as fOut:
          for oFile in outFiles:
            resFile = resultsFiles[self.getSizeFromFile(oFile)]
            oldId = self.getIdFromFile(oFile)
            rline = self.getResultsLine(resFile, oldId)

            fOut.write(rline.replace('Output #  {},'.format(oldId), 'Output #  {},'.format(i)))
            newFile = self._getPath(os.path.basename(oFile.replace('out{}.pdb'.format(oldId), 'out{}.pdb'.format(i))))
            os.rename(oFile, newFile)
            newOutFiles.append(newFile)
            i+=1
        return newOutFiles, newResFile


    def getResultsLine(self, resFile, lineId):
      i=1
      with open(resFile) as f:
        for line in f:
          if i==lineId:
            return line
          i+=1

    def moveDic2Path(self, dicFiles):
      newDic = {}
      for size in dicFiles:
        oFile = dicFiles[size]
        newDic[size] = self._getPath(os.path.basename(oFile))
        os.rename(oFile, newDic[size])
      return newDic

    def move2Path(self, files):
        newFiles = []
        for oFile in files:
          newFiles.append(self._getPath(os.path.basename(oFile)))
          os.rename(oFile, newFiles[-1])
        return newFiles

    def getTmpSizePath(self, pocketSize, file=''):
        return self._getTmpPath('Size_{}/{}'.format(pocketSize, file))

    def getPDBName(self):
        extraFiles = os.listdir(self._getExtraPath())
        for file in extraFiles:
            if '.gpf' in file:
                return file.split('/')[-1].split('.')[0]

    def getPDBFileName(self):
        extraFiles = os.listdir(self._getExtraPath())
        for file in extraFiles:
          if '.pdbqt' in file:
            return self._getExtraPath(file)

    def createPML(self, pocketFiles):
      pdbName = self.getPDBName()
      pdbFile = self._getExtraPath('{}.pdbqt'.format(pdbName))
      outFile = self._getExtraPath('{}_out.pdbqt'.format(pdbName))
      pmlFile = self._getExtraPath('{}.pml'.format(pdbName))
      with open(pdbFile) as fpdb:
        outStr = '\n'.join(fpdb.read().split('\n')[:-2])

      # Creates a pdb(qt) with the HETATOM corresponding to pocket points
      for file in pocketFiles:
        outStr += self.formatPocketStr(file)
      with open(outFile, 'w') as f:
        f.write(outStr)

      # Creates the pml for pymol visualization
      with open(pmlFile, 'w') as f:
        f.write(PML_STR.format(outFile.split('/')[-1]))

    def formatPocketStr(self, pocketFile):
      numId = pocketFile.split('out')[1].split('.')[0]
      outStr = ''
      with open(pocketFile) as f:
        for line in f:
          line = line.split()
          replacements = ['HETATM', line[1], 'APOL', 'STP', 'C', numId, *line[5:-1], 'Ve']
          pdbLine = self.writePDBLine(replacements)
          outStr += pdbLine
      return outStr

    def writePDBLine(self, j):
      '''j: elements to write in the pdb'''
      j[0] = j[0].ljust(6)  # atom#6s
      j[1] = j[1].rjust(5)  # aomnum#5d
      j[2] = j[2].center(4)  # atomname$#4s
      j[3] = j[3].ljust(3)  # resname#1s
      j[4] = j[4].rjust(1)  # Astring
      j[5] = j[5].rjust(4)  # resnum
      j[6] = str('%8.3f' % (float(j[6]))).rjust(8)  # x
      j[7] = str('%8.3f' % (float(j[7]))).rjust(8)  # y
      j[8] = str('%8.3f' % (float(j[8]))).rjust(8)  # z\
      j[9] = str('%6.2f' % (float(j[9]))).rjust(6)  # occ
      j[10] = str('%6.2f' % (float(j[10]))).ljust(6)  # temp
      j[11] = str('%8.3f' % (float(j[11]))).rjust(10)
      j[12] = j[12].rjust(3)  # elname
      return "\n%s%s %s %s %s%s    %s%s%s%s%s%s%s" % \
             (j[0], j[1], j[2], j[3], j[4], j[5], j[6], j[7], j[8], j[9], j[10], j[11], j[12])

    # --------------------------- INFO functions -----------------------------------

    def _citations(self):
        return []

    def _warnings(self):
      """ Try to find warnings on define params. """
      import re
      warnings = []
      inpFile = os.path.abspath(self.getPDBFileName())
      with open(inpFile) as f:
        fileStr = f.read()
      if re.search('\nHETATM', fileStr):
        warnings.append('The structure you are inputing has some *heteroatoms* (ligands).\n'
                        'This will affect the results as its volume is also taken as target.')

      return warnings