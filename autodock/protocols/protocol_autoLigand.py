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
protein structure using the AutoLigand software

"""

import os, shutil

from pwem.protocols import EMProtocol
from pwem.objects.data import AtomStruct
from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam, EnumParam, FloatParam, LEVEL_ADVANCED
from pyworkflow.utils.path import copyTree, makePath
from pyworkflow.protocol import params

from pwchem.objects import SetOfPockets
from pwchem.constants import *
from pwchem.utils import splitPDBLine, runOpenBabel, generate_gpf, calculate_centerMass
from pwchem import Plugin as pwchem_plugin

from autodock import Plugin as autodock_plugin
from autodock.objects import AutoLigandPocket


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
        group = form.addGroup('Grids generation')
        group.addParam('prevGrid', BooleanParam, label='Use a previously generated grid: ', default=False,
                       help='Use a previous ADT grid to run autoligand')

        group.addParam('inputGrid', PointerParam, pointerClass="GridADT", allowsNull=False,
                       label='Input ADT grid:', condition='prevGrid',
                       help="ADT grid previously generated with the ADT generate grid protocol on the protein"
                            " of interest")

        group.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct",
                      label='Input atomic structure:', condition='not prevGrid', allowsNull=False,
                      help="The atom structure to search pockets in")
        group.addParam('radius', FloatParam, label='Grid radius for whole protein: ',
                       allowsNull=False, condition='not prevGrid',
                       help='Radius of the Autodock grid for the whole protein.'
                            'The wizard will provide for an approximation')
        group.addParam('spacing', FloatParam, default=1, label='Step size (A)',
                      condition='not prevGrid',
                      help="Distance between each point in the electrostatic grid."
                           " This value is used to adjust the radius as number of "
                           "(x,y,z) points : radius/spacing = number of points along"
                           " 3 dimensions ")

        group = form.addGroup('Autoligand')
        group.addParam('fillType', EnumParam, default=NUMBER,
                       label='Mode of use',
                       choices=self.fillChoices,
                       help='List of parameters which will be written from the selected protocol')

        group.addParam('nFillPoints', IntParam, default=100, label='Number of fill points',
                       condition="fillType==0",
                       help="Number of fill points to use. The resulting grids will have this size")

        group.addParam('iniFillPoints', IntParam, default=50, label='Initial fill points',
                       condition="fillType==1",
                       help="Number of fill points to use. The resulting grids will have this size")
        group.addParam('endFillPoints', IntParam, default=100, label='Final fill points',
                       condition="fillType==1",
                       help="Number of fill points to use. The resulting grids will have this size")
        group.addParam('stepFillPoints', IntParam, default=10, label='Step of fill points',
                       condition="fillType==1",
                       help="Number of fill points to use. The resulting grids will have this size")

        group.addParam('propShared', FloatParam, default=0.5, label='Proportion of points for overlapping',
                      condition="fillType==1",
                      help="Min proportion of points (from the smaller) of two pockets to be considered overlapping")

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- Steps functions --------------------
    def _insertAllSteps(self):
        iniDep = [self._insertFunctionStep('convertInputStep')]
        gridDep = [self._insertFunctionStep('generateGridStep', prerequisites=iniDep)]

        predDeps, clustDeps = [], iniDep.copy()
        if self.fillType.get() == NUMBER:
            clustDeps += [self._insertFunctionStep('predictPocketStep', self.nFillPoints.get(), prerequisites=gridDep)]
        else:
            for pocketSize in range(self.iniFillPoints.get(), self.endFillPoints.get()+1, self.stepFillPoints.get()):
                predDeps += [self._insertFunctionStep('predictPocketStep', pocketSize, prerequisites=gridDep)]
                #Cluster pockets when each is predicted and when the previous clustering is done (iniDep for the first)
                clustDeps += [self._insertFunctionStep('clusterPocketsStep', pocketSize,
                                                       prerequisites=[predDeps[-1], clustDeps[-1]])]

        self._insertFunctionStep('createOutputStep', prerequisites=clustDeps)

    def convertInputStep(self):
        '''Moves necessary files to current extra path'''
        self.receptorFile = self.getStructFileName()
        if os.path.splitext(self.receptorFile)[1] != '.pdbqt':
            self.receptorFile = self.convertReceptor2PDBQT(self.receptorFile)

    def generateGridStep(self):
        outDir = self._getExtraPath()
        if self.prevGrid:
            shutil.copytree(self.getStructDir(), outDir, dirs_exist_ok=True)
        else:
            radius = self.radius.get()
            strFile = self.getStructFileName()
            if os.path.splitext(strFile)[1] == '.pdbqt':
                pdbFile = os.path.abspath(self._getTmpPath('pdbInput.pdb'))
                args = ' -ipdbqt {} -opdb -O {}'.format(os.path.abspath(strFile), pdbFile)
                runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())
            else:
                pdbFile = strFile

            structure, x_center, y_center, z_center = calculate_centerMass(pdbFile)
            npts = (radius * 2) / self.spacing.get()

            makePath(outDir)
            gpf_file = generate_gpf(self.receptorFile, spacing=self.spacing.get(),
                                    xc=x_center, yc=y_center, zc=z_center,
                                    npts=npts, outDir=outDir)

            args = "-p {} -l {}.glg".format(gpf_file, self.getStructName())
            self.runJob(autodock_plugin.getAutodockPath("autogrid4"), args, cwd=outDir)

    def predictPocketStep(self, pocketSize):
        pdbName = self.getStructName()
        program = "AutoLigand"
        args = ' -r {} -p {}'.format(pdbName, pocketSize)

        copyTree(self._getExtraPath(), self.getTmpSizePath(pocketSize))
        self.runJob(pwchem_plugin.getMGLPath('bin/pythonsh'),
                    autodock_plugin.getADTPath('%s.py' % program) + args,
                    cwd=self.getTmpSizePath(pocketSize))

    def clusterPocketsStep(self, pocketSize):
        outFiles, resFile = self.getOutFiles(pocketSize)
        for oFile in outFiles:
            overlaps = []
            curPockets = self.finalPockets.copy()
            newPock = AutoLigandPocket(oFile, None, resFile)
            for curPock in curPockets:
              propOverlap = self.calculateOverlap(newPock, curPock)
              if propOverlap > self.propShared.get():
                  overlaps.append(curPock)

            #New envelope is accepted if the average energy of the overlaps is higher (worse)
            if len(overlaps) > 0:
              if self.getMinEnergy(overlaps) > newPock.getEnergyPerVol():
                for overPock in overlaps:
                  self.finalPockets.remove(overPock)
                newPock.setNClusters(len(overlaps)+1)
                self.finalPockets += [newPock]
              else:
                for overPock in overlaps:
                  self.finalPockets[self.finalPockets.index(overPock)].incrNClusters()
            else:
              self.finalPockets += [newPock]


    def createOutputStep(self):
        outFiles, resultsFile = self.organizeOutput()

        outPockets = SetOfPockets(filename=self._getPath('pockets.sqlite'))
        for oFile in outFiles:
            pock = AutoLigandPocket(os.path.abspath(oFile), self.receptorFile, os.path.abspath(resultsFile))
            outPockets.append(pock)

        outHETMFile = outPockets.buildPDBhetatmFile()

        self._defineOutputs(outputPockets=outPockets)


    # --------------------------- Utils functions --------------------
    def getIdFromFile(self, file):
        return int(os.path.basename(file).split('out')[1].split('.')[0])

    def getSizeFromFile(self, file, results = False):
        if results:
            return int(os.path.basename(file).split('_')[1].split('Result')[0])
        else:
            return int(os.path.basename(file).split('_')[1].split('out')[0])

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
            pocketSize, outFiles = self.nFillPoints.get(), []
            oldFiles, oldResultsFile = self.getOutFiles(pocketSize)
            for oFile in oldFiles:
              outFiles.append(self._getPath(os.path.basename(oFile)))
              os.rename(oFile, outFiles[-1])
            
            resultsFile = self._getPath(os.path.basename(oldResultsFile))
            os.rename(oldResultsFile, resultsFile)
        else:
            resultsFiles = {}
            for pocketSize in self.getPocketsSizes(self.finalPockets):
                _, resFile = self.getOutFiles(pocketSize)
                resultsFiles[pocketSize] = resFile
            outFiles, resultsFile = self.manageIds(self.finalPockets, resultsFiles)

        return outFiles, resultsFile

    def manageIds(self, finalPockets, resultsFiles):
        newResFile = self._getPath('{}_Results.txt'.format(self.getStructName()))
        newOutFiles, i = [], 1
        with open(newResFile, 'w') as fOut:
          for pocket in finalPockets:
            oFile = pocket.getFileName()
            resFile = resultsFiles[self.getSizeFromFile(oFile)]
            oldId = self.getIdFromFile(oFile)

            rline = self.getResultsLine(resFile, oldId)
            rline = rline.replace('Output #  {},'.format(oldId), 'Output #  {},'.format(i))
            rline = rline.replace('\n', ', NumberOfClusters = {}\n'.format(pocket.getNClusters()))

            fOut.write(rline)
            newFile = self._getPath(os.path.basename(oFile.replace('out{}.pdb'.format(oldId), 'out{}.pdb'.format(i))))
            os.rename(oFile, newFile)
            newOutFiles.append(newFile)
            i+=1
        return newOutFiles, newResFile

    def convertReceptor2PDBQT(self, proteinFile):
      inName, inExt = os.path.splitext(os.path.basename(proteinFile))
      oFile = os.path.abspath(os.path.join(self._getExtraPath(inName + '.pdbqt')))

      args = ' -v -r %s -o %s' % (proteinFile, oFile)
      self.runJob(pwchem_plugin.getMGLPath('bin/pythonsh'),
                  autodock_plugin.getADTPath('Utilities24/prepare_receptor4.py') + args)
      return oFile

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

    def getStructName(self):
        return self.getStructFileName().split('/')[-1].split('.')[0]

    def getStructFileName(self):
        if self.prevGrid:
            return self.inputGrid.get().getProteinFile()
        else:
            return self.inputAtomStruct.get().getFileName()

    def getStructDir(self):
        atomStructFn = self.getStructFileName()
        return '/'.join(atomStructFn.split('/')[:-1])

    def getOutFileName(self):
      pdbName = self.getStructName()
      return self._getExtraPath('{}_out.pdbqt'.format(pdbName))

    # --------------------------- INFO functions -----------------------------------

    def _citations(self):
        return []

    def _warnings(self):
      """ Try to find warnings on define params. """
      import re
      warnings = []
      inpFile = os.path.abspath(self.getStructFileName())
      with open(inpFile) as f:
        fileStr = f.read()
      if re.search('\nHETATM', fileStr):
        warnings.append('The structure you are inputing has some *heteroatoms* (ligands).\n'
                        'This will affect the results as its volume is also taken as target.')

      return warnings