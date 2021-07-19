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
from pyworkflow.protocol.params import PointerParam, IntParam, EnumParam
from autodock import Plugin as autodock_plugin
from autodock.objects import AutoLigandPocket
from pwchem.objects import SetOfPockets
from pyworkflow.utils.path import copyTree
from ..constants import PML_STR

class ProtChemAutoLigand(EMProtocol):
    """Perform a pocket find experiment with autoLigand. See the help at
       http://autodock.scripps.edu/faqs-help/manual/autodock-4-2-user-guide/AutoDock4.2_UserGuide.pdf"""
    _label = 'autoLigand'
    fillChoices = ['Number', 'Range']

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputGrid', PointerParam, pointerClass="GridADT",
                       label='Input grid:', allowsNull=False,
                       help="The grid must be prepared for autodock")

        # form.addParam('fillType', EnumParam, default=0,
        #             label='Parameters to write',
        #            choices=self.fillChoices,
        #           help='List of parameters which will be written from the selected protocol')

        form.addParam('nFillPoints', IntParam, default=100, label='Number of fill points',
                      help="Number of fill points to use. The resulting grids will have this size")

        #form.addParam('iniFillPoints', IntParam, default=50, label='Initial fill points',
         #                     condition="fillType==1",
          #                    help="Number of fill points to use. The resulting grids will have this size")
           #     form.addParam('endFillPoints', IntParam, default=100, label='Final fill points',
            #                  condition="fillType==1",
             #                 help="Number of fill points to use. The resulting grids will have this size")
              #  form.addParam('stepFillPoints', IntParam, default=10, label='Step of fill points',
               #               condition="fillType==1",
                #              help="Number of fill points to use. The resulting grids will have this size")'''

        #form.addParallelSection(threads=1, mpi=1)

    # --------------------------- Steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('predictPocketStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
        '''Moves necessary files to current extra path'''
        gridPath = self.inputGrid.get().getFileName()
        prevExtraPath = '/'.join(gridPath.split('/')[:-1])
        copyTree(prevExtraPath, self._getExtraPath())

    def predictPocketStep(self):
        pdbName = self.getPDBName()

        program = "AutoLigand"
        args = ' -r {} -p {}'.format(pdbName, self.nFillPoints.get())

        self.runJob(autodock_plugin.getMGLPath('bin/pythonsh'),
                    autodock_plugin.getADTPath('%s.py' % program) + args,
                    cwd=self._getExtraPath())

    def createOutputStep(self):
        outFiles = []
        for file in os.listdir(self._getExtraPath()):
          if 'FILL' in file:
            os.rename(self._getExtraPath(file), self._getPath(file))
            outFiles.append(self._getPath(file))
          elif 'Result' in file:
            os.rename(self._getExtraPath(file), self._getPath(file))
            resultsFile = self._getPath(file)
        self.createPML(outFiles)

        outPockets = SetOfPockets(filename=self._getExtraPath('pockets.sqlite'))
        for oFile in outFiles:
            pock = AutoLigandPocket(oFile, resultsFile)
            outPockets.append(pock)
        self._defineOutputs(outputPockets=outPockets)

    # --------------------------- Utils functions --------------------

    def getPDBName(self):
        extraFiles = os.listdir(self._getExtraPath())
        for file in extraFiles:
            if '.gpf' in file:
                return file.split('/')[-1].split('.')[0]

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
      inpFile = os.path.abspath(self.inputAtomStruct.get().getFileName())
      with open(inpFile) as f:
        fileStr = f.read()
      if re.search('\nHETATM', fileStr):
        warnings.append('The structure you are inputing has some *heteroatoms* (ligands).\n'
                        'This will affect the results as its volume is also taken as target.')

      return warnings