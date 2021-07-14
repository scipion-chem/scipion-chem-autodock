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

import os

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, IntParam, EnumParam
from autodock import Plugin as autodock_plugin
from autodock.objects import autoLigandPocket
from pwchem.objects import SetOfPockets
from pyworkflow.utils.path import copyTree

class ProtChemAutoLigand(EMProtocol):
    """Perform a pocket find experiment with autoLigand. See the help at
       http://autodock.scripps.edu/faqs-help/manual/autodock-4-2-user-guide/AutoDock4.2_UserGuide.pdf"""
    _label = 'autoLigand'
    fillChoices = ['Number', 'Range']

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputGrid', PointerParam, pointerClass="GridAGD",
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

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('predictPocketStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
        '''Moves necessary files to current extra path'''
        gridPath = self.inputGrid.get().getFileName()
        prevExtraPath = '/'.join(gridPath.split('/')[:-1]) + '/extra/'
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

        outPockets = SetOfPockets(filename=self._getExtraPath('pockets.sqlite'))
        for oFile in outFiles:
            pock = autoLigandPocket(oFile, resultsFile)
            outPockets.append(pock)
        self._defineOutputs(outputPockets=outPockets)

    def _citations(self):
        return ['Morris2009']


    def getPDBName(self):
        extraFiles = os.listdir(self._getExtraPath())
        for file in extraFiles:
            if '.pdbqt' in file:
                return file.split('/')[-1].split('.')[0]