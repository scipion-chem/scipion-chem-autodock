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

import os, glob

from pyworkflow.protocol import params
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import getBaseFileName, performBatchThreading, replaceInFiles
from pwchem import Plugin as pwchemPlugin
from pwchem.constants import RDKIT_DIC

from autodock import Plugin as autodockPlugin
from autodock.protocols import ProtChemAutodockGPU

encoderOptions = ['GraphResNet', 'ChemProp', 'COATI']
preprocOptions = ['None', 'LigEff', 'PowerOf']


def average(l):
  return sum(l) / len(l)

class ProtEncoderDockScoring(ProtChemAutodockGPU):
  """Trains a prediction of the docking score over a set of small molecules using a encoder-regressor model.
  This model must be trained for each receptor separately with docking scores.
  """
  _label = 'Encoder-regressor dock scoring'
  _program = ""

  def _defineParams(self, form):
    form.addSection(label="Prediction")
    group = form.addGroup('Input')
    group.addParam('inputSmallMolecules', params.PointerParam, pointerClass="SetOfSmallMolecules",
                   label='Input small molecules: ', allowsNull=False,
                   help="Input small molecules to be scored with the model")

    group = form.addGroup('Training')
    group.addParam('loadModel', params.BooleanParam, label='Load pretrained model: ', default=False,
                   help='Whether to load a pretrained model or use some docked molecules to train a new one')
    group.addParam('pretrainedModel', params.PathParam, label='Pretrained model: ', default='', condition='loadModel',
                   help='Choose the pretrained model from in the scipion framework using the wizard or set the '
                        'path to the model directory if found somewhere else')

    group.addParam('dockedMols', params.PointerParam, pointerClass="SetOfSmallMolecules",
                   label='Input docked molecules for training: ', condition='not loadModel',
                   help='Input SetOfSmallMolecules that must be docked to the input receptor and whose scores will be '
                        'used to train the model')
    group.addParam('scoreName', params.StringParam, label='Docking score: ', default='', condition='not loadModel',
                   help='Docking score to use for the model training')
    group.addParam('scoreMerge', params.EnumParam, label='Merge strategy: ', choices=['Min', 'Max', 'Mean'], default=0,
                   expertLevel=params.LEVEL_ADVANCED, condition='not loadModel',
                   help='How to merge the scores if several values are found for the same molecule (because of '
                        'conformers or poses')

    group.addParam('batch', params.IntParam, label='Batch size: ', default=256,
                   expertLevel=params.LEVEL_ADVANCED, condition='not loadModel',
                   help='Batch size to use for training.')
    group.addParam('seed', params.IntParam, label='Random seed: ', default=44,
                   expertLevel=params.LEVEL_ADVANCED, condition='not loadModel',
                   help='Random seed.')

    group = form.addGroup('Encoder')
    group.addParam('encoder', params.EnumParam, label='Encoder model to use: ', choices=encoderOptions, default=1,
                   help='Encoder model to use in the regression.')
    group.addParam('encLayers', params.StringParam, label='Encoder layers: ', default='[2, 2, 2, 1]',
                   condition='encoder==0', expertLevel=params.LEVEL_ADVANCED,
                   help='Number of layers and filter size for the graph ResNet model.')
    group.addParam('filtSizes', params.StringParam, label='Filter sizes: ', default='[64, 128, 256, 256]',
                   condition='encoder==0', expertLevel=params.LEVEL_ADVANCED,
                   help='Filter sizes for each of the graph ResNet layers')
    group.addParam('encPool', params.EnumParam, label='Encoder pooling: ', choices=['Max', 'Mean'], default=0,
                   condition='encoder==0', expertLevel=params.LEVEL_ADVANCED,
                   help='Pooling applied to pass from atomic to molecular features (max, mean)')

    group = form.addGroup('Regression')
    group.addParam('preprocess', params.EnumParam, label='Score preprocessing: ', choices=preprocOptions, default=0,
                   expertLevel=params.LEVEL_ADVANCED,
                   help='Preprocessing to perform over the molecule scores: None, transform to ligand efficiency '
                        '(from energy) or pass an exponencial filter to enlarge high score importances.')
    group.addParam('regLayers', params.StringParam, label='Regression layers: ', default='[256, 256, 128]',
                   help='Number of layers and neurons in each of them to use for the regression model '
                        '(Fully connected neural network). In list of integers format e.g:\n'
                        '[256, 256, 128] defines the default regressor with 3 hidden layers of 256, 256 and 128 '
                        'neurons respectively.')

    form.addParallelSection(threads=4, mpi=1)

  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
      if not self.loadModel.get():
        self._insertFunctionStep('trainingStep')
      self._insertFunctionStep('createOutputStep')

  def trainingStep(self):
    dMols = self.dockedMols.get()
    inDockFile = self.writeInputDocksFile(dMols)

    args = f'{inDockFile} {self.getInputSMIFile()}'
    autodockPlugin.runScript(self, 'convertToSMIs.py', args, envDict=RDKIT_DIC)
    self.mergeSMIscores()

    self.runJob('conda activate dockRegressor && ')



  def getInputDockFile(self):
    return self._getExtraPath('inputDock.csv')

  def getInputSMIFile(self):
    return self._getExtraPath('inputSMIs.csv')

  def writeInputDocksFile(self, mols):
    inDocksFiles = self.getInputDockFile()
    with open(inDocksFiles, 'w') as f:
      for m in mols:
        score = getattr(m, self.scoreName.get())
        f.write(f'{m.getPoseFile()},{score}\n')
    return inDocksFiles

  def getMergeFunction(self):
    mergeDic = {0: min, 1: max, 2: average}
    return mergeDic[self.scoreMerge.get()]

  def mergeSMIscores(self):
    # parsing
    smiDic, smiFile = {}, self.getInputSMIFile()
    with open(smiFile) as f:
      for line in f:
        smi, score = line.strip().split(',')
        if smi not in smiDic:
          smiDic[smi] = []
        smiDic[smi].append(score)

    # merging
    mFunc = self.getMergeFunction()
    for smi, scores in smiDic.items():
      smiDic[smi] = mFunc(scores)

    # rewritting
    with open(smiFile, 'w') as f:
      for smi, score in smiDic.items():
        f.write(f'{smi},{score}\n')

    return smiDic







  def converToSMI(self, mols):
    smiDir = self.getInputSMIDir()
    if not os.path.exists(smiDir):
      os.makedirs(smiDir)

    molDir = self.copyInputMolsInDir(mols)
    args = ' --multiFiles -iD "{}" --pattern "{}" -of smi --outputDir "{}"'. \
      format(molDir, '*', smiDir)
    pwchemPlugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=smiDir)
    return smiDir

  def parseSMIs(self, smiDir):
    smiDic = {}
    for smiFile in os.listdir(smiDir):
      smiPath = os.path.join(smiDir, smiFile)
      with open(smiPath) as f:
        smi = f.read().strip()



  def copyInputMolsInDir(self, mols):
    oDir = os.path.abspath(self._getTmpPath('inMols'))
    if not os.path.exists(oDir):
      os.makedirs(oDir)

    for mol in mols:
      os.link(mol.getFileName(), os.path.join(oDir, os.path.split(mol.getFileName())[-1]))
    return oDir

  def getInputSMIDir(self):
    return os.path.abspath(self._getExtraPath('inputSMI'))



