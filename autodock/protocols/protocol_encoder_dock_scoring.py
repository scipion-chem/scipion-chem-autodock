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

import os, glob, shutil

from pyworkflow.protocol import params
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import getBaseFileName, performBatchThreading, replaceInFiles
from pwchem import Plugin as pwchemPlugin
from pwchem.constants import RDKIT_DIC

from autodock import Plugin as autodockPlugin
from autodock.protocols import ProtChemAutodockGPU
from autodock.constants import GCR_DIC

encoderOptions = ['ResNet', 'ChemProp', 'COATI']
preprocOptions = ['None', 'LigEff', 'PowerOf']

gonnaTrain = 'not loadModel or (loadModel and doTrain)'

def average(l):
  return sum(l) / len(l)

class ProtEncoderDockScoring(ProtChemAutodockGPU):
  """Trains a prediction of the docking score over a set of small molecules using a encoder-regressor model.
  This model must be trained for each receptor separately with docking scores.
  """
  _label = 'Encoder-regressor dock scoring'
  _program = ""

  def _defineParams(self, form):
    form.addHidden(params.USE_GPU, params.BooleanParam, default=True,
                   label="Use GPU for execution: ",
                   help="This protocol has both CPU and GPU implementation.\
                                             Select the one you want to use.")

    form.addHidden(params.GPU_LIST, params.StringParam, default='0', label="Choose GPU IDs",
                   help="Add a list of GPU devices that can be used")

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

    group.addParam('doTrain', params.BooleanParam, label='Train model: ', default=False, condition='loadModel',
                   help='Whether to load a pretrained model or use some docked molecules to train a new one')
    group.addParam('modelName', params.StringParam, label='Model name: ', default='', condition='not loadModel',
                   help='Set a name of the system model')
    group.addParam('dockedMols', params.PointerParam, pointerClass="SetOfSmallMolecules",
                   label='Input docked molecules for training: ', condition=gonnaTrain,
                   help='Input SetOfSmallMolecules that must be docked to the input receptor and whose scores will be '
                        'used to train the model')
    group.addParam('scoreName', params.StringParam, label='Docking score: ', default='', condition=gonnaTrain,
                   help='Docking score to use for the model training')
    group.addParam('scoreMerge', params.EnumParam, label='Merge strategy: ', choices=['Min', 'Max', 'Mean'], default=0,
                   expertLevel=params.LEVEL_ADVANCED, condition=gonnaTrain,
                   help='How to merge the scores if several values are found for the same molecule (because of '
                        'conformers or poses')

    group.addParam('batch', params.IntParam, label='Batch size: ', default=256,
                   expertLevel=params.LEVEL_ADVANCED, condition=gonnaTrain,
                   help='Batch size to use for training.')
    group.addParam('seed', params.IntParam, label='Random seed: ', default=44,
                   expertLevel=params.LEVEL_ADVANCED, condition=gonnaTrain,
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

  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
      if not self.loadModel.get() or (self.loadModel.get() and self.doTrain.get()):
        self._insertFunctionStep('trainingStep')
      self._insertFunctionStep('predictionStep')
      self._insertFunctionStep('createOutputStep')

  def trainingStep(self):
    encoderName = self.getEnumText('encoder').lower()
    sysName = self.getSystemName()
    confFile = self.writeConfFile()
    scriptName = self.getScriptPath()

    # todo: save model in autodock/models or similar (COATI LEFT)
    dMols = self.dockedMols.get()
    smisFile = os.path.abspath(self.buildSMIsFile(dMols, writeScores=True))
    args = f'--config {confFile} -n {sysName} -e {encoderName} -ef {autodockPlugin.getChemPropFile()} ' \
           f'-p {smisFile} --doTest --testBest 1 --doTrain '

    modelsPath = os.path.abspath(autodockPlugin.getPluginHome('models'))
    pwchemPlugin.runCondaCommand(self, args, GCR_DIC, f'python {scriptName}', cwd=self._getPath())

    shutil.copytree(os.path.abspath(self._getPath(sysName)), os.path.join(modelsPath, sysName), dirs_exist_ok=True)


  def predictionStep(self):
    confFile = self.writeConfFile()
    sysName = self.getSystemName()
    scriptName = self.getScriptPath()

    inMols = self.inputSmallMolecules.get()
    smisFile = os.path.abspath(self.buildSMIsFile(inMols, writeScores=False))
    args = f'--config {confFile} -n {sysName} -d {sysName} -p {smisFile} --doPredict ' \
           f'-ef {autodockPlugin.getChemPropFile()} '

    modelsPath = os.path.abspath(autodockPlugin.getPluginHome('models'))
    shutil.copytree(os.path.join(modelsPath, sysName), os.path.abspath(self._getPath(sysName)), dirs_exist_ok=True)

    pwchemPlugin.runCondaCommand(self, args, GCR_DIC, f'python {scriptName}', cwd=self._getPath())


  def createOutputStep(self):
    scoreDic = self.getScoreDic()
    outputSet = self.inputSmallMolecules.get().createCopy(self._getPath(), copyInfo=True)
    for mol in self.inputSmallMolecules.get():
      nMol = mol.clone()
      molFile = nMol.getFileName()
      setattr(nMol, '_gcrScore', params.Float(scoreDic[molFile]))
      outputSet.append(nMol)


    outputSet.updateMolClass()
    self._defineOutputs(outputSmallMolecules=outputSet)
    self._defineSourceRelation(self.inputSmallMolecules, outputSet)




  def getOutputCSV(self):
    sysName = self.getSystemName()
    return self._getPath(os.path.join(sysName, 'results/predictions.csv'))


  def getScoreDic(self):
    mapDic = self.parseCSVDic(self.getMapSMIFile())
    smiScoreDic = self.parseCSVDic(self.getOutputCSV())
    scoreDic = {molFile: float(eval(smiScoreDic[smi])[0]) for molFile, smi in mapDic.items()}
    return scoreDic

  def parseCSVDic(self, csvFile):
    smiDic = {}
    with open(csvFile) as f:
      for line in f:
        sline = line.strip().split(',')
        smiDic[sline[0]] = sline[1]
    return smiDic

  def getScriptPath(self):
    return autodockPlugin.getGCRPath(f'{autodockPlugin.getEnvName(GCR_DIC)}/main.py')

  def getConfFile(self):
    return os.path.abspath(self._getExtraPath('config.yaml'))

  def writeConfFile(self):
    confFile = self.getConfFile()
    gpuIdx = getattr(self, params.GPU_LIST).get().split(',')[0].strip()
    regLayers = eval(self.regLayers.get().strip()) + [1]
    with open(confFile, 'w') as f:
      f.write(f'frozenEncoder: True\n'
              f'layers: {self.encLayers.get()}\n'
              f'filtSizes: {self.filtSizes.get()}\n'
              f'pool: "{self.getEnumText("encPool").lower()}"\n\n'
              f'prePropFunc: {self.getEnumText("preprocess")}\n'
              f'regLayers: {regLayers}\n\n'
              f'seed: {self.seed.get()}\n'
              f'cuda: "cuda:{gpuIdx}"\n'
              f'batchSize: {self.batch.get()}\n')
    return confFile

  def buildSMIsFile(self, dMols, writeScores=True):
    getFileFunc = 'getPoseFile' if writeScores else 'getFileName'
    molFile = getattr(dMols.getFirstItem(), getFileFunc)()

    inp = 'train' if writeScores else 'predict'
    smiFile = self.getInputSMIFile(inp)
    if not self.checkHasSMI(molFile):
      inDockFile = self.writeInputDocksFile(dMols, writeScores)
      args = f'{inDockFile} {smiFile} {self.getMapSMIFile(inp)}'
      autodockPlugin.runScript(self, 'convertToSMIs.py', args, envDict=RDKIT_DIC)
    else:
      self.writeMeekoSMIs(dMols, writeScores)

    if writeScores:
      self.mergeSMIscores()
    return smiFile

  def getSystemName(self):
    if self.loadModel.get():
      sysName = self.pretrainedModel.split('/')[-1]
    else:
      if self.modelName.get().strip():
        sysName = self.modelName.get().strip()
      else:
        sysFile = self.dockedMols.get().getProteinFile()
        sysName = f"{sysFile.split('/')[-1]}_{self.getObjId()}"
    return sysName

  def writeMeekoSMIs(self, mols, writeScores=True):
    inp = 'train' if writeScores else 'predict'
    getFileFunc = 'getPoseFile' if writeScores else 'getFileName'
    with open(self.getMapSMIFile(inp), 'w') as fMap:
      with open(self.getInputSMIFile(), 'w') as f:
        for mol in mols:
          molFile = getattr(mol, getFileFunc)()

          smi = self.parseMeekoSMI(molFile)
          fMap.write(f'{molFile},{smi}\n')

          line = smi
          if writeScores:
            score = getattr(mol, self.scoreName.get())
            line += f',{score}'
          line += '\n'

          f.write(line)

  def parseMeekoSMI(self, molFile):
    with open(molFile) as f:
      for line in f:
        if 'REMARK SMILES' in line:
          smi = line.split('REMARK SMILES')[1].strip()
          return smi

  def checkHasSMI(self, molFile):
    with open(molFile) as f:
      txt = f.read()
    return 'REMARK SMILES' in txt

  def getInputDockFile(self):
    return self._getExtraPath('inputDock.csv')

  def getInputSMIFile(self, inp='train'):
    return self._getExtraPath(f'inputSMIs_{inp}.csv')

  def getMapSMIFile(self, inp='predict'):
    return self._getExtraPath(f'mapSMIs_{inp}.csv')

  def writeInputDocksFile(self, mols, writeScores=True):
    getFileFunc = 'getPoseFile' if writeScores else 'getFileName'
    inDocksFiles = self.getInputDockFile()
    with open(inDocksFiles, 'w') as f:
      for m in mols:
        line = getattr(m, getFileFunc)()
        if writeScores:
          score = getattr(m, self.scoreName.get())
          line += f",{score}"
        line += '\n'

        f.write(line)
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

