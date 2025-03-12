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

import os, shutil

from pyworkflow.protocol import params

from pwchem.utils import performBatchThreading, findThreadFiles, concatFiles
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

    group.addParam('predictInParts', params.BooleanParam, label='Predict in batches: ', default=False,
                   expertLevel=params.LEVEL_ADVANCED,
                   help='Whether to perform the predictions in input batches to avoid overloading the memory')
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

    form.addParallelSection(threads=4, mpi=1)

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
    smisFile = self.buildSMIsFile(dMols, writeScores=True)[0]
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
    smisFiles = self.buildSMIsFile(inMols, writeScores=False)

    modelsPath = os.path.abspath(autodockPlugin.getPluginHome('models'))
    shutil.copytree(os.path.join(modelsPath, sysName), os.path.abspath(self._getPath(sysName)), dirs_exist_ok=True)

    for i, smisFile in enumerate(smisFiles):
      args = f'--config {confFile} -n {sysName} -d {sysName} -p {smisFile} --doPredict ' \
             f'-ef {autodockPlugin.getChemPropFile()} -it {i}'

      pwchemPlugin.runCondaCommand(self, args, GCR_DIC, f'python {scriptName}', cwd=self._getPath())


  def createOutputStep(self):
    scoreDic = self.getScoreDic()
    outputSet = self.inputSmallMolecules.get().createCopy(self._getPath(), copyInfo=True)
    for mol in self.inputSmallMolecules.get():
      nMol = mol.clone()
      molFile = nMol.getFileName()
      if molFile in scoreDic:
        setattr(nMol, '_gcrScore', params.Float(scoreDic[molFile]))
        outputSet.append(nMol)


    outputSet.updateMolClass()
    self._defineOutputs(outputSmallMolecules=outputSet)
    self._defineSourceRelation(self.inputSmallMolecules, outputSet)




  def getOutputCSV(self):
    sysName = self.getSystemName()
    oFile = os.path.abspath(self._getPath(os.path.join(sysName, 'results/predictions.csv')))
    if not os.path.exists(oFile):
      threadFiles = findThreadFiles(oFile)
      concatFiles(threadFiles, oFile, remove=True)

    return oFile


  def getScoreDic(self):
    mapDic = self.parseCSVDic(self.getMapSMIFile(writeScores=False))
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
    nt = self.numberOfThreads.get()
    performBatchThreading(self.buildSMIsFileThread, dMols, nt, cloneItem=True, writeScores=writeScores)

    if not self.predictInParts.get() or writeScores:
      smiFiles = [self.mergeSMIFiles(writeScores)]
    else:
      smiFile = self.getInputSMIFile(writeScores)
      smiFiles = findThreadFiles(smiFile)

    self.mergeSMIFiles(writeScores, key='map')
    self.mergeSMIFiles(writeScores, key='dock')

    smiFiles = [os.path.abspath(smiFile) for smiFile in smiFiles]
    return smiFiles

  def mergeSMIFiles(self, writeScores, key='smi'):
    funcDic = {'smi': 'getInputSMIFile', 'map': 'getMapSMIFile', 'dock': 'getInputDockFile'}
    getFileFunc = funcDic[key]
    smiFile = getattr(self, getFileFunc)(writeScores)
    smiThreadFiles = findThreadFiles(smiFile)
    if len(smiThreadFiles) > 0:
      concatFiles(smiThreadFiles, smiFile, remove=True)

    return smiFile

  def buildSMIsFileThread(self, dMols, outLists, it, writeScores=True):
    getFileFunc = 'getPoseFile' if writeScores else 'getFileName'
    molFile = getattr(dMols[0], getFileFunc)()

    smiFile = self.getInputSMIFile(writeScores, it)
    if not self.checkHasSMI(molFile):
      inMolsFile = self.writeInputMolsFile(dMols, writeScores, it)
      args = f'{inMolsFile} {smiFile} {self.getMapSMIFile(writeScores, it)}'
      autodockPlugin.runScript(self, 'convertToSMIs.py', args, envDict=RDKIT_DIC, popen=True)
    else:
      smiFile, mapFile = self.writeMeekoSMIs(dMols, writeScores, it)

    if writeScores:
      self.mergeSMIscores(smiFile)




  def getSystemName(self):
    if self.loadModel.get():
      sysName = self.pretrainedModel.get().split('/')[-1]
    else:
      if self.modelName.get().strip():
        sysName = self.modelName.get().strip()
      else:
        sysFile = self.dockedMols.get().getProteinFile()
        sysName = f"{sysFile.split('/')[-1]}_{self.getObjId()}"
    return sysName

  def writeMeekoSMIs(self, mols, writeScores=True, it=None):
    getFileFunc = 'getPoseFile' if writeScores else 'getFileName'

    smiText, mapText = '', ''
    smiFile, mapFile = self.getInputSMIFile(writeScores, it), self.getMapSMIFile(writeScores, it)
    for mol in mols:
      molFile = getattr(mol, getFileFunc)()
      smi = self.parseMeekoSMI(molFile)
      mapText += f'{molFile},{smi}\n'

      line = smi
      if writeScores:
        score = getattr(mol, self.scoreName.get())
        line += f',{score}'
      smiText += f'{line}\n'

    with open(mapFile, 'w') as fMap:
      fMap.write(mapText)
    with open(smiFile, 'w') as f:
      f.write(smiText)

    return smiFile, mapFile

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

  def getInputDockFile(self, writeScores=True, it=None):
    inp = 'train' if writeScores else 'predict'
    csvFile = self._getExtraPath(f'inputDock_{inp}.csv')
    if it is not None:
      csvFile = csvFile.replace('.csv', f'_{it}.csv')
    return os.path.abspath(csvFile)

  def getInputSMIFile(self, writeScores=True, it=None):
    inp = 'train' if writeScores else 'predict'
    smiFile = self._getExtraPath(f'inputSMIs_{inp}.csv')
    if it is not None:
        smiFile = smiFile.replace('.csv', f'_{it}.csv')
    return os.path.abspath(smiFile)

  def getMapSMIFile(self, writeScores=True, it=None):
    inp = 'train' if writeScores else 'predict'
    smiFile = self._getExtraPath(f'mapSMIs_{inp}.csv')
    if it is not None:
        smiFile = smiFile.replace('.csv', f'_{it}.csv')
    return os.path.abspath(smiFile)

  def writeInputMolsFile(self, mols, writeScores=True, it=None):
    getFileFunc = 'getPoseFile' if writeScores else 'getFileName'
    inDocksFiles = self.getInputDockFile(writeScores, it)

    text = ''
    for m in mols:
      line = getattr(m, getFileFunc)()
      if writeScores:
        score = getattr(m, self.scoreName.get())
        line += f",{score}"
      text += f'{line}\n'

    with open(inDocksFiles, 'w') as f:
      f.write(text)
    return inDocksFiles

  def getMergeFunction(self):
    mergeDic = {0: min, 1: max, 2: average}
    return mergeDic[self.scoreMerge.get()]

  def mergeSMIscores(self, smiFile):
    # parsing
    smiDic = {}
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
    smiText = ''
    for smi, score in smiDic.items():
      smiText += f'{smi},{score}\n'

    with open(smiFile, 'w') as f:
      f.write(smiText)

    return smiFile

