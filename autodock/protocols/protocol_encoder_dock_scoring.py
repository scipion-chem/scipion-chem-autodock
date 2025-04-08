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

from pwchem.utils import performBatchThreading, findThreadFiles, concatFiles, splitFile, makeSubsets
from pwchem import Plugin as pwchemPlugin
from pwchem.constants import RDKIT_DIC
from pwchem.objects import SmallMoleculesLibrary

from autodock import Plugin as autodockPlugin
from autodock.protocols import ProtChemAutodockGPU
from autodock.constants import GCR_DIC

RESNET, CHEMPROP, COATI = 'ResNet', 'ChemProp', 'COATI'

encoderOptions = [RESNET, CHEMPROP, COATI]
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
    group.addParam('useLibrary', params.BooleanParam, label='Use library as input : ', default=False,
                   expertLevel=params.LEVEL_ADVANCED,
                   help='Whether to use a SMI library SmallMoleculesLibrary object as input')

    group.addParam('inputLibrary', params.PointerParam, pointerClass="SmallMoleculesLibrary",
                   label='Input library: ', condition='useLibrary',
                   help="Input Small molecules library to predict")
    group.addParam('inputSmallMolecules', params.PointerParam, pointerClass="SetOfSmallMolecules",
                   label='Input small molecules: ', condition='not useLibrary',
                   help="Input small molecules to be scored with the model")
    group.addParam('applyFilter', params.BooleanParam, label='Filter results: ', default=False,
                   help='Whether to filter the results by score')
    group.addParam('outThres', params.FloatParam, label='Score threshold: ', default=-7.0, condition='applyFilter',
                   help='Score threshold to use. Molecules with scores over this threshold will not be registered '
                        'in the output')

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
    group.addParam('nTrain', params.IntParam, label='Number of epochs: ', default=20, condition=gonnaTrain,
                   help='Number of epochs for training.')
    group.addParam('scoreName', params.StringParam, label='Docking score: ', default='', condition=gonnaTrain,
                   help='Docking score to use for the model training')
    group.addParam('scoreMerge', params.EnumParam, label='Merge strategy: ', choices=['Min', 'Max', 'Mean'], default=0,
                   expertLevel=params.LEVEL_ADVANCED, condition=gonnaTrain,
                   help='How to merge the scores if several values are found for the same molecule (because of '
                        'conformers or poses')

    group.addParam('predictInParts', params.BooleanParam, label='Predict in batches: ', default=True,
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
  def getGPUInputs(self, nGPUs):
    inMols, inSMIFiles = [None for i in range(nGPUs)], [None for i in range(nGPUs)]
    if not self.useLibrary.get():
      inMols = makeSubsets(self.inputSmallMolecules.get(), nGPUs, cloneItem=True)
    else:
      oDir = os.path.abspath(self._getTmpPath())
      libFile = os.path.abspath(self.inputLibrary.get().getFileName())
      inSMIFiles = splitFile(libFile, n=nGPUs, oDir=oDir, remove=False)
    return inMols, inSMIFiles

  def _insertAllSteps(self):
      tSteps = []
      if not self.loadModel.get() or (self.loadModel.get() and self.doTrain.get()):
        tSteps.append(self._insertFunctionStep('trainingStep'))

      gpuIdxs = getattr(self, params.GPU_LIST).get().split(',')
      inMols, inSMIFiles = self.getGPUInputs(len(gpuIdxs))

      pSteps = []
      for i, gId in enumerate(gpuIdxs):
        pSteps.append(self._insertFunctionStep('predictionStep', gId, inMols[i], inSMIFiles[i], prerequisites=tSteps))

      self._insertFunctionStep('createOutputStep', prerequisites=pSteps)

  def trainingStep(self):
    encoderName = self.getEnumText('encoder').lower()
    sysName = self.getSystemName()
    gpuIdx = getattr(self, params.GPU_LIST).get().split(',')[0].strip()
    confFile = self.writeConfFile(gpuIdx)
    scriptName = self.getScriptPath()

    dMols = self.dockedMols.get()
    smisFile = self.buildSMIsFile(dMols, writeScores=True)[0]
    args = f'--config {confFile} -n {sysName} -e {encoderName} ' \
           f'-p {smisFile} --doTest --testBest 1 --doTrain --trainN {self.nTrain.get()} '

    if encoderName == CHEMPROP.lower():
      args += f'-ef {autodockPlugin.getChemPropFile()} '

    modelsPath = os.path.abspath(autodockPlugin.getPluginHome('models'))
    pwchemPlugin.runCondaCommand(self, args, GCR_DIC, f'python {scriptName}', cwd=self._getPath())

    shutil.copytree(os.path.abspath(self._getPath(sysName)), os.path.join(modelsPath, sysName), dirs_exist_ok=True)


  def predictionStep(self, gpuIdx, inMols=None, inSMIFile=None):
    os.mkdir(self._getExtraPath(f'gpu_{gpuIdx}'))
    confFile = self.writeConfFile(gpuIdx)
    sysName = self.getSystemName()
    scriptName = self.getScriptPath()

    if inMols:
      smisFiles = self.buildSMIsFile(inMols, writeScores=False, gpuIdx=gpuIdx)
    elif inSMIFile:
      nt = self.numberOfThreads.get()
      smiFile = self.getInputSMIFile(writeScores=False, gpuIdx=gpuIdx)
      os.link(inSMIFile, smiFile)
      smisFiles = splitFile(smiFile, n=nt, remove=True, pref='inputSMIs_predict')

    modelsPath = os.path.abspath(autodockPlugin.getPluginHome('models'))
    shutil.copytree(os.path.join(modelsPath, sysName), os.path.abspath(self._getPath(sysName)), dirs_exist_ok=True)

    for i, smisFile in enumerate(smisFiles):
      it = f'{gpuIdx}{i}'
      args = f'--config {confFile} -n {sysName} -d {sysName} -p {smisFile} --doPredict -it {it} '
      if self.getEnumText('encoder') == CHEMPROP:
        args += f'-ef {autodockPlugin.getChemPropFile()} '

      pwchemPlugin.runCondaCommand(self, args, GCR_DIC, f'python {scriptName}', cwd=self._getPath())

  def createOutputStep(self):
    smiScoreDic = self.getScoreDic()

    if self.useLibrary.get():
        inLib, oLibFile = self.inputLibrary.get(), self._getPath('outputLibrary.smi')

        mapDic = inLib.getLibraryMap(fullLine=True)
        with open(oLibFile, 'w') as f:
          for smiName, score in smiScoreDic.items():
            if not self.applyFilter.get() or score < self.outThres.get():
              f.write(f'{mapDic[smiName]}\t{score}\n')

        prevHeaders = inLib.getHeaders()
        outputLib = inLib.clone()
        outputLib.setFileName(oLibFile)
        outputLib.setHeaders(prevHeaders + ['Conplex_score'])
        self._defineOutputs(outputLibrary=outputLib)

    else:
        scoreDic = self.mapMolScoreDic(smiScoreDic)
        outputSet = self.inputSmallMolecules.get().createCopy(self._getPath(), copyInfo=True)
        for mol in self.inputSmallMolecules.get():
          nMol = mol.clone()
          molFile = nMol.getFileName()
          if molFile in scoreDic:
            score = scoreDic[molFile]
            if not self.applyFilter.get() or score < self.outThres.get():
              setattr(nMol, '_gcrScore', params.Float(score))
              outputSet.append(nMol)
        outputSet.updateMolClass()
        self._defineOutputs(outputSmallMolecules=outputSet)

  ############# UTILS FUNCTIONS ###################

  def getOutputCSV(self):
    sysName = self.getSystemName()
    oFile = os.path.abspath(self._getPath(os.path.join(sysName, 'results/predictions.csv')))
    if not os.path.exists(oFile):
      threadFiles = findThreadFiles(oFile)
      concatFiles(threadFiles, oFile, remove=True, skipHead=1)

    return oFile

  def mapMolScoreDic(self, smiScoreDic):
    '''Maps the smi to the roiginal files and retuns: {molFile: score}'''
    mapDic = self.parseCSVDic(self.getMapSMIFile(writeScores=False))
    scoreDic = {molFile: float(smiScoreDic[smi]) for molFile, smi in mapDic.items() if smi in smiScoreDic}
    return scoreDic

  def writeSMIOutput(self, smi, smiName, oDir):
    oFile = os.path.join(oDir, f'{smiName}.smi')
    with open(oFile, 'w') as f:
      f.write(f'{smi} {smiName}\n')
    return oFile

  def getScoreDic(self):
    '''Return a dic as {smi: score}'''
    return self.parseCSVDic(self.getOutputCSV(), isScore=True)

  def parseCSVDic(self, csvFile, isScore=False):
    '''Returns a dic: {molFile: smi}'''
    smiDic = {}
    with open(csvFile) as f:
      for line in f:
        sline = line.strip().split(',')
        smiDic[sline[0]] = eval(sline[1])[0] if isScore else sline[1]
    return smiDic

  def getScriptPath(self):
    return autodockPlugin.getGCRPath(f'{autodockPlugin.getEnvName(GCR_DIC)}/main.py')

  def getConfFile(self, gpuIdx):
    return os.path.abspath(self._getExtraPath(f'config_{gpuIdx}.yaml'))

  def writeConfFile(self, gpuIdx):
    confFile = self.getConfFile(gpuIdx)
    regLayers = eval(self.regLayers.get().strip()) + [1]
    with open(confFile, 'w') as f:
      f.write(f'frozenEncoder: True\n'
              f'layers: {self.encLayers.get()}\n'
              f'filtSizes: {self.filtSizes.get()}\n'
              f'pool: "{self.getEnumText("encPool").lower()}"\n\n'
              f'prePropFunc: {self.getEnumText("preprocess")}\n'
              f'regLayers: {regLayers}\n\n'
              f'seed: {self.seed.get()}\n'
              f'device: "cuda:{gpuIdx}"\n'
              f'batchSize: {self.batch.get()}\n')
    return confFile

  def buildSMIsFile(self, dMols, writeScores=True, gpuIdx=0):
    nt = self.numberOfThreads.get()
    performBatchThreading(self.buildSMIsFileThread, dMols, nt, cloneItem=True, writeScores=writeScores, gpuIdx=gpuIdx)

    if not self.predictInParts.get() or writeScores:
      smiFiles = [self.mergeSMIFiles(writeScores, gpuIdx=gpuIdx)]
    else:
      smiFile = self.getInputSMIFile(writeScores, gpuIdx=gpuIdx)
      smiFiles = findThreadFiles(smiFile)

    self.mergeSMIFiles(writeScores, key='map', gpuIdx=gpuIdx)
    self.mergeSMIFiles(writeScores, key='dock', gpuIdx=gpuIdx)

    smiFiles = [os.path.abspath(smiFile) for smiFile in smiFiles]
    return smiFiles

  def mergeSMIFiles(self, writeScores, key='smi', gpuIdx=0):
    funcDic = {'smi': 'getInputSMIFile', 'map': 'getMapSMIFile', 'dock': 'getInputDockFile'}
    getFileFunc = funcDic[key]
    smiFile = getattr(self, getFileFunc)(writeScores, gpuIdx=gpuIdx)
    smiThreadFiles = findThreadFiles(smiFile)
    if len(smiThreadFiles) > 0:
      concatFiles(smiThreadFiles, smiFile, remove=True)

    return smiFile

  def buildSMIsFileThread(self, dMols, outLists, it, writeScores=True, gpuIdx=0):
    smiFile = self.getInputSMIFile(writeScores, it, gpuIdx=gpuIdx)

    fMol = dMols[0]
    molFile = fMol.getFileName()
    if molFile.endswith('.smi'):
      self.writeSMIs(dMols, writeScores, it, origin='file', gpuIdx=gpuIdx)
    else:
      getFileFunc = 'getPoseFile' if writeScores else 'getFileName'
      molFile = getattr(dMols[0], getFileFunc)()

      if self.checkHasSMI(molFile):
        smiFile, mapFile = self.writeSMIs(dMols, writeScores, it, origin='Meeko', gpuIdx=gpuIdx)
      else:
        inMolsFile = self.writeInputMolsFile(dMols, writeScores, it, gpuIdx=gpuIdx)
        mapFile = self.getMapSMIFile(writeScores, it, gpuIdx=gpuIdx)
        args = f'{inMolsFile} {smiFile} {mapFile}'
        autodockPlugin.runScript(self, 'convertToSMIs.py', args, envDict=RDKIT_DIC, popen=True)

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

  def writeSMIs(self, mols, writeScores=True, it=None, origin='Meeko', gpuIdx=0):
    getFileFunc = 'getPoseFile' if (writeScores and origin == 'Meeko') else 'getFileName'

    smiText, mapText = '', ''
    smiFile, mapFile = self.getInputSMIFile(writeScores, it, gpuIdx=gpuIdx), \
                       self.getMapSMIFile(writeScores, it, gpuIdx=gpuIdx)
    for mol in mols:
      molFile = getattr(mol, getFileFunc)()
      smi = self.parseMeekoSMI(molFile) if origin == 'Meeko' else self.parseFileSMI(molFile)
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

  def parseFileSMI(self, molFile):
    with open(molFile) as f:
      lines = f.read().strip().split('\n')
      smi = lines[-1].split()[0].strip()
    return smi

  def parseMeekoSMI(self, molFile):
    with open(molFile) as f:
      for line in f:
        if 'REMARK SMILES' in line:
          smi = line.split('REMARK SMILES')[1].strip()
          break
    return smi

  def checkHasSMI(self, molFile):
    with open(molFile) as f:
      txt = f.read()
    return 'REMARK SMILES' in txt

  def getInputDockFile(self, writeScores=True, it=None, gpuIdx=0):
    oDir = self._getExtraPath(f'gpu_{gpuIdx}')
    inp = 'train' if writeScores else 'predict'
    csvFile = os.path.join(oDir, f'inputDock_{inp}.csv')
    if it is not None:
      csvFile = csvFile.replace('.csv', f'_{it}.csv')
    return os.path.abspath(csvFile)

  def getInputSMIFile(self, writeScores=True, it=None, gpuIdx=0):
    oDir = self._getExtraPath(f'gpu_{gpuIdx}')
    inp = 'train' if writeScores else 'predict'
    smiFile = os.path.join(oDir, f'inputSMIs_{inp}.csv')
    if it is not None:
        smiFile = smiFile.replace('.csv', f'_{it}.csv')
    return os.path.abspath(smiFile)

  def getMapSMIFile(self, writeScores=True, it=None, gpuIdx=0):
    oDir = self._getExtraPath(f'gpu_{gpuIdx}')
    inp = 'train' if writeScores else 'predict'
    smiFile = os.path.join(oDir, f'mapSMIs_{inp}.csv')
    if it is not None:
        smiFile = smiFile.replace('.csv', f'_{it}.csv')
    return os.path.abspath(smiFile)

  def writeInputMolsFile(self, mols, writeScores=True, it=None, gpuIdx=0):
    getFileFunc = 'getPoseFile' if writeScores else 'getFileName'
    inDocksFiles = self.getInputDockFile(writeScores, it, gpuIdx=gpuIdx)

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

