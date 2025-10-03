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

import os, json

from pyworkflow.protocol import params
import pyworkflow.object as pwobj

from pwem.protocols import EMProtocol

from pwchem.utils import RESIDUES1TO3, getBaseName
from pwchem.objects import SetOfSmallMolecules, SmallMolecule

from autodock import Plugin as autodockPlugin
from autodock.objects import RingtailDatabase
from autodock.constants import VINA

scoreOptions = {'None': None, 'Energy': 'e', 'Ligand efficiency': 'le',
                'Energy percentile': 'pe', 'Ligand efficiency percentile': 'ple'}
clusterOptions = {'None': None, 'Morgan fingerprints': 'mfpc', 'Interactions fingerprints': 'ifpc'}


def splitSDF(sdfFile, oriName=None):
  '''Split multi molecule sdf files.
  Takes into account different conformations of the same molecule.
  '''

  def getMolName(molText, splitConf=True):
    molName, confId = molText.split('\n')[0].strip(), None
    if '_i' in molName and splitConf:
      molName, confId = molName.split('_i')
    return molName, confId

  molDic, confFile = {}, None
  oDir = os.path.dirname(sdfFile)

  with open(sdfFile) as f:
    sdfText = f.read().strip()
  molTexts = [molText.strip() for molText in sdfText.split('$$$$') if molText.strip()]

  if len(molTexts) > 1:
    if oriName is not None:
      confFile = sdfFile.replace('.sdf', f'_{oriName}.sdf')
      os.rename(sdfFile, confFile)

    for molText in molTexts:
      molName, confId = getMolName(molText)
      oFile = os.path.join(oDir, f'{molName}.sdf')
      if confId is not None:
        oFile = oFile.replace('.sdf', f'_{confId}.sdf')

      if molName not in molDic:
        molDic[molName] = []
      else:
        poseId = len(molDic[molName]) + 1
        oFile = oFile.replace('.sdf', f'_{poseId}.sdf')

      molDic[molName].append(oFile)
      with open(oFile, 'w') as fo:
        fo.write(f'{molText}\n\n$$$$')

  else:
    molName, confId = getMolName(molTexts[0])
    molDic = {molName: [sdfFile]}

  return molDic, confFile

class ProtRingtailFilter(EMProtocol):
  """Executes a serie of Ringtail filters over a docking database"""
  _label = 'Ringtail database filter'

  def _defineParams(self, form):
    form.addSection(label='Input')
    group = form.addGroup('Input')
    group.addParam('inputRingtail', params.PointerParam, pointerClass="RingtailDatabase",
                   label='Input Ringtail database: ', help="Input Ringtail database to filter")
    group.addParam('bookmark', params.StringParam, label='Bookmark name: ', default='passing_results',
                   help='Bookmark for the group of molecules passing the filter')

    group = form.addGroup('Score filter')
    group.addParam('scoreFilt', params.EnumParam, label='Score filtering: ', choices=list(scoreOptions.keys()), default=0,
                   help='Filter ligands by the selected score.\nLigand efficiency = energy / number of heavy atoms'
                        '\nPercentile expressed as percentage e.g. 1 for top 1 percent.')
    group.addParam('scoreValue', params.FloatParam, label='Highest energy: ', condition='scoreFilt!=0',
                   help='Specify the worst (highest) energy value accepted')

    group = form.addGroup('Structure filter')
    group.addParam('setMax', params.BooleanParam, label='Set maximum ligand size: ', default=False,
                   help='Whether to set the maximum ligand size in number of atoms to pass the filter.')
    group.addParam('maxAtoms', params.IntParam, label='Maximum number of atoms: ', condition='setMax', default=20,
                   help='Maximum number of heavy atoms (non-hydrogens) a ligand may have')
    group.addParam('smarts', params.StringParam, label='SMARTS substructure: ', default='',
                   help='SMARTS substructure the ligands need to fulfill to pass the filter.')

    group = form.addGroup('Clustering')
    group.addParam('cluster', params.EnumParam, label='Clustering method: ', choices=list(clusterOptions.keys()), default=0,
                   help='Cluster filtered ligands by Tanimoto distance of the selected fingerprints with Butina '
                        'clustering and output ligand with lowest ligand efficiency from each cluster. '
                        '\n\nMorgan: Useful for selecting chemically dissimilar ligands.'
                        '\nInteractions: Useful for enhancing selection of ligands with diverse interactions.')
    group.addParam('clustCut', params.FloatParam, label='Cluster cut-off: ', default=0.5, condition='cluster!=0',
                   help='Clustering cutoff to use.')

    group = form.addGroup('Output')
    group.addParam('outMols', params.BooleanParam, label='Output mols: ', default=True,
                   help='Generate a small molecules output')
    group.addParam('outBestMol', params.BooleanParam, label='Output only best pose: ', default=True,
                   help='Output only the best scoring pose for each molecule')

    form.addSection('Interactions')
    group = form.addGroup('Selection')
    group.addParam('selChain', params.StringParam, label='Select chain: ', default='',
                   help='Specify the receptor chain to define an interaction.')
    group.addParam('selResidue', params.StringParam, label='Select residue: ', default='',
                   help='Specify the residue position to define an interaction')
    group.addParam('selAtom', params.StringParam, label='Select atom: ', default='',
                   help='Specify the atom to define an interaction.\n'
                        ' It can be left empty, then the whole residue will be considered')

    group = form.addGroup('Interactions filter')
    group.addParam('minHB', params.IntParam, label='Minimum number of hydrogen bonds: ', default=0,
                   help='Minimum number of hydrogen bonds a ligand may have')
    group.addParam('hbInt', params.TextParam, label='Hydrogen bond interactions: ', default='',
                   help='Specify the receptor selection to filter for hydrogen bonds.\n'
                        'Use the wizard to append a receptor selection')
    group.addParam('vdwInt', params.TextParam, label='Van de Waals interactions: ', default='',
                   help='Specify the receptor selection to filter for Van der Waals interactions.\n'
                        'Use the wizard to append a receptor selection')


  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
    self._insertFunctionStep('filterStep')
    self._insertFunctionStep('createOutputStep')

  def filterStep(self):
    inDB = self.inputRingtail.get()
    kwargs = {'bookmark': self.bookmark.get()}

    scFilt = scoreOptions[self.getEnumText("scoreFilt")]
    kwargs['scoreDic'] = {scFilt: self.scoreValue.get()} if scFilt else {}
    clFilt = clusterOptions[self.getEnumText("cluster")]
    kwargs['clusterDic'] = {clFilt: self.clustCut.get()} if clFilt else {}

    kwargs['maxAtoms'] = self.maxAtoms.get() if self.setMax.get() else None
    kwargs['smarts'] = self.smarts.get().strip() if self.smarts.get().strip() else None
    kwargs['minHB'] = self.minHB.get()
    kwargs['vdwIntLines'] = self.vdwInt.get().strip()
    kwargs['hbIntLines'] = self.vdwInt.get().strip()

    kwargs['outDir'] = os.path.abspath(self._getPath()) if self.outMols.get() else None
    kwargs['outBest'] = self.outBestMol.get()

    args = inDB.buildFilterArgs(**kwargs)
    autodockPlugin.runRingtail(self, args, cwd=self._getPath())

  def createOutputStep(self):
    inDB = self.inputRingtail.get()
    if self.outMols.get():
      outDir = os.path.abspath(self._getPath())
      posesDic = self.getOutputDic()

      outputSet = SetOfSmallMolecules().create(outputPath=outDir)
      for molName, sdfFiles in posesDic.items():
        energies, effs = self.parsePosesFile(sdfFiles[0])
        for i, sFile in enumerate(sdfFiles):
          newSmallMol = SmallMolecule(smallMolFilename=sFile, type='AutoDock')
          newSmallMol.setMolName(molName)
          newSmallMol.setPoseFile(sFile)
          newSmallMol.setPoseId(i+1)
          newSmallMol._energy = pwobj.Float(energies[i])
          newSmallMol._ligandEfficiency = pwobj.Float(effs[i])

          outputSet.append(newSmallMol)

      outputSet.proteinFile.set(inDB.getReceptorFile())
      outputSet.setDocked(True)
      self._defineOutputs(outputSmallMolecules=outputSet)

    else:
      dbFile = os.path.abspath(inDB.getFileName())
      outputDB = RingtailDatabase(filename=dbFile)
      outputDB.createSumFile(self.getSumPath())
      self._defineOutputs(outputRingtail=outputDB)


  ################ UTILS functions ##############################

  def createElementLine(self):
    chainStr = json.loads(self.selChain.get())['chain']
    resDic = json.loads(self.selResidue.get())
    resName, resNum = RESIDUES1TO3[resDic['residues']], resDic['index'].split('-')[0]

    selStr = f'{chainStr}:{resName}:{resNum}:'
    if self.selAtom.get():
      atomDic = json.loads(self.selAtom.get())
      selStr += f'{atomDic["atom"]}'
    return selStr

  def getInteractionsArgs(self, vdw=True):
    args = ''
    if vdw and self.vdwInt.get().strip():
      vwdList = self.vdwInt.get().strip().split('\n')
      args += f'-vdw {"-vdw ".join(vwdList)} '
    elif not vdw and self.hbInt.get().strip():
      hbList = self.hbInt.get().strip().split('\n')
      args += f'-hb {"-hb ".join(hbList)} '
    return args

  def getOutputDic(self):
    '''Return a dic of the form: {molName: [sdfFileNames]}
    '''
    sdfFilesDic = {getBaseName(f): [self._getPath(f)] for f in os.listdir(self._getPath()) if '.sdf' in f}
    if not self.outBestMol.get():
      newSDFiles = {}
      for molName, sFiles in sdfFilesDic.items():
        molDic, _ = splitSDF(sFiles[0], oriName='poses')
        newSDFiles.update(molDic)
      sdfFilesDic = newSDFiles
    return sdfFilesDic

  def parsePosesFile(self, pFile):
    '''Return two lists with energies and ligand efficiencies stored in a sdf file
    '''
    energies, effs = [], []
    energy, eff = False, False

    with open(pFile) as f:
      for line in f:
        if energy:
          energies = eval(line.strip())
        if eff:
          effs = eval(line.strip())

        if 'Binding energies' in line:
          energy, eff = True, False
        elif 'Ligand effiencies' in line:
          energy, eff = False, True
        else:
          energy, eff = False, False

    return energies, effs


  def getSumPath(self):
    return self._getExtraPath('ringSum.txt')

  def _summary(self):
    s = []
    if os.path.exists(self.getSumPath()):
      with open(self.getSumPath()) as f:
        s.append(f.read())
    return s
