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

import os, glob

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, IntParam, FloatParam, STEPS_PARALLEL, BooleanParam, \
  LEVEL_ADVANCED, USE_GPU, GPU_LIST, StringParam
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath, createLink

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import runOpenBabel, generate_gpf, calculate_centerMass, getBaseFileName
from pwchem import Plugin as pwchem_plugin
from pwchem.constants import MGL_DIC

from autodock import Plugin as autodock_plugin


class ProtChemAutodock(EMProtocol):
  """Perform a docking experiment with autodock. Grid must be generated in this protocol in order to take into
       account ligands atom types. See the help at
       http://autodock.scripps.edu/faqs-help/manual/autodock-4-2-user-guide/AutoDock4.2_UserGuide.pdf"""
  _label = 'Autodock docking'
  _program = ""

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)
    self.stepsExecutionMode = STEPS_PARALLEL

  def _defineParams(self, form):
    form.addHidden(USE_GPU, BooleanParam, default=True,
                   label="Use GPU for execution: ",
                   help="This protocol has both CPU and GPU implementation.\
                                         Select the one you want to use.")

    form.addHidden(GPU_LIST, StringParam, default='0',
                   expertLevel=LEVEL_ADVANCED, label="Choose GPU IDs",
                   help="Add a list of GPU devices that can be used")

    form.addSection(label='Input')
    group = form.addGroup('Receptor specification')
    group.addParam('wholeProt', BooleanParam, label='Dock on whole protein: ', default=True,
                   help='Whether to dock on a whole protein surface or on specific regions')

    # Docking on whole protein
    group.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct",
                   label='Input atomic structure:', condition='wholeProt',
                   help="The atom structure to use as receptor in the docking")
    group.addParam('radius', FloatParam, label='Grid radius for whole protein: ',
                   condition='wholeProt', allowsNull=False,
                   help='Radius of the Autodock grid for the whole protein')

    # Docking on pockets
    group.addParam('inputStructROIs', PointerParam, pointerClass="SetOfStructROIs",
                   label='Input pockets:', condition='not wholeProt',
                   help="The protein structural ROIs to dock in")
    group.addParam('pocketRadiusN', FloatParam, label='Grid radius vs StructROI radius: ',
                   condition='not wholeProt', default=1.1, allowsNull=False,
                   help='The radius * n of each StructROI will be used as grid radius')

    group.addParam('spacing', FloatParam, label='Spacing of the grids: ', default=0.5, allowsNull=False,
                   condition='not wholeProt', help='Spacing of the generated Autodock grids')

    group = form.addGroup('Docking')
    group.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                   label='Input small molecules:', allowsNull=False,
                   help="Input small molecules to be docked with AutoDock")
    group.addParam('rmsTol', FloatParam, label='Cluster tolerance (A): ', default=2.0, expertLevel=LEVEL_ADVANCED,
                   help='Maximum RMSD for 2 docked structures to be in the same cluster')
    group.addParam('gaRun', IntParam, label='Number of positions per ligand', default=10)

    form.addSection(label="Genetic algorithm")
    form.addParam('gaPop', IntParam, label='Population size', default=150,
                  help='This is the number of individuals in the population. Each individual is a coupling of a '
                       'genotype and its associated phenotype')
    form.addParam('gaNumEvals', IntParam, label='Number of evaluations', default=2500000,
                  help='Set the maximum number of energy evaluations performed during each GA, LGA, or LS run')
    form.addParam('gaNumGens', IntParam, label='Number of generations', default=27000,
                  help='This is the maximum number of generations simulated during each GA or LGA run')
    form.addParam('gaElitism', IntParam, label='Elitism', default=1,
                  help='This is the number of top individuals that are guaranteed to survive into the next generation.')
    form.addParam('gaMutationRate', FloatParam, label='Mutation rate', default=0.02,
                  help='This is a floating point number from 0 to 1, representing the probability that a particular '
                       'gene is mutated')
    form.addParam('gaCrossOverRate', FloatParam, label='Crossover rate', default=0.8,
                  help='This is a floating point number from 0 to 1 denoting the crossover rate. Crossover rate '
                       'is the expected number of pairs in the population that will exchange genetic material')
    form.addParam('gaWindowSize', IntParam, label='Window size', default=10,
                  help='This is the number of preceding generations to take into consideration when deciding the '
                       'threshold for the worst individual in the current population')
    form.addParam('lsFreq', FloatParam, label='Local search frequency', default=0.06,
                  help='This is the probability of any particular phenotype being subjected to local search')

    form.addSection(label="Local search")
    form.addParam('swMaxIts', IntParam, label='Max. Number of iterations', default=300,
                  help='This is the maximum number of iterations that the local search procedure applies to the '
                       'phenotype of any given individual, per generation')
    form.addParam('swMaxSucc', IntParam, label='Successes in a raw', default=4,
                  help='Number of successes before changing rho in Solis & Wets algorithms')
    form.addParam('swMaxFail', IntParam, label='Failures in a raw', default=4,
                  help='Number of failures before changing rho in Solis & Wets algorithms.')
    form.addParam('swRho', FloatParam, label='Initial variance', default=1.0,
                  help='It defines the size of the local search, and specifies the size of the local space to sample')
    form.addParam('swLbRho', FloatParam, label='Variance lower bound', default=0.01,
                  help='This is the lower bound on rho, the variance for making changes to genes (i.e. translations, '
                       'orientation and torsions)')
    form.addParallelSection(threads=4, mpi=1)

  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
    cId = self._insertFunctionStep('convertStep', prerequisites=[])

    dockSteps = []
    if self.wholeProt:
      gridId = self._insertFunctionStep('generateGridsStep', prerequisites=[cId])
      for mol in self.inputSmallMolecules.get():
        dockId = self._insertFunctionStep('dockStep', mol.clone(), prerequisites=[gridId])
        dockSteps.append(dockId)
    else:
      for pocket in self.inputStructROIs.get():
        gridId = self._insertFunctionStep('generateGridsStep', pocket.clone(), prerequisites=[cId])
        for mol in self.inputSmallMolecules.get():
          dockId = self._insertFunctionStep('dockStep', mol.clone(), pocket.clone(), prerequisites=[gridId])
          dockSteps.append(dockId)

    self._insertFunctionStep('createOutputStep', prerequisites=dockSteps)

  def convertStep(self):
    self.ligandFileNames = []
    for mol in self.inputSmallMolecules.get():
      fnSmall, smallDir = self.convert2PDBQT(mol.clone(), self._getExtraPath('conformers'))
      self.ligandFileNames.append(fnSmall)

    self.receptorFile = self.getOriginalReceptorFile()
    if not '.pdbqt' in self.receptorFile:
      self.receptorFile = self.convertReceptor2PDBQT(self.receptorFile)

  def generateGridsStep(self, pocket=None):
    fnReceptor = self.receptorFile
    outDir = self.getOutputPocketDir(pocket)
    if self.wholeProt:
      radius = self.radius.get()
      # Use the original pdb for mass center
      pdbFile = self.getOriginalReceptorFile()
      if not os.path.splitext(pdbFile)[1] == '.pdb':
        pdbFile = self.convertReceptor2PDB(pdbFile)
      structure, x_center, y_center, z_center = calculate_centerMass(pdbFile)
    else:
      radius = (pocket.getDiameter() / 2) * self.pocketRadiusN.get()
      x_center, y_center, z_center = pocket.calculateMassCenter()

    makePath(outDir)
    npts = (radius * 2) / self.spacing.get()
    gpf_file = generate_gpf(fnReceptor, spacing=self.spacing.get(),
                            xc=x_center, yc=y_center, zc=z_center,
                            npts=npts, outDir=outDir, ligandFns=self.ligandFileNames)

    args = "-p {} -l {}.glg".format(gpf_file, self.getReceptorName())
    self.runJob(autodock_plugin.getAutodockPath("autogrid4"), args, cwd=outDir)

  def dockStep(self, mol, pocket=None):
    # Prepare grid
    fnReceptor = self.receptorFile
    nameReceptor, extReceptor = os.path.splitext(os.path.basename(fnReceptor))
    outDir = self.getOutputPocketDir(pocket)

    molFn = self.getMolLigandName(mol)
    molName, _ = os.path.splitext(os.path.basename(molFn))
    createLink(molFn, os.path.join(outDir, molName + '.pdbqt'))
    makePath(outDir)

    fnDPF = os.path.abspath(os.path.join(outDir, molName + ".dpf"))
    args = " -l %s -r %s -o %s - x no" % (molFn, fnReceptor, fnDPF)
    args += " -p ga_pop_size=%d" % self.gaPop.get()
    args += " -p ga_num_evals=%d" % self.gaNumEvals.get()
    args += " -p ga_num_generations=%d" % self.gaNumGens.get()
    args += " -p ga_elitism=%d" % self.gaElitism.get()
    args += " -p ga_mutation_rate=%f" % self.gaMutationRate.get()
    args += " -p ga_crossover_rate=%f" % self.gaCrossOverRate.get()
    args += " -p ga_window_size=%d" % self.gaWindowSize.get()
    args += " -p sw_max_its=%d" % self.swMaxIts.get()
    args += " -p sw_max_succ=%d" % self.swMaxSucc.get()
    args += " -p sw_max_fail=%d" % self.swMaxFail.get()
    args += " -p sw_rho=%f" % self.swRho.get()
    args += " -p sw_lb_rho=%d" % self.swLbRho.get()
    args += " -p ls_search_freq=%f" % self.lsFreq.get()
    args += " -p ga_run=%d" % self.gaRun.get()
    args += " -p rmstol=%f" % self.rmsTol.get()

    self.runJob(pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh'),
                autodock_plugin.getADTPath('Utilities24/prepare_dpf42.py') + args,
                cwd=outDir)

    fnDLG = fnDPF.replace('.dpf', '.dlg')
    if not getattr(self, USE_GPU):
      args = "-p %s -l %s" % (fnDPF, fnDLG)
      self.runJob(autodock_plugin.getAutodockPath("autodock4"), args, cwd=outDir)
    else:
      fnDPF = self.commentFirstLine(fnDPF)
      args = '-I {}'.format(fnDPF)
      autodock_plugin.runAutodockGPU(self, args, outDir)

  def createOutputStep(self):
    outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())
    for pocketDir in self.getPocketDirs():
      pocketDic = {}
      gridId = self.getGridId(pocketDir)
      for dlgFile in self.getDockedLigandsFiles(pocketDir):
        molName = getBaseFileName(dlgFile)
        pocketDic[molName] = self.parseDockedMolsDLG(dlgFile)
        for modelId in pocketDic[molName]:
          pdbqtFile = self._getPath('{}_{}.pdbqt'.format(molName, modelId))
          with open(pdbqtFile, 'w') as f:
            f.write(pocketDic[molName][modelId]['pdb'])
          pocketDic[molName][modelId]['file'] = pdbqtFile

      for smallMol in self.inputSmallMolecules.get():
        molFile = smallMol.getFileName()
        molName = getBaseFileName(molFile)
        molDic = pocketDic[molName]

        for posId in molDic:
          newSmallMol = SmallMolecule()
          newSmallMol.copy(smallMol, copyId=False)
          newSmallMol._energy = pwobj.Float(molDic[posId]['energy'])
          if 'ki' in molDic[posId]:
            newSmallMol._ligandEfficiency = pwobj.Float(molDic[posId]['ki'])
          else:
            newSmallMol._ligandEfficiency = pwobj.Float(None)
          if os.path.getsize(molDic[posId]['file']) > 0:
            newSmallMol.poseFile.set(molDic[posId]['file'])
            newSmallMol.setPoseId(posId)
            newSmallMol.gridId.set(gridId)
            newSmallMol.setMolClass('Autodock4')
            newSmallMol.setDockId(self.getObjId())

            outputSet.append(newSmallMol)

    outputSet.proteinFile.set(self.getOriginalReceptorFile())
    outputSet.setDocked(True)
    self._defineOutputs(outputSmallMolecules=outputSet)
    self._defineSourceRelation(self.inputSmallMolecules, outputSet)

  ########################### Utils functions ############################
  def getGridId(self, outDir):
    return outDir.split('_')[-1]

  def convert2PDBQT(self, smallMol, oDir):
    '''Convert ligand to pdbqt using prepare_ligand4 of ADT'''
    inFile = smallMol.getFileName()
    if os.path.splitext(inFile)[1] not in ['.pdb', '.mol2', '.pdbq']:
      # Convert to formats recognized by ADT
      outName, outDir = os.path.splitext(os.path.basename(inFile))[0], os.path.abspath(self._getTmpPath())
      args = ' -i "{}" -of mol2 --outputDir "{}" --outputName {}'.format(os.path.abspath(inFile),
                                                                         os.path.abspath(outDir), outName)
      pwchem_plugin.runScript(self, 'obabel_IO.py', args, env='plip', cwd=outDir)
      inFile = self._getTmpPath(outName + '.mol2')

    if not os.path.exists(oDir):
      os.mkdir(oDir)

    inName, inExt = os.path.splitext(os.path.basename(inFile))
    oFile = os.path.abspath(os.path.join(oDir, smallMol.getUniqueName() + '.pdbqt'))

    if inExt != '.pdbqt':
      args = ' -l {} -o {}'.format(inFile, oFile)
      self.runJob(pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh'),
                  autodock_plugin.getADTPath('Utilities24/prepare_ligand4.py') + args)
    else:
      createLink(inFile, oFile)
    return oFile, oDir

  def convertReceptor2PDB(self, proteinFile):
    inName, inExt = os.path.splitext(os.path.basename(proteinFile))
    oFile = os.path.abspath(os.path.join(self._getTmpPath(inName + '.pdb')))

    args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(proteinFile), oFile)
    runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())

    return oFile

  def convertReceptor2PDBQT(self, proteinFile):
    inName, inExt = os.path.splitext(os.path.basename(proteinFile))
    oFile = os.path.abspath(os.path.join(self._getExtraPath(inName + '.pdbqt')))

    args = ' -v -r %s -o %s' % (proteinFile, oFile)
    self.runJob(pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh'),
                autodock_plugin.getADTPath('Utilities24/prepare_receptor4.py') + args)

    return oFile

  def getReceptorName(self):
    if self.wholeProt:
      atomStructFn = self.getOriginalReceptorFile()
      return atomStructFn.split('/')[-1].split('.')[0]
    else:
      return self.inputStructROIs.get().getProteinName()

  def getOriginalReceptorFile(self):
    if self.wholeProt:
      return self.inputAtomStruct.get().getFileName()
    else:
      return self.inputStructROIs.get().getProteinFile()

  def getReceptorDir(self):
    atomStructFn = self.getOriginalReceptorFile()
    return '/'.join(atomStructFn.split('/')[:-1])

  def getLigandsFileNames(self):
    ligFns = []
    for mol in self.inputSmallMolecules.get():
      ligFns.append(mol.getFileName())
    return ligFns

  def getOutputPocketDir(self, pocket=None):
    if pocket == None:
      outDir = self._getExtraPath('pocket_1')
    else:
      outDir = self._getExtraPath('pocket_{}'.format(pocket.getObjId()))
    return outDir

  def getPocketDirs(self):
    dirs = []
    for file in os.listdir(self._getExtraPath()):
      d = self._getExtraPath(file)
      if os.path.isdir(d) and 'pocket' in file:
        dirs.append(d)
    dirs.sort()
    return dirs

  def getMolLigandName(self, mol):
    molName = mol.getUniqueName()
    for ligName in self.ligandFileNames:
      if molName + '.pdbqt' in ligName:
        return ligName

  def parseDockedMolsDLG(self, fnDlg):
    molDic = {}
    i = 1
    with open(fnDlg) as fRes:
      for line in fRes:
        if line.startswith('DOCKED: MODEL'):
          posId = line.split()[-1]
          molDic[posId] = {'pdb': ''}
        elif line.startswith('DOCKED: USER    Estimated Free Energy'):
          molDic[posId]['energy'] = line.split()[8]
        elif line.startswith('DOCKED: USER    Estimated Inhibition'):
          molDic[posId]['ki'] = line.split()[7]
        elif line.startswith('DOCKED: REMARK') or line.startswith('TER'):
          molDic[posId]['pdb'] += line[8:]
        elif line.startswith('DOCKED: ATOM'):
          molDic[posId]['pdb'] += line[8:]
          i += 1
    return molDic

  def _citations(self):
    return ['Morris2009']

  def commentFirstLine(self, fn):
    with open(fn) as f:
      all = '# ' + f.read()
    with open(fn, 'w') as f:
      f.write(all)
    return fn

  def getDockedLigandsFiles(self, outDir):
    return list(glob.glob(os.path.join(outDir, '*.dlg')))
