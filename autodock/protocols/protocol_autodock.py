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
  LEVEL_ADVANCED, USE_GPU, GPU_LIST, StringParam, EnumParam, LabelParam, TextParam
import pyworkflow.object as pwobj
from pyworkflow.utils.path import makePath, createLink

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import runOpenBabel, generate_gpf, calculate_centerMass, getBaseFileName
from pwchem import Plugin as pwchem_plugin
from pwchem.constants import MGL_DIC

from autodock import Plugin as autodock_plugin


LGA, GA, LS, SA = 0, 1, 2, 3
searchDic = {LGA: 'Lamarckian Genetic Algorithm', GA: 'Genetic Algorithm', LS: 'Local Search',
             SA: 'Simulated annealing'}

class ProtChemAutodockBase(EMProtocol):
    """Base class protocol for AutoDock docking protocols"""

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_PARALLEL

    def _defineFlexParams(self, form):
        group = form.addGroup('Flexible residues')
        group.addParam('doFlexRes', BooleanParam, label='Add flexible residues: ', default=False,
                       help='Whether to add residues of the receptor which will be treated as flexible by AutoDock')
        group.addParam('flexChain', StringParam, label='Residue chain: ', condition='doFlexRes',
                       help='Specify the protein chain')
        group.addParam('flexPosition', StringParam, label='Flexible residues: ', condition='doFlexRes',
                       help='Specify the residues to make flexible')
        group.addParam('addFlex', LabelParam, label='Add defined residues', condition='doFlexRes',
                       help='Here you can define flexible residues which will be added to the list below.'
                            'Be aware that it will be a range from the first to the last you choose')
        group.addParam('flexList', TextParam, width=70, default='', label='List of flexible residues: ',
                       condition='doFlexRes', help='List of chain | residues to make flexible. \n')

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Receptor specification')
        group.addParam('fromReceptor', EnumParam, label='Dock on : ', default=1,
                       choices=['Whole protein', 'SetOfStructROIs'],
                       help='Whether to dock on a whole protein surface or on specific regions')

        # Docking on whole protein
        group.addParam('inputAtomStruct', PointerParam, pointerClass="AtomStruct",
                       label='Input atomic structure: ', condition='fromReceptor == 0',
                       help="The atom structure to use as receptor in the docking")
        group.addParam('radius', FloatParam, label='Grid radius for whole protein: ',
                       condition='fromReceptor == 0', allowsNull=False,
                       help='Radius of the Autodock grid for the whole protein')

        # Docking on pockets
        group.addParam('inputStructROIs', PointerParam, pointerClass="SetOfStructROIs",
                       label='Input pockets: ', condition='not fromReceptor == 0',
                       help="The protein structural ROIs to dock in")
        group.addParam('pocketRadiusN', FloatParam, label='Grid radius vs StructROI radius: ',
                       condition='not fromReceptor == 0', default=1.1, allowsNull=False,
                       help='The radius * n of each StructROI will be used as grid radius')

        group.addParam('spacing', FloatParam, label='Spacing of the grids: ', default=0.5, allowsNull=False,
                       help='Spacing of the generated Autodock grids')

        self._defineFlexParams(form)

        group = form.addGroup('Docking')
        group.addParam('inputSmallMolecules', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Input small molecules: ', allowsNull=False,
                       help="Input small molecules to be docked with AutoDock")
        group.addParam('nRuns', IntParam, label='Number of docking runs: ', default=50,
                       help='Number of independent runs using the selected strategy. \nDifferent docking positions will be '
                            'found for each of them.')
        group.addParam('rmsTol', FloatParam, label='Cluster tolerance (A): ', default=2.0, expertLevel=LEVEL_ADVANCED,
                       help='Maximum RMSD for 2 docked structures to be in the same cluster')

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

    def getOriginalReceptorFile(self):
        if self.fromReceptor.get() == 0:
            return self.inputAtomStruct.get().getFileName()
        else:
            return self.inputStructROIs.get().getProteinFile()

    def buildFlexReceptor(self, receptorFn):
        molName = getBaseFileName(receptorFn)
        allFlexRes = self.parseFlexRes(molName)
        flexFn, rigFn = self.getFlexFiles()
        args = ' -r {} -s {} -g {} -x {}'.format(receptorFn, allFlexRes, rigFn, flexFn)
        self.runJob(pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh'),
                    autodock_plugin.getADTPath('Utilities24/prepare_flexreceptor4.py') + args, cwd=self._getExtraPath())
        return flexFn, rigFn

    def getFlexFiles(self):
        return os.path.abspath(self._getExtraPath('receptor_flex.pdbqt')), \
               os.path.abspath(self._getExtraPath('receptor_rig.pdbqt'))

    def parseFlexRes(self, molName):
        allFlexDic = {}
        for line in self.flexList.get().split('\n'):
            if line.strip():
                chain, res = line.strip().split(':')
                if chain in allFlexDic:
                    allFlexDic[chain] += res + '_'
                else:
                    allFlexDic[chain] = res + '_'

        allFlexStr = ''
        for chain in allFlexDic:
            allFlexStr += '{}:{}:{},'.format(molName, chain, allFlexDic[chain][:-1])
        return allFlexStr[:-1]

    def getOutputPocketDir(self, pocket=None):
        if pocket == None:
          outDir = os.path.abspath(self._getExtraPath('pocket_1'))
        else:
          outDir = os.path.abspath(self._getExtraPath('pocket_{}'.format(pocket.getObjId())))
        return outDir

    def getPocketDirs(self):
        dirs = []
        for file in os.listdir(self._getExtraPath()):
            d = self._getExtraPath(file)
            if os.path.isdir(d) and 'pocket' in file:
                dirs.append(d)
        dirs.sort()
        return dirs

    def getReceptorName(self):
        if self.fromReceptor.get() == 0:
            atomStructFn = self.getOriginalReceptorFile()
            return atomStructFn.split('/')[-1].split('.')[0]
        else:
            return self.inputStructROIs.get().getProteinName()

class ProtChemAutodock(ProtChemAutodockBase):
  """Perform a docking experiment with autodock. Grid must be generated in this protocol in order to take into
       account ligands atom types. See the help at
       http://autodock.scripps.edu/faqs-help/manual/autodock-4-2-user-guide/AutoDock4.2_UserGuide.pdf"""
  _label = 'Autodock docking'
  _program = ""


  def _defineParams(self, form):
    form.addHidden(USE_GPU, BooleanParam, default=True,
                   label="Use GPU for execution: ",
                   help="This protocol has both CPU and GPU implementation.\
                                         Select the one you want to use.")

    form.addHidden(GPU_LIST, StringParam, default='0', label="Choose GPU IDs",
                   help="Add a list of GPU devices that can be used")

    super()._defineParams(form)

    form.addSection(label="Search")
    form.addParam('searchType', EnumParam, label='Search type: ', default=LGA,
                  choices=list(searchDic.values()),
                  help='Type of search to perform in the docking. It defines how the positions of the molecules will '
                       'be searched.\n For more details check the page 55 of the AutoDock4 manual '
                       'https://autodock.scripps.edu/wp-content/uploads/sites/56/2022/04/AutoDock4.2.6_UserGuide.pdf')

    group = form.addGroup('Search parameters', condition='searchType in [{}, {}, {}]'.format(LGA, GA, LS))
    group.addParam('gaPop', IntParam, label='Population size', default=150,
                   help='This is the number of individuals in the population. Each individual is a coupling of a '
                        'genotype and its associated phenotype')
    group.addParam('gaNumEvals', IntParam, label='Number of evaluations: ', default=2500000,
                   help='Set the maximum number of energy evaluations performed during each GA, LGA, or LS run')

    group = form.addGroup('Genetic algorithm', condition='searchType in [{}, {}]'.format(LGA, GA))
    group.addParam('gaNumGens', IntParam, label='Number of generations: ', default=27000,
                  help='This is the maximum number of generations simulated during each GA or LGA run')
    group.addParam('gaElitism', IntParam, label='Elitism: ', default=1,
                  help='This is the number of top individuals that are guaranteed to survive into the next generation.')
    group.addParam('gaMutationRate', FloatParam, label='Mutation rate: ', default=0.02, expertLevel=LEVEL_ADVANCED,
                  help='This is a floating point number from 0 to 1, representing the probability that a particular '
                       'gene is mutated')
    group.addParam('gaCrossOverRate', FloatParam, label='Crossover rate: ', default=0.8, expertLevel=LEVEL_ADVANCED,
                  help='This is a floating point number from 0 to 1 denoting the crossover rate. Crossover rate '
                       'is the expected number of pairs in the population that will exchange genetic material')
    group.addParam('gaWindowSize', IntParam, label='Window size: ', default=10, expertLevel=LEVEL_ADVANCED,
                  help='This is the number of preceding generations to take into consideration when deciding the '
                       'threshold for the worst individual in the current population')

    group = form.addGroup("Local search", condition='searchType in [{}, {}]'.format(LGA, LS))
    group.addParam('lsFreq', FloatParam, label='Local search frequency: ', default=0.06,
                   condition='searchType == {}'.format(LGA),
                   help='This is the probability of any particular phenotype being subjected to local search')
    group.addParam('localType', EnumParam, label='Local search method: ',
                   choices=['Classical Solis and Wets', 'pseudo-Solis and Wets'], default=1,
                   help='Whether to use the classical Solis and Wets local searcher, using the method of uniform '
                        'variances for changes in translations, orientations, and torsions; or the pseudo-Solis and '
                        'Wets local searcher. This method maintains the relative proportions of variances for the '
                        'translations in Å and the rotations in radians')
    group.addParam('swMaxIts', IntParam, label='Number of iterations: ', default=300,
                   help='This is the maximum number of iterations that the local search procedure applies to the '
                        'phenotype of any given individual, per generation')
    group.addParam('swMaxSucc', IntParam, label='Successes in a raw: ', default=4, expertLevel=LEVEL_ADVANCED,
                   help='Number of successes before changing rho in Solis & Wets algorithms')
    group.addParam('swMaxFail', IntParam, label='Failures in a raw: ', default=4, expertLevel=LEVEL_ADVANCED,
                   help='Number of failures before changing rho in Solis & Wets algorithms.')
    group.addParam('swRho', FloatParam, label='Initial variance: ', default=1.0, expertLevel=LEVEL_ADVANCED,
                   help='It defines the size of the local search, and specifies the size of the local space to sample')
    group.addParam('swLbRho', FloatParam, label='Variance lower bound: ', default=0.01, expertLevel=LEVEL_ADVANCED,
                   help='This is the lower bound on rho, the variance for making changes to genes (i.e. translations, '
                        'orientation and torsions)')

    group = form.addGroup("Simulated annealing", condition='searchType == {}'.format(SA))
    line = group.addLine('Ligand initial state: ', expertLevel=LEVEL_ADVANCED,
                         help='[e0max] This keyword stipulates that the ligand’s initial state cannot have an energy '
                              'greater than the first value, nor can there be more than the second value’s number of '
                              'retries.')
    line.addParam('e0max1', FloatParam, label='Maximum initial energy: ', default=0.00)
    line.addParam('e0max2', IntParam, label='Maximum number of retries: ', default=10000)

    line = group.addLine('Ligand step sizes: ', expertLevel=LEVEL_ADVANCED,
                         help='[(t/q/d)step] Defines the maximum translation (A) / angular (º) / dihedral (º) jump '
                              'for the first cycle that the ligand may make in one simulated annealing step. When '
                              '“trnrf” is less than 1, the reduction factor is multiplied with the tstep at the end of '
                              'each cycle, to give the new value for the next cycle.')
    line.addParam('tstep', FloatParam, label='Translation: ', default=2.00)
    line.addParam('qstep', FloatParam, label='Angular: ', default=50.0)
    line.addParam('dstep', FloatParam, label='Dihedral: ', default=50.0)

    group.addParam('rt0_R', FloatParam, label='Initial annealing temperature: ', default=252,
                   help='[rt0/R] Initial absolute temperature. It will be multiplied by the gas constant '
                        'R (1.987 cal mol -1 K -1)for the input docking.')
    group.addParam('SA_schedule', EnumParam, label='SA schedule: ', default=0, choices=['Linear', 'Geometric'],
                   help='The default “linear_schedule” uses a linear or arithmetic temperature reduction schedule'
                        ' during Monte Carlo simulated annealing. The “geometric_schedule” keyword uses instead a '
                        'geometric reduction schedule, according to the rtrf parameter described next. At the end of '
                        'each cycle, the temperature is reduced by (rt0/cycles).')
    group.addParam('rtrf', FloatParam, label='Temperature reduction factor: ', default=0.90,
                   condition='SA_schedule == 1',
                   help='Annealing temperature reduction factor, g. At the end of each cycle, the '
                        'annealing temperature is multiplied by this factor, to give that of the next cycle.')
    group.addParam('cycles', IntParam, label='Number of temperature reduction cycles.: ', default=50,
                   help='Number of temperature reduction cycles.')

    group.addParam('accs', IntParam, label='Maximum number of accepted steps per cycle: ', default=30000,
                   help='Maximum number of accepted steps per cycle.', expertLevel=LEVEL_ADVANCED)
    group.addParam('rejs', IntParam, label='Maximum number of rejected steps per cycle: ', default=30000,
                   help='Maximum number of rejected steps per cycle', expertLevel=LEVEL_ADVANCED)

    group.addParam('select', EnumParam, label='Select state to following cycle: ',
                   choices=['Minimum', 'Last'], default=0, expertLevel=LEVEL_ADVANCED,
                   help='State selection flag. This character can be either m for the minimum state, or l for the last'
                        ' state found during each cycle, to begin the following cycle.')
    line = group.addLine('Per-cycle reduction factors: ', expertLevel=LEVEL_ADVANCED,
                         help='Per-cycle reduction factor for translation / orientation / torsinal dihedralsteps.')
    line.addParam('trnrf', FloatParam, label='Translation: ', default=1.0)
    line.addParam('quarf', FloatParam, label='Orientation: ', default=1.0)
    line.addParam('dihrf', FloatParam, label='Torsional dihedral: ', default=1.0)

    form.addParallelSection(threads=4, mpi=1)

  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
    cId = self._insertFunctionStep('convertStep', prerequisites=[])

    dockSteps = []
    if self.fromReceptor.get() == 0:
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
    if self.fromReceptor.get() == 0:
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

    if self.doFlexRes:
      flexFn, fnReceptor = self.buildFlexReceptor(fnReceptor)

    npts = (radius * 2) / self.spacing.get()
    gpf_file = generate_gpf(fnReceptor, spacing=self.spacing.get(),
                            xc=x_center, yc=y_center, zc=z_center,
                            npts=npts, outDir=outDir, ligandFns=self.ligandFileNames)

    args = "-p {} -l {}.glg".format(gpf_file, self.getReceptorName())
    self.runJob(autodock_plugin.getAutodockPath("autogrid4"), args, cwd=outDir)

  def dockStep(self, mol, pocket=None):
    # Prepare grid
    if self.doFlexRes:
        flexReceptorFn, receptorFn = self.getFlexFiles()
    else:
        flexReceptorFn, receptorFn = None, self.receptorFile
    outDir = self.getOutputPocketDir(pocket)

    molFn = self.getMolLigandName(mol)
    dpfFile = self.writeDPF(outDir, molFn, receptorFn, flexReceptorFn)

    fnDLG = dpfFile.replace('.dpf', '.dlg')
    if not getattr(self, USE_GPU):
      args = "-p %s -l %s" % (dpfFile, fnDLG)
      self.runJob(autodock_plugin.getAutodockPath("autodock4"), args, cwd=outDir)
    else:
      fnDPF = self.commentFirstLine(dpfFile)
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
          pdbqtFile = self._getPath('g{}_{}_{}.pdbqt'.format(gridId, molName, modelId))
          with open(pdbqtFile, 'w') as f:
            f.write(pocketDic[molName][modelId]['pdb'])
          pocketDic[molName][modelId]['file'] = pdbqtFile

      for smallMol in self.inputSmallMolecules.get():
        molFile = smallMol.getFileName()
        molName = getBaseFileName(molFile)
        if molName in pocketDic:
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

    # todo: manage several receptor conformations when flexible docking
    outputSet.proteinFile.set(self.getOriginalReceptorFile())
    outputSet.setDocked(True)
    self._defineOutputs(outputSmallMolecules=outputSet)
    self._defineSourceRelation(self.inputSmallMolecules, outputSet)

  ########################### Utils functions ############################

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

  def getReceptorDir(self):
    atomStructFn = self.getOriginalReceptorFile()
    return '/'.join(atomStructFn.split('/')[:-1])

  def writeDPF(self, outDir, molFn, receptorFn, flexFn=None):
      molName, _ = os.path.splitext(os.path.basename(molFn))
      createLink(molFn, os.path.join(outDir, molName + '.pdbqt'))
      makePath(outDir)

      fnDPF = os.path.abspath(os.path.join(outDir, molName + ".dpf"))
      baseDPF = os.path.abspath(os.path.join(outDir, molName + "_base.dpf"))
      # General parameters
      args = " -l %s -r %s -o %s" % (molFn, receptorFn, baseDPF)
      if flexFn:
          args += ' -x ' + flexFn

      self.runJob(pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh'),
                  autodock_plugin.getADTPath('Utilities24/prepare_dpf42.py') + args, cwd=outDir)

      with open(baseDPF) as f:
          myDPFstr = '\n'.join(f.read().split('\n')[:19])

      myDPFstr += '\n\n# Search parameters\n'
      if self.searchType.get() in [LGA, GA, LS]:
        myDPFstr += "ga_pop_size %d\n" % self.gaPop.get()
        myDPFstr += "ga_num_evals %d\n" % self.gaNumEvals.get()

      if self.searchType.get() in [LGA, GA]:
        myDPFstr += '\n# Genetic algorithm parameters\n'
        myDPFstr += "ga_num_generations %d\n" % self.gaNumGens.get()
        myDPFstr += "ga_elitism %d\n" % self.gaElitism.get()
        myDPFstr += "ga_mutation_rate %f\n" % self.gaMutationRate.get()
        myDPFstr += "ga_crossover_rate %f\n" % self.gaCrossOverRate.get()
        myDPFstr += "ga_window_size %d\n" % self.gaWindowSize.get()
        myDPFstr += "set_ga\n"

      if self.searchType.get() in [LGA, LS]:
        myDPFstr += '\n# Local search parameters\n'
        myDPFstr += "sw_max_its %d\n" % self.swMaxIts.get()
        myDPFstr += "sw_max_succ %d\n" % self.swMaxSucc.get()
        myDPFstr += "sw_max_fail %d\n" % self.swMaxFail.get()
        myDPFstr += "sw_rho %f\n" % self.swRho.get()
        myDPFstr += "sw_lb_rho %d\n" % self.swLbRho.get()
        myDPFstr += "ls_search_freq %f\n" % self.lsFreq.get()
        myDPFstr += 'set_psw1\n' if self.localType.get() == 1 else 'set_sw1\n'

      if self.searchType.get() == SA:
        myDPFstr += '\n# Simulated annealing parameters\n'
        myDPFstr += "e0max %f %d\n" % (self.e0max1.get(), self.e0max2.get())
        myDPFstr += "tstep %f\n" % self.tstep.get()
        myDPFstr += "qstep %f\n" % self.qstep.get()
        myDPFstr += "dstep %f\n" % self.dstep.get()
        myDPFstr += "rt0 %f\n" % (self.rt0_R.get() * 1.987)
        myDPFstr += 'linear_schedule\n' if self.SA_schedule.get() == 0 else 'geometric_schedule\n'
        if self.SA_schedule.get() == 1:
            myDPFstr += "rtrf %f\n" % self.rtrf.get()
        myDPFstr += "cycles %d\n" % self.cycles.get()
        myDPFstr += "accs %f\n" % self.accs.get()
        myDPFstr += "rejs %f\n" % self.rejs.get()
        myDPFstr += "select %s\n" % ('m' if self.select.get() == 0 else 'l')
        myDPFstr += "trnrf %f\n" % self.trnrf.get()
        myDPFstr += "quarf %f\n" % self.quarf.get()
        myDPFstr += "dihrf %f\n" % self.dihrf.get()

        myDPFstr += "simanneal %d\n" % self.nRuns.get()

      if self.searchType.get() == LGA:
        myDPFstr += "ga_run %d\n" % self.nRuns.get()
      elif self.searchType.get() == GA:
        myDPFstr += "do_global_only %d\n" % self.nRuns.get()
      elif self.searchType.get() == LS:
        myDPFstr += "do_local_only %d\n" % self.nRuns.get()

      myDPFstr += '\n# Analysis parameters\n'
      myDPFstr += "rmstol %f\nanalysis" % self.rmsTol.get()

      with open(fnDPF, 'w') as f:
          f.write(myDPFstr)
      return fnDPF


  def getLigandsFileNames(self):
    ligFns = []
    for mol in self.inputSmallMolecules.get():
      ligFns.append(mol.getFileName())
    return ligFns


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
          molDic[posId]['energy'] = line.split('=')[1].split('kcal/mol')[0]
        elif line.startswith('DOCKED: USER    Estimated Inhibition'):
          molDic[posId]['ki'] = line.split()[7]
        elif line.startswith('DOCKED: REMARK') or line.startswith('TER'):
          molDic[posId]['pdb'] += line[8:]
        elif line.startswith('DOCKED: ATOM'):
          molDic[posId]['pdb'] += line[8:]
          i += 1
    return molDic


  def _validate(self):
      vals = []
      if getattr(self, USE_GPU) and self.searchType.get() != LGA:
          vals.append('The use of a search type other than LGA for the GPU version is not possible. Either use the CPU '
                      'version or another search type.')
      return vals

  def commentFirstLine(self, fn):
    with open(fn) as f:
      all = '# ' + f.read()
    with open(fn, 'w') as f:
      f.write(all)
    return fn

  def getDockedLigandsFiles(self, outDir):
    return list(glob.glob(os.path.join(outDir, '*.dlg')))
