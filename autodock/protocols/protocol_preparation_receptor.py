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
import os, json

from subprocess import run

from pyworkflow.protocol import params
from pwem.objects.data import AtomStruct

from pwchem import Plugin as pwchem_plugin
from pwchem.utils import cleanPDB
from pwchem.constants import MGL_DIC, RDKIT_DIC
from pwchem.protocols import ProtChemPrepareReceptor

from autodock import Plugin as autodockPlugin
from autodock.protocols.protocol_autodock import ProtChemAutodockBase

class ProtChemADTPrepare(ProtChemPrepareReceptor, ProtChemAutodockBase):
    """
User IA Manual: PreparationReceptor Protocol

The PreparationReceptor protocol is responsible for converting a macromolecular
structure into a format suitable for docking, specifically the PDBQT format
required by AutoDock-based tools. It ensures that the receptor structure
includes all necessary chemical and structural features while preserving
biological relevance for binding site analysis.

To begin, the user must provide a receptor structure, typically in PDB format.
This structure should represent the rigid part of the docking system, such as a
protein or nucleic acid. The protocol processes this input by first validating
its geometry and ensuring that all atoms are properly defined. Hydrogens are
added where necessary, and special attention is given to standardizing atom
names and resolving common formatting issues.

A key step in this process involves assigning atomic charges, typically
Gasteiger charges, and detecting atom types compatible with the AutoDock force
field. The protocol ensures that non-polar hydrogens are merged appropriately
and that any alternate conformations or heteroatoms unrelated to the binding
site are removed unless explicitly preserved by the user.

The user can control whether to include or exclude water molecules, metal ions,
or cofactors, depending on the nature of the system. Additionally, it is
possible to define whether the receptor should be treated as a rigid entity or
prepared for limited flexibility in downstream protocols that support flexible
residues.

Once the structure has been processed, it is exported as a PDBQT file that
captures all required information, including torsion constraints (if any),
partial charges, and docking-specific atom types. This output can be passed
directly to grid generation or docking protocols in Scipion-Chem.

In summary, this protocol prepares the receptor by transforming a standard
biomolecular structure into a chemically complete and docking-compatible format.
It provides essential preconditions for accurate ligand docking and ensures
smooth integration into AutoDock workflows.
"""
    def _defineParamsBasic(self, form, condition='True'):
        choicesRepair = ['None', 'Bonds', 'Hydrogens', 'Bonds hydrogens']
        if self.typeRL=="target":
            choicesRepair.append('Check hydrogens')
        preparation = form.addGroup('Preparation', condition=condition)
        preparation.addParam('repair', params.EnumParam, choices=choicesRepair,
                      default=0, label='Repair action:',
                      help='Bonds: build a single bond from each atom with no bonds to its closest neighbor\n'
                           'Hydrogens: add hydrogens\n'
                           'Bonds hydrogens: build bonds and add hydrogens\n'
                           'Check hydrogens: add hydrogens only if there are none already')
        preparation.addParam('preserveCharges', params.EnumParam, choices=['Add gasteiger charges', 'Preserve input charges',
                                                                    'Preserve charges of specific atoms'],
                      default=0, label='Charge handling')
        preparation.addParam('chargeAtoms', params.StringParam, default="", condition='preserveCharges==2',
                      label='Atoms to preserve charge', help='Separated by commas: Zn, Fe, ...')
        preparation.addParam('nphs', params.BooleanParam, default=True,
                      label='Merge charges and remove non-polar hydrogens')
        preparation.addParam('lps', params.BooleanParam, default=True,
                      label='Merge charges and remove lone pairs')
        preparation.addParam('waters', params.BooleanParam, default=True,
                      label='Remove water residues')
        if self.typeRL=="target":
            preparation.addParam('nonstdres', params.BooleanParam, default=True,
                          label='Remove chains composed entirely of non-standard residues')
            preparation.addParam('nonstd', params.BooleanParam, default=False,
                          label='Remove non-standard residues from all chains')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('preparationStep')
        self._insertFunctionStep('createOutputStep')

    def callPrepare(self, prog, args, outDir, popen=False):
        if self.repair.get()==3:
            args+=' -A bonds_hydrogens'
        elif self.repair.get()==1:
            args += ' -A bonds'
        elif self.repair.get()==2:
            args += ' -A hydrogens'
        elif self.repair.get()==4:
            args += ' -A checkhydrogens'

        if self.preserveCharges.get()==1:
            args+=" -C"
        elif self.preserveCharges.get()==2:
            for atom in self.chargeAtoms.get().split(','):
                args+=" -p %s" % atom.strip()

        cleanup = " -U "
        first = True
        if self.nphs.get():
            cleanup += " nphs"
            first = False
        if self.lps.get():
            if not first:
                cleanup += "_"
            cleanup += "lps"
            first = False
        if self.waters.get():
            if not first:
                cleanup += "_"
            cleanup += "waters"
            first = False
        if self.typeRL == "target" and self.nonstdres.get():
            if not first:
                cleanup += "_"
            cleanup += "nonstdres"
        if cleanup != "-U ":
            args += cleanup

        if self.typeRL == "target" and self.nonstd.get():
            args += " -e"

        if not popen:
            self.runJob(pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh '),
                        autodockPlugin.getADTPath(f'Utilities24/{prog}.py ') + args, cwd=outDir)
        else:
            fullProgram = pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh ') + \
                          autodockPlugin.getADTPath(f'Utilities24/{prog}.py ')
            run(fullProgram + args, cwd=outDir, shell=True)

    def createOutputStep(self):
        fnOut = self._getExtraPath('atomStruct.pdbqt')
        if os.path.exists(fnOut):
            target = AtomStruct(filename=fnOut)
            self._defineOutputs(outputStructure=target)
            self._defineSourceRelation(self.inputAtomStruct, target)


MEEKO, MGL = 0, 1

class ProtChemADTPrepareReceptor(ProtChemADTPrepare):
    """Prepare receptor using Autodocking Tools from MGL"""
    _label = 'target preparation ADT'
    _program = ""

    def _defineParams(self, form):
        self.typeRL="target"
        form.addSection(label='Input')
        form.addParam('inputAtomStruct', params.PointerParam, pointerClass="AtomStruct",
                      label='Atomic Structure:', allowsNull=False,
                      help='Input Atomic structure to prepare for Autodock docking')

        self.defineCleanParams(form, w=False)

        form.addParam('prepProg', params.EnumParam, choices=['Meeko', 'MGLTools'], label='Preparation program:',
                      default=MEEKO, help='Software to use for the receptor preparation')

        ProtChemADTPrepare._defineParamsBasic(self, form, condition=f'prepProg=={MGL}')

        form.addParam('doZnDock', params.BooleanParam, label='Perform Zn metalloprotein preparation: ', default=False,
                       expertLevel=params.LEVEL_ADVANCED, condition=f'prepProg=={MGL}',
                       help='Whether to use the scripts for preparing the metalloprotein receptor containing Zn')

    def preparationStep(self):
        #Clean PDB
        pdbIni = self.inputAtomStruct.get().getFileName()
        filename = os.path.splitext(os.path.basename(pdbIni))[0]
        fnPdb = self._getExtraPath('%s_clean.pdb' % filename)

        chainIds = None
        if self.rchains.get():
            chainJson = json.loads(self.chain_name.get())  # From wizard dictionary
            if 'chain' in chainJson:
                chainIds = [chainJson["chain"].upper().strip()]
            elif 'model-chain' in chainJson:
                modelChains = chainJson["model-chain"].upper().strip()
                chainIds = [x.split('-')[1] for x in modelChains.split(',')]

        het2keep = self.het2keep.get().split(', ')
        inFile = self.getInpFilePDB()
        cleanedPDB = cleanPDB(inFile, fnPdb,False, self.HETATM.get(), chainIds, het2keep)

        fnOut = self.getReceptorPDBQT()
        if self.prepProg.get() == MGL:
            args = ' -r %s -o %s' % (os.path.abspath(cleanedPDB), fnOut)
            ProtChemADTPrepare.callPrepare(self, "prepare_receptor4", args, outDir=self._getExtraPath())

            if self.doZnDock.get():
                zincPrepPath = autodockPlugin.getVinaScriptsPath('zinc_pseudo.py')

                auxOut = fnOut.replace('.pdbqt', '_tz.pdbqt')
                args = ' -r {} -o {}'.format(fnOut, auxOut)
                fullProgram = '%s && %s %s' % (pwchem_plugin.getEnvActivationCommand(RDKIT_DIC), 'python', zincPrepPath)
                self.runJob(fullProgram, args, cwd=self._getExtraPath())
        else:
            outBase = os.path.splitext(fnOut)[0]
            args = f' -i {os.path.abspath(cleanedPDB)} -o {outBase} -a -p'
            autodockPlugin.runMeekoReceptor(self, args)

    def createOutputStep(self):
        fnOut = self.getReceptorPDBQT()
        fnOut = fnOut if not self.doZnDock.get() else fnOut.replace('.pdbqt', '_tz.pdbqt')
        if os.path.exists(fnOut):
            target = AtomStruct(filename=fnOut)
            self._defineOutputs(outputStructure=target)
            self._defineSourceRelation(self.inputAtomStruct, target)

    def _validate(self):
        errors = []
        if self.rchains.get() and not self.chain_name.get():
            errors.append('You must specify the chains to be maintained')
        return errors

    def getInpFilePDB(self):
        inFile = self.inputAtomStruct.get().getFileName()
        base, ext = os.path.splitext(inFile)

        pdbFile = f"{base}.pdb"
        if os.path.exists(pdbFile):
            return os.path.abspath(pdbFile)

        return os.path.abspath(self._getPath(f'{getBaseName(inFile)}{ext}'))