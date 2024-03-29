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

from autodock import Plugin as autodock_plugin
from autodock.protocols.protocol_autodock import ProtChemAutodockBase

class ProtChemADTPrepare(ProtChemPrepareReceptor, ProtChemAutodockBase):
    def _defineParamsBasic(self, form):
        choicesRepair = ['None', 'Bonds', 'Hydrogens', 'Bonds hydrogens']
        if self.typeRL=="target":
            choicesRepair.append('Check hydrogens')
        preparation = form.addGroup('Preparation')
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
        self._insertFunctionStep('createOutput')

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
            self.runJob(pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh'),
                        autodock_plugin.getADTPath('Utilities24/%s.py'%prog)+args, cwd=outDir)
        else:
            fullProgram = pwchem_plugin.getProgramHome(MGL_DIC, 'bin/pythonsh ') + \
                          autodock_plugin.getADTPath('Utilities24/%s.py' % prog)
            run(fullProgram + args, cwd=outDir, shell=True)

    def createOutput(self):
        fnOut = self._getExtraPath('atomStruct.pdbqt')
        if os.path.exists(fnOut):
            target = AtomStruct(filename=fnOut)
            self._defineOutputs(outputStructure=target)
            self._defineSourceRelation(self.inputAtomStruct, target)

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
        ProtChemADTPrepare._defineParamsBasic(self, form)

        form.addParam('doZnDock', params.BooleanParam, label='Perform Zn metalloprotein preparation: ', default=False,
                       expertLevel=params.LEVEL_ADVANCED,
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
        cleanedPDB = cleanPDB(self.inputAtomStruct.get().getFileName(), fnPdb,
                               False, self.HETATM.get(), chainIds, het2keep)

        fnOut = self.getReceptorPDBQT()
        args = ' -r %s -o %s' % (os.path.abspath(cleanedPDB), fnOut)
        ProtChemADTPrepare.callPrepare(self, "prepare_receptor4", args, outDir=self._getExtraPath())

        if self.doZnDock.get():
            zincPrepPath = autodock_plugin.getVinaScriptsPath('zinc_pseudo.py')

            auxOut = fnOut.replace('.pdbqt', '_tz.pdbqt')
            args = ' -r {} -o {}'.format(fnOut, auxOut)
            fullProgram = '%s && %s %s' % (pwchem_plugin.getEnvActivationCommand(RDKIT_DIC), 'python', zincPrepPath)
            self.runJob(fullProgram, args, cwd=self._getExtraPath())

    def createOutput(self):
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
