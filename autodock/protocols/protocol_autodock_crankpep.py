# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
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
import re
import shutil, os
import zipfile

from pwchem.utils import pdbqt2other
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from autodock import Plugin
from autodock.objects import GridADT


class ProtCrankPep(EMProtocol):
    """Perform a docking experiment with AutoDock-CrankPep https://github.com/ccsb-scripps/ADCP"""
    _label = 'AutoDock-CrankPep docking'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputGrids', params.PointerParam, pointerClass="SetOfGridADT",
                      label='Target grid file:', allowsNull=False,
                      help='Grid prepared with AGFR describing the receptor.')


        conformers = form.addGroup("Parameters")
        conformers.addParam('nRuns', params.IntParam, label='Number of docking runs: ', default=50,
                            help='Number of independent runs using the stochastic procedure. \n'
                            'Different docking positions will be found for each of them.')
        conformers.addParam('cyc', params.BooleanParam, default=False,
                            label='Cyclic peptide through backbone: ',
                            help='Choose whether there are cyclic peptides through the backbone.')
        conformers.addParam('cys', params.BooleanParam, default=False,
                            label='Cyclic peptide through CYS-S-S-CYS: ',
                            help='Choose whether there are cyclic peptides through CYS-S-S-CYS.')
        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
        self._insertFunctionStep('createFilesStep')
        #self._insertFunctionStep('createOutputStep')

    def createFilesStep(self):
        targetFile = os.path.abspath(self.inputGrids.get().getFirstItem().getFileName())
        for grid in self.inputGrids.get():
            protFile = (grid.getAttributeValue('_peptideFile'))
            protName = os.path.splitext(os.path.basename(protFile))[0]
            pdbFile = os.path.abspath(os.path.join(self._getExtraPath(protName + '.pdb')))
            pdbqt2other(self, protFile, pdbFile)

            seq = self.getSequenceFromPdb(pdbFile)

            args = [f'-t {targetFile} -s {seq} -N {self.nRuns.get()} -o {protName}_docking -ref {pdbFile}']

            if(self.cyc.get()):
                args.append(f'-cyc')
            if (self.cys.get()):
                args.append(f'-cys')

            resultsFolder = os.path.abspath(os.path.join(self._getPath(), protName))
            os.makedirs(resultsFolder, exist_ok=True)
            Plugin.runADCP(self, args, cwd=resultsFolder)

    def createOutputStep(self):
        recFile = os.path.abspath(self.inputAtomStruct.get().getFileName())
        for prot in self.inputSmallMolecules.get():
            protFile = os.path.abspath(prot.getFileName())
            protName = os.path.splitext(os.path.basename(protFile))[0]
            logFile = os.path.abspath(self._getExtraPath(f'{protName}.log'))

            data = self.getInfo(logFile)

            fileName = os.path.join(self._getExtraPath(), f"{protName}.trg")
            grid = GridADT(fileName, proteinFile=recFile, spacing=data['spacing'], massCX=data['center'][0], massCY=data['center'][1], massCZ=data['center'][2], tool='AGFR')
            grid._peptideFile = protFile
            grid._XLength = data['length'][0]
            grid._YLength = data['length'][1]
            grid._ZLength = data['length'][2]
            grid._XSize = data['size'][0]
            grid._YSize = data['size'][1]
            grid._ZSize = data['size'][2]
            grid._numPockets = data['numPockets']

            self._defineOutputs(outputGrid=grid)


# --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        return validations

    def _warnings(self):
        warnings = []
        return warnings

# --------------------------- UTILS functions -----------------------------------
    def getInfo(self, logFile):
        if not os.path.exists(logFile):
            raise FileNotFoundError(f"AGFR log file not found: {logFile}")

        with open(logFile, 'r') as f:
            log = f.read()

        data = {}

        try:
            data['center'] = tuple(
                map(float, re.findall(r'Box center:\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)', log)[0]))
            data['length'] = tuple(
                map(float, re.findall(r'Box length:\s+([\d\.\-]+)\s+([\d\.\-]+)\s+([\d\.\-]+)', log)[0]))
            data['size'] = tuple(map(int, re.findall(r'Box size\s+: +(\d+)\s+(\d+)\s+(\d+)', log)[0]))
        except IndexError:
            data['center'] = data['length'] = data['size'] = (None, None, None)

        spacing_match = re.search(r'spacing\s+: +([\d\.]+)', log)
        data['spacing'] = float(spacing_match.group(1)) if spacing_match else None

        pocket_match = re.search(r'found\s+(\d+)\s+pocket\(s\)', log)
        data['numPockets'] = int(pocket_match.group(1)) if pocket_match else 0

        return data

    def getSequenceFromPdb(self, pdbFile):
        """
        Extracts peptide sequence from a PDB file.
        Uppercase for helix (H), lowercase for coil/other.
        """
        sequence = []
        helix_residues = set()
        with open(pdbFile, 'r') as f:
            for line in f:
                if line.startswith("HELIX"):
                    start_res = int(line[21:25].strip())
                    end_res = int(line[33:37].strip())
                    chain = line[19]
                    helix_residues.update((chain, i) for i in range(start_res, end_res+1))
        seq_dict = {}
        for line in open(pdbFile, 'r'):
            if line.startswith("ATOM") and line[13:15].strip() == "CA":
                resname = line[17:20].strip()
                chain = line[21]
                resnum = int(line[22:26].strip())
                aa = self.threetToOne(resname)
                if aa:
                    if (chain, resnum) in helix_residues:
                        aa = aa.upper()
                    else:
                        aa = aa.lower()
                    seq_dict[(chain, resnum)] = aa

        for key in sorted(seq_dict):
            sequence.append(seq_dict[key])

        return ''.join(sequence)


    def threetToOne(self, resname):
        """Convert 3-letter amino acid code to 1-letter."""
        aa_dict = {
            'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G',
            'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L', 'MET':'M', 'ASN':'N',
            'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V',
            'TRP':'W', 'TYR':'Y'
        }
        return aa_dict.get(resname.upper())
