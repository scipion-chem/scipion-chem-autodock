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
import os

import pyworkflow
from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pwchem.utils import pdbqt2other
from pwem.protocols import EMProtocol
from pyworkflow.object import Float
from pyworkflow.protocol import params

from autodock import Plugin


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
        self._insertFunctionStep('createOutputStep')

    def createFilesStep(self):
        targetFile = os.path.abspath(self.inputGrids.get().getFirstItem().getFileName())
        for grid in self.inputGrids.get():
            protFile = grid.getAttributeValue('_peptideFile')
            protName = self.getBaseName(protFile)
            pdbFile = os.path.abspath(os.path.join(self._getExtraPath(f"{protName}.pdb")))

            pdbqt2other(self, protFile, pdbFile)

            seq = self.getSequenceFromPdb(pdbFile)
            args = [f'-t {targetFile} -s {seq} -N {self.nRuns.get()} -o {protName}_docking -ref {pdbFile}']

            if(self.cyc.get()):
                args.append('-cyc')
            if (self.cys.get()):
                args.append('-cys')

            resultsFolder = self.getResultsFolder(protName)
            Plugin.runADCP(self, args, cwd=resultsFolder)

    def createOutputStep(self):
        logFile = os.path.abspath(self._getPath('logs/run.stdout'))
        outputLogData = self.readOutputData(logFile)
        i = 0
        outputMols = SetOfSmallMolecules().create(outputPath=self._getPath())

        for grid in self.inputGrids.get():
            protFile = grid.getAttributeValue('_peptideFile')
            recFile = grid.getAttributeValue('_proteinFile')
            protName = self.getBaseName(protFile)
            resultsFolder = self.getResultsFolder(protName)
            rankedFiles = self.rankedFiles(resultsFolder)
            data = outputLogData[i]
            mappingFile = os.path.abspath(grid.getFileName())

            for file in rankedFiles:
                m = re.search(r'_ranked_(\d+)', file)
                if not m:
                    continue

                modeNum = int(m.group(1))
                if modeNum not in data:
                    continue

                affinity = data[modeNum]["affinity"]
                energy = data[modeNum]["energy"]
                bestRun = data[modeNum]["bestRun"]
                poseFile = os.path.abspath(os.path.join(resultsFolder, file))

                newMol = self.createDockedMolecule(
                    protFile, recFile, protName,
                    poseFile, mappingFile,
                    modeNum, energy, affinity, bestRun
                )

                outputMols.append(newMol)

            i += 1

        outputMols.setDocked(True)
        recFile = grid.getAttributeValue('_proteinFile')
        outputMols.setProteinFile(recFile)
        self._defineOutputs(outputSmallMolecules=outputMols)


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

        spacingMatch = re.search(r'spacing\s+: +([\d\.]+)', log)
        data['spacing'] = float(spacingMatch.group(1)) if spacingMatch else None

        pocketMatch = re.search(r'found\s+(\d+)\s+pocket\(s\)', log)
        data['numPockets'] = int(pocketMatch.group(1)) if pocketMatch else 0

        return data

    def getSequenceFromPdb(self, pdbFile):
        """
        Extracts peptide sequence from a PDB file.
        Uppercase for helix (H), lowercase for coil/other.
        """
        helixResidues = self.extractHelixResidues(pdbFile)
        seqDict = self.extractResidues(pdbFile, helixResidues)

        sequence = [seqDict[key] for key in sorted(seqDict)]
        return ''.join(sequence)

    def extractHelixResidues(self, pdbFile):
        helixResidues = set()
        with open(pdbFile, 'r') as f:
            for line in f:
                if line.startswith("HELIX"):
                    startRes = int(line[21:25].strip())
                    endRes = int(line[33:37].strip())
                    chain = line[19]
                    helixResidues.update((chain, i) for i in range(startRes, endRes + 1))
        return helixResidues

    def extractResidues(self, pdbFile, helixResidues):
        seqDict = {}
        with open(pdbFile, 'r') as f:
            for line in f:
                if not (line.startswith("ATOM") and line[13:15].strip() == "CA"):
                    continue

                resname = line[17:20].strip()
                chain = line[21]
                resnum = int(line[22:26].strip())
                aa = self.threetToOne(resname)
                if not aa:
                    continue

                # Helix uppercase, coil lowercase
                aa = aa.upper() if (chain, resnum) in helixResidues else aa.lower()
                seqDict[(chain, resnum)] = aa

        return seqDict


    def threetToOne(self, resname):
        """Convert 3-letter amino acid code to 1-letter."""
        aaDict = {
            'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G',
            'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L', 'MET':'M', 'ASN':'N',
            'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V',
            'TRP':'W', 'TYR':'Y'
        }
        return aaDict.get(resname.upper())

    def rankedFiles(self, directory):
        allFiles = os.listdir(directory)
        rankedFiles = [f for f in allFiles if '_ranked_' in f]
        return rankedFiles

    def readOutputData(self, logFile):
        allRuns = []
        dockingData = {}
        tableStarted = False

        with open(logFile, 'r') as f:
            for line in f:
                line = line.strip()
                if self.isStartSearch(line):
                    if dockingData:
                        allRuns.append(dockingData)
                        dockingData = {}
                    tableStarted = False
                    continue

                if self.isTableHeader(line):
                    tableStarted = True
                    continue

                if self.isEndSearch(line):
                    if dockingData:
                        allRuns.append(dockingData)
                        dockingData = {}
                    tableStarted = False
                    continue

                if not tableStarted or not self.isTableLine(line):
                    continue

                parsed = self.parseTableLine(line)
                if parsed:
                    mode, affinity, energy, bestRun = parsed
                    dockingData[mode] = {"affinity": affinity, "energy": energy, "bestRun": bestRun}

        if dockingData:
                allRuns.append(dockingData)

        return allRuns

    def isStartSearch(self, line):
        return line.startswith("Performing search")

    def isTableHeader(self, line):
        return line.startswith("mode |  affinity")

    def isTableLine(self, line):
        return line and not line.startswith("|") and not line.startswith("-----")

    def isEndSearch(self, line):
        return line.startswith("clean up")

    def parseTableLine(self, line):
        parts = line.split()
        if len(parts) >= 9 and parts[0].isdigit():
            return int(parts[0]), float(parts[1]), float(parts[7]), int(parts[8])
        return None

    def getBaseName(self, filePath):
        return os.path.splitext(os.path.basename(filePath))[0]

    def getResultsFolder(self, name):
        folder = os.path.abspath(os.path.join(self._getPath(), name))
        os.makedirs(folder, exist_ok=True)
        return folder

    def createDockedMolecule(self, protFile, recFile, protName, poseFile, mappingFile,
                              modeNum, energy, affinity, bestRun):

        mol = SmallMolecule(
            smallMolFilename=protFile,
            proteinFile=recFile,
            molName=protName,
            type='Autodock-CrankPep'
        )

        mol.setPoseFile(poseFile)
        mol.setMappingFile(mappingFile)
        mol.setConfId(modeNum)
        mol.setGridId(1)
        mol.setPoseId(modeNum)
        mol.setDockId(bestRun)
        mol.setEnergy(energy)
        mol._ligandAffinity = Float()
        mol.setAttributeValue('_ligandAffinity', affinity)

        return mol
