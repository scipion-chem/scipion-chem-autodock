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

from pwchem import MGL_DIC
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwem.convert import cifToPdb
from pwem.protocols import EMProtocol
from pyworkflow.object import String, Float, Integer
from pyworkflow.protocol import params

from autodock import Plugin

INPUT_TYPE = ['SetOfSmallMolecules', 'SetOfAtomStructs']


class ProtCrankPep(EMProtocol):
    """Perform a docking experiment with AutoDock-CrankPep https://github.com/ccsb-scripps/ADCP"""
    _label = 'AutoDock-CrankPep docking'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputAtomStruct', params.PointerParam, pointerClass="AtomStruct",
                      label='Receptor protein:', allowsNull=False,
                      help='It must be in pdbqt format.')

        form.addParam('inputPeptides', params.PointerParam, pointerClass="SetOfAtomStructs,SetOfSmallMolecules",
                      label='Input peptides:', allowsNull=True,
                      help='It must be in pdb or mol2 format, you may use Schrodinger convert to change it.')

        conformers = form.addGroup("Grid")
        conformers.addParam('padding', params.FloatParam, default=4.0,
                            label='Padding: ',
                            help='Amount of padding added to each side of the box.')
        conformers.addParam('flexRes', params.BooleanParam, default=False,
                            label='Flexible residues: ',
                            help='Are there flexible residues?.')
        conformers.addParam('flexibleString', params.StringParam, condition="flexRes", default='',
                            label='Flexible residues string: ',
                            help='Input flexible residues in the format "A:ILE10,VAL32;B:SER48".')

        conformers = form.addGroup("Docking")
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
        self._insertFunctionStep('gridGenerationStep')
        self._insertFunctionStep('extractFileStep')
        self._insertFunctionStep('dockingStep')
        self._insertFunctionStep('createOutputStep')

    def gridGenerationStep(self):
        recFile = os.path.abspath(self.inputAtomStruct.get().getFileName())
        peptides = self.inputPeptides.get()
        for prot in peptides:
            protFile, protName, _ = self.getProtInfo(prot)

            protFile = self.convertCifIfNeeded(protFile, self._getExtraPath())

            if not protFile.endswith('.pdbqt'):
                protFile = self.preparePeptidePDBQT(protFile, self._getExtraPath())

            args = [f'-r {recFile} -l {os.path.abspath(protFile)} -o {protName} -P {self.padding.get()}']

            if(self.flexRes.get()):
                args.append(f'-f {self.flexibleString.get()}')

            Plugin.runAGFR(self, args, cwd=self._getExtraPath())

    def extractFileStep(self):
        peptides = self.inputPeptides.get()
        for prot in peptides:
            _, protName, _ = self.getProtInfo(prot)

            extraDir = self._getExtraPath()
            outputDir = self._getPath(f'{protName}')
            os.makedirs(outputDir, exist_ok=True)

            zipFile = os.path.join(extraDir, f"{protName}.trg")
            if not os.path.exists(zipFile):
                self.error(f"AGFR output ZIP not found: {zipFile}")
                return

            zipOutput = os.path.join(outputDir, f"{protName}.trg")
            shutil.copy(zipFile, zipOutput)
            os.remove(zipFile)

            with zipfile.ZipFile(zipOutput, 'r') as zip_ref:
                zip_ref.extractall(outputDir)

    def createOutputStep(self):
        logFile = os.path.abspath(self._getPath('logs/run.stdout'))
        outputLogData = self.readOutputData(logFile)
        i = 0
        outputMols = SetOfSmallMolecules().create(outputPath=self._getPath())

        peptides = self.inputPeptides.get()
        recFile = self.inputAtomStruct.get().getFileName()
        for prot in peptides:
            protFile, protName, _ = self.getProtInfo(prot)
            resultsFolder = self.getResultsFolder(f'{protName}_docking')
            rankedFiles = self.rankedFiles(resultsFolder)
            data = outputLogData[i]
            targetFile = os.path.abspath(os.path.join(self._getPath(f"{protName}"), f"{protName}.trg"))

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
                    poseFile, targetFile,
                    modeNum, energy, affinity, bestRun
                )

                outputMols.append(newMol)

            i += 1

        outputMols.setDocked(True)
        outputMols.setProteinFile(recFile)
        self._defineOutputs(outputSmallMolecules=outputMols)

    def dockingStep(self):
        peptides = self.inputPeptides.get()

        for prot in peptides:
            protFile, protName, _ = self.getProtInfo(prot)
            targetFile = os.path.join(self._getPath(f"{protName}"), f"{protName}.trg")
            pdbFile = os.path.abspath(os.path.join(self._getExtraPath(f"{protName}.pdb")))

            cifToPdb(protFile, pdbFile)

            seq = self.getSequenceFromPdb(pdbFile)
            args = [f'-t {os.path.abspath(targetFile)} -s {seq} -N {self.nRuns.get()} -o {protName}_docking -ref {pdbFile}']

            if (self.cyc.get()):
                args.append('-cyc')
            if (self.cys.get()):
                args.append('-cys')

            resultsFolder = self.getResultsFolder(f'{protName}_docking')
            Plugin.runADCP(self, args, cwd=resultsFolder)



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
    def getResultsFolder(self, name):
        folder = os.path.abspath(os.path.join(self._getPath(), f'{name}_docking'))
        os.makedirs(folder, exist_ok=True)
        return folder

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
        """Read docking log and return a list of docking data per search."""

        with open(logFile, 'r') as f:
            allRuns = list(self.parseBlocks(f))

        return allRuns

    def parseBlocks(self, f):
        dockingData = {}
        tableStarted = False

        for rawLine in f:
            line = rawLine.strip()

            if self._isSearchBoundary(line):
                yield from self._flushData(dockingData)
                dockingData = {}
                tableStarted = False
                continue

            if self.isTableHeader(line):
                tableStarted = True
                continue

            if not tableStarted:
                continue

            parsed = self._processTableLine(line)
            if parsed:
                mode, affinity, energy, bestRun = parsed
                dockingData[mode] = {
                    "affinity": affinity,
                    "energy": energy,
                    "bestRun": bestRun,
                }

        yield from self._flushData(dockingData)

    def _isSearchBoundary(self, line):
        return self.isStartSearch(line) or self.isEndSearch(line)

    def _flushData(self, dockingData):
        if dockingData:
            yield dockingData

    def _processTableLine(self, line):
        if not self.isTableLine(line):
            return None
        return self.parseTableLine(line)

    def isStartSearch(self, line):
        return line.startswith("Performing search")

    def isTableHeader(self, line):
        return line.startswith("mode |  affinity")

    def isTableLine(self, line):
        return bool(line) and not line.startswith(("|", "-----"))

    def isEndSearch(self, line):
        return line.startswith("clean up")

    def parseTableLine(self, line):
        parts = line.split()
        if len(parts) >= 9 and parts[0].isdigit():
            return int(parts[0]), float(parts[1]), float(parts[7]), int(parts[8])
        return None

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

    def preparePeptidePDBQT(self, peptideFile, outDir):
        """
        Convert a peptide PDB/MOL2 file to PDBQT using prepare_ligand4.py via pythonsh.
        """
        ext = os.path.splitext(peptideFile)[1].lower()

        if ext == '.cif':
            pdbFile = os.path.join(outDir, os.path.splitext(os.path.basename(peptideFile))[0] + '.pdb')
            cifToPdb(peptideFile, pdbFile)
            peptideFile = pdbFile

        pdbqtFile = os.path.join(outDir, os.path.splitext(os.path.basename(peptideFile))[0] + '.pdbqt')
        prog = 'prepare_ligand4'
        pythonsh = Plugin.getProgramHome(MGL_DIC, 'bin/pythonsh ')
        scriptPath = Plugin.getADTPath(f'Utilities24/{prog}.py ')
        program = pythonsh + scriptPath

        arguments = f"-l {os.path.abspath(peptideFile)} -o {os.path.abspath(pdbqtFile)}"

        self.runJob(program, arguments, cwd=outDir)

        return pdbqtFile

    def getProtInfo(self, prot):
        filePath = (prot.getFileName())
        name = os.path.splitext(os.path.basename(filePath))[0]
        ext = os.path.splitext(filePath)[1].lower()
        return filePath, name, ext

    def convertCifIfNeeded(self, filePath, outDir):
        if filePath.lower().endswith('.cif'):
            pdbFile = os.path.join(outDir, os.path.splitext(os.path.basename(filePath))[0] + '.pdb')
            if not os.path.exists(pdbFile):
                cifToPdb(filePath, pdbFile)
            return pdbFile
        return filePath