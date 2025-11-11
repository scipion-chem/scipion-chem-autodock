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

import pyworkflow
from pwchem.objects import SetOfSmallMolecules
from pwem.objects import SetOfAtomStructs
from pwem.protocols import EMProtocol
from pyworkflow.object import String, Float, Integer
from pyworkflow.protocol import params

from autodock import Plugin
from autodock.objects import GridADT, SetOfGridADT

INPUT_TYPE = ['SetOfSmallMolecules', 'SetOfAtomStructs']


class ProtGenerateTargetFile(EMProtocol):
    """Prepare target file for AGFR docking."""
    _label = 'Grid generation with AGFR'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputAtomStruct', params.PointerParam, pointerClass="AtomStruct",
                      label='Receptor protein:', allowsNull=False,
                      help='It must be in pdbqt format.')
        #todo input smallMols or atomStructs
        form.addParam('inputType', params.EnumParam, label='Input type: ', default=INPUT_TYPE[0],
                      choices=INPUT_TYPE, allowsNull=False)

        form.addParam('inputSmallMolecules', params.PointerParam, pointerClass="SetOfSmallMolecules",
                      #condition=f'inputType=={INPUT_TYPE[0]}',
                      label='Set of small molecules:', allowsNull=True,
                      help='It must be in pdb or mol2 format, you may use Schrodinger convert to change it.')
        form.addParam('inputAtomStructs', params.PointerParam, pointerClass="SetOfAtomStructs",
                      #condition=f'inputType=={INPUT_TYPE[1]}',
                      label='Set of atom structures:', allowsNull=True,
                      help='It must be in pdb or mol2 format, you may use Schrodinger convert to change it.')

        conformers = form.addGroup("Parameters")
        conformers.addParam('padding', params.FloatParam, default=4.0,
                            label='Padding: ',
                            help='Amount of padding added to each side of the box.')
        conformers.addParam('flexRes', params.BooleanParam, default=False,
                            label='Flexible residues: ',
                            help='Are there flexible residues?.')
        conformers.addParam('flexibleString', params.StringParam, condition="flexRes", default='',
                            label='Flexible residues string: ',
                            help='Input flexible residues in the format "A:ILE10,VAL32;B:SER48".')
        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
        self._insertFunctionStep('createFileStep')
        self._insertFunctionStep('extractFileStep')
        self._insertFunctionStep('createOutputStep')

    def createFileStep(self):
        recFile = os.path.abspath(self.inputAtomStruct.get().getFileName())
        if self.inputType.get() == 'SetOfSmallMolecules':
            peptides = self.inputSmallMolecules.get()
        else:
            peptides = self.inputAtomStructs.get()
        for prot in peptides:
            protFile = os.path.abspath(prot.getFileName())
            protName = os.path.splitext(os.path.basename(protFile))[0]
            args = [f'-r {recFile} -l {protFile} -o {protName} -P {self.padding.get()}']

            if(self.flexRes.get()):
                args.append(f'-f {self.flexibleString.get()}')

            Plugin.runAGFR(self, args, cwd=self._getExtraPath())

    def extractFileStep(self):
        if self.inputType.get() == 'SetOfSmallMolecules':
            peptides = self.inputSmallMolecules.get()
        else:
            peptides = self.inputAtomStructs.get()
        for prot in peptides:
            protFile = os.path.abspath(prot.getFileName())
            protName = os.path.splitext(os.path.basename(protFile))[0]
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
        recFile = os.path.abspath(self.inputAtomStruct.get().getFileName())
        grids = SetOfGridADT(filename=self._getPath('setOfGrids.sqlite'))
        if self.inputType.get() == 'SetOfSmallMolecules':
            peptides = self.inputSmallMolecules.get()
        else:
            peptides = self.inputAtomStructs.get()
        for prot in peptides:
            protFile = os.path.abspath(prot.getFileName())
            protName = os.path.splitext(os.path.basename(protFile))[0]
            logFile = os.path.abspath(self._getExtraPath(f'{protName}.log'))

            data = self.getInfo(logFile)

            fileName = os.path.join(self._getPath(f"{protName}"), f"{protName}.trg")
            grid = GridADT(fileName, proteinFile=recFile, spacing=data['spacing'], massCX=data['center'][0], massCY=data['center'][1], massCZ=data['center'][2], tool='AGFR')
            grid._peptideFile = pyworkflow.object.String()
            grid._XLength = pyworkflow.object.Float()
            grid._YLength = pyworkflow.object.Float()
            grid._ZLength = pyworkflow.object.Float()
            grid._XSize = pyworkflow.object.Float()
            grid._YSize = pyworkflow.object.Float()
            grid._ZSize = pyworkflow.object.Float()
            grid._numPockets = pyworkflow.object.Integer()

            grid.setAttributeValue('_peptideFile', protFile)
            grid.setAttributeValue('_XLength' , data['length'][0])
            grid.setAttributeValue('_YLength' , data['length'][1])
            grid.setAttributeValue('_ZLength' , data['length'][2])
            grid.setAttributeValue('_XSize', data['size'][0])
            grid.setAttributeValue('_YSize', data['size'][1])
            grid.setAttributeValue('_ZSize', data['size'][2])
            grid.setAttributeValue('_numPockets', data['numPockets'])
            grids.append(grid)

        self._defineOutputs(outputGrids=grids)


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

