# # -*- coding: utf-8 -*-
# # # **************************************************************************
# # # *
# # # * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# # # *
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************
import sys

from rdkit import Chem
import openbabel

from utils import parseMoleculeFile

def pdbqt_to_smi(pdbqt_file):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdbqt", "can")

    mol = openbabel.OBMol()
    success = obConversion.ReadFile(mol, pdbqt_file)

    if success:
        smiles = obConversion.WriteString(mol).strip()
        return smiles.split()[0]
    else:
        return "Unable to load pdbqt"

if __name__ == "__main__":
    '''Use: python <scriptName> <listFile> 
    Convert the molecule files listed in the input file to SMI. 
    Input file must have each mol path in each line as the first column of a csv. 
    Further line info will be kept in the putput file
    '''
    inFile = sys.argv[1]
    outFile = sys.argv[2]


#####################################################################
    with open(outFile, 'w') as fo:
        with open(inFile) as f:
            for line in f:
                molFile = line.split(',')[0].strip()
                if molFile.endswith('.pdbqt'):
                    smi = pdbqt_to_smi(molFile)
                    mol = Chem.MolFromSmiles(smi)
                else:
                    mol = parseMoleculeFile(molFile)

                canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
                fo.write(line.replace(molFile, canonical_smiles))

