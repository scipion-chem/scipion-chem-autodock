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
import os, sys

from rdkit import Chem
from meeko import MoleculePreparation


def getBaseFileName(filename):
    return os.path.splitext(os.path.basename(filename))[0]

def readLigands(ligandsFiles, removeHs=False):
    mols_dict = {}

    for molFile in ligandsFiles:
        if molFile.endswith('.mol2'):
            m = Chem.MolFromMol2File(molFile, removeHs=removeHs)
            mols_dict[m] = molFile

        elif molFile.endswith('.mol'):
            m = Chem.MolFromMolFile(molFile, removeHs=removeHs)
            mols_dict[m] = molFile

        elif molFile.endswith('.pdb'):
            m = Chem.MolFromPDBFile(molFile, removeHs=removeHs)
            mols_dict[m] = molFile

        elif molFile.endswith('.smi'):
            f = open(molFile, "r")
            firstline = next(f)
            m = Chem.MolFromSmiles(str(firstline), removeHs=removeHs)
            mols_dict[m] = molFile

        elif molFile.endswith('.sdf'):
            suppl = Chem.SDMolSupplier(molFile, removeHs=removeHs)
            for mol in suppl:
                mols_dict[mol] = molFile

    mols = list(mols_dict.keys())

    return mols_dict, mols

def parseParams(paramsFile):
    paramsDic = {}
    with open(paramsFile) as f:
        for line in f:
            key, value = line.strip().split('::')
            if key == 'ligandFiles' or key == 'moleculesFiles':
                paramsDic[key] = value.strip().split()
            else:
                paramsDic[key] = value.strip()
    return paramsDic

if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile> 
    '''
    paramsDic = parseParams(sys.argv[1])
    ligandFiles = paramsDic['ligandFiles']
    keepHs = eval(paramsDic['keepNonPolar'])
    hydra = eval(paramsDic['hydrate'])

    outDir = paramsDic['outDir']


#####################################################################
    molFileDic, mols = readLigands(ligandFiles)
    outFiles = []
    if len(mols) > 0:
        preparator = MoleculePreparation(keep_nonpolar_hydrogens=keepHs, hydrate=hydra)
        for mol in mols:
            preparator.prepare(mol)
            preparator.show_setup()

            inFile = molFileDic[mol]
            outFile = os.path.join(outDir, getBaseFileName(inFile)) + '.pdbqt'
            preparator.write_pdbqt_file(outFile)
            outFiles.append(outFile)

        with open(os.path.join(outDir, 'meeko_files.txt'), 'w') as f:
            f.write('\n'.join(outFiles))

    else:
        print('None of the input molecules could be read by RDKit.\n'
              'Preparing the molecules with the RDKit ligand preparation protocol might solve this issue')


