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

from utils import getMolFilesDic, parseParams, getBaseName

def fixLigand(mol):
    """
    Add explicit hydrogens to carbon atoms and ensure tetravalent nitrogens are assigned a +1 charge
    """
    chemProblems = Chem.DetectChemistryProblems(mol)
    for chemProblem in chemProblems:
        if chemProblem.GetType() == 'AtomValenceException':
            atom = mol.GetAtomWithIdx(chemProblem.GetAtomIdx())
            if atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 0 and atom.GetExplicitValence() == 4:
                atom.SetFormalCharge(1)
    Chem.SanitizeMol(mol)
    for a in mol.GetAtoms():
        rad = a.GetNumRadicalElectrons()
        if rad:
            a.SetNumExplicitHs(rad)
            a.SetNumRadicalElectrons(0)
    return Chem.AddHs(Chem.RemoveHs(mol), addCoords=True)

def getBiggestFrag(mol):
    frags = Chem.GetMolFrags(mol, asMols=True)
    maxi, bMol = 0, None
    for mol in frags:
        if mol.GetNumAtoms() > maxi:
            maxi, bMol = mol.GetNumAtoms(), mol
    return bMol

if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile> 
    '''
    paramsDic = parseParams(sys.argv[1], listParams=['ligandFiles', 'moleculesFiles'], sep='::')
    ligandFiles = paramsDic['ligandFiles']
    hydra = eval(paramsDic['hydrate'])

    outDir = paramsDic['outDir']
    writeOut = paramsDic['writeOut'] if 'writeOut' in paramsDic else False


#####################################################################
    molFileDic, _ = getMolFilesDic(ligandFiles)
    outFiles = []
    if len(molFileDic) > 0:
        preparator = MoleculePreparation(hydrate=hydra)
        for mol, molFile in molFileDic.items():
            mol = fixLigand(mol)
            mol = getBiggestFrag(mol)
            preparator.prepare(mol)

            outFile = os.path.join(outDir, getBaseName(molFile)) + '.pdbqt'
            preparator.write_pdbqt_file(outFile)
            outFiles.append(outFile)

        if writeOut:
            with open(os.path.join(outDir, 'meeko_files.txt'), 'w') as f:
                f.write('\n'.join(outFiles))

    else:
        print('None of the input molecules could be read by RDKit.\n'
              'Preparing the molecules with the RDKit ligand preparation protocol might solve this issue')

