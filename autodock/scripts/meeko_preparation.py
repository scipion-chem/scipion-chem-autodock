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

from meeko import MoleculePreparation

from utils import getMolFilesDic, parseParams, getBaseName

if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile> 
    '''
    paramsDic = parseParams(sys.argv[1], listParams=['ligandFiles', 'moleculesFiles'], sep='::')
    ligandFiles = paramsDic['ligandFiles']
    hydra = eval(paramsDic['hydrate'])

    outDir = paramsDic['outDir']


#####################################################################
    molFileDic, mols = getMolFilesDic(ligandFiles)
    outFiles = []
    if len(mols) > 0:
        preparator = MoleculePreparation(hydrate=hydra)
        for mol in mols:
            preparator.prepare(mol)

            inFile = molFileDic[mol]
            outFile = os.path.join(outDir, getBaseName(inFile)) + '.pdbqt'
            preparator.write_pdbqt_file(outFile)
            outFiles.append(outFile)

        with open(os.path.join(outDir, 'meeko_files.txt'), 'w') as f:
            f.write('\n'.join(outFiles))

    else:
        print('None of the input molecules could be read by RDKit.\n'
              'Preparing the molecules with the RDKit ligand preparation protocol might solve this issue')

