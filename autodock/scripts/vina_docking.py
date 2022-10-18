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

from vina import Vina

def getBaseFileName(filename):
    return os.path.splitext(os.path.basename(filename))[0]

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
    '''Use: python <scriptName> <paramsFile> <outputDir>
    '''
    pDic = parseParams(sys.argv[1])
    ligandFiles = pDic['ligandFiles']
    receptorFile = pDic['receptorFile']
    outDir = pDic['outDir']


#####################################################################

    v = Vina(sf_name=pDic['scoreName'], cpu=int(pDic['nCPUs']))
    v.set_receptor(rigid_pdbqt_filename=receptorFile)

    outFiles = []
    for molFile in ligandFiles:
        v.set_ligand_from_file(molFile)
        v.compute_vina_maps(center=eval(pDic['boxCenter']), box_size=eval(pDic['boxSize']))
        v.dock(exhaustiveness=int(pDic['exhaust']), n_poses=int(pDic['nPoses']),
               min_rmsd=float(pDic['minRMSD']), max_evals=int(pDic['maxEvals']))

        outFile = os.path.join(outDir, getBaseFileName(molFile)) + '.pdbqt'
        v.write_poses(pdbqt_filename=outFile)
        outFiles.append(outFile)

    with open(os.path.join(outDir, 'docked_files.txt'), 'w') as f:
        f.write('\n'.join(outFiles))


