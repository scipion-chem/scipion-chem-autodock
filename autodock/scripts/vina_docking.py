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
from multiprocessing import Process, Queue, Pipe

from utils import getBaseName, parseParams

def runVina(q, kwargs):
    try:
        v, outDir = kwargs['v'], kwargs['outDir']
        molFile, pDic = kwargs['molFile'], kwargs['pDic']
        outFile = getOutFile(molFile, outDir)
        v.set_ligand_from_file(molFile)
        v.dock(exhaustiveness=int(pDic['exhaust']), n_poses=int(pDic['nPoses']),
               min_rmsd=float(pDic['minRMSD']), max_evals=int(pDic['maxEvals']))

        v.write_poses(pdbqt_filename=outFile)
        q.put(("SUCCESS", outFile))
    except Exception as e:
        q.put(("ERROR", e))
def safeRunVina(kwargs):
    q = Queue()
    p = Process(target=runVina, args=(q, kwargs))
    p.start()
    p.join()  # Add timeout to prevent hanging
    if not q.empty():  # Check if the subprocess put a result
        status, payload = q.get()
        if status == "SUCCESS":
            return payload  # Return the result
        else:
            raise payload  # Re-raise the exception
    else:
        raise RuntimeError("C++ function did not return any result")

def getOutFile(molFile, outDir):
    return os.path.join(outDir, getBaseName(molFile)) + '.pdbqt'

if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile> <outputDir>
    '''
    pDic = parseParams(sys.argv[1], listParams=['ligandFiles', 'moleculesFiles'], sep='::')
    ligandFiles = pDic['ligandFiles']
    receptorFile = pDic['receptorFile']
    flexRecFile = pDic['flexRecFile'] if 'flexRecFile' in pDic else None

    mapsName = pDic['mapsName']
    outDir = pDic['outDir']

#####################################################################

    v = Vina(sf_name=pDic['scoreName'], cpu=int(pDic['nCPUs']))
    v.set_receptor(rigid_pdbqt_filename=receptorFile, flex_pdbqt_filename=flexRecFile)

    if pDic['scoreName'] == 'Vina':
        v.compute_vina_maps(center=eval(pDic['boxCenter']), box_size=eval(pDic['boxSize']))
    else:
        v.write_maps(mapsName, pDic['gpfFile'])
        v.load_maps(mapsName)

    outFiles, failed = [], []
    for molFile in ligandFiles:
        try:
            kwargs = {'v': v, 'molFile': molFile, 'outDir': outDir, 'pDic': pDic}
            result = safeRunVina(kwargs)
            if isinstance(result, Exception):
                raise result
            outFiles.append(result)
        except:
            failed.append(molFile)

    with open(os.path.join(outDir, f'docked_files_{pDic["it"]}.txt'), 'w') as f:
        f.write('\n'.join(outFiles))

    if failed:
        with open(os.path.join(outDir, f'failed_docked_files_{pDic["it"]}.txt'), 'w') as f:
            f.write('\n'.join(failed))


