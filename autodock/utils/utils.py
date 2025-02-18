# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

import os

def splitSDF(sdfFile, oriName='conformers'):
  '''Split sdf conformer files'''
  with open(sdfFile) as f:
    sdfText = f.read()

  mols = sdfText.split('$$$$')[:-1]
  if len(mols) > 1:
    oFiles = []
    for i, molText in enumerate(mols):
      oFiles += [sdfFile.replace('.sdf', f'_{i + 1}.sdf')]
      with open(oFiles[-1], 'w') as fo:
        fo.write(f'{molText.strip()}\n\n$$$$')
    confFile = sdfFile.replace('.sdf', f'_{oriName}.sdf')
    os.rename(sdfFile, confFile)
  else:
    oFiles = [sdfFile]
    confFile = None

  return oFiles, confFile