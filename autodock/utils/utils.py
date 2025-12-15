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

def splitSDF(sdfFile, oriName=None, oDir=None):
  '''Split multi molecule sdf files.
  Takes into account different conformations of the same molecule.
  '''
  
  def getMolName(molText):
    molName, confId = molText.split('\n')[0].strip(), None
    if '_i' in molName:
      molName, confId = molName.split('_i')
      confId = int(confId)
    return molName, confId

  molDic, confFile = {}, None
  oDir = os.path.dirname(sdfFile) if oDir is None else oDir

  with open(sdfFile) as f:
    sdfText = f.read().strip()
  molTexts = [molText.strip() for molText in sdfText.split('$$$$') if molText.strip()]

  if len(molTexts) > 1:
    if oriName is not None:
      confFile = sdfFile.replace('.sdf', f'_{oriName}.sdf')
      os.rename(sdfFile, confFile)

    for molText in molTexts:
      molName, confId = getMolName(molText)
      if confId is None:
          confId = 0
      confId += 1
      oFile = os.path.join(oDir, f'{molName}.sdf')
      if confId is not None:
        oFile = oFile.replace('.sdf', f'-{confId}.sdf')

      if molName not in molDic:
        molDic[molName] = []
      with open(oFile, 'w') as fo:
        fo.write(f'{molText}\n\n$$$$')
      molDic[molName].append((oFile, confFile))

  elif len(molTexts) == 1:
    molName, confId = getMolName(molTexts[0])
    molDic = {molName: [(sdfFile, confFile)]}

  return molDic
