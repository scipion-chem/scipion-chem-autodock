# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

from pyworkflow.protocol import params

from pwem.protocols import EMProtocol

from autodock import Plugin as autodockPlugin
from autodock.objects import RingtailDatabase

inputOptions = ['AutoDock-GPU', 'Vina']

class ProtRingtailCreation(EMProtocol):
  """Create a Ringtail database with the results of a docking experiment with AutoDock-GPU or Vina"""
  _label = 'Ringtail database creation'
  _program = ""

  def _defineParams(self, form):
    form.addSection(label='Input')
    group = form.addGroup('Input')
    group.addParam('inputSmallMolecules', params.PointerParam, pointerClass="SetOfSmallMolecules",
                   label='Input small molecules: ',
                   help="Input small molecules that will be processed by Ringtail")
    group.addParam('inOption', params.EnumParam, label='Molecules docked by: ', choices=inputOptions, default=0,
                   help='Define the origin of the input molecules')

  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
    self._insertFunctionStep('createOutputStep')

  def getDLGDir(self):
    molFile = self.inputSmallMolecules.get().getFirstItem().getPoseFile()
    molDir = os.path.dirname(molFile)
    return os.path.abspath(os.path.join('/'.join(os.path.split(molDir)[:-1]), 'extra'))

  def getSumPath(self):
    return self._getExtraPath('ringSum.txt')

  def createOutputStep(self):
    dlgDir = self.getDLGDir()
    args = f'write --file_path {dlgDir} --recursive -o ringtail.db '
    if self.inOption.get() == 1:
      args += '-m vina '
    autodockPlugin.runRingtail(self, args, cwd=self._getPath())

    outputDB = RingtailDatabase(filename=self._getPath('ringtail.db'))
    outputDB.createSumFile(self.getSumPath())
    self._defineOutputs(outputRingtail=outputDB)


  def _summary(self):
    s = []
    if os.path.exists(self.getSumPath()):
      with open(self.getSumPath()) as f:
        s.append(f.read())
    return s
