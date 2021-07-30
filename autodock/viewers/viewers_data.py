# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

import pwem.viewers.views as views
import pwem.viewers.showj as showj
import pyworkflow.viewer as pwviewer
from ..protocols import ProtChemAutoLigand
from pwchem.viewers import PyMolViewer
import os

class SetOfDatabaseIDView(views.ObjectView):
    """ Customized ObjectView for SetOfDatabaseID. """
    def __init__(self, project, inputid, path, other='',
                 viewParams={}, **kwargs):
        defaultViewParams = {showj.MODE: 'metadata',
                             showj.RENDER: '_PDBLigandImage'}
        defaultViewParams.update(viewParams)
        views.ObjectView.__init__(self, project, inputid, path, other,
                                  defaultViewParams, **kwargs)

class AutoLigandViewer(pwviewer.Viewer):
  _label = 'Viewer pockets'
  _environments = [pwviewer.DESKTOP_TKINTER]
  _targets = [ProtChemAutoLigand]

  def _validate(self):
    return []

  # =========================================================================
  # ShowAtomStructs
  # =========================================================================

  def getInputAtomStructFile(self):
    return os.path.abspath(self.protocol.inputAtomStruct.get().getFileName())

  def _visualize(self, obj, **kwargs):
    pdbFileName = self.protocol.getPDBName()
    outDir = os.path.abspath(self.protocol._getExtraPath())
    pymolFile = outDir + '/' + pdbFileName + '.pml'

    pymolV = PyMolViewer(project=self.getProject())
    pymolV.visualize(pymolFile, cwd=outDir)
