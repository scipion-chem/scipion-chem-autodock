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

import pyworkflow.viewer as pwviewer
from pyworkflow.protocol import params

import pwem.viewers.views as views
import pwem.viewers.showj as showj

from autodock.objects import RingtailDatabase

class SetOfDatabaseIDView(views.ObjectView):
    """ Customized ObjectView for SetOfDatabaseID. """
    def __init__(self, project, inputid, path, other='',
                 viewParams={}, **kwargs):
        defaultViewParams = {showj.MODE: 'metadata',
                             showj.RENDER: '_PDBLigandImage'}
        defaultViewParams.update(viewParams)
        views.ObjectView.__init__(self, project, inputid, path, other,
                                  defaultViewParams, **kwargs)

class ViewerRingtail(pwviewer.ProtocolViewer):
  _label = 'Viewer ringtail database'
  _targets = [RingtailDatabase]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Visualization of ringtail database')
    group = form.addGroup('Select bookmark')
    group.addParam('bookmark', params.StringParam, default='passing_results', label='Database bookmark: ',
                  help='Select the bookmark stored in the database to be displayed')

    group = form.addGroup('Displays')
    group.addParam('displayPymol', params.LabelParam, label='Display molecules: ',
                   help='Display molecules in the bookmark using PyMol')
    group.addParam('displayPlot', params.LabelParam, label='Display energies: ',
                   help='Display a plot with the energies vs ligand efficiencies of the molecules in the bookmark')

  def getBookmarks(self):
    return self.protocol.getBookmarks()

  def _getVisualizeDict(self):
    return {
      'displayPymol': self._viewMols,
      'displayPlot': self._viewEnergies,
    }

  def _viewMols(self, e=None):
    molDB = self.getInputDB()
    molDB.displayPlot(bookmark=self.bookmark.get().strip(), pymol=True)

  def _viewEnergies(self, e=None):
    molDB = self.getInputDB()
    molDB.displayPlot(bookmark=self.bookmark.get().strip(), pymol=False)

  def getInputDB(self):
    molDB = None
    if isinstance(self.protocol, RingtailDatabase):
      molDB = self.protocol
    elif hasattr(self.protocol, 'outputRingtail'):
      molDB = getattr(self.protocol, 'outputRingtail')
    else:
      print('Cannot find outputRingtail')
    return molDB


