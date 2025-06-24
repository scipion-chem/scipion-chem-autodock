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

"""
This protocol is used to perform generate a pharmacophore from a structural ROI from AutoSite, which contains
atomtype information

"""


from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam
from pyworkflow.protocol import params

from pwchem.objects import PharmacophoreChem, PharmFeature

class ProtChemAutoSiteGenPharmacophore(EMProtocol):
    """Generate a Pharmacophore from a AutoSite structural ROI 
    
    User IA Manual: AutoSitePharmacophore Protocol in Scipion-Chem

The AutoSitePharmacophore protocol is intended for the identification of
pharmacophoric features within predicted binding pockets on a target receptor.
It functions by analyzing previously generated AutoSite map outputs, typically
derived from grid-based site prediction, and extracting chemical features such
as hydrogen bond donors and acceptors, hydrophobic regions, and aromatic
centers that are relevant for ligand binding.

To use this protocol effectively, the user must first supply a valid AutoSite
output directory, which contains the spatial information of predicted ligand
binding sites in a receptor. This data is used to extract pharmacophoric points.
The protocol optionally allows the inclusion of a reference ligand in PDBQT
format. When such a ligand is provided, the pharmacophore model can be filtered
or refined to better reflect features relevant to known interactions or to
constrain the model to regions near the ligand?s pose. This improves the
biological relevance of the predicted features and aids downstream applications
such as virtual screening, scaffold hopping, or de novo design.

The user has control over the filtering behavior through several parameters. A
distance cutoff determines how close a pharmacophoric feature must be to the
reference ligand to be retained. This spatial constraint helps focus the model
on relevant portions of the pocket. Additionally, one can specify which types of
pharmacophoric features to include in the final model, such as hydrophobic,
hydrogen-bond donor, or acceptor features, depending on the goals of the
analysis. Further, the output can be limited to features found only in selected
binding sites or cavities, especially when AutoSite has predicted multiple
regions of interest.

The protocol generates as output a pharmacophore description in standard formats
that can be visualized within Scipion or exported for use in external tools.
These features can then guide ligand design strategies or be used in
pharmacophore-based screening to identify new potential hits. Ultimately, this
protocol serves as a bridge between structure-based binding site prediction and
ligand-centric design approaches, enabling a more interpretable and
chemically-relevant description of receptor interaction potential."""
    _label = 'AutoSite pharmacophore'

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Input')
        group.addParam('inputStructROIs', PointerParam, pointerClass="SetOfStructROIs",
                      label='Input AutoSite structural ROIs:', allowsNull=False,
                      help="Input the set of sutoSite structural ROIs you want to generate a pharmacophore from")

        form.addParam('inputStructROISelect', params.StringParam, label="Select structural ROI: ", important=True,
                      help='Select the specific ROI you want to use for the pharmacophore generation')


    # --------------------------- Steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    def createOutputStep(self):
        typeDic = {'C': 'Hydrophobe', 'H': 'Donor', 'O': 'Acceptor'}
        roiObj = self.getSelectedROI()
        clusts = self.getClustersInfo(roiObj)

        outPharm = PharmacophoreChem().create(outputPath=self._getPath())
        outPharm.setProteinFile(roiObj.getProteinFile())
        for point in clusts:
            pharmFeat = PharmFeature(type=typeDic[point['type']], radius=point['rad'],
                                     x=point['x'], y=point['y'], z=point['z'])
            outPharm.append(pharmFeat)

        self._defineOutputs(outputPharmacophore=outPharm)

    # --------------------------- Utils functions --------------------
    def getSelectedROI(self):
        for roi in self.inputStructROIs.get():
            if roi.__str__() == self.inputStructROISelect.get():
                myROI = roi
                break

        return myROI

    def getClustersInfo(self, asROI):
        clusts = []
        with open(asROI._extraFile.get()) as f:
            for line in f:
                sl = line.split()
                clusts.append({'x': sl[5], 'y': sl[6], 'z': sl[7], 'rad': sl[-3], 'type': sl[2]})
        return clusts

    # --------------------------- INFO functions -----------------------------------

    def _citations(self):
        return []

    def _validate(self):
        vals = []
        if self.getSelectedROI().getPocketClass() != 'AutoSite':
            vals.append('Unfortunately, only AutoSite structural ROIs can be used to generate pharmacophores '
                        'using this protocol')
        return vals
