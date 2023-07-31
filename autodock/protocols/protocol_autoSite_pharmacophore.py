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

import os

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, BooleanParam, IntParam, EnumParam, FloatParam, LEVEL_ADVANCED
from pyworkflow.protocol import params

from pwchem.objects import PharmacophoreChem, PharmFeature
from pwchem.constants import *
from pwchem import Plugin as pwchem_plugin

from autodock import Plugin as autodock_plugin



class ProtChemAutoSiteGenPharmacophore(EMProtocol):
    """Generate a Pharmacophore from a AutoSite structural ROI """
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
