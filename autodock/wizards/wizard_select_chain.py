# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Alberto Manuel Parra Pérez (amparraperez@gmail.com)
# *
# * Biocomputing Unit, CNB-CSIC
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
This wizard will extract the chains from a atomic structure (pdb) file in
order to select it in the protocol.
Then, it will load the structure and will take all chain related
information such as name and number of residues.
"""

from pwchem.wizards import *
from pwchem.utils import RESIDUES1TO3

from autodock.protocols import *

SelectChainWizard().addTarget(protocol=ProtChemADTPrepareReceptor,
                              targets=['chain_name'],
                              inputs=['inputAtomStruct'],
                              outputs=['chain_name'])


SelectChainWizardQT().addTarget(protocol=ProtChemAutodock,
                                targets=['flexChain'],
                                inputs=[{'fromReceptor': ['inputAtomStruct', 'inputStructROIs']}],
                                outputs=['flexChain'])

SelectChainWizardQT().addTarget(protocol=ProtChemVina,
                                targets=['flexChain'],
                                inputs=[{'fromReceptor': ['inputAtomStruct', 'inputStructROIs']}],
                                outputs=['flexChain'])

SelectMultiChainWizard().addTarget(protocol=ProtChemADTPrepareReceptor,
                                   targets=['chain_name'],
                                   inputs=['inputAtomStruct'],
                                   outputs=['chain_name'])

SelectElementWizard().addTarget(protocol=ProtChemAutoSiteGenPharmacophore,
                                targets=['inputStructROISelect'],
                                inputs=['inputStructROIs'],
                                outputs=['inputStructROISelect'])

SelectResidueWizardQT().addTarget(protocol=ProtChemAutodock,
                                  targets=['flexPosition'],
                                  inputs=[{'fromReceptor': ['inputAtomStruct', 'inputStructROIs']}, 'flexChain'],
                                  outputs=['flexPosition'])

SelectResidueWizardQT().addTarget(protocol=ProtChemVina,
                                  targets=['flexPosition'],
                                  inputs=[{'fromReceptor': ['inputAtomStruct', 'inputStructROIs']}, 'flexChain'],
                                  outputs=['flexPosition'])

class AddFlexibleWizard(EmWizard):
  _targets = [(ProtChemAutodock, ['addFlex']), (ProtChemVina, ['addFlex'])]

  def show(self, form, *params):
    protocol = form.protocol
    chainStr, resStr = protocol.flexChain.get(),  protocol.flexPosition.get()
    chainId = json.loads(chainStr)['chain']
    resIds, posIdxs = json.loads(resStr)['residues'], json.loads(resStr)['index'].split('-')

    allResStr = []
    for i, rId in enumerate(resIds):
        allResStr.append(RESIDUES1TO3[rId] + str(int(posIdxs[0]) + i))

    form.setVar('flexList', protocol.flexList.get() +
                '{}:{}\n'.format(chainId, '_'.join(allResStr)))


