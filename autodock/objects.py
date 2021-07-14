# -*- coding: utf-8 -*-
#  **************************************************************************
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

import pyworkflow.object as pwobj
import pwem.objects.data as data
from pwchem.objects import proteinPocket


class AutodockGrid(data.EMFile):
    """A search grid in the file format of Autodock"""
    def __init__(self, **kwargs):
        data.EMFile.__init__(self, **kwargs)


class autoLigandPocket(proteinPocket):
    """ Represent a pocket file from autoLigand"""

    def __init__(self, filename=None, resultsFile=None, **kwargs):
        proteinPocket.__init__(self, filename, **kwargs)
        print(filename)
        print(filename.split('out'))
        print(filename.split('out')[1].split('.'))
        self.pocketId = int(filename.split('out')[1].split('.')[0])
        self.properties = self.parseFile(resultsFile)
        self.setObjId(self.pocketId)

    def __str__(self):
        s = 'Fpocket pocket {}\nFile: {}'.format(self.pocketId, self.getFileName())
        return s

    def getVolume(self):
        return self.properties['Total Volume']

    def getEnergyPerVol(self):
        return self.properties['Total Energy per Vol']

    def parseFile(self, filename):
        props, i = {}, 1
        with open(filename) as f:
            for line in f:
                if i==self.pocketId:
                    line = line.split(',')
                    print(line)
                    #Volume
                    props[line[1].split('=')[0].strip()] = float(line[1].split('=')[1].strip())
                    #Energy/vol
                    props[line[2].strip()] = float(line[3].split('=')[1].strip())

        return props

