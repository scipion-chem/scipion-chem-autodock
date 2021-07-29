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

import pwem.objects.data as data
from .constants import ATTRIBUTES_MAPPING as AM
from pwchem.objects import ProteinPocket
from pyworkflow.object import (Object, Float, Integer, String,
                               OrderedDict, CsvList, Boolean, Set, Pointer,
                               Scalar, List)



class AutodockGrid(data.EMFile):
    """A search grid in the file format of Autodock"""
    def __init__(self, **kwargs):
        data.EMFile.__init__(self, **kwargs)


class AutoLigandPocket(ProteinPocket):
    """ Represent a pocket file from autoLigand"""

    def __init__(self, filename=None, proteinFile=None, resultsFile=None, **kwargs):
        self.pocketId = int(filename.split('out')[1].split('.')[0])
        self.properties = self.parseFile(resultsFile, filename)
        kwargs.update(self.getKwargs(self.properties, AM))
        ProteinPocket.__init__(self, filename, proteinFile, **kwargs)
        self.setObjId(self.pocketId)

        #Build contact atoms
        cAtoms = self.buildContactAtoms(calculate=True)
        self.setContactAtoms(self.encodeIds(self.getAtomsIds(cAtoms)))
        cResidues = self.getResiduesFromAtoms(cAtoms)
        self.setContactResidues(self.encodeIds(self.getResiduesIds(cResidues)))

    def __str__(self):
        s = 'Fpocket pocket {}\nFile: {}'.format(self.pocketId, self.getFileName())
        return s

    def getVolume(self):
        return self.properties['Total Volume']

    def getEnergyPerVol(self):
        return self.properties['Total Energy per Vol']

    def parseFile(self, resultsFile, filename):
        props, i = {}, 1
        with open(resultsFile) as f:
            for line in f:
                if i==self.pocketId:
                    line = line.split(',')
                    #Volume
                    props[line[1].split('=')[0].strip()] = float(line[1].split('=')[1].strip())
                    #Energy/vol
                    props[line[2].strip()] = float(line[3].split('=')[1].strip())
                i+=1

        with open(filename) as f:
            npts, points = 0, []
            for line in f:
                line = line.split()
                points += [tuple(map(float, line[5:8]))]
                npts += 1
            props['nPoints'] = npts
            self._points = points
        props['class'] = 'AutoLigand'
        return props

    def getPoints(self):
        return self._points


class GridADT(data.EMFile):
    """ Represent a grid file in map (ASCIII) format generated with ADT"""
    def __init__(self, filename=None, **kwargs):
        data.EMFile.__init__(self, filename, **kwargs)
        self._radius = Float(kwargs.get('radius', None))
        self._spacing = Float(kwargs.get('spacing', None))
        self._massCenter = List(kwargs.get('massCenter', None))
        self._npts = Integer(kwargs.get('npts', None))

    def __str__(self):
        return '{} (Radius={}, Spacing={})'.format(self.__class__.__name__, self.getRadius(), self.getSpacing())

    def getRadius(self):
        return self._radius

    def setRadius(self, value):
        self._radius.set(value)

    def getSpacing(self):
        return self._spacing

    def setSpacing(self, value):
        self._spacing.set(value)

    def getMassCenter(self):
        return self._massCenter

    def setMassCenter(self, values):
        self._massCenter.set(values)

    def getNumberOfPoints(self):
        return self._npts

    def setNumberOfPoints(self, value):
        self._npts.set(value)