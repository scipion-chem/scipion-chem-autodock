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

import os

import pwem.objects.data as data
from pyworkflow.object import Float, Integer, String

from autodock import Plugin
from autodock.constants import ADGPU


class AutodockGrid(data.EMFile):
    """A search grid in the file format of Autodock"""
    def __init__(self, **kwargs):
        data.EMFile.__init__(self, **kwargs)


class GridADT(data.EMFile):
    """ Represent a grid file in map (ASCIII) format generated with ADT"""
    def __init__(self, filename=None, proteinFile=None, **kwargs):
        data.EMFile.__init__(self, filename, **kwargs)
        self._proteinFile = String(proteinFile)
        self._radius = Float(kwargs.get('radius', None))
        self._spacing = Float(kwargs.get('spacing', None))
        self._massCX = Float(kwargs.get('massCX', None))
        self._massCY = Float(kwargs.get('massCY', None))
        self._massCZ = Float(kwargs.get('massCZ', None))
        self._npts = Integer(kwargs.get('npts', None))

    def __str__(self):
        return '{} (Radius={}, Spacing={})'.format(self.__class__.__name__, self.getRadius(), self.getSpacing())

    def getRadius(self):
        return self._radius.get()

    def setRadius(self, value):
        self._radius.set(value)

    def getSpacing(self):
        return self._spacing.get()

    def setSpacing(self, value):
        self._spacing.set(value)

    def getMassCenter(self):
        return [self._massCX.get(), self._massCY.get(), self._massCZ.get()]

    def setMassCenter(self, values):
        self._massCX.set(values[0])
        self._massCY.set(values[1])
        self._massCZ.set(values[2])

    def getNumberOfPoints(self):
        return self._npts.get()

    def setNumberOfPoints(self, value):
        self._npts.set(value)

    def getProteinFile(self):
        return self._proteinFile.get()

    def getFilesDirectory(self):
        return '/'.join(self.getProteinFile().split('/')[:-1])


class RingtailDatabase(data.EMFile):
    """A Scipion object to refer to a RingTail virtual screening database"""
    def __init__(self, receptorFile=None, dbType=ADGPU, **kwargs):
        data.EMFile.__init__(self, **kwargs)
        self._receptorFile = self._radius = String(receptorFile)
        self._type = String(dbType)

    def getDBDir(self):
        inDBFile = os.path.abspath(self.getFileName())
        return os.path.dirname(inDBFile)

    def getSummary(self):
        args = f'read --input_db {self.getFileName()} -su'
        ringSum = Plugin.runRingtail(None, args, popen=True, getOutput=True)
        return ringSum.decode("utf-8")

    def createSumFile(self, path):
        with open(path, 'w') as f:
            f.write(self.getSummary())

    def getReceptorFile(self):
        return self._receptorFile.get()

    def setReceptorFile(self, file):
        self._receptorFile.set(String(file))

    def getType(self):
        return self._type.get()

    def setType(self, value):
        self._type.set(String(value))

    def getBookmarks(self):
        from autodock import Plugin
        from pwchem.constants import RDKIT_DIC

        inFile = os.path.abspath(self.getFileName())
        oFile = os.path.join(os.path.dirname(inFile), 'bookmarks.txt')
        Plugin.runScript(None, 'getRTBookmarks.py', f'{inFile} {oFile}', envDict=RDKIT_DIC, popen=True)

        with open(oFile) as f:
            bookMarks = eval(f.read().strip())
        return bookMarks

    def displayPlot(self, bookmark=None, pymol=True):
        from autodock import Plugin
        inDB = os.path.abspath(self.getFileName())
        flag = 'pymol' if pymol else "plot"

        args = f'read -i {inDB} --{flag} '
        if bookmark:
            args += f'-s {bookmark} '
        cwd = self.getDBDir()
        Plugin.runRingtail(None, args, popen=True, cwd=cwd)

    def buildFilterArgs(self, bookmark='base', scoreDic=None, clusterDic=None,
                        maxAtoms=None, smarts=None, minHB=0, vdwIntLines='', hbIntLines='',
                        outDir=None, outBest=True):
        dbFile = os.path.abspath(self.getFileName())
        args = f'read -i {dbFile} -s {bookmark} '
        if self.getType() == 'Vina':
            args += '-m vina '

        if scoreDic:
            for scoreName, scoreValue in scoreDic.items():
                args += f'-{scoreName} {scoreValue} '

        if clusterDic:
            for clustName, clustValue in clusterDic.items():
                args += f'-{clustName} {clustValue} '

        if maxAtoms:
            args += f'-mna {maxAtoms} '

        if smarts:
            args += f'--ligand_substruct {smarts.strip()} '

        if minHB > 0:
            args += f'-hc {minHB} '

        if vdwIntLines:
            vdwIntList = vdwIntLines.split("\n")
            args += f'-vdw {"-vdw ".join(vdwIntList)} '

        if hbIntLines:
            hbIntList = hbIntLines.split("\n")
            args += f'-hb {"-hb ".join(hbIntList)} '

        if outDir:
            args += f'-sdf {outDir} --individual_sdf_files '
            if not outBest:
                args += '-oap '

        return args

    def performBaseFilter(self, bookmark='base'):
        args = self.buildFilterArgs(bookmark, scoreDic={'e': 1000})
        cwd = self.getDBDir()
        Plugin.runRingtail(None, args, cwd=cwd, popen=True)

