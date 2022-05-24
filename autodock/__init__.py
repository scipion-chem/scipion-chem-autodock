# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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
This package contains the protocols for
manipulation of atomic struct objects
"""

import os
from .bibtex import _bibtexStr

import pwem

from pwchem import Plugin as pwchemPlugin
from pwchem.constants import MGL_DIC

_logo = 'autodock.png'
AUTODOCK, AUTODOCK_DEFAULT_VERSION = 'autodock', '4.2.6'
VINA, VINA_DEFAULT_VERSION = 'vina', '1.1.2'

class Plugin(pwem.Plugin):
    @classmethod
    def defineBinaries(cls, env):
        ADT_INSTALLED = 'adt_installed'
        adtCommands = 'wget {} -O {} --no-check-certificate && '.format(cls.getADTSuiteUrl(), cls.getADTTar())
        adtCommands += 'tar -xf {} --strip-components 1 && '.format(cls.getADTTar())
        adtCommands += 'rm {} && touch {}'.format(cls.getADTTar(), ADT_INSTALLED)
        adtCommands = [(adtCommands, ADT_INSTALLED)]

        env.addPackage(AUTODOCK, version=AUTODOCK_DEFAULT_VERSION,
                       tar='void.tgz',
                       commands=adtCommands,
                       default=True)

        VINA_INSTALLED = 'vina_installed'
        vinaCommands = 'wget {} -O {} --no-check-certificate && '.format(cls.getVinaURL(), cls.getVinaTgz())
        vinaCommands += 'tar -zxf {} --strip-components 1 && '.format(cls.getVinaTgz())
        vinaCommands += 'rm {} && touch {}'.format(cls.getVinaTgz(), VINA_INSTALLED)
        vinaCommands = [(vinaCommands, VINA_INSTALLED)]

        env.addPackage(VINA, version=VINA_DEFAULT_VERSION,
                       tar='void.tgz',
                       commands=vinaCommands,
                       default=True)


    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar('AUTODOCK_HOME', AUTODOCK + '-' + AUTODOCK_DEFAULT_VERSION)
        cls._defineEmVar('VINA_HOME', VINA + '-' + VINA_DEFAULT_VERSION)

    @classmethod
    def getPluginHome(cls, path=""):
        import pwchem
        fnDir = os.path.split(pwchem.__file__)[0]
        return os.path.join(fnDir,path)


    @classmethod
    def getADTPath(cls, path=''):
        return pwchemPlugin.getProgramHome(MGL_DIC, os.path.join('MGLToolsPckgs', 'AutoDockTools', path))

    @classmethod
    def getAutodockPath(cls, path=''):
        return os.path.join(cls.getVar('AUTODOCK_HOME'),path)

    @classmethod
    def getVinaPath(cls, path=''):
      return os.path.join(cls.getVar('VINA_HOME'), path)

    @classmethod
    def getADTSuiteUrl(cls):
      return 'https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/autodocksuite-4.2.6-x86_64Linux2.tar'

    @classmethod
    def getVinaURL(cls):
      return 'https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz'

    @classmethod
    def getADTTar(cls):
      return AUTODOCK + '-' + AUTODOCK_DEFAULT_VERSION + '.tar'

    @classmethod
    def getVinaTgz(cls):
      return VINA + '-' + VINA_DEFAULT_VERSION + '.tgz'

    @classmethod
    def runVina(cls, protocol, program="vina", args=None, cwd=None):
        if program == 'vina':
            program = cls.getVinaPath('bin/vina')
        protocol.runJob(program, args, env=cls.getEnviron(), cwd=cwd)


