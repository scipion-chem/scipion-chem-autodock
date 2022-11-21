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
AUTODOCK_DIC = {'name': 'autodock', 'version': '4.2.6', 'home': 'AUTODOCK_HOME'}
VINA_DIC = {'name': 'vina', 'version': '1.1.2', 'home': 'VINA_HOME'}
ASITE_DIC = {'name': 'AutoSite', 'version': '1.0', 'home': 'AUTOSITE_HOME'}

class Plugin(pwem.Plugin):
    @classmethod
    def defineBinaries(cls, env):
        cls.addADTPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addAutoSitePackage(env, default=bool(cls.getCondaActivationCmd()))

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(AUTODOCK_DIC['home'], AUTODOCK_DIC['name'] + '-' + AUTODOCK_DIC['version'])
        cls._defineEmVar(VINA_DIC['home'], VINA_DIC['name'] + '-' + VINA_DIC['version'])
        cls._defineEmVar(ASITE_DIC['home'], ASITE_DIC['name'] + '-' + ASITE_DIC['version'])

    @classmethod
    def addADTPackage(cls, env, default=False):
        ADT_INSTALLED = 'adt_installed'
        adtCommands = 'wget {} -O {} --no-check-certificate && '.format(cls.getADTSuiteUrl(), cls.getADTTar())
        adtCommands += 'tar -xf {} --strip-components 1 && '.format(cls.getADTTar())
        adtCommands += 'rm {} && touch {}'.format(cls.getADTTar(), ADT_INSTALLED)
        adtCommands = [(adtCommands, ADT_INSTALLED)]

        env.addPackage(AUTODOCK_DIC['name'], version=AUTODOCK_DIC['version'],
                       tar='void.tgz',
                       commands=adtCommands,
                       default=True)

        VINA_INSTALLED = 'vina_installed'
        vinaCommands = 'wget {} -O {} --no-check-certificate && '.format(cls.getVinaURL(), cls.getVinaTgz())
        vinaCommands += 'tar -zxf {} --strip-components 1 && '.format(cls.getVinaTgz())
        vinaCommands += 'rm {} && touch {}'.format(cls.getVinaTgz(), VINA_INSTALLED)
        vinaCommands = [(vinaCommands, VINA_INSTALLED)]

        env.addPackage(VINA_DIC['name'], version=VINA_DIC['version'],
                       tar='void.tgz',
                       commands=vinaCommands,
                       default=True)

    @classmethod
    def addAutoSitePackage(cls, env, default=False):
      # todo: check how only adtools or adfr needed
      ASITE_INSTALLED = 'asite_installed'
      asiteCommands = 'wget {} -O {} --no-check-certificate && '.format(cls.getADFRSuiteUrl(), cls.getASITETar())
      asiteCommands += 'tar -zxf {} --strip-components 1 && '.format(cls.getASITETar())
      asiteCommands += './install.sh -d . -c 0 -l && '.format(cls.getASITETar())
      asiteCommands += 'rm {} && touch {}'.format(cls.getASITETar(), ASITE_INSTALLED)
      asiteCommands = [(asiteCommands, ASITE_INSTALLED)]

      env.addPackage(ASITE_DIC['name'], version=ASITE_DIC['version'],
                     tar='void.tgz',
                     commands=asiteCommands,
                     default=True)

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
        return os.path.join(cls.getVar('AUTODOCK_HOME'), path)

    @classmethod
    def getVinaPath(cls, path=''):
      return os.path.join(cls.getVar('VINA_HOME'), path)

    @classmethod
    def getASitePath(cls, path=''):
      return os.path.join(cls.getVar('AUTOSITE_HOME'), path)

    @classmethod
    def getADTSuiteUrl(cls):
      return 'https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/autodocksuite-{}-x86_64Linux2.tar'.\
        format(AUTODOCK_DIC['version'])

    @classmethod
    def getADFRSuiteUrl(cls):
      return 'https://ccsb.scripps.edu/adfr/download/1038/ADFRsuite_x86_64Linux_{}.tar.gz'.\
        format(ASITE_DIC['version'])

    @classmethod
    def getVinaURL(cls):
      return 'https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_{}_linux_x86.tgz'.\
        format(VINA_DIC['version'].replace('.', '_'))

    @classmethod
    def getADTTar(cls):
      return AUTODOCK_DIC['name'] + '-' + AUTODOCK_DIC['version'] + '.tar'

    @classmethod
    def getASITETar(cls):
      return ASITE_DIC['name'] + '-' + ASITE_DIC['version'] + '.tgz'

    @classmethod
    def getVinaTgz(cls):
      return VINA_DIC['name'] + '-' + VINA_DIC['version'] + '.tgz'

    @classmethod
    def runVina(cls, protocol, program="vina", args=None, cwd=None):
        if program == 'vina':
            program = cls.getVinaPath('bin/vina')
        protocol.runJob(program, args, env=cls.getEnviron(), cwd=cwd)


