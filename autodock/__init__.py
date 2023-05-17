# **************************************************************************
# *
# * Authors:	Carlos Oscar Sorzano (coss@cnb.csic.es)
# *			 	Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *			 	Martín Salinas Antón (martin.salinas@cnb.csic.es)
# *
# * Unidad deBioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307USA
# *
# *All comments concerning this program package may be sent to the
# *e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This package contains the protocols for the manipulation of atomic struct objects
"""

import os, subprocess, json
from .bibtex import _bibtexStr
from .constants import *
from .install_helper import InstallHelper

import pwem
import pyworkflow.utils as pwutils

from pwchem import Plugin as pwchemPlugin
from pwchem.constants import MGL_DIC

# Pluging variables
_logo = 'autodock_logo.png'
__version__ = ALPHA_VERSION
AUTODOCK_DIC = {'name': 'autoDock', 'version': '4.2.6', 'home': 'AUTODOCK_HOME'}
VINA_DIC = {'name': 'vina', 'version': '1.2', 'home': 'VINA_HOME'}
ASITE_DIC = {'name': 'autoSite', 'version': '1.0', 'home': 'AUTOSITE_HOME'}
MEEKO_DIC = {'name': 'meeko', 'version': '0.3.3', 'home': 'MEEKO_HOME'}
ADGPU_DIC = {'name': 'autoDockGPU', 'version': DEFAULT_VERSION, 'home': 'ADGPU_HOME'}

enVars = {'GPU_INCLUDE_PATH': pwem.Config.CUDA_BIN.replace('bin', 'include'), 'GPU_LIBRARY_PATH': pwem.Config.CUDA_LIB}

class Plugin(pwem.Plugin):
	@classmethod
	def _defineVariables(cls):
		cls._defineEmVar(AUTODOCK_DIC['home'], os.path.join(pwem.Config.EM_ROOT, AUTODOCK_DIC['name'] + '-' + AUTODOCK_DIC['version']))
		cls._defineEmVar(VINA_DIC['home'], os.path.join(pwem.Config.EM_ROOT, VINA_DIC['name'] + '-' + VINA_DIC['version']))
		cls._defineEmVar(ASITE_DIC['home'], os.path.join(pwem.Config.EM_ROOT, ASITE_DIC['name'] + '-' + ASITE_DIC['version']))
		cls._defineEmVar(ADGPU_DIC['home'], os.path.join(pwem.Config.EM_ROOT, ADGPU_DIC['name'] + '-' + ADGPU_DIC['version']))
	
	@classmethod
	def defineBinaries(cls, env):
		"""
        This function defines the binaries for each package.
        """
		cls.addADTPackage(env)
		cls.addVinaPackage(env)
		cls.addAutoSitePackage(env)
		cls.addMeekoPackage(env)
		cls.addAutoDockGPUPackage(env)

	@classmethod
	def getEnviron(cls):
		""" Create the needed environment for AutoDock-GPU installation. """
		environ = pwutils.Environ(os.environ)
		environ.update({
			'PATH': pwem.Config.CUDA_BIN,
			'LD_LIBRARY_PATH': pwem.Config.CUDA_LIB,
	 	}, position=pwutils.Environ.END)
		environ.update(enVars)

		return environ

	@classmethod
	def addADTPackage(cls, env, default=True):
		""" This function provides the neccessary commands for installing AutoDock. """
		# Instantiating the install helper
		installer = InstallHelper()

		# Installing protocol
		installer.getExtraFile(cls.getADTSuiteUrl(), cls.getADTTar(), 'ADT_DOWNLOADED')\
			.addCommand(f'tar -xf {cls.getADTTar()} --strip-components 1 && rm {cls.getADTTar()}', 'ADT_EXTRACTED')\
			.addProtocolPackage(env, AUTODOCK_DIC['name'], version=AUTODOCK_DIC['version'], dependencies=['wget'], default=default)

	@classmethod
	def addVinaPackage(cls, env, default=True):
		""" This function provides the neccessary commands for installing Vina. """
		# Instantiating the install helper
		installer = InstallHelper()

		# Installing protocol
		installer.getCondaEnvCommand(VINA_DIC['name'], cls.getVar(VINA_DIC['home']), binaryVersion=VINA_DIC['version'], pythonVersion='3.10', requirementsFile=False)\
			.addCondaPackages(VINA_DIC['name'], ['numpy=1.23.5', 'swig', 'boost-cpp', 'sphinx_rtd_theme', 'vina'], binaryVersion=VINA_DIC['version'], channel='conda-forge')\
			.addProtocolPackage(env, VINA_DIC['name'], VINA_DIC['version'], ['conda'], default=default)

	@classmethod
	def addAutoSitePackage(cls, env, default=True):
		""" This function provides the neccessary commands for installing AutoSite. """
		# TODO: Check how only adtools or adfr needed
		# Instantiating the install helper
		installer = InstallHelper()

		# Installing protocol
		installer.getExtraFile(cls.getADFRSuiteUrl(), ASITE_DIC['home'], 'ASITE_DOWNLOADED', fileName=cls.getASITETar())\
			.addCommand(f'tar -zxf {cls.getASITETar()} --strip-components 1 && rm {cls.getASITETar()}', 'ASITE_EXTRACTED')\
			.addCommand('./install.sh -d . -c 0 -l', 'ASITE_INSTALLED')\
			.addProtocolPackage(env, ASITE_DIC['name'], ASITE_DIC['version'], ['wget'], default=default)

	@classmethod
	def addMeekoPackage(cls, env, default=True):
		MEEKO_INSTALLED = 'meeko_installed'
		meekoCommands = '%s %s && ' % (cls.getCondaActivationCmd(), cls.getVar("RDKIT_ENV_ACTIVATION"))
		meekoCommands += 'pip install meeko && touch {}'.format(MEEKO_INSTALLED)
		meekoCommands = [(meekoCommands, MEEKO_INSTALLED)]

		env.addPackage(MEEKO_DIC['name'], version=MEEKO_DIC['version'],
					 tar='void.tgz', commands=meekoCommands, default=default)

	@classmethod
	def addAutoDockGPUPackage(cls, env, default=True):
		ADGPU_INSTALLED = 'ADGPU_INSTALLED'

		compCap = None
		compCapDic = cls.getNVIDIACompCapDic()

		nvidiaName = cls.getNVIDIAName()
		if nvidiaName in compCapDic:
			compCap = compCapDic[nvidiaName]

		adGPUCommands = 'cd .. && rm -r %s && git clone %s && cd %s && ' % \
						(ADGPU_DIC['name'], cls.getAutoDockGPUGithub(), ADGPU_DIC['name'])
		adGPUCommands += 'make DEVICE=GPU OVERLAP=ON '
		if compCap:
			adGPUCommands += 'TARGETS={} '.format(compCap)
		adGPUCommands += '&& touch {}'.format(ADGPU_INSTALLED)
		adGPUCommands = [(adGPUCommands, ADGPU_INSTALLED)]


		env.addPackage(ADGPU_DIC['name'], tar='void.tgz', commands=adGPUCommands, default=default,
					 vars=enVars, updateCuda=True)

	# ---------------------------------- Protocol functions-----------------------
	@classmethod
	def runAutodockGPU(cls, protocol, args, cwd=None):
		""" Run autodock gpu command from a given protocol """
		program = ''
		progDir = pwchemPlugin.getProgramHome(ADGPU_DIC, path='bin')
		binPaths = os.listdir(progDir)
		for binName in binPaths:
			if 'autodock_gpu' in binName:
				program = os.path.join(progDir, binName)
				break

		if program:
			protocol.runJob(program, args, env=cls.getEnviron(), cwd=cwd)
		else:
			print('No autodock_gpu binary was found in {}'.format(progDir))

	@classmethod
	def runVina(cls, protocol, program="vina", args=None, cwd=None):
		if program == 'vina':
			program = cls.getVinaPath('bin/vina')
		protocol.runJob(program, args, env=cls.getEnviron(), cwd=cwd)

	@classmethod
	def runScript(cls, protocol, scriptName, args, env, cwd=None, popen=False):
		""" Run rdkit command from a given protocol. """
		scriptName = cls.getScriptsDir(scriptName)
		fullProgram = '%s %s && %s %s' % (cls.getCondaActivationCmd(), cls.getEnvActivation(env), 'python', scriptName)
		if not popen:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
		else:
			subprocess.check_call(fullProgram + args, cwd=cwd, shell=True)
		
	# ---------------------------------- Utils functions-----------------------
	@classmethod
	def getNVIDIACompCapDic(cls):
		with open(cls.getPluginHome('utils/NVIDIA_ComputeCapabilities.json')) as f:
			jDic = json.load(f)
		return jDic

	@classmethod
	def getNVIDIAName(cls):
		return subprocess.check_output("nvidia-smi -L", shell=True).decode().strip().split(':')[1].split('(')[0]\
			.lower().replace('nvidia', '').strip()

	@classmethod
	def getPluginHome(cls, path=""):
		import autodock
		fnDir = os.path.split(autodock.__file__)[0]
		return os.path.join(fnDir, path)

	@classmethod
	def getEnvActivation(cls, env):
		activation = cls.getVar("{}_ENV_ACTIVATION".format(env.upper()))
		return activation

	@classmethod
	def getScriptsDir(cls, scriptName):
		return cls.getPluginHome('scripts/%s' % scriptName)

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
	def getAutoDockGPUGithub(cls):
		return 'https://github.com/ccsb-scripps/AutoDock-GPU.git'

	@classmethod
	def getADTTar(cls):
		return AUTODOCK_DIC['name'] + '-' + AUTODOCK_DIC['version'] + '.tar'

	@classmethod
	def getASITETar(cls):
		return ASITE_DIC['name'] + '-' + ASITE_DIC['version'] + '.tgz'