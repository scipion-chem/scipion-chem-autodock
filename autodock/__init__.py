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
ADGPU_DIC = {'name': 'autoDockGPU', 'version': DEFAULT_VERSION, 'home': 'ADGPU_HOME'}
VINA_DIC = {'name': 'vina', 'version': '1.2.3', 'home': 'VINA_HOME'}
VINAGPU_DIC = {'name': 'vinaGPU', 'version': DEFAULT_VERSION, 'home': 'VINAGPU_HOME'}
ASITE_DIC = {'name': 'autoSite', 'version': DEFAULT_VERSION, 'home': 'AUTOSITE_HOME'}
MEEKO_DIC = {'name': 'meeko', 'version': '0.3.3'}


# Installation variables
enVars = {'GPU_INCLUDE_PATH': pwem.Config.CUDA_BIN.replace('bin', 'include'), 'GPU_LIBRARY_PATH': pwem.Config.CUDA_LIB}

class Plugin(pwem.Plugin):
	"""
    Definition of class variables. For each package, a variable will be created.
    _<packageNameInLowercase>Home will contain the full path of the package, ending with a folder whose name will be <packageNameFirstLetterLowercase>-<defaultPackageVersion> variable.
        For example: _atdHome = "~/Documents/scipion/software/em/autoDock-4.2.6"
    
    Inside that package, for each binary, there will also be another variable.
    _<binaryNameInLowercase>Binary will be a folder inside _<packageNameInLowercase>Home and its name will be <binaryName>.
        For example: _atdBinary = "~/Documents/scipion/software/em/autoDock-4.2.6/AutoDock"
    """
	# AutoDock
	_atdHome = os.path.join(pwem.Config.EM_ROOT, AUTODOCK_DIC['name'] + '-' + AUTODOCK_DIC['version'])
	_atdBinary = os.path.join(_atdHome, 'AutoDock')

	# AutoDockGPU
	_atdgpuHome = os.path.join(pwem.Config.EM_ROOT, ADGPU_DIC['name'] + '-' + ADGPU_DIC['version'])
	_atdgpuBinary = os.path.join(_atdgpuHome, 'AutoDockGPU')

	# Vina
	_vinaHome = os.path.join(pwem.Config.EM_ROOT, VINA_DIC['name'] + '-' + VINA_DIC['version'])
	_vinaBinary = os.path.join(_vinaHome, 'AutoDock-Vina')

	# VinaGPU
	_vinagpuHome = os.path.join(pwem.Config.EM_ROOT, VINAGPU_DIC['name'] + '-' + VINAGPU_DIC['version'])
	_vinagpuBinary = os.path.join(_vinagpuHome, 'AutoDock-VinaGPU')

	# AutoSite
	_asiteHome = os.path.join(pwem.Config.EM_ROOT, ASITE_DIC['name'] + '-' + ASITE_DIC['version'])
	_asiteBinary = _asiteHome

	@classmethod
	def _defineVariables(cls):
		cls._defineEmVar(AUTODOCK_DIC['home'], cls._atdHome)
		cls._defineEmVar(ADGPU_DIC['home'], cls._atdgpuHome)
		cls._defineEmVar(VINA_DIC['home'], cls._vinaHome)
		cls._defineEmVar(VINAGPU_DIC['home'], cls._vinagpuHome)
		cls._defineEmVar(ASITE_DIC['home'], cls._asiteHome)
	
	@classmethod
	def defineBinaries(cls, env):
		"""
        This function defines the binaries for each package.
        """
		#cls.addADTPackage(env)
		#cls.addAutoDockGPUPackage(env)
		#cls.addVinaPackage(env)
		cls.addVinaGPUPackage(env)
		#cls.addAutoSitePackage(env)

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
		installer = InstallHelper(AUTODOCK_DIC['name'], packageHome=cls.getVar(AUTODOCK_DIC['home']), packageVersion=AUTODOCK_DIC['version'])

		# Installing package
		installer.getExtraFile(cls.getADTSuiteUrl(), 'ADT_DOWNLOADED', fileName=cls.getADTTar())\
			.addCommand(f'tar -xf {cls.getADTTar()} --strip-components 1 && rm {cls.getADTTar()}', 'ADT_EXTRACTED')\
			.addPackage(env, dependencies=['wget'], default=default)

	@classmethod
	def addAutoDockGPUPackage(cls, env, default=True):
		""" This function provides the neccessary commands for installing AutoDock-GPU. """
		# Instantiating the install helper
		installer = InstallHelper(ADGPU_DIC['name'], packageHome=cls.getVar(ADGPU_DIC['home']), packageVersion=ADGPU_DIC['version'])

		# Getting Nvidia card data
		compCapDic = cls.getNVIDIACompCapDic()
		nvidiaName = cls.getNVIDIAName()
		compCap = compCapDic[nvidiaName] if nvidiaName in compCapDic else None
		targetsFlag = f' TARGETS={compCap}' if compCap else ''

		# Installing package
		installer.getCloneCommand(cls.getAutoDockGPUGithub(), binaryFolderName=cls._atdBinary, targeName='ATDGPU_CLONED')\
			.addCommand(f'cd {cls._atdBinary} && make DEVICE=GPU OVERLAP=ON{targetsFlag}', 'ATDGPU_COMPILED')\
			.addPackage(env, dependencies=['git', 'make'], default=default, vars=enVars, updateCuda=True)

	@classmethod
	def addVinaPackage(cls, env, default=True):
		""" This function provides the neccessary commands for installing AutoDock-Vina. """
		# Instantiating the install helper
		installer = InstallHelper(VINA_DIC['name'], packageHome=cls.getVar(VINA_DIC['home']), packageVersion=VINA_DIC['version'])

		# Installing package
		installer.getCloneCommand('https://github.com/ccsb-scripps/AutoDock-Vina.git', targeName='ADT_VINA_CLONED')\
			.getCondaEnvCommand(pythonVersion='3.10', requirementsFile=False)\
			.addCondaPackages(['numpy=1.23.5', 'swig', 'boost-cpp', 'sphinx_rtd_theme', f'vina={VINA_DIC["version"]}'], channel='conda-forge')\
			.addPackage(env, ['git', 'conda'], default=default)
	
	@classmethod
	def addVinaGPUPackage(cls, env, default=True):
		""" This function procvides the neccessary commands for installing AutoDock-VinaGPU. """
		# Instantiating the install helper
		installer = InstallHelper(VINAGPU_DIC['name'], packageHome=cls.getVar(VINAGPU_DIC['home']), packageVersion=VINAGPU_DIC['version'])

		# Defining variables
		boostFoldername = 'boost'
		boostFilename = boostFoldername + '.tar.gz'

		# Installing package
		installer.getCloneCommand('https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0.git', binaryFolderName=cls._vinagpuBinary, targeName='VINA_GPU_CLONED')\
			.getExtraFile('https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/boost_1_82_0.tar.gz', 'BOOST_DOWNLOADED', fileName=boostFilename)\
			.addCommand(f'mkdir -p {boostFoldername} && tar -xf {boostFilename} --strip-components 1 -C {boostFoldername} && rm {boostFilename}', 'BOOST_EXTRACTED')\
			.getCondaEnvCommand(requirementsFile=False)\
			.addCondaPackages(['cuda'], channel="\"nvidia/label/cuda-11.5.0\"")\
			.addPackage(env, dependencies=['wget', 'tar', 'conda', 'make'], default=default)

	@classmethod
	def addAutoSitePackage(cls, env, default=True):
		""" This function provides the neccessary commands for installing AutoSite. """
		# TODO: Check how only adtools or adfr needed
		# Instantiating the install helper
		installer = InstallHelper(ASITE_DIC['name'], packageHome=cls.getVar(ASITE_DIC['home']), packageVersion=ASITE_DIC['version'])

		# Generating AutoSite installation commands
		installer.getExtraFile(cls.getADFRSuiteUrl(), 'ASITE_DOWNLOADED', fileName=cls.getASITETar())\
			.addCommand(f'tar -zxf {cls.getASITETar()} --strip-components 1 && rm {cls.getASITETar()}', 'ASITE_EXTRACTED')\
			.addCommand('./install.sh -d . -c 0 -l', 'ASITE_INSTALLED')
		
		# Generating meeko installation commands
		installer.addCommand(f'{cls.getCondaActivationCmd()} {cls.getVar("RDKIT_ENV_ACTIVATION")} && pip install {MEEKO_DIC["name"]}=={MEEKO_DIC["version"]}', 'MEEKO_INSTALLED')
		
		# Adding package
		installer.addPackage(env, ['wget', 'conda'], default=default)

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
		return cls.getVar("{}_ENV_ACTIVATION".format(env.upper()))

	@classmethod
	def getScriptsDir(cls, scriptName):
		return cls.getPluginHome('scripts/%s' % scriptName)

	@classmethod
	def getVinaPath(cls, path=''):
		return os.path.join(cls.getVar('VINA_HOME'), path)

	@classmethod
	def getPackagePath(cls, package='AUTODOCK', path=''):
		return os.path.join(cls.getVar('{}_HOME'.format(package.upper())), path)

	@classmethod
	def getVinaScriptsPath(cls, scriptName=''):
		return cls.getPackagePath(package='VINA', path='AutoDock-Vina/example/autodock_scripts/{}'.format(scriptName))

	@classmethod
	def getADTPath(cls, path=''):
		return pwchemPlugin.getProgramHome(MGL_DIC, os.path.join('MGLToolsPckgs', 'AutoDockTools', path))

	@classmethod
	def getADTSuiteUrl(cls):
		return 'https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/autodocksuite-{}-x86_64Linux2.tar'.\
			format(AUTODOCK_DIC['version'])

	@classmethod
	def getADFRSuiteUrl(cls):
		return 'https://ccsb.scripps.edu/adfr/download/1038/ADFRsuite_x86_64Linux_{}.tar.gz'.\
			format(ASITE_DIC['version'])

	@classmethod
	def getAutoDockGPUGithub(cls):
		return 'https://github.com/ccsb-scripps/AutoDock-GPU.git'

	@classmethod
	def getVinaGithub(cls):
		return 'https://github.com/ccsb-scripps/AutoDock-Vina.git'

	@classmethod
	def getADTTar(cls):
		return AUTODOCK_DIC['name'] + '-' + AUTODOCK_DIC['version'] + '.tar'

	@classmethod
	def getASITETar(cls):
		return ASITE_DIC['name'] + '-' + ASITE_DIC['version'] + '.tgz'

	@classmethod
	def getASITETar(cls):
		return ASITE_DIC['name'] + '-' + ASITE_DIC['version'] + '.tgz'