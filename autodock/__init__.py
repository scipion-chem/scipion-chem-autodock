# **************************************************************************
# *
# * Authors:	Carlos Oscar Sorzano (coss@cnb.csic.es)
# *			 	Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *			 	Martín Salinas Antón (martin.salinas@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This package contains the protocols for the manipulation of atomic struct objects
"""

# General imports
import os, subprocess, json

# Scipion em imports
import pwem
import pyworkflow.utils as pwutils
from scipion.install.funcs import InstallHelper

# Plugin imports
from pwchem import Plugin as pwchemPlugin
from pwchem.constants import MGL_DIC, RDKIT_DIC
from pwchem.utils import insistentRun

from .bibtex import _bibtexStr
from .constants import *

# Pluging variables
_logo = 'autodock_logo.png'
__version__ = ALPHA_VERSION
chemPropFile = 'example_model_v2_regression_mol.ckpt'

# Installation variables
enVars = {'GPU_INCLUDE_PATH': pwem.Config.CUDA_BIN.replace('bin', 'include'), 'GPU_LIBRARY_PATH': pwem.Config.CUDA_LIB}

class Plugin(pwchemPlugin):
    """
    Definition of class variables. For each package, a variable will be created.
    _<packageNameInLowercase>Home will contain the full path of the package, ending with a folder whose name will be <packageNameFirstLetterLowercase>-<defaultPackageVersion> variable.
        For example: _atdHome = "~/Documents/scipion/software/em/autoDock-4.2.6"

    Inside that package, for each binary, there will also be another variable.
    _<binaryNameInLowercase>Binary will be a folder inside _<packageNameInLowercase>Home and its name will be <binaryName>.
        For example: _atdBinary = "~/Documents/scipion/software/em/autoDock-4.2.6/AutoDock"
    """

    # AutoDockGPU
    _atdgpuHome = pwchemPlugin.getDefPath(ADGPU_DIC)
    _atdgpuBinary = os.path.join(_atdgpuHome, 'AutoDockGPU')

    # Vina
    _vinaHome = pwchemPlugin.getDefPath(VINA_DIC)
    _vinaBinary = os.path.join(_vinaHome, 'AutoDock-Vina')

    # VinaGPU
    _vinagpuHome = pwchemPlugin.getDefPath(VINAGPU_DIC)
    _vinagpuBinary = os.path.join(_vinagpuHome, 'AutoDock-VinaGPU')

    # AutoSite
    _asiteHome = pwchemPlugin.getDefPath(ASITE_DIC)
    _asiteBinary = _asiteHome

    # Ringtail
    _ringtailHome = pwchemPlugin.getDefPath(RINGTAIL_DIC)

    # Scrubber
    _scrubberHome = pwchemPlugin.getDefPath(SCRUBBER_DIC)

    # GCR models
    _gcrHome = pwchemPlugin.getDefPath(GCR_DIC)


    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(ADGPU_DIC['home'], cls._atdgpuHome)
        cls._defineEmVar(VINA_DIC['home'], cls._vinaHome)
        cls._defineEmVar(VINAGPU_DIC['home'], cls._vinagpuHome)
        cls._defineEmVar(ASITE_DIC['home'], cls._asiteHome)
        cls._defineEmVar(SCRUBBER_DIC['home'], cls._scrubberHome)
        cls._defineEmVar(RINGTAIL_DIC['home'], cls._ringtailHome)
        cls._defineEmVar(GCR_DIC['home'], cls._gcrHome)

    @classmethod
    def defineBinaries(cls, env):
        """
        This function defines the binaries for each package.
        """
        cls.addADTPackage(env)
        cls.addAutoDockGPUPackage(env)
        cls.addVinaPackage(env)
        cls.addVinaGPUPackage(env)
        cls.addAutoSitePackage(env)
        cls.addRingtailPackage(env)
        cls.addScrubberPackage(env)
        cls.addGCRPackage(env)

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
        installer = InstallHelper(AUTODOCK_DIC['name'], packageHome=cls.getVar(AUTODOCK_DIC['home']),
                                                            packageVersion=AUTODOCK_DIC['version'])

        installer.getCondaEnvCommand() \
            .addCondaPackages([f'autodock={AUTODOCK_DIC["version"]}'], channel='bioconda', targetName='ADT_CONDA') \
            .addPackage(env, ['conda'], default=default)

    @classmethod
    def addAutoDockGPUPackage(cls, env, default=True):
        """ This function provides the neccessary commands for installing AutoDock-GPU. """
        # Instantiating the install helper
        installer = InstallHelper(ADGPU_DIC['name'], packageHome=cls.getVar(ADGPU_DIC['home']), packageVersion=ADGPU_DIC['version'])

        # Only install AutoDockGPU on systems with Nvidia GPUs
        if cls.getGPUPlatform() == 'nvidia':
            # Getting Nvidia card data
            compCapDic = cls.getNVIDIACompCapDic()
            nvidiaName = cls.getNVIDIAName()
            compCap = compCapDic[nvidiaName] if nvidiaName in compCapDic else None
            targetsFlag = f' TARGETS={compCap}' if compCap else ''

            # Installing package
            make_cmd = f'cd {cls._atdgpuBinary} && CC=/usr/bin/gcc-12 CXX=/usr/bin/g++-12 make DEVICE=GPU OVERLAP=ON{targetsFlag}'

            installer.getCloneCommand(cls.getAutoDockGPUGithub(),
                                      binaryFolderName=cls._atdgpuBinary,
                                      targeName='ATDGPU_CLONED') \
                .addCommand(make_cmd, 'ATDGPU_COMPILED') \
                .addPackage(env, dependencies=['git', 'make'], default=default)

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
    def addVinaGPUPackage(cls, env, default=True, testSoft=False):
        """ This function procvides the neccessary commands for installing AutoDock-VinaGPU. """
        # Instantiating the install helper
        installer = InstallHelper(VINAGPU_DIC['name'], packageHome=cls.getVar(VINAGPU_DIC['home']),
            packageVersion=VINAGPU_DIC['version'])

        # Defining GPU platform and OpenCL version
        gpuPlatform = '-DNVIDIA_PLATFORM' if cls.getGPUPlatform() == 'nvidia' else '-DAMD_PLATFORM'
        openCLVersion = '-DOPENCL_2_0' if cls.getOpenCLVersion() == '2.0' else '-DOPENCL_3_0'

        # Cloning AutoDock-VinaGPU
        installer.getCloneCommand(f'https://github.com/DeltaGroupNJUPT/Vina-GPU-{VINAGPU_DIC["version"]}.git',
            binaryFolderName=cls._vinagpuBinary, targeName='VINA_GPU_CLONED')

        # Downloading and extracting Boost library
        installer.getExtraFile('https://archives.boost.io/release/1.74.0/source/boost_1_74_0.tar.gz',
            'BOOST_DOWNLOADED', fileName=boostFilename) \
            .addCommand(f'mkdir -p {boostFoldername} && tar -xf {boostFilename} --strip-components 1 -C '
            f'{boostFoldername} && rm {boostFilename}', 'BOOST_EXTRACTED') \
            .addCommand(f'cd {boostFoldername} && ./bootstrap.sh --with-libraries=program_options,system,filesystem && ./b2',
            'BOOST_INSTALLED').getCondaEnvCommand(requirementsFile=False)

        # Installing CUDA in a Conda enviroment
        installer.addCondaPackages(['cuda'], channel="\"nvidia/label/cuda-11.7.0\"")

        # Defining AutoDock-VinaGPU makefile location
        softwares = [f'AutoDock-Vina-GPU-{VINAGPU_DIC["version"]}',
                     f'QuickVina2-GPU-{VINAGPU_DIC["version"]}',
                     f'QuickVina-W-GPU-{VINAGPU_DIC["version"]}']

        for i, soft in enumerate(softwares):
            softDir = os.path.join(cls._vinagpuBinary, soft)
            makefile = os.path.join(softDir, 'Makefile')
            softBin = f'{soft[:-2]}-{soft[-1]}'
            oldStrConfig = f"/home/shidi/Vina-GPU-{VINAGPU_DIC['version']}"

            # Modifying makefile and compiling
            installer.addCommand(
                f"{cls.getEnvActivationCommand(VINAGPU_DIC)} && python3 {makefileModifier} "
                f"{makefile} {boostPath} $CONDA_PREFIX {openCLVersion} {gpuPlatform} && "
                f"sed -i 's|{oldStrConfig}|{cls._vinagpuBinary}|g' "
                f"{softDir}/input_file_example/2bm2_config.txt",
                f"MAKEFILE_{i}_MODIFIED"
            )

            if testSoft:
                installer.addCommand(f'{cls.getEnvActivationCommand(VINAGPU_DIC)} && '
                    f'make source && ./{softBin} --config ./input_file_example/2bm2_config.txt',
                    f'VINAGPU_{i}_TESTED', workDir=softDir)

            installer.addCommand(f'{cls.getEnvActivationCommand(VINAGPU_DIC)} && make clean && make',
                f'VINAGPU_{i}_COMPILED', workDir=softDir)

        # Adding package
        installer.addPackage(env,dependencies=['wget', 'tar', 'conda', 'make', 'clinfo'], default=default)

        # Cloning AutoDock-VinaGPU
        installer.getCloneCommand('https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0.git',
                                  binaryFolderName=cls._vinagpuBinary,targeName='VINA_GPU_CLONED')

        # Downloading and extracting Boost library
        (installer.getExtraFile('https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/boost_1_82_0.tar.gz',
                               'BOOST_DOWNLOADED', fileName=boostFilename) \
            .addCommand(f'mkdir -p {boostFoldername} && tar -xf {boostFilename} --strip-components 1 -C '
            f'{boostFoldername} && rm {boostFilename}','BOOST_EXTRACTED') \
            .getCondaEnvCommand(requirementsFile=False))

        # Installing CUDA in a Conda enviroment
        installer.addCondaPackages(['cuda'], channel="\"nvidia/label/cuda-11.5.0\"")

    @classmethod
    def addAutoSitePackage(cls, env, default=True):
        """ This function provides the neccessary commands for installing AutoSite. """
        # Instantiating the install helper
        installer = InstallHelper(ASITE_DIC['name'], packageHome=cls.getVar(ASITE_DIC['home']), packageVersion=ASITE_DIC['version'])

        # Generating AutoSite installation commands
        installer.getExtraFile(cls.getADFRSuiteUrl(), 'ASITE_DOWNLOADED', fileName=cls.getASITETar())\
            .addCommand(f'tar -zxf {cls.getASITETar()} --strip-components 1 && rm {cls.getASITETar()}', 'ASITE_EXTRACTED')\
            .addCommand('./install.sh -d . -c 0 -l', 'ASITE_INSTALLED')

        # Adding package
        installer.addPackage(env, dependencies=['wget', 'conda'], default=default)

    @classmethod
    def addRingtailPackage(cls, env, default=True):
        """ This function provides the necessary commands for installing Ringtail and Meeko. """
        # Instantiating the install helper
        installer = InstallHelper(RINGTAIL_DIC['name'], packageHome=cls.getVar(RINGTAIL_DIC['home']),
                                                            packageVersion=RINGTAIL_DIC['version'])

        # Installing package
        installer.addCommand(f'{cls.getEnvActivationCommand(RDKIT_DIC)} && '
                                                 f'conda install -y ringtail={RINGTAIL_DIC["version"]} pymol-open-source -c conda-forge',
                                                 'RINGTAIL_INSTALLED'). \
            addCommand(f'{cls.getEnvActivationCommand(RDKIT_DIC)} && pip install --no-deps '
                                 f'{MEEKO_DIC["name"]}=={MEEKO_DIC["version"]} prody==2.4', 'MEEKO_INSTALLED'). \
            addCommand(f'{cls.getEnvActivationCommand(RDKIT_DIC)} && conda install -c conda-forge -y gemmi=0.7.3', 'MEEKO_DEPS_INSTALLED'). \
            addPackage(env, dependencies=['conda'], default=default)

    @classmethod
    def addScrubberPackage(cls, env, default=True):
        """ This function provides the necessary commands for installing Scrubber. """
        # Instantiating the install helper
        installer = InstallHelper(SCRUBBER_DIC['name'], packageHome=cls.getVar(SCRUBBER_DIC['home']),
                                                            packageVersion=SCRUBBER_DIC['version'])

        # Installing package
        installer.getCloneCommand(cls.getScrubberGithub(), targeName='SCRUBBER_CLONED'). \
            addCommand(f'conda create --name {cls.getEnvName(SCRUBBER_DIC)} python=3.10 -y'). \
            addCommand(f'{cls.getEnvActivationCommand(SCRUBBER_DIC)} && cd molscrub && pip install -e .',
                                 'SCRUBBER_INSTALLED'). \
            addCommand(f'{cls.getEnvActivationCommand(SCRUBBER_DIC)} && conda install -y -c conda-forge autogrid=4.2.9',
                                 'AUTOGRID_INSTALLED'). \
            addPackage(env, dependencies=['git', 'conda', 'pip'], default=default)

    @classmethod
    def addGCRPackage(cls, env, default=False):
        """ This function provides the necessary commands for installing GEDS. """
        # Instantiating the install helper
        installer = InstallHelper(GCR_DIC['name'], packageHome=cls.getVar(GCR_DIC['home']),
                                                            packageVersion=GCR_DIC['version'])

        # Installing package
        installer.getCloneCommand(cls.getGCRGithub(), binaryFolderName=cls.getEnvName(GCR_DIC), targeName='GCR_CLONED'). \
            addCommand(f'cd {cls.getEnvName(GCR_DIC)} && conda env create -f environment.yml', 'GCR_INSTALLED'). \
            addCommand(f'wget {cls.getChempropModelLink()}', 'CHEMPROP_MODEL'). \
            addCommand(f'{cls.getEnvActivationCommand(GCR_DIC)} && git clone {cls.getCOATIGithub()} && '
                                 f'cd COATI && pip install .', 'COATI_INSTALLED'). \
            addPackage(env, dependencies=['git', 'conda', 'pip'], default=default)

    # ---------------------------------- Protocol functions-----------------------
    @classmethod
    def runAutoDock4(cls, protocol, args, cwd=None, popen=False):
        fullProgram = f'{cls.getEnvActivationCommand(AUTODOCK_DIC)} && autodock4 '
        if not popen:
            protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
        else:
            subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

    @classmethod
    def runAutodockGPU(cls, protocol, args, cwd=None):
        """ Run autodock gpu command from a given protocol """
        program = ''
        progDir = pwchemPlugin.getProgramHome(ADGPU_DIC, path='AutoDockGPU/bin')
        binPaths = os.listdir(progDir)
        for binName in binPaths:
            if 'autodock_gpu' in binName:
                program = os.path.join(progDir, binName)
                break

        if program:
            kwargs = {"cwd": cwd}
            insistentRun(protocol, program, args, **kwargs)
        else:
            print('No autodock_gpu binary was found in {}'.format(progDir))

    @classmethod
    def runVina(cls, protocol, program="vina", args=None, cwd=None):
        if program == 'vina':
            program = cls.getVinaPath('bin/vina')
        protocol.runJob(program, args, env=cls.getEnviron(), cwd=cwd)

    @classmethod
    def runVinaGPU(cls, protocol, program, args):
        """ Run Vina GPU command from a given protocol """
        progPath = cls.getVinaGPUBinary(program)
        progDir = os.path.dirname(progPath)

        if os.path.exists(progPath):
            kwargs = {"cwd": progDir}
            insistentRun(protocol, progPath, args, **kwargs)
        else:

            print('No Vina GPU binary was found in {}'.format(progDir))

    @classmethod
    def getVinaGPUBinary(cls, program):
        progBin = f'{program}-GPU-{VINAGPU_DIC["version"]}'
        progDir = pwchemPlugin.getProgramHome(VINAGPU_DIC, path=f'AutoDock-VinaGPU/{progBin}')
        progPath = os.path.join(progDir, f'{progBin[:-2]}-{progBin[-1]}')
        return progPath

    @classmethod
    def runMeekoLigand(cls, protocol, args, cwd=None, popen=False):
        fullProgram = f'{cls.getEnvActivationCommand(RDKIT_DIC)} && mk_prepare_ligand.py '
        if not popen:
            protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
        else:
            subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

    @classmethod
    def runMeekoReceptor(cls, protocol, args, cwd=None, popen=False):
        fullProgram = f'{cls.getEnvActivationCommand(RDKIT_DIC)} && mk_prepare_receptor.py '
        if not popen:
            protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
        else:
            subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

    @classmethod
    def runScrubber(cls, protocol, args, cwd=None, popen=False):
        fullProgram = f'{cls.getEnvActivationCommand(SCRUBBER_DIC)} && scrub.py '
        if not popen:
            protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
        else:
            subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

    @classmethod
    def runAutogrid(cls, protocol, args, cwd=None, popen=False):
        fullProgram = f'{cls.getEnvActivationCommand(SCRUBBER_DIC)} && scrub.py '
        if not popen:
            protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
        else:
            subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

    @classmethod
    def runRingtail(cls, protocol, args, cwd=None, popen=False, getOutput=False):
        fullProgram = f'{cls.getEnvActivationCommand(RDKIT_DIC)} && rt_process_vs '
        if not popen:
            protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
        else:
            if getOutput:
                return subprocess.check_output(f'{fullProgram} {args}', cwd=cwd, shell=True)
            else:
                subprocess.Popen(f'{fullProgram} {args}', cwd=cwd, shell=True)

    @classmethod
    def runScript(cls, protocol, scriptName, args, envDict, cwd=None, popen=False):
        """ Run rdkit command from a given protocol. """
        scriptName = cls.getScriptsDir(scriptName)
        fullProgram = '%s && %s %s' % (cls.getEnvActivationCommand(envDict), 'python', scriptName)
        if not popen:
            protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
        else:
            subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

    @classmethod
    def runADCP(cls, protocol, args, cwd=None, popen=False):
        fullProgram = pwchemPlugin.getProgramHome(ASITE_DIC, path='bin/adcp')
        if not popen:
            protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
        else:
            subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

    @classmethod
    def runAGFR(cls, protocol, args, cwd=None, popen=False):
        fullProgram = pwchemPlugin.getProgramHome(ASITE_DIC, path='bin/agfr')
        if not popen:
            protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
        else:
            subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)

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
    def getScriptsDir(cls, scriptName):
        return cls.getPluginHome('scripts/%s' % scriptName)

    @classmethod
    def getModelsDir(cls, modelName=''):
        return os.path.join(cls.getPluginHome('models'), modelName)

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
    def getGCRPath(cls, path=''):
        return pwchemPlugin.getProgramHome(GCR_DIC, path)

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
    def getScrubberGithub(cls):
        return 'https://github.com/forlilab/molscrub.git'

    @classmethod
    def getGCRGithub(cls):
        return 'https://github.com/DaniDelHoyo/GCR_Regression_ForliLab.git'

    @classmethod
    def getChempropModelLink(cls, modelFile='example_model_v2_regression_mol.ckpt'):
        return f'https://github.com/chemprop/chemprop/raw/refs/heads/main/tests/data/{modelFile}'

    @classmethod
    def getChemPropFile(cls):
        return os.path.abspath(pwchemPlugin.getProgramHome(GCR_DIC, chemPropFile))

    @classmethod
    def getCOATIGithub(cls):
        return 'https://github.com/terraytherapeutics/COATI.git'

    @classmethod
    def getASITETar(cls):
        return ASITE_DIC['name'] + '-' + ASITE_DIC['version'] + '.tgz'


    @classmethod
    def getOpenCLVersion(cls):
        """ This function returns the OpenCL version available in the current device. """
        try:
            # If OpenCL is not installed, the following command will return an error
            version = subprocess.run(["clinfo | grep \"Platform Version\" | awk \'{print $4}\'"], capture_output=True, shell=True)
            return version.stdout.decode('utf-8').strip()
        except Exception:
            return None

    @classmethod
    def getGPUPlatform(cls):
        """ This function returns the GPU platform of the current device. """
        try:
            # If no Nvidia drivers are present, the following command will reuturn an error
            subprocess.run(["nvidia-smi"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            return 'nvidia'
        except Exception:
            return 'amd'
