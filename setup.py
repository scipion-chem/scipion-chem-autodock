"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
import os, re

here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

VERSION_FILE = os.path.join(here, "autodock", "version.py")

def get_version():
    with open(VERSION_FILE, "r") as f:
        content = f.read()
    match = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]', content, re.M)
    if not match:
        raise RuntimeError(f"Unable to find version string in {VERSION_FILE}")
    return match.group(1)

setup(
    name='scipion-chem-autodock',
    version=get_version(),
    description='Scipion framework plugin for the use of tools provided by Autodock software suite',
    long_description=long_description,
    url='https://github.com/scipion-chem/scipion-chem-autodock',
    author='Carlos Oscar Sorzano',
    author_email='scipion@cnb.csic.es',
    keywords='scipion pdb scipion-3',
    packages=find_packages(),
    install_requires=[requirements],
    include_package_data=True,
    package_data={
       'autodock': ['autodock_logo.png'],
    },
    entry_points={
        'pyworkflow.plugin': 'autodock = autodock'
    }
)
