"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

from dynamo import __version__

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='scipion-em-dynamo',  # Required
    version=__version__,  # Required
    description='Plugin to use DYNAMO within the Scipion framework',  # Required
    long_description=long_description,  # Optional
    url='https://github.com/scipion-em/scipion-em-dynamo',  # Optional
    author='I2PC, Scipion Team',  # Optional
    author_email='scipion@cnb.csic.es',  # Optional
    keywords='scipion-3.0 electron-microscopy cryo-em structural-biology image-processing',  # Optional
    packages=find_packages(),
    install_requires=[requirements],
    entry_points={
        'pyworkflow.plugin': 'dynamo = dynamo'
    },
    package_data={  # Optional
       'dynamo': ['icon.png', 'protocols.conf'],
    },
    classifiers=[  # Optional
                      # How mature is this project? Common values are
                      #   3 - Alpha
                      #   4 - Beta
                      #   5 - Production/Stable
                      'Development Status :: 4 - Beta',

                      # Indicate who your project is intended for
                      #   'Intended Audience :: Users',

                      # Pick your license as you wish
                      'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

                      # Specify the Python versions you support here. In particular, ensure
                      # that you indicate whether you support Python 2, Python 3 or both.
                      'Programming Language :: Python :: 3'

                      # Specify Framework
                      'Framework :: Scipion :: 3'
                  ],
)
