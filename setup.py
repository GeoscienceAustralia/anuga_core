#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Cournapeau David <cournape@gmail.com>
#               2010 Fabian Pedregosa <fabian.pedregosa@inria.fr>
# License: 3-clause BSD
#
# Setup.py taken from scikit learn

descr = """A set of python modules for modelling the effect of tsunamis and flooding"""

import sys
import os
import shutil
from distutils.command.clean import clean as Clean


#==============================================================================
DISTNAME = 'anuga'
DESCRIPTION = 'A set of python modules for tsunami and flood modelling'
with open('README.rst') as f:
    LONG_DESCRIPTION = f.read()
MAINTAINER = 'Stephen Roberts'
MAINTAINER_EMAIL = 'stephen.roberts@anu.edu.au'
URL = "https://github.com/anuga-community/anuga_core"
LICENSE = 'GPL'
DOWNLOAD_URL = "http://sourceforge.net/projects/anuga/"
VERSION = '3.1.9'
#===============================================================================



###############################################################################
# Optional setuptools features
# We need to import setuptools early, if we want setuptools features,
# as it monkey-patches the 'setup' function

# For some commands, use setuptools
SETUPTOOLS_COMMANDS = set([
    'develop', 'release', 'bdist_egg', 'bdist_rpm',
    'bdist_wininst', 'install_egg_info', 'build_sphinx',
    'egg_info', 'easy_install', 'upload', 'bdist_wheel',
    '--single-version-externally-managed',
])


if len(SETUPTOOLS_COMMANDS.intersection(sys.argv)) > 0:
    import setuptools
    extra_setuptools_args = dict(
        zip_safe=False,  # the package can run out of an .egg file
        include_package_data=True,
      	install_requires=['cython',
                         'numpy',
                         'dill',
                         'future',
                         'gdal >= 3.0.4'
                         'gitpython',
                         'matplotlib',
                         'meshpy',
                         'netcdf4',
                         'Pmw',
                         'pymetis',
                         'pytest',
                         'pytz',
                         'scipy',
                         'utm']
    )
else:
    extra_setuptools_args = dict()

###############################################################################

class CleanCommand(Clean):
    description = "Remove build artifacts from the source tree"

    def run(self):
        Clean.run(self)
        if os.path.exists('build'):
            shutil.rmtree('build')
        for dirpath, dirnames, filenames in os.walk('anuga'):
            for filename in filenames:
                if (filename.endswith('.so') or filename.endswith('.pyd')
                        or filename.endswith('.pyc')):
                    os.unlink(os.path.join(dirpath, filename))
            for dirname in dirnames:
                if dirname == '__pycache__':
                    shutil.rmtree(os.path.join(dirpath, dirname))


###############################################################################
def configuration(parent_package='', top_path=None):
    if os.path.exists('MANIFEST'):
        os.remove('MANIFEST')

    from numpy.distutils.misc_util import Configuration
    config = Configuration(None, parent_package, top_path)

    # Avoid non-useful msg:
    # "Ignoring attempt to set 'name' (from ... "
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('anuga')

    return config



def setup_package():

    metadata = dict(name=DISTNAME,
                    maintainer=MAINTAINER,
                    maintainer_email=MAINTAINER_EMAIL,
                    description=DESCRIPTION,
                    license=LICENSE,
                    url=URL,
                    version=VERSION,
                    download_url=DOWNLOAD_URL,
                    long_description=LONG_DESCRIPTION,
                    classifiers=['Intended Audience :: Science/Research',
                                 'Intended Audience :: Developers',
                                 'License :: OSI Approved',
                                 'Programming Language :: C',
                                 'Programming Language :: C++',
                                 'Programming Language :: Python',
                                 'Topic :: Software Development',
                                 'Topic :: Scientific/Engineering',
                                 'Operating System :: Microsoft :: Windows',
                                 'Operating System :: POSIX',
                                 'Operating System :: Unix',
                                 'Operating System :: MacOS',
                                 'Programming Language :: Python :: 3.7',
                                 'Programming Language :: Python :: 3.8',
                                 'Programming Language :: Python :: 3.9',
                                 'Programming Language :: Python :: 3.10',
                                 ],
                    cmdclass={'clean': CleanCommand},
                    **extra_setuptools_args)



    metadata['version'] = VERSION
    metadata['configuration'] = configuration

    from numpy.distutils.core import setup

    setup(**metadata)


if __name__ == "__main__":
    setup_package()
