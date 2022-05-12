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
URL = "http://anuga.anu.edu.au"
LICENSE = 'GPL'
DOWNLOAD_URL = "http://sourceforge.net/projects/anuga/"
VERSION = '3.0.2'
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
      	install_requires=['pytest',
                          'numpy',
                          'scipy',
                          'netcdf4',
                          'matplotlib',
                          'gdal',
                          'dill',
                          'cython',
                          'future',   # FIXME(Ole): Deprecate
                          'mpi4py',
                          'Pmw',
                          'triangle',
                          'pymetis',
                          'gitpython'
                          'pytz',
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

def store_revision_info(destination_path='anuga', version=VERSION, verbose=False):
    """Obtain current revision from git and store it.

    Author: Stephen Roberts (stephen.roberts@anu.edu.au)

    CreationDate: August 2020

    Description:
        This function obtains current version from Git and stores it
        is a Python file named 'revision.py' for use with
        get_version_info()

        If git is not available on the system PATH, an Exception is thrown
    """
   

    # Git revision information (relies on the gitpython package)
    # https://stackoverflow.com/questions/14989858/get-the-current-git-hash-in-a-python-script
    try:
        import git
        repo = git.Repo(search_parent_directories=True)
    except:
        # Create dummy values for git revision info
        __git_sha__ = 'No git sha available'
        __git_committed_datetime__ = 'No git date available'
        __version__ = version

        msg = ('Could not import git module. ANUGA will still work, but will not store '
            'revision information in output file. You may need to install python git '
            'e.g. as pip install gitpython')
        #raise Warning(msg)  # I can't remember why does this cause ANUGA to stop instead of just issuing the warning (Ole)?
        print('WARNING', msg)
    else:
        __git_sha__ = repo.head.object.hexsha
        __git_committed_datetime__ = repo.head.object.committed_datetime
        __version__ = version

 
    if verbose: print('git_sha: ',__git_sha__)
    if verbose: print('git_committed_datetime: ',__git_committed_datetime__)
    if verbose: print('version: ',__version__ )

    # Determine absolute filename
    if destination_path[-1] != os.sep:
        destination_path += os.sep
        
    filename = destination_path + 'revision.py'

    fid = open(filename, 'w')

    docstring = 'Stored git revision info.\n\n'
    docstring += 'This file provides the git sha id and commit date for the installed '
    docstring += 'revision of ANUGA.\n'
    docstring += 'The file is automatically generated and should not '
    docstring += 'be modified manually.\n'
    fid.write('"""%s"""\n\n' %docstring)
    
    fid.write('__git_sha__ = "%s"\n'%__git_sha__)
    fid.write('__git_committed_datetime__ = "%s"\n'%__git_committed_datetime__)
    fid.write('__version__ = "%s"\n'%__version__)
    fid.close()


    if verbose is True:
        print('Revision info stored to %s' % filename)



def setup_package():

    # Update anuga/revision.py file
    store_revision_info()

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
                                 ],
                    cmdclass={'clean': CleanCommand},
                    **extra_setuptools_args)



    metadata['version'] = VERSION
    metadata['configuration'] = configuration

    from numpy.distutils.core import setup

    setup(**metadata)


if __name__ == "__main__":


    setup_package()
