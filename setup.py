#!/usr/bin/env python
""" ANUGA models the effect of tsunamis and flooding upon a terrain mesh.
    In typical usage, a Domain class is created for a particular piece of
    terrain. Boundary conditions are specified for the domain, such as inflow
    and outflow, and then the simulation is run.

    This is the public API to ANUGA. It provides a toolkit of often-used
    modules, which can be used directly by including the following line in
    the user's code:

    import anuga
        
    This usage pattern abstracts away the internal heirarchy of the ANUGA
    system, allowing the user to concentrate on writing simulations without
    searching through the ANUGA source tree for the functions that they need.
    
    Also, it isolates the user from "under-the-hood" refactorings.
"""

from __future__ import division, print_function

DOCLINES = __doc__.split("\n")

import os
import sys
import subprocess


if sys.version_info[:2] < (2, 6) or (3, 0) <= sys.version_info[0:2] :
    raise RuntimeError("Python version 2.6, 2.7. ")

# if sys.version_info[0] >= 3:
#     import builtins
# else:
#     import __builtin__ as builtins
import __builtin__ as builtins


CLASSIFIERS = """\
Development Status :: 2 - Development
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: C
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

MAJOR               = 1
MINOR               = 3
MICRO               = 1
ISRELEASED          = False
VERSION             = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


# Return the svn revision as a string
def svn_revision():

    return filter(str.isdigit, "$Revision$")


# Return the git revision as a string
def git_revision():

    #return "Unknown"

    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION

# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'): os.remove('MANIFEST')


# This is the numpy/scipy hack: Set a global variable so that the main
# anuga __init__ can detect if it is being loaded by the setup routine, to
# avoid attempting to load components that aren't built yet.
builtins.__ANUGA_SETUP__ = True


def get_version_info():
    # Adding the git rev number needs to be done inside write_version_py(),
    # otherwise the import of anuga.version messes up the build under Python 3.
    FULLVERSION = VERSION
    SVN_REVISION = svn_revision()
    
    if os.path.exists('.git'):
        GIT_REVISION = git_revision()
    elif os.path.exists('anuga/version.py'):
        # must be a source distribution, use existing version file
        try:
            from anuga.version import git_revision as GIT_REVISION
        except ImportError:
            raise ImportError("Unable to import git_revision. Try removing " \
                              "anuga/version.py and the build directory " \
                              "before building.")
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev+' + SVN_REVISION
        #FULLVERSION += '.dev+' + GIT_REVISION[:7]
         

    return FULLVERSION, GIT_REVISION, SVN_REVISION


def write_version_py(filename='anuga/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM ANUGA SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
svn_revision = '%(svn_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    FULLVERSION, GIT_REVISION, SVN_REVISION = get_version_info()

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version' : FULLVERSION,
                       'git_revision' : GIT_REVISION,
                       'svn_revision' : SVN_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()


def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)

    config.add_subpackage('anuga')

    config.get_version('anuga/version.py') # sets config.version

    return config



    
def setup_package():
    src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    # Rewrite the version file everytime
    write_version_py()

    metadata = dict(
        name = 'anuga',
        maintainer = "Anuga Developers",
        maintainer_email = "anuga-user@lists.sourceforge.net",
        description = DOCLINES[0],
        long_description = "\n".join(DOCLINES[2:]),
        url = "http://anuga.anu.edu.au",
        author = "Stephen Roberts, Ole Nielsen et al.",
        download_url = "http://sourceforge.net/projects/anuga/",
        license = 'GPL',
        classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
        platforms = ["Windows", "Linux", "Mac OS-X", "Unix"]
    )

    # Run build
    if len(sys.argv) >= 2 and ('--help' in sys.argv[1:] or
            sys.argv[1] in ('--help-commands', 'egg_info', '--version',
                            'clean')):
        # Use setuptools for these commands (they don't work well or at all
        # with distutils).  For normal builds use distutils.
        try:
            from setuptools import setup
        except ImportError:
            from distutils.core import setup

        FULLVERSION, GIT_REVISION = get_version_info()
        metadata['version'] = FULLVERSION
    else:
        if len(sys.argv) >= 2 and sys.argv[1] == 'bdist_wheel':
            # bdist_wheel needs setuptools
            import setuptools
        from numpy.distutils.core import setup
        cwd = os.path.abspath(os.path.dirname(__file__))
        print(cwd)
        if not os.path.exists(os.path.join(cwd, 'PKG-INFO')):
            # Generate Cython sources, unless building from source release
            #generate_cython()
            pass
        metadata['configuration'] = configuration

    try:
        setup(**metadata)
    finally:
        del sys.path[0]
        os.chdir(old_path)
    return


if __name__ == '__main__':
    setup_package()
