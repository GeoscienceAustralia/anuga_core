"""
Obtain current revision from git and store it.

Author: Stephen Roberts (stephen.roberts@anu.edu.au)

CreationDate: May 2022

Description:
        
This script obtains current version from setup.py and Git commit info
and stores it in a Python file named 'revision.py' for use with get_version_info()
"""

import re
import os


# ===================================================
# Read VERSION from setup.py file
# ===================================================
with open('setup.py') as infile:
    for line in infile:
        match = re.match(r'VERSION =', line)
        if match != None:
            VERSION = re.findall('\d.\d.\ddev|\d.\d.\d',line)[0]


destination_path='anuga'
version=VERSION
verbose=True
    
   
    

# Git revision information (relies on the gitpython package)
# https://stackoverflow.com/questions/14989858/get-the-current-git-hash-in-a-python-script
try:
    # Test if we are using a git repository
    import git
    repo = git.Repo(search_parent_directories=True)
except:
    msg = ('No git revision data available')
    #raise Warning(msg)  # I can't remember why does this cause ANUGA to stop instead of just issuing the warning (Ole)?
    raise(msg)
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

