"""Wrapper around data_audit.

This wrapper is specific to ANUGA.

This wrapper will run the

Specify what extensions, directories and files should be ignored by
the data_audit process.

These will generally be specific to each software project.
"""

from data_audit import IP_verified as IP_engine

# Ignore source code files
extensions_to_ignore = ['.py','.c', '.h', '.cpp', '.f', '.bat', '.m','.sh','.awk']

# Ignore LaTeX documents
extensions_to_ignore += ['.tex', '.sty', '.cls', '.bib', '.def']

# Ignore pdf and doc documents
extensions_to_ignore += ['.pdf', '.doc']

# Ignore generated stuff 
extensions_to_ignore += ['.pyc', '.o', '.so', '~']
extensions_to_ignore += ['.aux', '.log', '.idx', 'ilg', '.ind',
                         '.bbl', '.blg', '.syn', '.toc']

# Ignore license files themselves
extensions_to_ignore += ['.lic']    


# Ignore certain other files,
files_to_ignore = ['README.txt', 'LICENSE.txt', 'Makefile',
                   '.temp', 'SConstruct', 'SConscript', 'log.ini']

# Ignore directories
directories_to_ignore = ['.svn', 'misc', '.metadata']
#directories_to_ignore += ['old_pyvolution_documentation']



def IP_verified(directory, verbose=False):

    result = IP_engine(directory,
                       extensions_to_ignore,
                       directories_to_ignore,
                       files_to_ignore,
                       verbose=verbose)    
    return result
    
