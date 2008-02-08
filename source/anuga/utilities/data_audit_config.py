"""Specify what extensions, directories and files should be ignored by
the data_audit process.

These will generally be specific to the software project.
"""


# Ignore source code files
extensions_to_ignore = ['.py','.c','.h', '.f'] #, '.gif', '.jpg', '.png']

# Ignore generated stuff 
extensions_to_ignore += ['.pyc', '.o', '.so', '~']
extensions_to_ignore += ['.aux', '.log', '.idx', 'ilg', '.ind',
                         '.bbl', '.blg']

# Ignore license files themselves
extensions_to_ignore += ['.lic']    


# Ignore certain other files
files_to_ignore = ['README.txt']

# Ignore directories
directories_to_ignore = ['anuga_work', 'pymetis', 'obsolete_code',
                         'anuga_parallel', 'anuga_viewer',
                         'planning', 'coding_standards',
                         'experimentation',
                         '.svn', 'misc', '.metadata']

