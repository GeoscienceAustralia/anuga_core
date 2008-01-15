"""Remove computer generated garbage such as

   *.py~
   *.pyc
   *.o   
   *.so
   *.dll

Note: Recompile ANUGA after running this script
"""

import os


extensions_to_delete = ['~',
                        '.pyc',       # Python
                        '.o', '.so', '.dll',  # C
                        '.aux', '.ps',        # LaTeX
                        '.sww']

filenames_to_delete = []    
for dirpath, dirnames, filenames in os.walk('.'):

    print 'Searching dir', dirpath
   
    if '.svn' in dirnames:
        dirnames.remove('.svn')  # don't visit SVN directories


    for filename in filenames:
        for ext in extensions_to_delete:
            if filename.endswith(ext):
                absname = os.path.join(dirpath, filename)
                print '  Flagged for deletion', absname
                filenames_to_delete.append(absname)    


print 
N = len(filenames_to_delete)             
if N > 0:
    msg = '%d files flagged for deletion. Proceed? (Y/N)[N]' %N
    answer = raw_input(msg)
    
    if answer.lower() == 'y':
        for filename in filenames_to_delete:
            print 'Deleting', filename
            os.remove(filename)
    else:
        print 'Nothing deleted'
else:
    print 'No files flagged for deletion'            



