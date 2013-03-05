#!/usr/bin/env python

"""Update, compile and create PDF and HTML from LaTeX file

Usage:
    python update_anuga_user_manual.py <options>
    
Options:
    --no_html: Skip automatic generation of html version
    


This script can for example be run from a cronjob:

  crontab -e

with content like 

SHELL=/bin/sh
PATH=/usr/local/sbin:/usr/local/bin:/sbin:/bin:/usr/sbin:/usr/bin:~/bin
PYTHONPATH=.:/home/ole/inundation/anuga_core/source:/home/ole/lib/python/site-packages:

# m h dom mon dow command
  32 6,10,14,18,22  *   *   *  ~/inundation/anuga_core/documentation/user_manual/update_anuga_user_manual.py > ~/inundation/anuga_core/documentation/user_manual/update_anuga.log
#
   
or 
 
# m h dom mon dow command
  42 *  *   *   *  ~/inundation/anuga_core/documentation/user_manual/update_anuga_user_manual.py > ~/inundation/anuga_core/documentation/user_manual/update_anuga.log
#   
   
Check function of crontab by reading mail using e.g. mutt   
   
   
Note UNIX only
 
"""

from os import system, chdir
from os.path import expanduser, split, join
from anuga.utilities.system_tools import get_revision_number, get_pathname_from_package
from anuga.config import major_revision
from sys import argv

# Determine absolute path for user manual

# Path for ANUGA
#anugapath = get_pathname_from_package('anuga')

# Strip trailing source/anuga of path
#basepath = split(split(anugapath)[0])[0]

# Add local path to user_manual
#docpath = join(join(basepath, 'documentation'), 'user_manual')
texfiles = ['anuga_user_manual', 
            'anuga_installation_guide',
            'anuga_whats_new',
            'anuga_internal_tools']

#print 'Moving to', docpath
#chdir(docpath) # Move to location of LaTeX files 
system('svn update') # Update from svn


do_html = False       
if len(argv) > 1:
    if argv[1] == '--no_html':
        do_html = False
    else:
        msg = 'Unknown option: %s' %argv[1]
        raise Exception(msg)

# Update version info
fid = open('version.tex')
lines = []
for line in fid.readlines():
    if line.startswith('\\release'):
        line = '\\release{%s}\n' %(major_revision)
            
    lines.append(line)
fid.close()

fid = open('version.tex', 'w')
fid.writelines(lines)
fid.close()

print 'Updated version info:'
for line in lines:    
    print line.strip()

for texfile in texfiles:
    print 'Processing %s' %texfile
    # Compile with LaTeX, makeindex etc
    for i in range(3):
        #system('latex --interaction=nonstopmode %s.tex' %texfile)
        system('pdflatex --interaction=nonstopmode -file-line-error %s.tex' %texfile)
        system('makeindex %s.idx' %texfile)
        system('makeindex mod%s.idx' %texfile)
        system('bibtex %s' %texfile)    

    # Create pdf file
    #system('dvips %s -o %s.ps' %((texfile,)*2))    
    #system('ps2pdf %s.ps' %texfile)    
     
    # Create html pages
    if do_html is True:
        system('latex2html %s' %texfile)
    else:
        print 'Skipping html version for %s as requested' %texfile

# Clean-up
#system('/bin/rm version.tex')
system('svn update') # Restore version file


print 'Cleanup aux tex files'
system('rm *.aux *.bbl *.blg *.idx *.ilg *.ind *.log *.out *.toc *.syn ')


# Print
print 'User manual compiled'

print 'Anuga version:', major_revision
print 'Build:', get_revision_number()
system('date')

