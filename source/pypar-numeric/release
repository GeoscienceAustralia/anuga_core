#!/usr/bin/env python
#
# Release version as tagged by
#   cvs tag release-x-y-z

import sys, os, string, os.path, time

import pypar
try:
  date = pypar.__date__
except:
  print 'WARNING: Could not obtain __date__ from pypar.py'
  date = time.ctime()

rev = pypar.__version__
  

package = 'pypar'
destination = os.path.expanduser('~/public_html/pypar')

# FIXME: temporarily commented out. Need to rewrite for Subversion.
#
# Get tag to release
#
#if len(sys.argv) > 1:
#  release = sys.argv[1]
#else:
#  from popen2 import popen2 
#  output, input = popen2('cvs log 2>/dev/null')
#  input.close()
#  lines = output.readlines()
#
#  found = 0  
#  for line in lines:
#    #print line,
#    
#    if found:
#      tag, revision = line.split(':')
#      tag = tag.strip()
#      print 'Found tag %s corresponding to CVS revision %s\n'\
#            %(tag, revision.strip()) 
#      break
#
#    if string.find(line, 'symbolic names:') == 0:
#      found = 1
#
#      
#i = string.find(tag,'-')
#rev = string.replace( tag[i:], '-', '_' )



curdir = os.getcwd()
os.chdir('/tmp')
     
      
# Export a clean directory
#s = 'cvs export -r %s %s' %(tag, package)
s = 'svn export %s %s' %(curdir, package)
print s
os.system(s)

#import sys; sys.exit()


release_name  = package + "_" + rev
s = 'mv %s %s' %(package, release_name)
print s
os.system(s)

# Create installation tree

s = 'mkdir %s/lib %s/lib/pypar %s/examples' %((release_name,)*3)
print s
os.system(s)

#cleanup dev stuff that shouldn't go into package
s = '/bin/rm %s/install.py %s/compile.py %s/Makefile' %((release_name,)*3)
print s
os.system(s)

s = 'mv %s/pypar.py %s/lib/pypar' %((release_name,)*2)
print s
os.system(s)

s = 'mv %s/__init__.py %s/lib/pypar' %((release_name,)*2)
print s
os.system(s)

s = 'mv %s/mpiext.c %s/lib/pypar' %((release_name,)*2)
print s
os.system(s)

s = 'mv %s/*.py %s/examples' %((release_name,)*2)
print s
os.system(s)

s = 'mv %s/examples/setup.py %s' %((release_name,)*2)
print s
os.system(s)

s = 'mv %s/pytiming %s/examples' %((release_name,)*2)
print s
os.system(s)

s = 'mv %s/ring_example.py %s/examples' %((release_name,)*2)
print s
os.system(s)

s = 'mv %s/runpytiming %s/examples' %((release_name,)*2)
print s
os.system(s)

s = 'mv %s/ctiming.c %s/examples' %((release_name,)*2)
print s
os.system(s)

# Make tarball and copy to destination
s = 'tar cvfz %s.tgz %s' %(release_name, release_name)
print s
os.system(s)

s = 'cp %s.tgz %s' %(release_name, destination)
print s
os.system(s)

s = 'cp %s/lib/pypar/pypar.py %s' %(release_name, destination)
print s
os.system(s)

s = 'cp %s/lib/pypar/mpiext.c %s' %(release_name, destination)
print s
os.system(s)

s = 'cp %s/examples/pytiming %s' %(release_name, destination)
print s
os.system(s)

s = 'cp %s/examples/ctiming.c %s' %(release_name, destination)
print s
os.system(s)

s = 'cp %s/examples/ring_example.py %s' %(release_name, destination)
print s
os.system(s)

s = 'cp %s/README %s' %(release_name, destination)
print s
os.system(s)

s = 'cp %s/DOC %s' %(release_name, destination)
print s
os.system(s)

# Update web page
#
print 'Updating WEB page ' + destination
input = open(destination + '/' + 'index.src', 'r')
output = open(destination + '/' + 'index.php', 'w')

output.write('<!-- AUTOMATICALLY GENERATED - EDIT index.src instead -->\n')
for line in input.readlines():
  line = string.replace(line, '<date>', date) 
  output.write(string.replace(line,'<filename>',release_name+'.tgz'))
output.close()    

#os.system('mv %s.tgz /home/web/dm_web/software/%s' %(release_name, package))
#os.system('cp %s/README /home/web/dm_web/software/%s' %(release_name, package))
#os.system('cp %s/DOC /home/web/dm_web/software/%s' %(release_name, package))

# Make soft link
#
s = 'cd %s; rm pypar.tgz; ln -s %s.tgz pypar.tgz' %(destination, release_name)
print s
os.system(s)

#Cleanup
s = '/bin/rm -f -r %s' %release_name
print s
os.system(s)


