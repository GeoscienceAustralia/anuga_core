"""Obtain the latest revision number and date from Subversion. 

Using the tool SubWCRev.exe on the download page on the SVN website or in the TortoiseSVN Folder under bin. 

To create run 

#SubWCRev.exe path\to\working\copy version.in version.h
SubWCRev.exe . ver.txt ver.py

FIXME: Works only under Win32

Ex:
version = 1798
status = 'Modified'
date = '2005/09/07 13:22:59'

"""



template = 'version.txt'
version = 'version.py'

#Write out template
txt = """version = $WCREV$
status = '$WCMODS?Modified:Not modified$'
date = '$WCDATE$'
"""

fid = open(template, 'w')
fid.write(txt)
fid.close()

#Run conversion
import os
cmd = 'SubWCRev.exe . %s %s' %(template, version)
#print cmd
err = os.system(cmd)
if err != 0:
    msg = 'Command %s could not execute.'
    msg += 'Make sure the program SubWCRev.exe is available on your path'
    raise Exception(msg)



#Obtain version
cmd = 'from %s import version, status, date' %version[:-3]
#print cmd
exec(cmd)

print 'Version: %d' %version
print 'Date: %s' %date
print 'Status: %s' %status 


