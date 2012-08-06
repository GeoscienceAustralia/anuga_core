"""Create a release of pypar

The release name will be compose of a major release code obtained from
pypar.py as well as the subversion version number.
For example

pypar-2.0.2_39

refers to pypar version 2.0.2 with subversion number 39.

The release files will be store locally in ~ (home) and also
copied to sourceforge where the release manager must log in
and create the Release using the File Releases interface.

This script assumes a Linux platform

"""

from os import sep, system, remove, popen, chdir, getcwd, listdir
from os.path import join
from tempfile import mktemp
from sys import platform, stdout



def get_revision_number():
    """Get the version number of the SVN
    NOTE: This requires that the command svn is on the system PATH
    (simply aliasing svn to the binary will not work)
    """

    # Error msg
    msg = 'Command "svn" is not '
    msg += 'recognised on the system PATH.\n\n'
    msg += 'Try to obtain the version info '
    msg += 'by using the command: "svn info".\n'
    msg += 'In this case, make sure svn is accessible on the system path. '
    msg += 'Simply aliasing svn to the binary will not work. '
      
    try:
        # The null stuff is so this section fails quitly.
        # This could cause the svn info command to fail due to
        # the redirection being bad on some platforms.
        # If that occurs then change this code.
        if platform[0:3] == 'win':
            system('svn up')
            fid = popen('svn info 2> null')
        else:
            system('svn up')            
            fid = popen('svn info 2>/dev/null')
	
    except:
        raise Exception(msg)
    else:
        #print 'Got version from svn'            
        version_info = fid.read()
      
        if version_info == '':
          raise Exception(msg)    
        else:
          pass
          print 'Got version from file'

            
    for line in version_info.split('\n'):
        if line.startswith('Revision:'):
            break

    fields = line.split(':')
    msg = 'Keyword "Revision" was not found anywhere in text: %s'\
          %version_info
    assert fields[0].startswith('Revision'), msg            

    try:
        revision_number = int(fields[1])
    except:
        msg = 'Revision number must be an integer. I got %s' %fields[1]
        msg += 'Check that the command svn is on the system path' 
        raise Exception(msg)                
        
    return revision_number



if __name__ == '__main__':
  
    if platform == 'win32':
        msg = 'This script is not written for Windows.'+\
              'Please run it on a Unix platform'
        raise Exception, msg

    # Get version number from file __metadat__
    # This can't be imported because source is a package
    # and mpiext needs to be compiled before anything can
    # be imported.
    # Instead we read it manually (until a better solution is found)

    fid = open('source/__metadata__.py')
    major_revision = None
    for line in fid.readlines():
        if line.startswith('__version__'):
            i = line.find('=')
            major_revision = str(line[i+1:].strip())[1:-1]
            
    if major_revision is None:
        raise 'No version was found'

    

    # line separator 
    lsep = '----------------------------------------------------------------------'

    # Get svn revision number and create
    # file with version info for release.
    # This will mean that the version currently checked out is
    # the one which will be released.

    svn_revision = get_revision_number()
    revision = '%s_%s' %(major_revision, svn_revision)
    print 'Creating pypar revision %s' %revision

    distro_filename = 'pypar-%s.tgz' %revision

    
    
    # Create area directory
    release_name = 'pypar_%s' %revision 
    release_dir = '~/%s' %release_name 
    s = 'mkdir %s' %release_dir
    try:
        print s    
        system(s)
    except:
        pass


    # Export a clean directory tree from the working copy to a temporary dir
    tmp_dir = mktemp()
    s = 'mkdir %s' %tmp_dir
    print s    
    system(s)
    
    distro_dir = join(tmp_dir, release_name) 
    s = 'mkdir %s' %distro_dir
    print s    
    system(s)    

    
    

    #-----------------------------
    # Get pypar source
    s = 'svn export -r %d source %s/source' %(svn_revision,
                                              distro_dir) 
    print s
    system(s)
    
    #-----------------------------
    # Copy license file to top dir   
    s = 'cp %s/source/LICENSE %s' %(distro_dir, distro_dir)     
    print s
    system(s)    
   
    
    #-----------------------------
    # Get demos
    s = 'svn export -r %d demos %s/demos' %(svn_revision,
                                            distro_dir) 
    print s
    system(s)

    
    #-----------------------------
    # Get documentation
    s = 'svn export -r %d documentation %s/documentation' %(svn_revision,
                                                            distro_dir) 
    print s
    system(s)
    


    # Zip everything up
    s = 'cd %s; tar cvfz %s *' %(tmp_dir, distro_filename)
    print s
    system(s)

    # Move distro to release area
    s = '/bin/mv %s/*.tgz %s' %(tmp_dir, release_dir) 
    print s
    system(s)

    # Clean up
    s = '/bin/rm -rf %s/pypar' %(tmp_dir) 
    print s
    system(s)


    #-----------------------------

    print 'Done'
    print
    print
    print lsep
    print 'The release files are in %s:' %release_dir
    system('ls -la %s' %release_dir)
    print lsep
    print
    print

    answer = raw_input('Do you want to upload this to sourceforge? Y/N [Y]')
    if answer.lower() != 'n':
        
        #------------------------------
        print 'Uploading to sourceforge'


        import os, os.path
        release_dir = os.path.expanduser(release_dir)
        os.chdir(release_dir)
        print 'Reading from', os.getcwd()

        s = 'rsync -avP -e ssh *.tgz uniomni@frs.sourceforge.net:uploads/'
        print s
        os.system(s)

        #from ftplib import FTP
        #ftp = FTP('upload.sourceforge.net')
        #print ftp.login() # Anonymous
        #print ftp.cwd('incoming')
        #
        #for filename in os.listdir('.'):
        #    print 'Uploading %s... ' %filename,
        #    stdout.flush()
        #
        #    fid=open(filename, 'rb')
        #    print ftp.storbinary('STOR %s' %filename, fid)
        #    fid.close()
        #
        #print 'FTP done'
        #print ftp.quit()

        print
        print lsep
        print '    ********************* NOTE *************************'
        print lsep
        print 'To complete this release you must log into'
        print 'http://sourceforge.net/projects/pypar as admin'
        print 'and complete the process by selecting File Releases '
        print 'in the admin menu there.'
        print lsep
        print
        print
        




