"""compile.py - compile Python C-extension

   Commandline usage: 
     python compile.py <filename>

   Usage from within Python:
     import compile
     compile.compile(<filename>,..)

   Ole Nielsen, Oct 2001      
"""   

def compile(FNs=None, CC=None, LD = None, SFLAG = None, verbose = 1):
  """compile(FNs=None, CC=None, LD = None, SFLAG = None):
  
     Compile FN(s) using compiler CC (e.g. mpicc), 
     Loader LD and shared flag SFLAG.
     If CC is absent use default compiler dependent on platform
     if LD is absent CC is used.
     if SFLAG is absent platform default is used
     FNs can be either one filename or a list of filenames
     In the latter case, the first will be used to name so file.
  """
  import os, string, sys, types
  
  # Input check
  #
  assert not FNs is None, "No filename provided"

  if not type(FNs) == types.ListType:
    FNs = [FNs]

  # Determine platform and compiler
  #
  if sys.platform == 'sunos5':  #Solaris
    if CC:
      compiler = CC
    else:  
      compiler = 'gcc'
    if LD:
      loader = LD
    else:  
      loader = compiler
    if SFLAG:
      sharedflag = SFLAG
    else:  
      sharedflag = 'G'
      
  elif sys.platform == 'osf1V5':  #Compaq AlphaServer
    if CC:
      compiler = CC
    else:  
      compiler = 'cc'
    if LD:
      loader = LD
    else:  
      loader = compiler
    if SFLAG:
      sharedflag = SFLAG
    else:  
      sharedflag = 'shared'    
      
  elif sys.platform == 'linux2':  #Linux
    if CC:
      compiler = CC
    else:  
      compiler = 'gcc'
    if LD:
      loader = LD
    else:  
      loader = compiler
    if SFLAG:
      sharedflag = SFLAG
    else:  
      sharedflag = 'shared'    
      
  elif sys.platform == 'darwin':  #Mac OS X:
    if CC:
      compiler = CC
    else:  
      compiler = 'cc'
    if LD:
      loader = LD
    else:  
      loader = compiler
    if SFLAG:
      sharedflag = SFLAG
    else:  
      sharedflag = 'bundle -flat_namespace -undefined suppress'
      
  else:
    if verbose: print "Unrecognised platform %s - revert to default"\
                %sys.platform
    
    if CC:
      compiler = CC
    else:  
      compiler = 'cc'
    if LD:
      loader = LD
    else:  
      loader = 'ld'
    if SFLAG:
      sharedflag = SFLAG
    else:  
      sharedflag = 'G'

            
  # Find location of include files
  #
  include = os.path.join(os.path.join(sys.exec_prefix, 'include'),
                     'python'+sys.version[:3])
  
  
  # Check existence of Python.h
  #
  headerfile = include + '/Python.h'
  try:
    open(headerfile, 'r')
  except:
    raise """Did not find Python header file %s.
    Make sure files for Python C-extensions are installed. 
    In debian linux, for example, you need to install a
    package called something like python2.1-dev""" %headerfile
  
  
  # Check filename(s)
  #
  object_files = ''
  for FN in FNs:        
    root, ext = os.path.splitext(FN)
    if ext == '':
      FN = FN + '.c'
    elif ext.lower() != '.c':
      raise Exception, "Unrecognised extension: " + FN
    
    try:
      open(FN,'r')
    except:    
      raise Exception, "Could not open: " + FN

    if not object_files: root1 = root  #Remember first filename        
    object_files += root + '.o '  
  
  
    # Compile
    #
    s = "%s -c %s -I%s -o %s.o" %(compiler, FN, include, root)
    if os.name == 'posix' and os.uname()[4] == 'x86_64':
      #Extra flags for 64 bit architectures
      #s += ' -fPIC -m64' #gcc
      s += ' -fPIC -tp amd64' #pgcc
      
    
    if verbose:
      print s
    else:
      s = s + ' 2> /dev/null' #Suppress errors
  
    try:
      os.system(s)
    except:
      raise 'Could not compile %s - please try manually' %FN  

  
  # Make shared library (*.so)
  s = "%s -%s %s -o %s.so" %(loader, sharedflag, object_files, root1)

  if os.name == 'posix' and os.uname()[4] == 'x86_64':
      pass
      #Extra flags for 64 bit architectures using Portland compilers
      s += ' -mcmodel=medium'
  
  if verbose:
    print s
  else:
    s = s + ' 2> /dev/null' #Suppress warnings
  
  try:  
    err=os.system(s)
  except:
    raise 'Could not link %s - please try manually' %root1
    



if __name__ == '__main__':

  import sys
  
  # Get file to compile
  #
  if len(sys.argv) < 2:
    print "Usage:"
    print "  python compile.py <filename>"
    sys.exit()

  compile(sys.argv[1], verbose=1)

