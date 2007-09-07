"""compile.py - compile Python C-extension

   Commandline usage: 
     python compile.py <filename>

   Usage from within Python:
     import compile
     compile.compile(<filename>,..)

   Ole Nielsen, Duncan Gray Oct 2001      
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


  libext = 'so' #Default extension (Unix)
  libs = ''
  version = sys.version[:3]
  
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

  elif sys.platform == 'cygwin':  #Cygwin (compilation same as linux)
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
    libext = 'dll'
    libs = '/lib/python%s/config/libpython%s.dll.a' %(version,version)
      
  elif sys.platform == 'win32':  #Windows
    if CC:
      compiler = CC
    else:  
      compiler = 'gcc.exe' #Some systems require this (a security measure?) 
    if LD:
      loader = LD
    else:  
      loader = compiler
    if SFLAG:
      sharedflag = SFLAG
    else:  
      sharedflag = 'shared'
     
    libext = 'dll'

    v = version.replace('.','')
    dllfilename = 'python%s.dll' %(v)
    libs = os.path.join(sys.exec_prefix,dllfilename)

      
      
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
  if sys.platform == 'win32':  #Windows
    python_include = os.path.join(sys.exec_prefix, 'include')    
  else:  
    python_include = os.path.join(os.path.join(sys.exec_prefix, 'include'),
                                  'python' + version)

  # Check existence of Python.h
  #
  headerfile = python_include + os.sep + 'Python.h'
  try:
    open(headerfile, 'r')
  except:
    raise """Did not find Python header file %s.
    Make sure files for Python C-extensions are installed. 
    In debian linux, for example, you need to install a
    package called something like python2.3-dev""" %headerfile



  #Add Python path + utilities to includelist (see ticket:31)
  #Assume there is only one 'utilities' dir under path dirs
  
  utilities_include_dir = None
  for pathdir in sys.path:

    utilities_include_dir = pathdir + os.sep + 'utilities'
    #print pathdir
    #print utilities_include_dir
    try:
      os.stat(utilities_include_dir)
    except OSError:
      pass
    else:
      #print 'Found %s to be used as include dir' %utilities_include_dir
      break

  # This is hacky since it
  # assumes the location of the compile_all that determines buildroot
  try:
    utilities_include_dir = buildroot + os.sep + "source" + os.sep + "anuga" \
                            + os.sep + 'utilities'
  except:
    # This will make compile work locally
    utilities_include_dir = '.'
    utilities_include_dir = '../utilities'


    
  try:
    os.stat(utilities_include_dir)
  except OSError: 
    utilities_include_dir = buildroot + os.sep + 'utilities'
  
  
  
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
      open(FN, 'r')
    except:
      raise Exception, "Could not open: " + FN

    if not object_files: root1 = root  #Remember first filename        
    object_files += root + '.o '  
  
  
    # Compile
    #
    if utilities_include_dir is None:    
      s = '%s -c %s -I"%s" -o "%s.o" -Wall -O3'\
          %(compiler, FN, python_include, root)
    else:
      if FN == "triangle.c" or FN == "mesh_engine_c_layer.c":
        s = '%s -c %s -I"%s" -I"%s" -o "%s.o" -O3 -DTRILIBRARY=1 -DNO_TIMER=1'\
            %(compiler, FN, python_include, utilities_include_dir, root)
      else:
        s = '%s -c %s -I"%s" -I"%s" -o "%s.o" -Wall -O3'\
            %(compiler, FN, python_include, utilities_include_dir, root)

    if os.name == 'posix' and os.uname()[4] == 'x86_64':
      #Extra flags for 64 bit architectures
      #Second clause will always fail on Win32 because uname is UNIX specific
      #but won't get past first clause

      #FIXME: Which one?
      #s += ' -fPIC'
      s += ' -fPIC -m64' 
      
      
    if verbose:
      print s
    else:
      s = s + ' 2> /dev/null' #Suppress errors
  
    try:
      err = os.system(s)
      if err != 0:
          raise 'Attempting to compile %s failed - please try manually' %FN 
    except:
      raise 'Could not compile %s - please try manually' %FN  

  
  # Make shared library (*.so or *.dll)
  if libs is "":
    s = '%s -%s %s -o %s.%s -lm' %(loader, sharedflag, object_files, root1, libext)
  else:
    s = '%s -%s %s -o %s.%s "%s" -lm' %(loader, sharedflag, object_files, root1, libext, libs)
  if verbose:
    print s
  else:
    s = s + ' 2> /dev/null' #Suppress warnings
  
  try:  
    err=os.system(s)
    if err != 0:        
        raise 'Atempting to link %s failed - please try manually' %root1     
  except:
    raise 'Could not link %s - please try manually' %root1
    

def can_use_C_extension(filename):
    """Determine whether specified C-extension
    can and should be used.
    """

    from anuga.config import use_extensions

    from os.path import splitext

    root, ext = splitext(filename)
    
    C=False
    if use_extensions:
        try:
            s = 'import %s' %root
            #print s
            exec(s)
        except:
            try:
                open(filename)
            except:
                msg = 'C extension %s cannot be opened' %filename
                print msg                
            else:    
                print '------- Trying to compile c-extension %s' %filename
            
                try:
                    compile(filename)
                except:
                    print 'WARNING: Could not compile C-extension %s'\
                          %filename
                else:
                    try:
                        exec('import %s' %root)
                    except:
                        msg = 'C extension %s seems to compile OK, '
                        msg += 'but it can still not be imported.'
                        raise msg
                    else:
                        C=True
        else:
            C=True
            
    if not C:
        pass
        print 'NOTICE: C-extension %s not used' %filename

    return C


if __name__ == '__main__':


  import sys, os
  from os.path import splitext
  
  if len(sys.argv) > 1:
      files = sys.argv[1:]
      for filename in files:
          root, ext = splitext(filename)

          if ext <> '.c':
              print 'WARNING (compile.py): Skipping %s. I only compile C-files.' %filename
      
  else:  
      #path = os.path.split(sys.argv[0])[0] or os.getcwd()
      path = '.'
      files = os.listdir(path)
      
      

  for filename in files:
      root, ext = splitext(filename)

      if ext == '.c':
          for x in ['.dll', '.so']:
              try:
                  os.remove(root + x)
              except:
                  pass

          print '---------------------------------------'      
          print 'Trying to compile c-extension %s in %s'\
                %(filename, os.getcwd())
          try:
            if filename == 'triang.c': 
              compile(['triang.c','triangle.c'])
            elif  filename == 'mesh_engine_c_layer.c': 
              compile(['mesh_engine_c_layer.c','triangle.c'])
            else:
              compile(filename)
          except Exception, e:
              print 'Could not compile C extension %s' %filename
          else:
              print 'C extension %s OK' %filename
          print    
        

