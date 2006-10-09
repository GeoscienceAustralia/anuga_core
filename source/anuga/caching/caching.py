# =============================================================================
# caching.py - Supervised caching of function results.
# Copyright (C) 1999, 2000, 2001, 2002 Ole Moller Nielsen
# Australian National University (1999-2003)
# Geoscience Australia (2003-present)
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License (http://www.gnu.org/copyleft/gpl.html)
#    for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
#
#
# Contact address: Ole.Nielsen@ga.gov.au
#
# Version 1.5.6 February 2002
# =============================================================================
 
"""Module caching.py - Supervised caching of function results.

Public functions:

cache(func,args) -- Cache values returned from func given args.
cachestat() --      Reports statistics about cache hits and time saved.
test() --       Conducts a basic test of the caching functionality.

See doc strings of individual functions for detailed documentation.
"""

# -----------------------------------------------------------------------------
# Initialisation code

# Determine platform
#
import os
if os.name in ['nt', 'dos', 'win32', 'what else?']:
  unix = 0
else:
  unix = 1

# Make default caching directory name
#
if unix:
  homedir = '~'
  CR = '\n'
else:
  homedir = 'c:'
  CR = '\r\n'  #FIXME: Not tested under windows
  
cachedir = homedir + os.sep + '.python_cache' + os.sep

# -----------------------------------------------------------------------------
# Options directory with default values - to be set by user
#

options = { 
  'cachedir': cachedir,  # Default cache directory 
  'maxfiles': 1000000,   # Maximum number of cached files
  'savestat': 1,         # Log caching info to stats file
  'verbose': 1,          # Write messages to standard output
  'bin': 1,              # Use binary format (more efficient)
  'compression': 1,      # Use zlib compression
  'bytecode': 0,         # Recompute if bytecode has changed
  'expire': 0            # Automatically remove files that have been accessed
                         # least recently
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def set_option(key, value):
  """Function to set values in the options directory.

  USAGE:
    set_option(key, value)

  ARGUMENTS:
    key --   Key in options dictionary. (Required)
    value -- New value for key. (Required)

  DESCRIPTION:
    Function to set values in the options directory.
    Raises an exception if key is not in options.
  """

  if options.has_key(key):
    options[key] = value
  else:
    raise KeyError(key)  # Key not found, raise an exception

# -----------------------------------------------------------------------------
# Function cache - the main routine

def cache(func, args=(), kwargs = {}, dependencies=None , cachedir=None,
          verbose=None, compression=None, evaluate=0, test=0, clear=0,
          return_filename=0):
  """Supervised caching of function results.

  USAGE:
    result = cache(func, args, kwargs, dependencies, cachedir, verbose,
                   compression, evaluate, test, return_filename)

  ARGUMENTS:
    func --            Function object (Required)
    args --            Arguments to func (Default: ())
    kwargs --          Keyword arguments to func (Default: {})    
    dependencies --    Filenames that func depends on (Default: None)
    cachedir --        Directory for cache files (Default: options['cachedir'])
    verbose --         Flag verbose output to stdout
                       (Default: options['verbose'])
    compression --     Flag zlib compression (Default: options['compression'])
    evaluate --        Flag forced evaluation of func (Default: 0)
    test --            Flag test for cached results (Default: 0)
    clear --           Flag delete cached results (Default: 0)    
    return_filename -- Flag return of cache filename (Default: 0)    

  DESCRIPTION:
    A Python function call of the form

      result = func(arg1,...,argn)

    can be replaced by

      from caching import cache
      result = cache(func,(arg1,...,argn))

  The latter form returns the same output as the former but reuses cached
  results if the function has been computed previously in the same context.
  'result' and the arguments can be simple types, tuples, list, dictionaries or
  objects, but not unhashable types such as functions or open file objects. 
  The function 'func' may be a member function of an object or a module.

  This type of caching is particularly useful for computationally intensive
  functions with few frequently used combinations of input arguments. Note that
  if the inputs or output are very large caching might not save time because
  disc access may dominate the execution time.

  If the function definition changes after a result has been cached it will be
  detected by examining the functions bytecode (co_code, co_consts,
  func_defualts, co_argcount) and it will be recomputed.

  LIMITATIONS:
    1 Caching uses the apply function and will work with anything that can be
      pickled, so any limitation in apply or pickle extends to caching. 
    2 A function to be cached should not depend on global variables
      as wrong results may occur if globals are changed after a result has
      been cached.

  -----------------------------------------------------------------------------
  Additional functionality:

  Keyword args
    Keyword arguments (kwargs) can be added as a dictionary of keyword: value
    pairs, following the syntax of the built-in function apply(). 
    A Python function call of the form
    
      result = func(arg1,...,argn, kwarg1=val1,...,kwargm=valm)    

    is then cached as follows

      from caching import cache
      result = cache(func,(arg1,...,argn), {kwarg1:val1,...,kwargm:valm})
    
    The default value of kwargs is {}  

  Explicit dependencies:
    The call
      cache(func,(arg1,...,argn),dependencies = <list of filenames>)
    Checks the size, creation time and modification time of each listed file.
    If any file has changed the function is recomputed and the results stored
    again.

  Specify caching directory:
    The call
      cache(func,(arg1,...,argn), cachedir = <cachedir>)
    designates <cachedir> where cached data are stored. Use ~ to indicate users
    home directory - not $HOME. The default is ~/.python_cache on a UNIX
    platform and c:/.python_cache on a Win platform.

  Silent operation:
    The call
      cache(func,(arg1,...,argn),verbose=0)
    suppresses messages to standard output.

  Compression:
    The call
      cache(func,(arg1,...,argn),compression=0)
    disables compression. (Default: compression=1). If the requested compressed
    or uncompressed file is not there, it'll try the other version.

  Forced evaluation:
    The call
      cache(func,(arg1,...,argn),evaluate=1)
    forces the function to evaluate even though cached data may exist.

  Testing for presence of cached result:
    The call
      cache(func,(arg1,...,argn),test=1)
    retrieves cached result if it exists, otherwise None. The function will not
    be evaluated. If both evaluate and test are switched on, evaluate takes
    precedence.
    
  Obtain cache filenames:
    The call    
      cache(func,(arg1,...,argn),return_filename=1)
    returns the hashed base filename under which this function and its
    arguments would be cached

  Clearing cached results:
    The call
      cache(func,'clear')
    clears all cached data for 'func' and
      cache('clear')
    clears all cached data.
 
    NOTE: The string 'clear' can be passed an *argument* to func using
      cache(func,('clear',)) or cache(func,tuple(['clear'])).

    New form of clear:
      cache(func,(arg1,...,argn),clear=1)
    clears cached data for particular combination func and args 
      
  """

  # Imports and input checks
  #
  import types, time, string

  if not cachedir:
    cachedir = options['cachedir']

  if verbose == None:  # Do NOT write 'if not verbose:', it could be zero.
    verbose = options['verbose']

  if compression == None:  # Do NOT write 'if not compression:',
                           # it could be zero.
    compression = options['compression']

  # Create cache directory if needed
  #
  CD = checkdir(cachedir,verbose)

  # Handle the case cache('clear')
  #
  if type(func) == types.StringType:
    if string.lower(func) == 'clear':
      clear_cache(CD,verbose=verbose)
      return

  # Handle the case cache(func, 'clear')
  #
  if type(args) == types.StringType:
    if string.lower(args) == 'clear':
      clear_cache(CD,func,verbose=verbose)
      return

  # Force singleton arg into a tuple
  #
  if type(args) != types.TupleType:
    args = tuple([args])
  
  # Check that kwargs is a dictionary
  #
  if type(kwargs) != types.DictType:
    raise TypeError    
    
  #print 'hashing' #FIXME: make faster hashing function
   
  # Hash arguments (and keyword args) to integer
  #
  arghash = myhash((args,kwargs))

  # Get sizes and timestamps for files listed in dependencies.
  # Force singletons into a tuple.
  #
  if dependencies and type(dependencies) != types.TupleType \
                  and type(dependencies) != types.ListType:
    dependencies = tuple([dependencies])
  deps = get_depstats(dependencies)

  # Extract function name from func object
  #
  funcname = get_funcname(func)

  # Create cache filename
  #
  FN = funcname+'['+`arghash`+']'  # The symbol '(' does not work under unix

  if return_filename:
    return(FN)

  if clear:
    for file_type in file_types:
      file_name = CD+FN+'_'+file_type
      for fn in [file_name, file_name + '.z']:
        if os.access(fn, os.F_OK):              
          if unix:
            os.remove(fn)
          else:
            # FIXME: os.remove doesn't work under windows        
            os.system('del '+fn)
          if verbose is True:
            print 'MESSAGE (caching): File %s deleted' %fn
        ##else:
        ##  print '%s was not accessed' %fn
    return None


  #-------------------------------------------------------------------        
  
  # Check if previous computation has been cached
  #
  if evaluate:
    Retrieved = None  # Force evaluation of func regardless of caching status.
    reason = 4
  else:
    (T, FN, Retrieved, reason, comptime, loadtime, compressed) = \
      CacheLookup(CD, FN, func, args, kwargs, deps, verbose, compression, \
                  dependencies)

  if not Retrieved:
    if test:  # Do not attempt to evaluate function
      T = None
    else:  # Evaluate function and save to cache
      if verbose:
        msg1(funcname, args, kwargs,reason)

      # Remove expired files automatically
      #
      if options['expire']:
        DeleteOldFiles(CD,verbose)
        
      # Save args before function is evaluated in case
      # they are modified by function
      #
      save_args_to_cache(CD,FN,args,kwargs,compression)

      # Execute and time function with supplied arguments
      #
      t0 = time.time()
      T = apply(func,args,kwargs)
      #comptime = round(time.time()-t0)
      comptime = time.time()-t0

      if verbose:
        msg2(funcname,args,kwargs,comptime,reason)

      # Save results and estimated loading time to cache
      #
      loadtime = save_results_to_cache(T, CD, FN, func, deps, comptime, \
                                       funcname, dependencies, compression)
      if verbose:
        msg3(loadtime, CD, FN, deps, compression)
      compressed = compression

  if options['savestat'] and (not test or Retrieved):
  ##if options['savestat']:
    addstatsline(CD,funcname,FN,Retrieved,reason,comptime,loadtime,compressed)

  return(T)  # Return results in all cases

# -----------------------------------------------------------------------------

def cachestat(sortidx=4, period=-1, showuser=None, cachedir=None):
  """Generate statistics of caching efficiency.

  USAGE:
    cachestat(sortidx, period, showuser, cachedir)

  ARGUMENTS:
    sortidx --  Index of field by which lists are (default: 4)
                Legal values are
                 0: 'Name'
                 1: 'Hits'
                 2: 'CPU'
                 3: 'Time Saved'
                 4: 'Gain(%)'
                 5: 'Size'
    period --   If set to -1 all available caching history is used.
                If set 0 only the current month is used (default -1).
    showuser -- Flag for additional table showing user statistics
                (default: None).
    cachedir -- Directory for cache files (default: options['cachedir']).

  DESCRIPTION:
    Logged caching statistics is converted into summaries of the form
    --------------------------------------------------------------------------
    Function Name   Hits   Exec(s)  Cache(s)  Saved(s)   Gain(%)      Size
    --------------------------------------------------------------------------
  """

  __cachestat(sortidx, period, showuser, cachedir)
  return

# -----------------------------------------------------------------------------

#Has mostly been moved to proper unit test
def test(cachedir=None,verbose=0,compression=None):
  """Test the functionality of caching.

  USAGE:
    test(verbose)

  ARGUMENTS:
    verbose --     Flag whether caching will output its statistics (default=0)
    cachedir --    Directory for cache files (Default: options['cachedir'])
    compression -- Flag zlib compression (Default: options['compression'])
  """
   
  import string, time

  # Initialise
  #
  import caching
  reload(caching)

  if not cachedir:
    cachedir = options['cachedir']

  if verbose is None:  # Do NOT write 'if not verbose:', it could be zero.
    verbose = options['verbose']
  
  if compression == None:  # Do NOT write 'if not compression:',
                           # it could be zero.
    compression = options['compression']
  else:
    try:
      set_option('compression', compression)
    except:
      test_error('Set option failed')      

  try:
    import zlib
  except:
    print
    print '*** Could not find zlib, default to no-compression      ***'
    print '*** Installing zlib will improve performance of caching ***'
    print
    compression = 0        
    set_option('compression', compression)    
  
  print  
  print_header_box('Testing caching module - please stand by')
  print    

  # Define a test function to be cached
  #
  def f(a,b,c,N,x=0,y='abcdefg'):
    """f(a,b,c,N)
       Do something time consuming and produce a complex result.
    """

    import string

    B = []
    for n in range(N):
      s = str(n+2.0/(n + 4.0))+'.a'*10
      B.append((a,b,c,s,n,x,y))
    return(B)
    
  # Check that default cachedir is OK
  #      
  CD = checkdir(cachedir,verbose)    
    
    
  # Make a dependency file
  #    
  try:
    DepFN = CD + 'testfile.tmp'
    DepFN_wildcard = CD + 'test*.tmp'
    Depfile = open(DepFN,'w')
    Depfile.write('We are the knights who say NI!')
    Depfile.close()
    test_OK('Wrote file %s' %DepFN)
  except:
    test_error('Could not open file %s for writing - check your environment' \
               % DepFN)

  # Check set_option (and switch stats off
  #    
  try:
    set_option('savestat',0)
    assert(options['savestat'] == 0)
    test_OK('Set option')
  except:
    test_error('Set option failed')    
    
  # Make some test input arguments
  #
  N = 5000  #Make N fairly small here

  a = [1,2]
  b = ('Thou shalt count the number three',4)
  c = {'Five is right out': 6, (7,8): 9}
  x = 3
  y = 'holy hand granate'

  # Test caching
  #
  if compression:
    comprange = 2
  else:
    comprange = 1

  for comp in range(comprange):
  
    # Evaluate and store
    #
    try:
      T1 = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, evaluate=1, \
                         verbose=verbose, compression=comp)
      if comp:                   
        test_OK('Caching evaluation with compression')
      else:     
        test_OK('Caching evaluation without compression')      
    except:
      if comp:
        test_error('Caching evaluation with compression failed - try caching.test(compression=0)')
      else:
        test_error('Caching evaluation failed - try caching.test(verbose=1)')

    # Retrieve
    #                           
    try:                         
      T2 = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
                         compression=comp) 

      if comp:                   
        test_OK('Caching retrieval with compression')
      else:     
        test_OK('Caching retrieval without compression')      
    except:
      if comp:
        test_error('Caching retrieval with compression failed - try caching.test(compression=0)')
      else:                                      
        test_error('Caching retrieval failed - try caching.test(verbose=1)')

    # Reference result
    #   
    T3 = f(a,b,c,N,x=x,y=y)  # Compute without caching
    
    if T1 == T2 and T2 == T3:
      if comp:
        test_OK('Basic caching functionality (with compression)')
      else:
        test_OK('Basic caching functionality (without compression)')
    else:
      test_error('Cached result does not match computed result')


  # Test return_filename
  #    
  try:
    FN = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
                       return_filename=1)    
    assert(FN[:2] == 'f[')
    test_OK('Return of cache filename')
  except:
    test_error('Return of cache filename failed')

  # Test existence of cachefiles
  #  
  try:
    (datafile,compressed0) = myopen(CD+FN+'_'+file_types[0],"rb",compression)
    (argsfile,compressed1) = myopen(CD+FN+'_'+file_types[1],"rb",compression)
    (admfile,compressed2) =  myopen(CD+FN+'_'+file_types[2],"rb",compression)
    test_OK('Presence of cache files')
    datafile.close()
    argsfile.close()
    admfile.close()
  except:
    test_error('Expected cache files did not exist') 
              
  # Test 'test' function when cache is present
  #      
  try:
    #T1 = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
    #                   evaluate=1)  
    T4 = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, test=1)
    assert(T1 == T4)

    test_OK("Option 'test' when cache file present")
  except:
    test_error("Option 'test' when cache file present failed")      

  # Test that 'clear' works
  #
  #try:
  #  caching.cache(f,'clear',verbose=verbose)
  #  test_OK('Clearing of cache files')
  #except:
  #  test_error('Clear does not work')
  try:
    caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, clear=1)    
    test_OK('Clearing of cache files')
  except:
    test_error('Clear does not work')  

  

  # Test 'test' function when cache is absent
  #      
  try:
    T4 = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, test=1)
    assert(T4 is None)
    test_OK("Option 'test' when cache absent")
  except:
    test_error("Option 'test' when cache absent failed")      
          
  # Test dependencies
  #
  T1 = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
                       dependencies=DepFN)  
  T2 = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
                       dependencies=DepFN)                     
                       
  if T1 == T2:
    test_OK('Basic dependencies functionality')
  else:
    test_error('Dependencies do not work')

  # Test basic wildcard dependency
  #
  T3 = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
                       dependencies=DepFN_wildcard)                     
    
  if T1 == T3:
    test_OK('Basic dependencies with wildcard functionality')
  else:
    test_error('Dependencies with wildcards do not work')


  # Test that changed timestamp in dependencies triggers recomputation
  
  # Modify dependency file
  Depfile = open(DepFN,'a')
  Depfile.write('You must cut down the mightiest tree in the forest with a Herring')
  Depfile.close()
  
  T3 = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
                       dependencies=DepFN, test = 1)                     
  
  if T3 is None:
    test_OK('Changed dependencies recognised')
  else:
    test_error('Changed dependencies not recognised')    
  
  # Test recomputation when dependencies have changed
  #
  T3 = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
                       dependencies=DepFN)                       
  if T1 == T3:
    test_OK('Recomputed value with changed dependencies')
  else:
    test_error('Recomputed value with changed dependencies failed')

  # Performance test (with statistics)
  # Don't really rely on this as it will depend on specific computer. 
  #

  set_option('savestat',1)

  N = 20*N   #Should be large on fast computers...
  tt = time.time()
  T1 = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose)
  t1 = time.time() - tt
  
  tt = time.time()
  T2 = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose)
  t2 = time.time() - tt
  
  if T1 == T2:
    if t1 > t2:
      test_OK('Performance test: relative time saved = %s pct' \
              %str(round((t1-t2)*100/t1,2)))
    #else:
    #  print 'WARNING: Performance a bit low - this could be specific to current platform'
  else:       
    test_error('Basic caching failed for new problem')
            
  # Test presence of statistics file
  #
  try: 
    DIRLIST = os.listdir(CD)
    SF = []
    for FN in DIRLIST:
      if string.find(FN,statsfile) >= 0:
        fid = open(CD+FN,'r')
        fid.close()
    test_OK('Statistics files present') 
  except:
    test_OK('Statistics files cannot be opened')          
      
  print_header_box('Show sample output of the caching function:')
  
  T2 = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=0)
  T2 = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=0)
  T2 = caching.cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=1)
  
  print_header_box('Show sample output of cachestat():')
  if unix:
    cachestat()    
  else:
    try:
      import time
      t = time.strptime('2030','%Y')
      cachestat()
    except:  
      print 'caching.cachestat() does not work here, because it'
      print 'relies on time.strptime() which is unavailable in Windows'
      
  print
  test_OK('Caching self test completed')    
      
            
  # Test setoption (not yet implemented)
  #

  
#==============================================================================
# Auxiliary functions
#==============================================================================

# Import pickler
# cPickle is used by functions mysave, myload, and compare
#
import cPickle  # 10 to 100 times faster than pickle
pickler = cPickle 

# Local immutable constants
#
comp_level = 1              # Compression level for zlib.
                            # comp_level = 1 works well.
textwidth1 = 16             # Text width of key fields in report forms.
textwidth2 = 132            # Maximal width of textual representation of
                            # arguments.
textwidth3 = 16             # Initial width of separation lines. Is modified.
textwidth4 = 50             # Text width in test_OK()
statsfile  = '.cache_stat'  # Basefilename for cached statistics.
                            # It will reside in the chosen cache directory.

file_types = ['Result',     # File name extension for cached function results.
              'Args',       # File name extension for stored function args.
              'Admin']      # File name extension for administrative info.

Reason_msg = ['OK',         # Verbose reasons for recomputation
              'No cached result', 
              'Dependencies have changed', 
              'Byte code or arguments have changed',
              'Recomputation was requested by caller',
	      'Cached file was unreadable']              
              
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def CacheLookup(CD, FN, func, args, kwargs, deps, verbose, compression, 
                dependencies):
  """Determine whether cached result exists and return info.

  USAGE:
    (T, FN, Retrieved, reason, comptime, loadtime, compressed) = \  
    CacheLookup(CD, FN, func, args, kwargs, deps, verbose, compression, \
                dependencies)

  INPUT ARGUMENTS:
    CD --            Cache Directory
    FN --            Suggested cache file name
    func --          Function object
    args --          Tuple of arguments
    kwargs --        Dictionary of keyword arguments    
    deps --          Dependencies time stamps
    verbose --       Flag text output
    compression --   Flag zlib compression
    dependencies --  Given list of dependencies
    
  OUTPUT ARGUMENTS:
    T --             Cached result if present otherwise None
    FN --            File name under which new results must be saved
    Retrieved --     True if a valid cached result was found
    reason --        0: OK (if Retrieved), 
                     1: No cached result, 
                     2: Dependencies have changed, 
                     3: Arguments or Bytecode have changed
                     4: Recomputation was forced
    comptime --      Number of seconds it took to computed cachged result
    loadtime --      Number of seconds it took to load cached result
    compressed --    Flag (0,1) if cached results were compressed or not 

  DESCRIPTION:
    Determine if cached result exists as follows:
    Load in saved arguments and bytecode stored under hashed filename.
    If they are identical to current arguments and bytecode and if dependencies
    have not changed their time stamp, then return cached result.

    Otherwise return filename under which new results should be cached.
    Hash collisions are handled recursively by calling CacheLookup again with a
    modified filename.
  """

  import time, string, types

  # Assess whether cached result exists - compressed or not.
  #
  if verbose:
    print 'Caching: looking for cached files %s_{%s,%s,%s}.z'\
           %(CD+FN, file_types[0], file_types[1], file_types[2])
  (datafile,compressed0) = myopen(CD+FN+'_'+file_types[0],"rb",compression)
  (argsfile,compressed1) = myopen(CD+FN+'_'+file_types[1],"rb",compression)
  (admfile,compressed2) =  myopen(CD+FN+'_'+file_types[2],"rb",compression)

  if not (argsfile and datafile and admfile) or \
     not (compressed0 == compressed1 and compressed0 == compressed2):
    # Cached result does not exist or files were compressed differently
    #
    # This will ensure that evaluation will take place unless all files are
    # present.

    reason = 1
    return(None,FN,None,reason,None,None,None) #Recompute using same filename

  compressed = compressed0  # Remember if compressed files were actually used
  datafile.close()

  # Retrieve arguments and adm. info
  #
  R, reason = myload(argsfile,compressed)  # The original arguments
  argsfile.close()
    
  ##if R == None and reason > 0:
  if reason > 0:
    return(None,FN,None,reason,None,None,None) #Recompute using same filename 
  else:   
    (argsref, kwargsref) = R

  R, reason = myload(admfile,compressed)
  admfile.close()  
  ##if R == None and reason > 0:
  if reason > 0:
    return(None,FN,None,reason,None,None,None) #Recompute using same filename 

  
  depsref  = R[0]  # Dependency statistics
  comptime = R[1]  # The computation time
  coderef  = R[2]  # The byte code
  funcname = R[3]  # The function name

  # Check if dependencies have changed
  #
  if dependencies and not compare(depsref,deps):
    if verbose:
      print 'MESSAGE (caching.py): Dependencies', dependencies, \
            'have changed - recomputing'
    # Don't use cached file - recompute
    reason = 2
    return(None,FN,None,reason,None,None,None)

  # Get bytecode from func
  #
  bytecode = get_bytecode(func)

  #print compare(argsref,args), 
  #print compare(kwargsref,kwargs),
  #print compare(bytecode,coderef)

  # Check if arguments or bytecode have changed
  #
  if compare(argsref,args) and compare(kwargsref,kwargs) and \
     (not options['bytecode'] or compare(bytecode,coderef)):

    # Arguments and dependencies match. Get cached results
    #
    T, loadtime, compressed, reason = load_from_cache(CD,FN,compressed)
    ###if T == None and reason > 0:  #This doesn't work if T is a numeric array
    if reason > 0:
      return(None,FN,None,reason,None,None,None) #Recompute using same FN 

    Retrieved = 1
    reason = 0

    if verbose:
      msg4(funcname,args,kwargs,deps,comptime,loadtime,CD,FN,compressed)

      if loadtime >= comptime:
        print 'WARNING (caching.py): Caching did not yield any gain.'
        print '                      Consider executing function ',
        print '('+funcname+') without caching.'
  else:

    # Non matching arguments or bytecodes signify a hash-collision.
    # This is resolved by recursive search of cache filenames
    # until either a matching or an unused filename is found.
    #
    (T,FN,Retrieved,reason,comptime,loadtime,compressed) = \
       CacheLookup(CD,FN+'x',func,args,kwargs,deps,verbose,compression, \
                   dependencies)

    # DEBUGGING
    # if not Retrieved:
    #   print 'Arguments did not match'
    # else:
    #   print 'Match found !'
    if not Retrieved:
      reason = 3     #The real reason is that args or bytecodes have changed.
                     #Not that the recursive seach has found an unused filename
    
  return((T, FN, Retrieved, reason, comptime, loadtime, compressed))

# -----------------------------------------------------------------------------

def clear_cache(CD,func=None, verbose=None):
  """Clear cache for func.

  USAGE:
     clear(CD, func, verbose)

  ARGUMENTS:
     CD --       Caching directory (required)
     func --     Function object (default: None)
     verbose --  Flag verbose output (default: None)

  DESCRIPTION:

    If func == None, clear everything,
    otherwise clear only files pertaining to func.
  """

  import os, re
   
  if CD[-1] != os.sep:
    CD = CD+os.sep
  
  if verbose == None:
    verbose = options['verbose']

  # FIXME: Windows version needs to be tested

  if func:
    funcname = get_funcname(func)
    if verbose:
      print 'MESSAGE (caching.py): Clearing', CD+funcname+'*'

    file_names = os.listdir(CD)
    for file_name in file_names:
      #RE = re.search('^' + funcname,file_name)  #Inefficient
      #if RE:
      if file_name[:len(funcname)] == funcname:
        if unix:
          os.remove(CD+file_name)
        else:
          os.system('del '+CD+file_name)
          # FIXME: os.remove doesn't work under windows
  else:
    file_names = os.listdir(CD)
    if len(file_names) > 0:
      if verbose:
        print 'MESSAGE (caching.py): Remove the following files:'
        for file_name in file_names:
            print file_name

        A = raw_input('Delete (Y/N)[N] ?')
      else:
        A = 'Y'	
	
      if A == 'Y' or A == 'y':
        for file_name in file_names:
          if unix:
            os.remove(CD+file_name)
          else:
            os.system('del '+CD+file_name)
            # FIXME: os.remove doesn't work under windows 
          #exitcode=os.system('/bin/rm '+CD+'* 2> /dev/null')

# -----------------------------------------------------------------------------

def DeleteOldFiles(CD,verbose=None):
  """Remove expired files

  USAGE:
    DeleteOldFiles(CD,verbose=None)
  """

  if verbose == None:
    verbose = options['verbose']

  maxfiles = options['maxfiles']

  # FIXME: Windows version

  import os
  block = 1000  # How many files to delete per invokation
  Files = os.listdir(CD)
  numfiles = len(Files)
  if not unix: return  # FIXME: Windows case ?

  if numfiles > maxfiles:
    delfiles = numfiles-maxfiles+block
    if verbose:
      print 'Deleting '+`delfiles`+' expired files:'
      os.system('ls -lur '+CD+'* | head -' + `delfiles`)            # List them
    os.system('ls -ur '+CD+'* | head -' + `delfiles` + ' | xargs /bin/rm')
                                                                  # Delete them
    # FIXME: Replace this with os.listdir and os.remove

# -----------------------------------------------------------------------------

def save_args_to_cache(CD,FN,args,kwargs,compression):
  """Save arguments to cache

  USAGE:
    save_args_to_cache(CD,FN,args,kwargs,compression)
  """

  import time, os, sys, types

  (argsfile, compressed) = myopen(CD+FN+'_'+file_types[1], 'wb', compression)

  if not argsfile:
    if verbose:
      print 'ERROR (caching): Could not open %s' %argsfile.name
    raise IOError

  mysave((args,kwargs),argsfile,compression)  # Save args and kwargs to cache
  argsfile.close()

  # Change access rights if possible
  #
  #if unix:
  #  try:
  #    exitcode=os.system('chmod 666 '+argsfile.name)
  #  except:
  #    pass
  #else:
  #  pass  # FIXME: Take care of access rights under Windows

  return

# -----------------------------------------------------------------------------

def save_results_to_cache(T, CD, FN, func, deps, comptime, funcname,
                          dependencies, compression):
  """Save computed results T and admin info to cache

  USAGE:
    save_results_to_cache(T, CD, FN, func, deps, comptime, funcname,
                          dependencies, compression)
  """

  import time, os, sys, types

  (datafile, compressed1) = myopen(CD+FN+'_'+file_types[0],'wb',compression)
  (admfile, compressed2) = myopen(CD+FN+'_'+file_types[2],'wb',compression)

  if not datafile:
    if verbose:
      print 'ERROR (caching): Could not open %s' %datafile.name
    raise IOError

  if not admfile:
    if verbose:
      print 'ERROR (caching): Could not open %s' %admfile.name
    raise IOError

  t0 = time.time()

  mysave(T,datafile,compression)  # Save data to cache
  datafile.close()
  #savetime = round(time.time()-t0,2)
  savetime = time.time()-t0  

  bytecode = get_bytecode(func)  # Get bytecode from function object
  admtup = (deps, comptime, bytecode, funcname)  # Gather admin info

  mysave(admtup,admfile,compression)  # Save admin info to cache
  admfile.close()

  # Change access rights if possible
  #
  #if unix:
  #  try:
  #    exitcode=os.system('chmod 666 '+datafile.name)
  #    exitcode=os.system('chmod 666 '+admfile.name)
  #  except:
  #    pass
  #else:
  #  pass  # FIXME: Take care of access rights under Windows

  return(savetime)

# -----------------------------------------------------------------------------

def load_from_cache(CD,FN,compression):
  """Load previously cached data from file FN

  USAGE:
    load_from_cache(CD,FN,compression)
  """

  import time

  (datafile, compressed) = myopen(CD+FN+'_'+file_types[0],"rb",compression)
  t0 = time.time()
  T, reason = myload(datafile,compressed)
  #loadtime = round(time.time()-t0,2)
  loadtime = time.time()-t0
  datafile.close() 

  return T, loadtime, compressed, reason

# -----------------------------------------------------------------------------

def myopen(FN,mode,compression=1):
  """Open file FN using given mode

  USAGE:
    myopen(FN,mode,compression=1)

  ARGUMENTS:
    FN --           File name to be opened
    mode --         Open mode (as in open)
    compression --  Flag zlib compression

  DESCRIPTION:
     if compression
       Attempt first to open FN + '.z'
       If this fails try to open FN
     else do the opposite
     Return file handle plus info about whether it was compressed or not.
  """

  import string

  # Determine if file exists already (if writing was requested)
  # This info is only used to determine if access modes should be set
  #
  if 'w' in mode or 'a' in mode:
    try:
      file = open(FN+'.z','r')
      file.close()
      new_file = 0
    except:
      try:
        file = open(FN,'r') 
        file.close()
        new_file = 0
      except:
        new_file = 1
  else:
    new_file = 0 #Assume it exists if mode was not 'w'
  

  compressed = 0
  if compression:
    try:
      file = open(FN+'.z',mode)
      compressed = 1
    except:
      try:
        file = open(FN,mode)
      except:
        file = None
  else:
    try:
      file = open(FN,mode)
    except:
      try:
        file = open(FN+'.z',mode)
        compressed = 1
      except:
        file = None

  # Now set access rights if it is a new file
  #
  if file and new_file:
    if unix:
      exitcode=os.system('chmod 666 '+file.name)
    else:
      pass  # FIXME: Take care of access rights under Windows

  return(file,compressed)

# -----------------------------------------------------------------------------

def myload(file, compressed):
  """Load data from file

  USAGE:
    myload(file, compressed)
  """

  reason = 0
  try:
    if compressed:
      import zlib

      RsC = file.read()
      try:
        Rs  = zlib.decompress(RsC)
      except:
        #  File "./caching.py", line 1032, in load_from_cache
        #  T = myload(datafile,compressed)
        #  File "./caching.py", line 1124, in myload
        #  Rs  = zlib.decompress(RsC)
        #  zlib.error: Error -5 while decompressing data
        #print 'ERROR (caching): Could not decompress ', file.name
        #raise Exception
        reason = 5  #(Unreadable file)
	return None, reason  
      
      
      del RsC  # Free up some space
      R   = pickler.loads(Rs)
    else:
      try:
        R = pickler.load(file)
      #except EOFError, e:
      except:
        #Catch e.g., file with 0 length or corrupted
        reason = 5  #(Unreadable file)
	return None, reason
      
  except MemoryError:
    import sys
    if options['verbose']:
      print 'ERROR (caching): Out of memory while loading %s, aborting' \
            %(file.name)

    # Raise the error again for now
    #
    raise MemoryError

  return R, reason

# -----------------------------------------------------------------------------

def mysave(T,file,compression):
  """Save data T to file

  USAGE:
    mysave(T,file,compression)

  """

  bin = options['bin']

  if compression:
    try:
      import zlib
    except:
      print
      print '*** Could not find zlib ***'
      print '*** Try to run caching with compression off ***'
      print "*** caching.set_option('compression', 0) ***"
      raise Exception
      

    try:
      Ts  = pickler.dumps(T, bin)
    except MemoryError:
      msg = '****WARNING (caching.py): Could not pickle data for compression.'
      msg += ' Try using compression = False'
      raise MemoryError, msg
    else:  
      #Compressed pickling      
      TsC = zlib.compress(Ts, comp_level)
      file.write(TsC)
  else:
      #Uncompressed pickling
      pickler.dump(T, file, bin)

      # FIXME: This may not work on Windoze network drives.
      # The error msg is IOError: [Errno 22] Invalid argument
      # Testing with small files was OK, though.
      # I think this is an OS problem.

      # Excerpt from http://www.ultraseek.com/support/faqs/4173.html
      
# The error is caused when there is a problem with server disk access (I/0). This happens at the OS level, and there is no controlling these errors through the Ultraseek application.
#
#Ultraseek contains an embedded Python interpreter. The exception "exceptions.IOError: [Errno 22] Invalid argument" is generated by the Python interpreter. The exception is thrown when a disk access operation fails due to an I/O-related reason.
#
#The following extract is taken from the site http://www.python.org:
#
#---------------------------------------------------------------------------------------------
#exception IOError
#Raised when an I/O operation (such as a print statement, the built-in open() function or a method of a file object) fails for an I/O-related reason, e.g., ``file not found'' or ``disk full''.
#This class is derived from EnvironmentError. See the discussion above for more information on exception instance attributes.
#---------------------------------------------------------------------------------------------
#
#The error code(s) that accompany exceptions are described at:
#http://www.python.org/dev/doc/devel//lib/module-errno.html
#
#You can view several postings on this error message by going to http://www.python.org, and typing the below into the search box:
#
#exceptions.IOError invalid argument Errno 22
        
      #try:
      #  pickler.dump(T,file,bin)
      #except IOError, e:
      #  print e
      #  msg = 'Could not store to %s, bin=%s' %(file, bin)
      #  raise msg 
      

# -----------------------------------------------------------------------------

def myhash(T):
  """Compute hashed integer from hashable values of tuple T

  USAGE:
    myhash(T)

  ARGUMENTS:
    T -- Tuple
  """

  import types

  # Get hash vals for hashable entries
  #
  if type(T) == types.TupleType or type(T) == types.ListType:
    hvals = []
    for k in range(len(T)):
      h = myhash(T[k])
      hvals.append(h)
    val = hash(tuple(hvals))
  elif type(T) == types.DictType:
    val = dicthash(T)
  else:
    try:
      val = hash(T)
    except:
      val = 1
      try:
        import Numeric
        if type(T) == Numeric.ArrayType:
          hvals = []        
          for e in T:
            h = myhash(e)
            hvals.append(h)          
          val = hash(tuple(hvals))
        else:
          val = 1  #Could implement other Numeric types here 
      except:    
        pass

  return(val)

# -----------------------------------------------------------------------------

def dicthash(D):
  """Compute hashed integer from hashable values of dictionary D

  USAGE:
    dicthash(D)
  """

  keys = D.keys()

  # Get hash values for hashable entries
  #
  hvals = []
  for k in range(len(keys)):
    try:
      h = hash(D[keys[k]])
      hvals.append(h)
    except:
      pass

  # Hash obtained values into one value
  #
  return(hash(tuple(hvals)))

# -----------------------------------------------------------------------------

def compare(A,B):
  """Safe comparison of general objects

  USAGE:
    compare(A,B)

  DESCRIPTION:
    Return 1 if A and B they are identical, 0 otherwise
  """

  try:
    identical = (A == B)
  except:
    try:
      identical = (pickler.dumps(A) == pickler.dumps(B))
    except:
      identical = 0

  return(identical)

# -----------------------------------------------------------------------------

def nospace(s):
  """Replace spaces in string s with underscores

  USAGE:
    nospace(s)

  ARGUMENTS:
    s -- string
  """

  import string

  newstr = ''
  for i in range(len(s)):
    if s[i] == ' ':
      newstr = newstr+'_'
    else:
      newstr = newstr+s[i]

  return(newstr)

# -----------------------------------------------------------------------------

def get_funcname(func):
  """Retrieve name of function object func (depending on its type)

  USAGE:
    get_funcname(func)
  """

  import types, string

  if type(func) == types.FunctionType:
    funcname = func.func_name
  elif type(func) == types.BuiltinFunctionType:
    funcname = func.__name__
  else:
    tab = string.maketrans("<>'","   ")
    tmp = string.translate(`func`,tab)
    tmp = string.split(tmp)
    funcname = string.join(tmp)

  funcname = nospace(funcname)
  return(funcname)

# -----------------------------------------------------------------------------

def get_bytecode(func):
  """ Get bytecode from function object.

  USAGE:
    get_bytecode(func)
  """

  import types

  if type(func) == types.FunctionType:
    bytecode = func.func_code.co_code
    consts = func.func_code.co_consts
    argcount = func.func_code.co_argcount    
    defaults = func.func_defaults     
  elif type(func) == types.MethodType:
    bytecode = func.im_func.func_code.co_code
    consts =  func.im_func.func_code.co_consts
    argcount =  func.im_func.func_code.co_argcount    
    defaults = func.im_func.func_defaults         
  else:
    #raise Exception  #Test only
    bytecode = None   #Built-in functions are assumed not to change
    consts = 0
    argcount = 0
    defaults = 0

  return (bytecode, consts, argcount, defaults)

# -----------------------------------------------------------------------------

def get_depstats(dependencies):
  """ Build dictionary of dependency files and their size, mod. time and ctime.

  USAGE:
    get_depstats(dependencies):
  """

  import types

  d = {}
  if dependencies:

    #Expand any wildcards
    import glob
    expanded_dependencies = []
    for FN in dependencies:
      expanded_FN = glob.glob(FN)
      
      expanded_dependencies += expanded_FN

    
    for FN in expanded_dependencies:
      if not type(FN) == types.StringType:
        errmsg = 'ERROR (caching.py): Dependency must be a string.\n'
        errmsg += '                    Dependency given: %s' %FN
        raise Exception, errmsg      
      if not os.access(FN,os.F_OK):
        errmsg = 'ERROR (caching.py): Dependency '+FN+' does not exist.'
        raise Exception, errmsg
      (size,atime,mtime,ctime) = filestat(FN)

      # We don't use atime because that would cause recomputation every time.
      # We don't use ctime because that is irrelevant and confusing for users.
      d.update({FN : (size,mtime)})

  return(d)

# -----------------------------------------------------------------------------

def filestat(FN):
  """A safe wrapper using os.stat to get basic file statistics
     The built-in os.stat breaks down if file sizes are too large (> 2GB ?)

  USAGE:
    filestat(FN)

  DESCRIPTION:
     Must compile Python with
     CFLAGS="`getconf LFS_CFLAGS`" OPT="-g -O2 $CFLAGS" \
              configure
     as given in section 8.1.1 Large File Support in the Libray Reference
  """

  import os, time

  try:
    stats = os.stat(FN)
    size  = stats[6]
    atime = stats[7]
    mtime = stats[8]
    ctime = stats[9]
  except:

    # Hack to get the results anyway (works only on Unix at the moment)
    #
    print 'Hack to get os.stat when files are too large'

    if unix:
      tmp = '/tmp/cach.tmp.'+`time.time()`+`os.getpid()`
      # Unique filename, FIXME: Use random number

      # Get size and access time (atime)
      #
      exitcode=os.system('ls -l --full-time --time=atime '+FN+' > '+tmp)
      (size,atime) = get_lsline(tmp)

      # Get size and modification time (mtime)
      #
      exitcode=os.system('ls -l --full-time '+FN+' > '+tmp)
      (size,mtime) = get_lsline(tmp)

      # Get size and ctime
      #
      exitcode=os.system('ls -l --full-time --time=ctime '+FN+' > '+tmp)
      (size,ctime) = get_lsline(tmp)

      try:
        exitcode=os.system('rm '+tmp)
        # FIXME: Gives error if file doesn't exist
      except:
        pass
    else:
      pass
      raise Exception  # FIXME: Windows case

  return(long(size),atime,mtime,ctime)

# -----------------------------------------------------------------------------

def get_lsline(FN):
  """get size and time for filename

  USAGE:
    get_lsline(file_name)

  DESCRIPTION:
    Read in one line 'ls -la' item from file (generated by filestat) and 
    convert time to seconds since epoch. Return file size and time.
  """

  import string, time

  f = open(FN,'r')
  info = f.read()
  info = string.split(info)

  size = info[4]
  week = info[5]
  mon  = info[6]
  day  = info[7]
  hour = info[8]
  year = info[9]

  str = week+' '+mon+' '+day+' '+hour+' '+year
  timetup = time.strptime(str)
  t = time.mktime(timetup)
  return(size, t)

# -----------------------------------------------------------------------------

def checkdir(CD,verbose=None):
  """Check or create caching directory

  USAGE:
    checkdir(CD,verbose):

  ARGUMENTS:
    CD -- Directory
    verbose -- Flag verbose output (default: None)

  DESCRIPTION:
    If CD does not exist it will be created if possible
  """

  import os
  import os.path

  if CD[-1] != os.sep: 
    CD = CD + os.sep  # Add separator for directories

  CD = os.path.expanduser(CD) # Expand ~ or ~user in pathname
  if not (os.access(CD,os.R_OK and os.W_OK) or CD == ''):
    try:
      exitcode=os.mkdir(CD)

      # Change access rights if possible
      #
      if unix:
        exitcode=os.system('chmod 777 '+CD)
      else:
        pass  # FIXME: What about acces rights under Windows?
      if verbose: print 'MESSAGE: Directory', CD, 'created.'
    except:
      print 'WARNING: Directory', CD, 'could not be created.'
      if unix:
        CD = '/tmp/'
      else:
        CD = 'C:'  
      print 'Using directory %s instead' %CD

  return(CD)

#==============================================================================
# Statistics
#==============================================================================

def addstatsline(CD,funcname,FN,Retrieved,reason,comptime,loadtime,
                 compression):
  """Add stats entry

  USAGE:
    addstatsline(CD,funcname,FN,Retrieved,reason,comptime,loadtime,compression)

  DESCRIPTION:
    Make one entry in the stats file about one cache hit recording time saved
    and other statistics. The data are used by the function cachestat.
  """

  import os, time

  try:
    TimeTuple = time.localtime(time.time())
    extension = time.strftime('%b%Y',TimeTuple)
    SFN = CD+statsfile+'.'+extension
    #statfile = open(SFN,'a')
    (statfile, dummy) = myopen(SFN,'a',compression=0)

    # Change access rights if possible
    #
    #if unix:
    #  try:
    #    exitcode=os.system('chmod 666 '+SFN)
    #  except:
    #    pass
  except:
    print 'Warning: Stat file could not be opened'

  try:
    if os.environ.has_key('USER'):
      user = os.environ['USER']
    else:
      user = 'Nobody'

    date = time.asctime(TimeTuple)

    if Retrieved:
      hit = '1'
    else:
      hit = '0'

    # Get size of result file
    #    
    if compression:
      stats = os.stat(CD+FN+'_'+file_types[0]+'.z')
    else:
      stats = os.stat(CD+FN+'_'+file_types[0])
  
    if stats: 
      size = stats[6]
    else:
      size = -1  # Error condition, but don't crash. This is just statistics  

    # Build entry
    #  
    entry = date             + ',' +\
            user             + ',' +\
            FN               + ',' +\
            str(int(size))   + ',' +\
            str(compression) + ',' +\
            hit              + ',' +\
            str(reason)      + ',' +\
            str(round(comptime,4)) + ',' +\
            str(round(loadtime,4)) +\
            CR
	    
    statfile.write(entry)
    statfile.close()
  except:
    print 'Warning: Writing of stat file failed'

# -----------------------------------------------------------------------------

# FIXME: should take cachedir as an optional arg
#
def __cachestat(sortidx=4,period=-1,showuser=None,cachedir=None):
  """  List caching statistics.

  USAGE:
    __cachestat(sortidx=4,period=-1,showuser=None,cachedir=None):

      Generate statistics of caching efficiency.
      The parameter sortidx determines by what field lists are sorted.
      If the optional keyword period is set to -1,
      all available caching history is used.
      If it is 0 only the current month is used.
      Future versions will include more than one month....
      OMN 20/8/2000
  """

  import os
  import os.path
  from string import split, rstrip, find, atof, atoi
  from time import strptime, localtime, strftime, mktime, ctime

  # sortidx = 4    # Index into Fields[1:]. What to sort by.

  Fields = ['Name', 'Hits', 'Exec(s)', \
            'Cache(s)', 'Saved(s)', 'Gain(%)', 'Size']
  Widths = [25,7,9,9,9,9,13]
  #Types = ['s','d','d','d','d','.2f','d']
  Types = ['s','d','.2f','.2f','.2f','.2f','d']  

  Dictnames = ['Function', 'User']

  if not cachedir:
    cachedir = checkdir(options['cachedir'])

  SD = os.path.expanduser(cachedir)  # Expand ~ or ~user in pathname

  if period == -1:  # Take all available stats
    SFILENAME = statsfile
  else:  # Only stats from current month  
       # MAKE THIS MORE GENERAL SO period > 0 counts several months backwards!
    TimeTuple = localtime(time())
    extension = strftime('%b%Y',TimeTuple)
    SFILENAME = statsfile+'.'+extension

  DIRLIST = os.listdir(SD)
  SF = []
  for FN in DIRLIST:
    if find(FN,SFILENAME) >= 0:
      SF.append(FN)

  blocksize = 15000000
  total_read = 0
  total_hits = 0
  total_discarded = 0
  firstday = mktime(strptime('2030','%Y'))
             # FIXME: strptime don't exist in WINDOWS ?
  lastday = 0

  FuncDict = {}
  UserDict = {}
  for FN in SF:
    input = open(SD+FN,'r')
    print 'Reading file ', SD+FN

    while 1:
      A = input.readlines(blocksize)
      if len(A) == 0: break
      total_read = total_read + len(A)
      for record in A:
        record = tuple(split(rstrip(record),','))
        #print record

        if len(record) in [8,9]:
          n = 0
          timestamp = record[n]; n=n+1
	
	  try:
            t = mktime(strptime(timestamp))
	  except:
            total_discarded = total_discarded + 1	  
	    continue    
	     
          if t > lastday:
            lastday = t
          if t < firstday:
            firstday = t

          user     = record[n]; n=n+1
          func     = record[n]; n=n+1

          # Strip hash-stamp off 
          #
          i = find(func,'[')
          func = func[:i]

          size        = atof(record[n]); n=n+1
          compression = atoi(record[n]); n=n+1
          hit         = atoi(record[n]); n=n+1
          reason      = atoi(record[n]); n=n+1   # Not used here    
          cputime     = atof(record[n]); n=n+1
          loadtime    = atof(record[n]); n=n+1

          if hit:
            total_hits = total_hits + 1
            saving = cputime-loadtime

            if cputime != 0:
              rel_saving = round(100.0*saving/cputime,2)
            else:
              #rel_saving = round(1.0*saving,2)
              rel_saving = 100.0 - round(1.0*saving,2)  # A bit of a hack

            info = [1,cputime,loadtime,saving,rel_saving,size]

            UpdateDict(UserDict,user,info)
            UpdateDict(FuncDict,func,info)
          else:
            pass #Stats on recomputations and their reasons could go in here
              
        else:
          #print 'Record discarded'
          #print record
          total_discarded = total_discarded + 1

    input.close()

  # Compute averages of all sums and write list
  #

  if total_read == 0:
    printline(Widths,'=')
    print 'CACHING STATISTICS: No valid records read'
    printline(Widths,'=')
    return

  print
  printline(Widths,'=')
  print 'CACHING STATISTICS: '+ctime(firstday)+' to '+ctime(lastday)
  printline(Widths,'=')
  #print '  Period:', ctime(firstday), 'to', ctime(lastday)
  print '  Total number of valid records', total_read
  print '  Total number of discarded records', total_discarded
  print '  Total number of hits', total_hits
  print

  print '  Fields', Fields[2:], 'are averaged over number of hits'
  print '  Time is measured in seconds and size in bytes'
  print '  Tables are sorted by', Fields[1:][sortidx]

  # printline(Widths,'-')

  if showuser:
    Dictionaries = [FuncDict, UserDict]
  else:
    Dictionaries = [FuncDict]

  i = 0
  for Dict in Dictionaries:
    for key in Dict.keys():
      rec = Dict[key]
      for n in range(len(rec)):
        if n > 0:
          rec[n] = round(1.0*rec[n]/rec[0],2)
      Dict[key] = rec

    # Sort and output
    #
    keylist = SortDict(Dict,sortidx)

    # Write Header
    #
    print
    #print Dictnames[i], 'statistics:'; i=i+1
    printline(Widths,'-')
    n = 0
    for s in Fields:
      if s == Fields[0]:  # Left justify
        s = Dictnames[i] + ' ' + s; i=i+1
        exec "print '%-" + str(Widths[n]) + "s'%s,"; n=n+1
      else:
        exec "print '%" + str(Widths[n]) + "s'%s,"; n=n+1
    print
    printline(Widths,'-')

    # Output Values
    #
    for key in keylist:
      rec = Dict[key]
      n = 0
      if len(key) > Widths[n]: key = key[:Widths[n]-3] + '...'
      exec "print '%-" + str(Widths[n]) + Types[n]+"'%key,";n=n+1
      for val in rec:
        exec "print '%" + str(Widths[n]) + Types[n]+"'%val,"; n=n+1
      print
    print

#==============================================================================
# Auxiliary stats functions
#==============================================================================

def UpdateDict(Dict,key,info):
  """Update dictionary by adding new values to existing.

  USAGE:
    UpdateDict(Dict,key,info)
  """

  if Dict.has_key(key):
    dinfo = Dict[key]
    for n in range(len(dinfo)):
      dinfo[n] = info[n] + dinfo[n]
  else:
    dinfo = info[:]  # Make a copy of info list

  Dict[key] = dinfo
  return Dict

# -----------------------------------------------------------------------------

def SortDict(Dict,sortidx=0):
  """Sort dictionary

  USAGE:
    SortDict(Dict,sortidx):

  DESCRIPTION:
    Sort dictionary of lists according field number 'sortidx'
  """

  import types

  sortlist  = []
  keylist = Dict.keys()
  for key in keylist:
    rec = Dict[key]
    if not type(rec) in [types.ListType, types.TupleType]:
      rec = [rec]

    if sortidx > len(rec)-1:
      if options['verbose']:
        print 'ERROR: Sorting index to large, sortidx = ', sortidx
      raise IndexError

    val = rec[sortidx]
    sortlist.append(val)

  A = map(None,sortlist,keylist)
  A.sort()
  keylist = map(lambda x: x[1], A)  # keylist sorted by sortidx

  return(keylist)

# -----------------------------------------------------------------------------

def printline(Widths,char):
  """Print textline in fixed field.

  USAGE:
    printline(Widths,char)
  """

  s = ''
  for n in range(len(Widths)):
    s = s+Widths[n]*char
    if n > 0:
      s = s+char

  print s

#==============================================================================
# Messages
#==============================================================================

def msg1(funcname,args,kwargs,reason):
  """Message 1

  USAGE:
    msg1(funcname,args,kwargs,reason):
  """

  import string
  #print 'MESSAGE (caching.py): Evaluating function', funcname,

  print_header_box('Evaluating function %s' %funcname)
  
  msg7(args,kwargs)
  msg8(reason)  
  
  print_footer()
  
  #
  # Old message
  #
  #args_present = 0
  #if args:
  #  if len(args) == 1:
  #    print 'with argument', mkargstr(args[0], textwidth2),
  #  else:
  #    print 'with arguments', mkargstr(args, textwidth2),
  #  args_present = 1      
  #    
  #if kwargs:
  #  if args_present:
  #    word = 'and'
  #  else:
  #    word = 'with'
  #      
  #  if len(kwargs) == 1:
  #    print word + ' keyword argument', mkargstr(kwargs, textwidth2)
  #  else:
  #    print word + ' keyword arguments', mkargstr(kwargs, textwidth2)
  #  args_present = 1            
  #else:
  #  print    # Newline when no keyword args present
  #        
  #if not args_present:   
  #  print '',  # Default if no args or kwargs present
    
    

# -----------------------------------------------------------------------------

def msg2(funcname,args,kwargs,comptime,reason):
  """Message 2

  USAGE:
    msg2(funcname,args,kwargs,comptime,reason)
  """

  import string

  #try:
  #  R = Reason_msg[reason]
  #except:
  #  R = 'Unknown reason'  
  
  #print_header_box('Caching statistics (storing) - %s' %R)
  print_header_box('Caching statistics (storing)')  
  
  msg6(funcname,args,kwargs)
  msg8(reason)

  print string.ljust('| CPU time:', textwidth1) + str(round(comptime,2)) + ' seconds'

# -----------------------------------------------------------------------------

def msg3(savetime, CD, FN, deps,compression):
  """Message 3

  USAGE:
    msg3(savetime, CD, FN, deps,compression)
  """

  import string
  print string.ljust('| Loading time:', textwidth1) + str(round(savetime,2)) + \
                     ' seconds (estimated)'
  msg5(CD,FN,deps,compression)

# -----------------------------------------------------------------------------

def msg4(funcname,args,kwargs,deps,comptime,loadtime,CD,FN,compression):
  """Message 4

  USAGE:
    msg4(funcname,args,kwargs,deps,comptime,loadtime,CD,FN,compression)
  """

  import string

  print_header_box('Caching statistics (retrieving)')
  
  msg6(funcname,args,kwargs)
  print string.ljust('| CPU time:', textwidth1) + str(round(comptime,2)) + ' seconds'
  print string.ljust('| Loading time:', textwidth1) + str(round(loadtime,2)) + ' seconds'
  print string.ljust('| Time saved:', textwidth1) + str(round(comptime-loadtime,2)) + \
        ' seconds'
  msg5(CD,FN,deps,compression)

# -----------------------------------------------------------------------------

def msg5(CD,FN,deps,compression):
  """Message 5

  USAGE:
    msg5(CD,FN,deps,compression)

  DESCRIPTION:
   Print dependency stats. Used by msg3 and msg4
  """

  import os, time, string

  print '|'
  print string.ljust('| Caching dir: ', textwidth1) + CD

  if compression:
    suffix = '.z'
    bytetext = 'bytes, compressed'
  else:
    suffix = ''
    bytetext = 'bytes'

  for file_type in file_types:
    file_name = FN + '_' + file_type + suffix
    print string.ljust('| ' + file_type + ' file: ', textwidth1) + file_name,
    stats = os.stat(CD+file_name)
    print '('+ str(stats[6]) + ' ' + bytetext + ')'

  print '|'
  if len(deps) > 0:
    print '| Dependencies:  '
    dependencies  = deps.keys()
    dlist = []; maxd = 0
    tlist = []; maxt = 0
    slist = []; maxs = 0
    for d in dependencies:
      stats = deps[d]
      t = time.ctime(stats[1])
      s = str(stats[0])
      #if s[-1] == 'L':
      #  s = s[:-1]  # Strip rightmost 'long integer' L off. 
      #              # FIXME: Unnecessary in versions later than 1.5.2

      if len(d) > maxd: maxd = len(d)
      if len(t) > maxt: maxt = len(t)
      if len(s) > maxs: maxs = len(s)
      dlist.append(d)
      tlist.append(t)
      slist.append(s)

    for n in range(len(dlist)):
      d = string.ljust(dlist[n]+':', maxd+1)
      t = string.ljust(tlist[n], maxt)
      s = string.rjust(slist[n], maxs)

      print '| ', d, t, ' ', s, 'bytes'
  else:
    print '| No dependencies'
  print_footer()

# -----------------------------------------------------------------------------

def msg6(funcname,args,kwargs):
  """Message 6

  USAGE:
    msg6(funcname,args,kwargs)
  """

  import string
  print string.ljust('| Function:', textwidth1) + funcname

  msg7(args,kwargs)
  
# -----------------------------------------------------------------------------    

def msg7(args,kwargs):
  """Message 7
  
  USAGE:
    msg7(args,kwargs):
  """
  
  import string
  
  args_present = 0  
  if args:
    if len(args) == 1:
      print string.ljust('| Argument:', textwidth1) + mkargstr(args[0], \
                         textwidth2)
    else:
      print string.ljust('| Arguments:', textwidth1) + \
            mkargstr(args, textwidth2)
    args_present = 1
            
  if kwargs:
    if len(kwargs) == 1:
      print string.ljust('| Keyword Arg:', textwidth1) + mkargstr(kwargs, \
                         textwidth2)
    else:
      print string.ljust('| Keyword Args:', textwidth1) + \
            mkargstr(kwargs, textwidth2)
    args_present = 1

  if not args_present:                
    print '| No arguments' # Default if no args or kwargs present

# -----------------------------------------------------------------------------

def msg8(reason):
  """Message 8
  
  USAGE:
    msg8(reason):
  """
  
  import string
    
  try:
    R = Reason_msg[reason]
  except:
    R = 'Unknown'  
  
  print string.ljust('| Reason:', textwidth1) + R
    
# -----------------------------------------------------------------------------

def print_header_box(line):
  """Print line in a nice box.
  
  USAGE:
    print_header_box(line)

  """
  global textwidth3

  import time

  time_stamp = time.ctime(time.time())
  line = time_stamp + '. ' + line
    
  N = len(line) + 1
  s = '+' + '-'*N + CR

  print s + '| ' + line + CR + s,

  textwidth3 = N

# -----------------------------------------------------------------------------
    
def print_footer():
  """Print line same width as that of print_header_box.
  """
  
  N = textwidth3
  s = '+' + '-'*N + CR    
      
  print s      
      
# -----------------------------------------------------------------------------

def mkargstr(args, textwidth, argstr = ''):
  """ Generate a string containing first textwidth characters of arguments.

  USAGE:
    mkargstr(args, textwidth, argstr = '')

  DESCRIPTION:
    Exactly the same as str(args) possibly followed by truncation,
    but faster if args is huge.
  """

  import types

  WasTruncated = 0

  if not type(args) in [types.TupleType, types.ListType, types.DictType]:
    if type(args) == types.StringType:
      argstr = argstr + "'"+str(args)+"'"
    else:
      #Truncate large Numeric arrays before using str()
      import Numeric
      if type(args) == Numeric.ArrayType:
        if len(args.flat) > textwidth:
          args = 'Array: ' + str(args.shape)

      argstr = argstr + str(args)
  else:
    if type(args) == types.DictType:
      argstr = argstr + "{"
      for key in args.keys():
        argstr = argstr + mkargstr(key, textwidth) + ": " + \
                 mkargstr(args[key], textwidth) + ", "
        if len(argstr) > textwidth:
          WasTruncated = 1
          break
      argstr = argstr[:-2]  # Strip off trailing comma      
      argstr = argstr + "}"

    else:
      if type(args) == types.TupleType:
        lc = '('
        rc = ')'
      else:
        lc = '['
        rc = ']'
      argstr = argstr + lc
      for arg in args:
        argstr = argstr + mkargstr(arg, textwidth) + ', '
        if len(argstr) > textwidth:
          WasTruncated = 1
          break

      # Strip off trailing comma and space unless singleton tuple
      #
      if type(args) == types.TupleType and len(args) == 1:
        argstr = argstr[:-1]    
      else:
        argstr = argstr[:-2]
      argstr = argstr + rc

  if len(argstr) > textwidth:
    WasTruncated = 1

  if WasTruncated:
    argstr = argstr[:textwidth]+'...'
  return(argstr)

# -----------------------------------------------------------------------------

def test_OK(msg):
  """Print OK msg if test is OK.
  
  USAGE
    test_OK(message)
  """

  import string
    
  print string.ljust(msg, textwidth4) + ' - OK' 
  
  #raise StandardError
  
# -----------------------------------------------------------------------------

def test_error(msg):
  """Print error if test fails.
  
  USAGE
    test_error(message)
  """
  
  print 'ERROR (caching.test): %s' %msg 
  print 'Please send this code example and output to '
  print 'Ole.Nielsen@anu.edu.au'
  print
  print
  
  #import sys
  #sys.exit()
  raise StandardError
