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

cache(my_F,args) -- Cache values returned from callable object my_F given args.
cachestat() --      Reports statistics about cache hits and time saved.
test() --       Conducts a basic test of the caching functionality.

See doc strings of individual functions for detailed documentation.
"""


# -----------------------------------------------------------------------------
# Initialisation code

# Determine platform
#



from os import getenv
import collections
import inspect
import types
import time
import sys

import os
if os.name in ['nt', 'dos', 'win32', 'what else?']:
  unix = False
else:
  unix = True

import anuga.utilities.log as log
from anuga.utilities import system_tools

import numpy as num

cache_dir = '.python_cache'

# Make default caching directory name
# We are changing the 'data directory' environment variable from
# INUNDATIONHOME to ANUGADATA - this gives a changeover.
if unix:
    homedir = getenv('ANUGADATA')
    if not homedir:
        homedir = getenv('INUNDATIONHOME')

    if not homedir:
        homedir = '~'
    else:
        # Since homedir will be a group area, individually label the caches
        user = getenv('LOGNAME')
        if not user:
            cache_dir += '_' + user
    
    CR = '\n'
else:
    homedir = 'c:'
    CR = '\r\n'  #FIXME: Not tested under windows
  
cachedir = os.path.join(homedir, cache_dir)

# It turns out hashes are no longer stable under Python3 (grr).
# https://stackoverflow.com/questions/27522626/hash-function-in-python-3-3-returns-different-results-between-sessions
# https://stackoverflow.com/questions/30585108/disable-hash-randomization-from-within-python-program
#
# The fix is to use another hashing library.
if system_tools.major_version == 3:
    import hashlib
    def hash(x):
        res = hashlib.sha256(str(x).encode()).hexdigest()
        #print('MY:', x, res)
        
        return res

# -----------------------------------------------------------------------------
# Options directory with default values - to be set by user
#

options = { 
  'cachedir': cachedir,  # Default cache directory 
  'maxfiles': 1000000,   # Maximum number of cached files
  'savestat': True,      # Log caching info to stats file
  'verbose': True,       # Write messages to standard output
  'bin': True,           # Use binary format (more efficient)
  'compression': True,   # Use zlib compression
  'bytecode': True,      # Recompute if bytecode has changed
  'expire': False        # Automatically remove files that have been accessed
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

  if key in options:
    options[key] = value
  else:
    raise KeyError(key)  # Key not found, raise an exception

# -----------------------------------------------------------------------------
# Function cache - the main routine

def cache(my_F, 
          args=(), 
          kwargs={}, 
          dependencies=None, 
          cachedir=None,
          verbose=None, 
          compression=None, 
          evaluate=False, 
          test=False, 
          clear=False,
          return_filename=False):
  """Supervised caching of function results. Also known as memoization.

  USAGE:
    result = cache(my_F, args, kwargs, dependencies, cachedir, verbose,
                   compression, evaluate, test, return_filename)

  ARGUMENTS:
    my_F --            Callable object (Required)
    args --            Arguments to my_F (Default: ())
    kwargs --          Keyword arguments to my_F (Default: {})    
    dependencies --    Filenames that my_F depends on (Default: None)
    cachedir --        Directory for cache files (Default: options['cachedir'])
    verbose --         Flag verbose output to stdout
                       (Default: options['verbose'])
    compression --     Flag zlib compression (Default: options['compression'])
    evaluate --        Flag forced evaluation of my_F (Default: False)
    test --            Flag test for cached results (Default: False)
    clear --           Flag delete cached results (Default: False)    
    return_filename -- Flag return of cache filename (Default: False)    

  DESCRIPTION:
    A Python function call of the form

      result = my_F(arg1,...,argn)

    can be replaced by

      from caching import cache
      result = cache(my_F,(arg1,...,argn))

  The latter form returns the same output as the former but reuses cached
  results if the function has been computed previously in the same context.
  'result' and the arguments can be simple types, tuples, list, dictionaries or
  objects, but not unhashable types such as functions or open file objects. 
  The function 'my_F' may be a member function of an object or a module.

  This type of caching is particularly useful for computationally intensive
  functions with few frequently used combinations of input arguments. Note that
  if the inputs or output are very large caching might not save time because
  disc access may dominate the execution time.

  If the function definition changes after a result has been cached it will be
  detected by examining the functions bytecode (co_code, co_consts,
  func_defaults, co_argcount) and it will be recomputed.

  LIMITATIONS:
    1 Caching uses function(*args, **kwargs) to evaluate and will work 
      with anything that can be pickled, so any limitation in function(,)
      or pickle extends to caching. 
    2 A function to be cached should not depend on global variables
      as wrong results may occur if globals are changed after a result has
      been cached.

  -----------------------------------------------------------------------------
  Additional functionality:

  Keyword args
    Keyword arguments (kwargs) can be added as a dictionary of keyword: value
    pairs, following Python's 'extended call syntax'. 
    
    A Python function call of the form
    
      result = my_F(arg1,...,argn, kwarg1=val1,...,kwargm=valm)    

    is then cached as follows

      from caching import cache
      result = cache(my_F,(arg1,...,argn), {kwarg1:val1,...,kwargm:valm})
    
    The default value of kwargs is {}  

  Explicit dependencies:
    The call
      cache(my_F,(arg1,...,argn), dependencies = <list of filenames>)
    Checks the size, creation time and modification time of each listed file.
    If any file has changed the function is recomputed and the results stored
    again.

  Specify caching directory:
    The call
      cache(my_F,(arg1,...,argn), cachedir = <cachedir>)
    designates <cachedir> where cached data are stored. Use ~ to indicate users
    home directory - not $HOME. The default is ~/.python_cache on a UNIX
    platform and c:/.python_cache on a Win platform.

  Silent operation:
    The call
      cache(my_F,(arg1,...,argn), verbose=False)
    suppresses messages to standard output.

  Compression:
    The call
      cache(my_F,(arg1,...,argn), compression=False)
    disables compression. (Default: compression=True). If the requested compressed
    or uncompressed file is not there, it'll try the other version.

  Forced evaluation:
    The call
      cache(my_F,(arg1,...,argn), evaluate=True)
    forces the function to evaluate even though cached data may exist.

  Testing for presence of cached result:
    The call
      cache(my_F,(arg1,...,argn), test=True)
    retrieves cached result if it exists, otherwise None. The function will not
    be evaluated. If both evaluate and test are switched on, evaluate takes
    precedence.
    ??NOTE: In case of hash collisions, this may return the wrong result as
    ??it only checks if *a* cached result is present. 
    # I think this was due to the bytecode option being False for some reason. (23/1/2009).
    
  Obtain cache filenames:
    The call    
      cache(my_F,(arg1,...,argn), return_filename=True)
    returns the hashed base filename under which this function and its
    arguments would be cached

  Clearing cached results:
    The call
      cache(my_F,'clear')
    clears all cached data for 'my_F' and
      cache('clear')
    clears all cached data.
 
    NOTE: The string 'clear' can be passed an *argument* to my_F using
      cache(my_F,('clear',)) or cache(my_F,tuple(['clear'])).

    New form of clear:
      cache(my_F,(arg1,...,argn), clear=True)
    clears cached data for particular combination my_F and args 
      
  """

  # Imports and input checks
  #
  import time, string

  if not cachedir:
    cachedir = options['cachedir']

  if verbose is None:  # Do NOT write 'if not verbose:', it could be zero.
    verbose = options['verbose']

  if compression is None:  # Do NOT write 'if not compression:',
                           # it could be zero.
    compression = options['compression']

  # Create cache directory if needed
  CD = checkdir(cachedir,verbose)

  # Handle the case cache('clear')
  if isinstance(my_F, str):
    if my_F.lower() == 'clear':
      clear_cache(CD,verbose=verbose)
      return

  # Handle the case cache(my_F, 'clear')
  if isinstance(args, str):
    if args.lower() == 'clear':
      clear_cache(CD,my_F,verbose=verbose)
      return

  # Force singleton arg into a tuple
  if not isinstance(args, tuple):
    args = tuple([args])
  
  # Check that kwargs is a dictionary
  if not isinstance(kwargs, dict):
    raise TypeError    
    
  # Hash arguments (and keyword args) to integer
  arghash = myhash((args, kwargs))
  
  # Get sizes and timestamps for files listed in dependencies.
  # Force singletons into a tuple.
  if dependencies and not isinstance(dependencies, (tuple, list)):
    dependencies = tuple([dependencies])
  deps = get_depstats(dependencies)

  # Extract function name from my_F object
  funcname = get_funcname(my_F)

  # Create cache filename
  FN = funcname+'_'+str(arghash)  
  #print()
  #print('FN', FN)
  #print('repr(arghash)', repr(arghash))
  #print('arghash', arghash)  
  #print()
  
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
            log.critical('MESSAGE (caching): File %s deleted' % fn)
        ##else:
        ##  log.critical('%s was not accessed' % fn)
    return None


  #-------------------------------------------------------------------        
  
  # Check if previous computation has been cached
  if evaluate is True:
    Retrieved = None  # Force evaluation of my_F regardless of caching status.
    reason = 5
  else:
    T, FN, Retrieved, reason, comptime, loadtime, compressed = \
        CacheLookup(CD, FN, my_F, 
                    args, kwargs, 
                    deps, 
                    verbose, 
                    compression,
                    dependencies)

  if not Retrieved:
    if test:  # Do not attempt to evaluate function
      T = None
    else:  # Evaluate function and save to cache
      if verbose is True:
        
        msg1(funcname, args, kwargs,reason)

      # Remove expired files automatically
      if options['expire']:
        DeleteOldFiles(CD,verbose)
        
      # Save args before function is evaluated in case
      # they are modified by function
      save_args_to_cache(CD,FN,args,kwargs,compression)

      # Execute and time function with supplied arguments
      t0 = time.time()

      T = my_F(*args, **kwargs) # Built-in 'apply' deprecated in Py3K    
      
      #comptime = round(time.time()-t0)
      comptime = time.time()-t0

      if verbose is True:
        msg2(funcname,args,kwargs,comptime,reason)

      # Save results and estimated loading time to cache
      loadtime = save_results_to_cache(T, CD, FN, my_F, deps, comptime, \
                                       funcname, dependencies, compression)
      if verbose is True:
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

# Has mostly been moved to proper unit test.
# What remains here includes example of the 
# cache statistics form.
def test(cachedir=None, verbose=False, compression=None):
  """Test the functionality of caching.

  USAGE:
    test(verbose)

  ARGUMENTS:
    verbose --     Flag whether caching will output its statistics (default=False)
    cachedir --    Directory for cache files (Default: options['cachedir'])
    compression -- Flag zlib compression (Default: options['compression'])
  """
   
  import string, time

  # Initialise
  #
  #import caching
  #reload(caching)

  if not cachedir:
    cachedir = options['cachedir']

  if verbose is None:  # Do NOT write 'if not verbose:', it could be zero.
    verbose = options['verbose']
  
  if compression is None:  # Do NOT write 'if not compression:',
                           # it could be zero.
    compression = options['compression']
  else:
    try:
      set_option('compression', compression)
    except:
      logtesterror('Set option failed')      

  try:
    import zlib
  except:
    log.critical()
    log.critical('*** Could not find zlib, default to no-compression      ***')
    log.critical('*** Installing zlib will improve performance of caching ***')
    log.critical()
    compression = 0        
    set_option('compression', compression)    
  
  log.critical('\nTesting caching module - please stand by\n')

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
    logtestOK('Wrote file %s' %DepFN)
  except:
    logtesterror('Could not open file %s for writing - check your environment' \
               % DepFN)

  # Check set_option (and switch stats off
  #    
  try:
    set_option('savestat',0)
    assert(options['savestat'] == 0)
    logtestOK('Set option')
  except:
    logtesterror('Set option failed')    
    
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
      T1 = cache(f,(a,b,c,N), {'x':x, 'y':y}, evaluate=1, \
                   verbose=verbose, compression=comp)
      if comp:                   
        logtestOK('Caching evaluation with compression')
      else:     
        logtestOK('Caching evaluation without compression')      
    except:
      if comp:
        logtesterror('Caching evaluation with compression failed - try caching.test(compression=0)')
      else:
        logtesterror('Caching evaluation failed - try caching.test(verbose=1)')

    # Retrieve
    #                           
    try:                         
      T2 = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
                   compression=comp) 

      if comp:                   
        logtestOK('Caching retrieval with compression')
      else:     
        logtestOK('Caching retrieval without compression')      
    except:
      if comp:
        logtesterror('Caching retrieval with compression failed - try caching.test(compression=0)')
      else:                                      
        logtesterror('Caching retrieval failed - try caching.test(verbose=1)')

    # Reference result
    #   
    T3 = f(a,b,c,N,x=x,y=y)  # Compute without caching
    
    if T1 == T2 and T2 == T3:
      if comp:
        logtestOK('Basic caching functionality (with compression)')
      else:
        logtestOK('Basic caching functionality (without compression)')
    else:
      logtesterror('Cached result does not match computed result')


  # Test return_filename
  #    
  try:
    FN = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
                 return_filename=1)    
    assert(FN[:2] == 'f[')
    logtestOK('Return of cache filename')
  except:
    logtesterror('Return of cache filename failed')

  # Test existence of cachefiles
  #  
  try:
    (datafile,compressed0) = myopen(CD+FN+'_'+file_types[0],"rb",compression)
    (argsfile,compressed1) = myopen(CD+FN+'_'+file_types[1],"rb",compression)
    (admfile,compressed2) =  myopen(CD+FN+'_'+file_types[2],"rb",compression)
    logtestOK('Presence of cache files')
    datafile.close()
    argsfile.close()
    admfile.close()
  except:
    logtesterror('Expected cache files did not exist') 
              
  # Test 'test' function when cache is present
  #      
  try:
    #T1 = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
    #                   evaluate=1)  
    T4 = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, test=1)
    assert(T1 == T4)

    logtestOK("Option 'test' when cache file present")
  except:
    logtesterror("Option 'test' when cache file present failed")      

  # Test that 'clear' works
  #
  #try:
  #  cache(f,'clear',verbose=verbose)
  #  logtestOK('Clearing of cache files')
  #except:
  #  logtesterror('Clear does not work')
  try:
    cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, clear=1)    
    logtestOK('Clearing of cache files')
  except:
    logtesterror('Clear does not work')  

  

  # Test 'test' function when cache is absent
  #      
  try:
    T4 = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, test=1)
    assert(T4 is None)
    logtestOK("Option 'test' when cache absent")
  except:
    logtesterror("Option 'test' when cache absent failed")      
          
  # Test dependencies
  #
  T1 = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
               dependencies=DepFN)  
  T2 = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
               dependencies=DepFN)                     
                       
  if T1 == T2:
    logtestOK('Basic dependencies functionality')
  else:
    logtesterror('Dependencies do not work')

  # Test basic wildcard dependency
  #
  T3 = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
               dependencies=DepFN_wildcard)                     
    
  if T1 == T3:
    logtestOK('Basic dependencies with wildcard functionality')
  else:
    logtesterror('Dependencies with wildcards do not work')


  # Test that changed timestamp in dependencies triggers recomputation
  
  # Modify dependency file
  Depfile = open(DepFN,'a')
  Depfile.write('You must cut down the mightiest tree in the forest with a Herring')
  Depfile.close()
  
  T3 = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
               dependencies=DepFN, test = 1)                     
  
  if T3 is None:
    logtestOK('Changed dependencies recognised')
  else:
    logtesterror('Changed dependencies not recognised')    
  
  # Test recomputation when dependencies have changed
  #
  T3 = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose, \
               dependencies=DepFN)                       
  if T1 == T3:
    logtestOK('Recomputed value with changed dependencies')
  else:
    logtesterror('Recomputed value with changed dependencies failed')

  # Performance test (with statistics)
  # Don't really rely on this as it will depend on specific computer. 
  #

  set_option('savestat',1)

  N = 20*N   #Should be large on fast computers...
  tt = time.time()
  T1 = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose)
  t1 = time.time() - tt
  
  tt = time.time()
  T2 = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=verbose)
  t2 = time.time() - tt
  
  if T1 == T2:
    if t1 > t2:
      logtestOK('Performance test: relative time saved = %s pct' \
              %str(round((t1-t2)*100/t1,2)))
  else:       
    logtesterror('Basic caching failed for new problem')
            
  # Test presence of statistics file
  #
  try: 
    DIRLIST = os.listdir(CD)
    SF = []
    for FN in DIRLIST:
      if string.find(FN,statsfile) >= 0:
        fid = open(CD+FN,'r')
        fid.close()
    logtestOK('Statistics files present') 
  except:
    logtestOK('Statistics files cannot be opened')          
      
  print_header_box('Show sample output of the caching function:')
  
  T2 = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=0)
  T2 = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=0)
  T2 = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=1)
  
  print_header_box('Show sample output of cachestat():')
  if unix:
    cachestat()    
  else:
    try:
      import time
      t = time.strptime('2030','%Y')
      cachestat()
    except:  
      log.critical('cachestat() does not work here, because it relies on '
                   'time.strptime() which is unavailable in Windows')
      
  logtestOK('Caching self test completed')    
      
            
  # Test setoption (not yet implemented)
  #

  
#==============================================================================
# Auxiliary functions
#==============================================================================

# Import pickler
# cPickle is used by functions mysave, myload, and compare
#
#import cPickle  # 10 to 100 times faster than pickle
#import pickle as pickler
import dill as pickler
#pickler = cPickle 

# Local immutable constants
#
comp_level = 1              # Compression level for zlib.
                            # comp_level = 1 works well.
textwidth1 = 16             # Text width of key fields in report forms.
#textwidth2 = 132            # Maximal width of textual representation of
textwidth2 = 300            # Maximal width of textual representation of
                            # arguments.
textwidth3 = 16             # Initial width of separation lines. Is modified.
textwidth4 = 50             # Text width in logtestOK()
statsfile  = '.cache_stat'  # Basefilename for cached statistics.
                            # It will reside in the chosen cache directory.

file_types = ['Result',     # File name extension for cached function results.
              'Args',       # File name extension for stored function args.
              'Admin']      # File name extension for administrative info.

Reason_msg = ['OK',         # Verbose reasons for recomputation
              'No cached result', 
              'Dependencies have changed', 
              'Arguments have changed',
              'Bytecode has changed',
              'Recomputation was requested by caller',
              'Cached file was unreadable']              
              
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

def CacheLookup(CD, FN, my_F, args, kwargs, deps, verbose, compression, 
                dependencies):
  r"""Determine whether cached result exists and return info.

  USAGE:
    (T, FN, Retrieved, reason, comptime, loadtime, compressed) = \  
    CacheLookup(CD, FN, my_F, args, kwargs, deps, verbose, compression, \
                dependencies)

  INPUT ARGUMENTS:
    CD --            Cache Directory
    FN --            Suggested cache file name
    my_F --          Callable object
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
                     3: Arguments have changed
                     4: Bytecode has changed
                     5: Recomputation was forced
                     6: Unreadable file
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

  import time, string

  # Assess whether cached result exists - compressed or not.
  #
  if verbose:
    log.critical('Caching: looking for cached files %s_{%s,%s,%s}.z'
                 % (CD+FN, file_types[0], file_types[1], file_types[2]))
  (datafile,compressed0) = myopen(CD+FN+'_'+file_types[0],"rb",compression)
  (argsfile,compressed1) = myopen(CD+FN+'_'+file_types[1],"rb",compression)
  (admfile,compressed2) =  myopen(CD+FN+'_'+file_types[2],"rb",compression)

  if verbose is True and deps is not None:
    log.critical('Caching: Dependencies are %s' % list(deps.keys()))

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
  R, reason = myload(argsfile, compressed)  # The original arguments
  argsfile.close()
  
  if reason > 0:
      # Recompute using same filename   
      return(None, FN, None, reason, None, None, None)
  else:   
      (argsref, kwargsref) = R

  R, reason = myload(admfile, compressed)
  admfile.close()  

  if reason > 0:
    return(None,FN,None,reason,None,None,None) # Recompute using same filename 

  
  depsref  = R[0]  # Dependency statistics
  comptime = R[1]  # The computation time
  coderef  = R[2]  # The byte code
  funcname = R[3]  # The function name

  # Check if dependencies have changed
  #
  if dependencies and not compare(depsref, deps):
    if verbose:
      log.critical('Dependencies %s have changed - recomputing' % dependencies)

    # Don't use cached file - recompute
    reason = 2
    return(None, FN, None, reason, None, None, None)

  # Get bytecode from my_F
  #
  bytecode = get_bytecode(my_F)

  # Check if arguments or bytecode have changed
  if compare(argsref, args) and compare(kwargsref, kwargs) and \
     (not options['bytecode'] or compare(bytecode, coderef)):

    # Arguments and dependencies match. Get cached results
    T, loadtime, compressed, reason = load_from_cache(CD, FN, compressed)
    if reason > 0:
        # Recompute using same FN     
        return(None, FN, None, reason, None, None, None)

    Retrieved = 1
    reason = 0

    if verbose:
      msg4(funcname,args,kwargs,deps,comptime,loadtime,CD,FN,compressed)

      if loadtime >= comptime:
        log.critical('Caching did not yield any gain.')
        log.critical('Consider executing function %s without caching.'
                     % funcname)
  else:

    # Non matching arguments or bytecodes signify a hash-collision.
    # This is resolved by recursive search of cache filenames
    # until either a matching or an unused filename is found.
    #
    (T, FN, Retrieved, reason, comptime, loadtime, compressed) = \
        CacheLookup(CD, FN+'x', my_F, args, kwargs, deps, 
                    verbose, compression, dependencies)

    # The real reason is that args or bytecodes have changed.
    # Not that the recursive seach has found an unused filename
    if not Retrieved:
      if not compare(bytecode, coderef):
        reason = 4 # Bytecode has changed
      else:   
        reason = 3 # Arguments have changed 
        
  # PADARN NOTE 17/12/12: Adding a special case to handle the existence of a 
  # FitInterpolate object. C Structures are serialised so they can be pickled.
  #---------------------------------------------------------------------------
  from anuga.fit_interpolate.general_fit_interpolate import FitInterpolate
  
  # Setup for quad_tree extension
  #from anuga.utilities import compile
  #if compile.can_use_C_extension('quad_tree_ext.c'):
  #import quad_tree_ext
  #else:
  #    msg = "C implementation of quad tree extension not avaliable"
  #    raise Exception(msg)

  # Setup for sparse_matrix extension
  #from anuga.utilities import compile
  #if compile.can_use_C_extension('sparse_matrix_ext.c'):

  #else:
  #    msg = "C implementation of sparse_matrix extension not avaliable"
  #    raise Exception(msg)

  import anuga.utilities.sparse_matrix_ext as sparse_matrix_ext
  import anuga.utilities.quad_tree_ext as quad_tree_ext
  from anuga.geometry.aabb import AABB

  if isinstance(T, FitInterpolate):
    if hasattr(T,"D"):
        T.D=sparse_matrix_ext.deserialise_dok(T.D)
    if hasattr(T,"AtA"):
        T.AtA=sparse_matrix_ext.deserialise_dok(T.AtA)
    if hasattr(T,"root"):
        T.build_quad_tree(verbose=verbose)
  #---------------------------------------------------------------------------

  return((T, FN, Retrieved, reason, comptime, loadtime, compressed))

# -----------------------------------------------------------------------------

def clear_cache(CD, my_F=None, verbose=None):
  """Clear cache for my_F.

  USAGE:
     clear(CD, my_F, verbose)

  ARGUMENTS:
     CD --       Caching directory (required)
     my_F --     Function object (default: None)
     verbose --  Flag verbose output (default: None)

  DESCRIPTION:

    If my_F is None, clear everything,
    otherwise clear only files pertaining to my_F.
  """

  import os, re
   
  if CD[-1] != os.sep:
    CD = CD+os.sep
  
  if verbose is None:
    verbose = options['verbose']

  # FIXME: Windows version needs to be tested

  if my_F:
    funcname = get_funcname(my_F)
    if verbose:
      log.critical('Clearing %s' % CD+funcname+'*')

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
        log.critical('Remove the following files:')
        for file_name in file_names:
            log.critical('     ' + file_name)

        A = input('Delete (Y/N)[N] ?')
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

  if verbose is None:
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
      log.critical('Deleting %d expired files:' % delfiles)
      os.system('ls -lur '+CD+'* | head -' + repr(delfiles))            # List them
    os.system('ls -ur '+CD+'* | head -' + repr(delfiles) + ' | xargs /bin/rm')
                                                                  # Delete them
    # FIXME: Replace this with os.listdir and os.remove

# -----------------------------------------------------------------------------

def save_args_to_cache(CD, FN, args, kwargs, compression):
  """Save arguments to cache

  USAGE:
    save_args_to_cache(CD,FN,args,kwargs,compression)
  """

  import time, os, sys

  (argsfile, compressed) = myopen(CD+FN+'_'+file_types[1], 'wb', compression)
  if argsfile is None:
    msg = 'ERROR (caching): Could not open argsfile for writing: %s' %FN
    raise IOError(msg)

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

def save_results_to_cache(T, CD, FN, my_F, deps, comptime, funcname,
                          dependencies, compression):
  """Save computed results T and admin info to cache

  USAGE:
    save_results_to_cache(T, CD, FN, my_F, deps, comptime, funcname,
                          dependencies, compression)
  """

  import time, os, sys
  verbose = False

  # PADARN NOTE 17/12/12: Adding a special case to handle the existence of a 
  # FitInterpolate object. C Structures are serialised so they can be pickled.
  #---------------------------------------------------------------------------
  from anuga.fit_interpolate.general_fit_interpolate import FitInterpolate
  
  import anuga.utilities.quad_tree_ext as quad_tree_ext
  import anuga.utilities.sparse_matrix_ext as sparse_matrix_ext
  from anuga.geometry.aabb import AABB

  if isinstance(T, FitInterpolate):
    if hasattr(T,"D"):
        T.D=sparse_matrix_ext.serialise_dok(T.D)
    if hasattr(T,"AtA"):
        T.AtA=sparse_matrix_ext.serialise_dok(T.AtA)
    if hasattr(T,"root"):
        T.root.root=None


  #---------------------------------------------------------------------------

  (datafile, compressed1) = myopen(CD+FN+'_'+file_types[0],'wb',compression)
  (admfile, compressed2) = myopen(CD+FN+'_'+file_types[2],'wb',compression)

  if not datafile:
    if verbose:
        log.critical('ERROR: Could not open %s' % datafile.name)
    raise IOError

  if not admfile:
    if verbose:
        log.critical('ERROR: Could not open %s' % admfile.name)
    raise IOError

  t0 = time.time()

  mysave(T,datafile,compression)  # Save data to cache
  datafile.close()
  #savetime = round(time.time()-t0,2)
  savetime = time.time()-t0  

  bytecode = get_bytecode(my_F)  # Get bytecode from function object
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

def load_from_cache(CD, FN, compression):
  """Load previously cached data from file FN

  USAGE:
    load_from_cache(CD,FN,compression)
  """

  import time

  (datafile, compressed) = myopen(CD+FN+'_'+file_types[0],"rb",compression)
  t0 = time.time()
  T, reason = myload(datafile,compressed)

  loadtime = time.time()-t0
  datafile.close() 

  return T, loadtime, compressed, reason

# -----------------------------------------------------------------------------

def myopen(FN, mode, compression=True):
  """Open file FN using given mode

  USAGE:
    myopen(FN, mode, compression=True)

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

  return(file, compressed)

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
        #raise Exception
        reason = 6  # Unreadable file
        return None, reason  
      
      
      del RsC  # Free up some space
      R = pickler.loads(Rs)
    else:
      try:
        R = pickler.load(file)
      #except EOFError, e:
      except:
        #Catch e.g., file with 0 length or corrupted
        reason = 6  # Unreadable file
        return None, reason
      
  except MemoryError:
    if options['verbose']:
      log.critical('ERROR: Out of memory while loading %s, aborting'
                   % file.name)

    # Raise the error again for now
    #
    raise MemoryError

  return R, reason

# -----------------------------------------------------------------------------

def mysave(T, file, compression):
  """Save data T to file

  USAGE:
    mysave(T, file, compression)

  """

  bin = options['bin']

  if compression:
    try:
      import zlib
    except:
      log.critical()
      log.critical('*** Could not find zlib ***')
      log.critical('*** Try to run caching with compression off ***')
      log.critical("*** caching.set_option('compression', 0) ***")
      raise Exception
      

    try:
      Ts  = pickler.dumps(T, bin)
    except MemoryError:
      msg = '****WARNING (caching.py): Could not pickle data for compression.'
      msg += ' Try using compression = False'
      raise MemoryError(msg)
    else:  
      # Compressed pickling      
      TsC = zlib.compress(Ts, comp_level)
      file.write(TsC)
  else:
      # Uncompressed pickling
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
      #  raise Exception(msg)
      

# -----------------------------------------------------------------------------

    
def myhash(T, ids=None):
  """Compute hashed integer from a range of inputs.
  If T is not hashable being e.g. a tuple T, myhash will recursively 
  hash the values individually

  USAGE:
    myhash(T)

  ARGUMENTS:
    T -- Anything
  """

  
  # Replacing Python2: if type(T) in [TupleType, ListType, DictType, InstanceType]:
  if isinstance(T, (tuple, list, dict)) or type(T) is type:
      # Keep track of unique id's to protect against infinite recursion
      if ids is None: ids = []

      # Check if T has already been encountered
      i = id(T) 

      if i in ids:
          return 0 # T has been hashed already      
      else:
          ids.append(i)
    
    
  # Start hashing

  # On some architectures None, False and True gets different hash values
  if T is None:
      return(-1)
  if T is False:
      return(0)
  if T is True:
      return(1)

  # Get hash values for hashable entries
  if isinstance(T, (tuple, list)):
      #print('LIST or TUPLE', T)
      hvals = ''
      for t in T:
          h = myhash(t, ids)
          hvals += str(h)
          
      val = hash(hvals)
  elif isinstance(T, dict):
      #print('DICT')
      
      I = list(T.items())    
      if system_tools.major_version == 2:
          # Make dictionary ordering unique
          I.sort()
      else:
          # As of Python 3.7 they now are ordered: https://mail.python.org/pipermail/python-dev/2017-December/151283.html

          pass
      val = myhash(I, ids)
  elif isinstance(T, num.ndarray):
      #print('NUM')
      T = num.array(T) # Ensure array is contiguous

      # Use mean value for efficiency
      val = hash(num.average(T.flat))
  elif callable(T):
      #print('CALLABLE')

      I = myhash(T.__dict__, ids)                
      val = myhash(I, ids)      
  elif type(T) is type: #isinstance(T, object):  # This is instead of the old InstanceType:
  #elif isinstance(T, object):  # This is instead of the old InstanceType:
      #print('OBJECT', T, dir(T), type(T)) 
      # Use the attribute values 
      val = myhash(T.__dict__, ids)
  else:
      # This must be a simple Python type that should hash without problems
      #print('ALL', T, type(T))
      val = hash(str(T))


  #print(ids, val)
  return(val)



def compare(A, B, ids=None):
    """Safe comparison of general objects

    USAGE:
      compare(A,B)

    DESCRIPTION:
      Return 1 if A and B they are identical, 0 otherwise
    """

    # Keep track of unique id's to protect against infinite recursion
    if ids is None: ids = {}

    # Check if T has already been encountered
    iA = id(A) 
    iB = id(B)     
    
    if (iA, iB) in ids:
        # A and B have been compared already
        return ids[(iA, iB)]
    else:
        ids[(iA, iB)] = True
    
    
    # Check if arguments are of same type
    if type(A) != type(B):
        return False
        
  
    # Compare recursively based on argument type
    if isinstance(A, (tuple, list)):
        N = len(A)
        if len(B) != N: 
            identical = False
        else:
            identical = True
            for i in range(N):
                if not compare(A[i], B[i], ids): 
                    identical = False; break
                    
    elif isinstance(A, dict):
        if len(A) != len(B):
            identical = False
        else:
            # Dictionaries are now ordered as of Python 3.7
            # Make dictionary ordering unique 
            #a = list(A.items()); a.sort()    
            #b = list(B.items()); b.sort()
            
            identical = compare(A, B, ids)
            
    elif isinstance(A, num.ndarray):
        # Use element by element comparison
        identical = num.all(A==B)

    #elif type(A) == types.InstanceType:
    elif type(A) is type:
        # Take care of special case where elements are instances            
        # Base comparison on attributes
        identical = compare(A.__dict__, 
                            B.__dict__, 
                            ids)
    else:       
        # Fall back to general code
        try:
            identical = (A == B)
        except:
            import pickle
            # Use pickle to compare data
            # The native pickler must be used
            # since the faster cPickle does not 
            # guarantee a unique translation
            # FIXME (Ole): Try to fall back on the dill pickler
            try:
                identical = (pickle.dumps(A,0) == pickle.dumps(B,0))
            except:
                identical = False

    # Record result of comparison and return            
    ids[(iA, iB)] = identical
    
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

def get_funcname(my_F):
  """Retrieve name of function object func (depending on its type)

  USAGE:
    get_funcname(my_F)
  """

  import string

  if type(my_F) == types.FunctionType:
    funcname = my_F.__name__
  elif type(my_F) == types.BuiltinFunctionType:
    funcname = my_F.__name__
  else:
    if system_tools.major_version == 3:
      tab = str.maketrans("<>'","   ")
      tmp = str.translate(repr(my_F), tab)
      tmp = str.split(tmp)
    elif system_tools.major_version == 2:
      tab = string.maketrans("<>'","   ")
      tmp = string.translate(repr(my_F), tab)
      tmp = string.split(tmp)
    else:
      raise Exception('Unsupported version: %' % system_tools.version)
      
    funcname = ' '.join(tmp)
    
    # Truncate memory address as in
    # class __main__.Dummy at 0x00A915D0
    index = funcname.find('at 0x')
    if index >= 0:
      funcname = funcname[:index+5] # Keep info that there is an address 

  funcname = nospace(funcname)
  return(funcname)

# -----------------------------------------------------------------------------

def get_bytecode(my_F):
    """ Get bytecode and associated values from function object.

    It is assumed that my_F is callable and there either 
    a function
    a class
    a method
    a callable object
    or a builtin function

    USAGE:
      get_bytecode(my_F)
    """

    if type(my_F) == types.FunctionType:
        # Function
        return get_func_code_details(my_F)
    elif type(my_F) == types.MethodType:
        # Method
        return get_func_code_details(my_F.__func__)
    elif type(my_F) in [types.BuiltinFunctionType, types.BuiltinMethodType]:      
        # Built-in functions are assumed not to change  
        return None, 0, 0, 0
    elif inspect.isclass(my_F):
        return get_func_code_details(my_F.__init__)
    elif hasattr(my_F, '__call__'):
        bytecode = get_func_code_details(my_F.__call__.__func__)
       
        # Add hash value of object to detect attribute changes
        return bytecode + (myhash(my_F),) 
    else:
        msg = 'Unknown function type: %s' % type(my_F)
        raise Exception(msg)          

    
  
def get_func_code_details(my_F):
  """Extract co_code, co_consts, co_argcount, func_defaults
  """
  
  bytecode = my_F.__code__.co_code
  consts = my_F.__code__.co_consts
  argcount = my_F.__code__.co_argcount    
  defaults = my_F.__defaults__       

  return bytecode, consts, argcount, defaults  

# -----------------------------------------------------------------------------

def get_depstats(dependencies):
  """ Build dictionary of dependency files and their size, mod. time and ctime.

  USAGE:
    get_depstats(dependencies):
  """

  d = {}
  if dependencies:

    #Expand any wildcards
    import glob
    expanded_dependencies = []
    for FN in dependencies:
      expanded_FN = glob.glob(FN)

      if expanded_FN == []:
        errmsg = 'ERROR (caching.py): Dependency '+FN+' does not exist.'
        raise Exception(errmsg)

      expanded_dependencies += expanded_FN

    
    for FN in expanded_dependencies:
      if not isinstance(FN, str):
        errmsg = 'ERROR (caching.py): Dependency must be a string.\n'
        errmsg += '                   Dependency given: %s' %FN
        raise Exception(errmsg)
      if not os.access(FN,os.F_OK):
        errmsg = 'ERROR (caching.py): Dependency '+FN+' does not exist.'
        raise Exception(errmsg)
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
    log.critical('Hack to get os.stat when files are too large')

    if unix:
      tmp = '/tmp/cach.tmp.'+repr(time.time())+repr(os.getpid())
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

  return(int(size),atime,mtime,ctime)

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

def checkdir(CD, verbose=None, warn=False):
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
      if verbose: log.critical('MESSAGE: Directory %s created.' % CD)
    except:
      if warn is True:
        log.critical('WARNING: Directory %s could not be created.' % CD)
      if unix:
        CD = '/tmp/'
      else:
        CD = 'C:'  
      if warn is True:
        log.critical('Using directory %s instead' % CD)

  return(CD)

checkdir(cachedir, warn=True)

#==============================================================================
# Statistics
#==============================================================================

def addstatsline(CD, funcname, FN, Retrieved, reason, comptime, loadtime,
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
    log.critical('Warning: Stat file could not be opened')

  try:
    if 'USER' in os.environ:
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
    log.critical('Warning: Writing of stat file failed')

# -----------------------------------------------------------------------------

# FIXME: should take cachedir as an optional arg
#
def __cachestat(sortidx=4, period=-1, showuser=None, cachedir=None):
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
  from string import split, rstrip, find
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
    log.critical('Reading file %s' % SD+FN)

    while True:
      A = input.readlines(blocksize)
      if len(A) == 0: break
      total_read = total_read + len(A)
      for record in A:
        record = tuple(split(rstrip(record),','))

        if len(record) == 9:
          timestamp = record[0]
        
          try:
            t = mktime(strptime(timestamp))
          except:
            total_discarded = total_discarded + 1         
            continue    
             
          if t > lastday:
            lastday = t
          if t < firstday:
            firstday = t

          user     = record[1]
          my_F     = record[2]

          # Strip hash-stamp off 
          #
          i = find(my_F,'[')
          my_F = my_F[:i]

          size        = float(record[3])

          # Compression kepword can be Boolean
          if record[4] in ['True', '1']:
            compression = 1
          elif record[4] in ['False', '0']:  
            compression = 0
          else:
            log.critical('Unknown value of compression %s' % str(record[4]))
            log.critical(str(record))
            total_discarded = total_discarded + 1            
            continue

          #compression = int(record[4]) # Can be Boolean
          hit         = int(record[5])
          reason      = int(record[6])   # Not used here    
          cputime     = float(record[7])
          loadtime    = float(record[8])

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
            UpdateDict(FuncDict,my_F,info)
          else:
            pass #Stats on recomputations and their reasons could go in here
              
        else:
          total_discarded = total_discarded + 1

    input.close()

  # Compute averages of all sums and write list
  #

  if total_read == 0:
    printline(Widths,'=')
    log.critical('CACHING STATISTICS: No valid records read')
    printline(Widths,'=')
    return

  log.critical()
  printline(Widths,'=')
  log.critical('CACHING STATISTICS: '+ctime(firstday)+' to '+ctime(lastday))
  printline(Widths,'=')
  log.critical('  Total number of valid records %d' % total_read)
  log.critical('  Total number of discarded records %d' % total_discarded)
  log.critical('  Total number of hits %d' % total_hits)
  log.critical()

  log.critical('  Fields %s are averaged over number of hits' % Fields[2:])
  log.critical('  Time is measured in seconds and size in bytes')
  log.critical('  Tables are sorted by %s' % Fields[1:][sortidx])

  if showuser:
    Dictionaries = [FuncDict, UserDict]
  else:
    Dictionaries = [FuncDict]

  i = 0
  for Dict in Dictionaries:
    for key in list(Dict.keys()):
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
    log.critical()
    printline(Widths,'-')
    n = 0
    for s in Fields:
      if s == Fields[0]:  # Left justify
        s = Dictnames[i] + ' ' + s; i=i+1
        #exec "print '%-" + str(Widths[n]) + "s'%s,"; n=n+1
        log.critical('%-*s' % (Widths[n], s))
        n += 1
      else:
        #exec "print '%" + str(Widths[n]) + "s'%s,"; n=n+1
        log.critical('%*s' % (Widths[n], s))
        n += 1
    log.critical()
    printline(Widths,'-')

    # Output Values
    #
    for key in keylist:
      rec = Dict[key]
      n = 0
      if len(key) > Widths[n]: key = key[:Widths[n]-3] + '...'
      #exec "print '%-" + str(Widths[n]) + Types[n]+"'%key,";n=n+1
      log.critical('%-*s' % (Widths[n], str(key)))
      n += 1
      for val in rec:
        #exec "print '%" + str(Widths[n]) + Types[n]+"'%val,"; n=n+1
        log.critical('%*s' % (Widths[n], str(key)))
        n += 1
      log.critical()
    log.critical()

#==============================================================================
# Auxiliary stats functions
#==============================================================================

def UpdateDict(Dict, key, info):
  """Update dictionary by adding new values to existing.

  USAGE:
    UpdateDict(Dict,key,info)
  """

  if key in Dict:
    dinfo = Dict[key]
    for n in range(len(dinfo)):
      dinfo[n] = info[n] + dinfo[n]
  else:
    dinfo = info[:]  # Make a copy of info list

  Dict[key] = dinfo
  return Dict

# -----------------------------------------------------------------------------

def SortDict(Dict, sortidx=0):
  """Sort dictionary

  USAGE:
    SortDict(Dict,sortidx):

  DESCRIPTION:
    Sort dictionary of lists according field number 'sortidx'
  """

  sortlist  = []
  keylist = list(Dict.keys())
  for key in keylist:
    rec = Dict[key]
    if not isinstance(rec, (list, tuple)):
      rec = [rec]

    if sortidx > len(rec)-1:
      msg = 'ERROR: Sorting index too large, sortidx = %s' % str(sortidx)
      raise IndexError(msg)

    val = rec[sortidx]
    sortlist.append(val)

  A = list(zip(sortlist,keylist))
  A.sort()
  keylist = [x[1] for x in A]  # keylist sorted by sortidx

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

  log.critical(s)

#==============================================================================
# Messages
#==============================================================================

def msg1(funcname, args, kwargs, reason):
  """Message 1

  USAGE:
    msg1(funcname, args, kwargs, reason):
  """

  import string

  print_header_box('Evaluating function %s' %funcname)
  
  msg7(args, kwargs)
  msg8(reason)  
  
  print_footer()
  
# -----------------------------------------------------------------------------

def msg2(funcname, args, kwargs, comptime, reason):
  """Message 2

  USAGE:
    msg2(funcname, args, kwargs, comptime, reason)
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

  log.critical(str.ljust('| CPU time:', textwidth1) +
               str(round(comptime,2)) + ' seconds')

# -----------------------------------------------------------------------------

def msg3(savetime, CD, FN, deps, compression):
  """Message 3

  USAGE:
    msg3(savetime, CD, FN, deps, compression)
  """

  import string
  log.critical(str.ljust('| Loading time:', textwidth1) + 
               str(round(savetime,2)) + ' seconds (estimated)')
  msg5(CD,FN,deps,compression)

# -----------------------------------------------------------------------------

def msg4(funcname, args, kwargs, deps, comptime, loadtime, CD, FN, compression):
  """Message 4

  USAGE:
    msg4(funcname, args, kwargs, deps, comptime, loadtime, CD, FN, compression)
  """

  import string

  print_header_box('Caching statistics (retrieving)')
  
  msg6(funcname,args,kwargs)
  log.critical(str.ljust('| CPU time:', textwidth1) +
               str(round(comptime,2)) + ' seconds')
  log.critical(str.ljust('| Loading time:', textwidth1) +
               str(round(loadtime,2)) + ' seconds')
  log.critical(str.ljust('| Time saved:', textwidth1) +
               str(round(comptime-loadtime,2)) + ' seconds')
  msg5(CD,FN,deps,compression)

# -----------------------------------------------------------------------------

def msg5(CD, FN, deps, compression):
  """Message 5

  USAGE:
    msg5(CD, FN, deps, compression)

  DESCRIPTION:
   Print dependency stats. Used by msg3 and msg4
  """

  import os, time, string

  log.critical('|')
  log.critical(str.ljust('| Caching dir: ', textwidth1) + CD)

  if compression:
    suffix = '.z'
    bytetext = 'bytes, compressed'
  else:
    suffix = ''
    bytetext = 'bytes'

  for file_type in file_types:
    file_name = FN + '_' + file_type + suffix
    stats = os.stat(CD+file_name)
    log.critical(str.ljust('| ' + file_type + ' file: ', textwidth1) +
                 file_name + '('+ str(stats[6]) + ' ' + bytetext + ')')

  log.critical('|')
  if len(deps) > 0:
    log.critical('| Dependencies:  ')
    dependencies  = list(deps.keys())
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
      d = str.ljust(dlist[n]+':', maxd+1)
      t = str.ljust(tlist[n], maxt)
      s = str.rjust(slist[n], maxs)

      log.critical('| %s %s %s bytes' % (d, t, s))
  else:
    log.critical('| No dependencies')
  print_footer()

# -----------------------------------------------------------------------------

def msg6(funcname, args, kwargs):
  """Message 6

  USAGE:
    msg6(funcname, args, kwargs)
  """

  import string
  log.critical(str.ljust('| Function:', textwidth1) + funcname)

  msg7(args, kwargs)
  
# -----------------------------------------------------------------------------    

def msg7(args, kwargs):
  """Message 7
  
  USAGE:
    msg7(args, kwargs):
  """
  
  import string
  
  args_present = 0  
  if args:
    if len(args) == 1:
      log.critical(str.ljust('| Argument:', textwidth1) +
                   mkargstr(args[0], textwidth2))
    else:
      log.critical(str.ljust('| Arguments:', textwidth1) +
                   mkargstr(args, textwidth2))
    args_present = 1
            
  if kwargs:
    if len(kwargs) == 1:
      log.critical(str.ljust('| Keyword Arg:', textwidth1) +
                   mkargstr(kwargs, textwidth2))
    else:
      log.critical(str.ljust('| Keyword Args:', textwidth1) +
                   mkargstr(kwargs, textwidth2))
    args_present = 1

  if not args_present:                
    log.critical('| No arguments')    # Default if no args or kwargs present

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
  
  log.critical(str.ljust('| Reason:', textwidth1) + R)
    
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

  log.critical(s + '| ' + line + CR + s)

  textwidth3 = N

# -----------------------------------------------------------------------------
    
def print_footer():
  """Print line same width as that of print_header_box.
  """
  
  N = textwidth3
  s = '+' + '-'*N + CR    
      
  log.critical(s)
      
# -----------------------------------------------------------------------------

def mkargstr(args, textwidth, argstr = '', level=0):
  """ Generate a string containing first textwidth characters of arguments.

  USAGE:
    mkargstr(args, textwidth, argstr = '', level=0)

  DESCRIPTION:
    Exactly the same as str(args) possibly followed by truncation,
    but faster if args is huge.
  """

  if level > 10:
      # Protect against circular structures
      return '...'
  
  WasTruncated = 0

  if not isinstance(args, (tuple, list, dict)):
    if isinstance(args, str):
      argstr = argstr + "'"+str(args)+"'"
    else:
      # Truncate large numeric arrays before using str()
      if isinstance(args, num.ndarray):
#        if len(args.flat) > textwidth:  
#        Changed by Duncan and Nick 21/2/07 .flat has problems with 
#        non-contigous arrays and ravel is equal to .flat except it 
#        can work with non-contiguous  arrays
        if len(num.ravel(args)) > textwidth:
          args = 'Array: ' + str(args.shape)

      argstr = argstr + str(args)
  else:
    if isinstance(args, dict):
      argstr = argstr + "{"
      for key in list(args.keys()):
        argstr = argstr + mkargstr(key, textwidth, level=level+1) + ": " + \
                 mkargstr(args[key], textwidth, level=level+1) + ", "
        if len(argstr) > textwidth:
          WasTruncated = 1
          break
      argstr = argstr[:-2]  # Strip off trailing comma      
      argstr = argstr + "}"

    else:
      if isinstance(args, tuple):
        lc = '('
        rc = ')'
      else:
        lc = '['
        rc = ']'
      argstr = argstr + lc
      for arg in args:
        argstr = argstr + mkargstr(arg, textwidth, level=level+1) + ', '
        if len(argstr) > textwidth:
          WasTruncated = 1
          break

      # Strip off trailing comma and space unless singleton tuple
      #
      if isinstance(args, tuple) and len(args) == 1:
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

def logtestOK(msg):
  """Print OK msg if test is OK.
  
  USAGE
    logtestOK(message)
  """

  import string
    
  log.critical(str.ljust(msg, textwidth4) + ' - OK' )
  
  #raise StandardError
  
# -----------------------------------------------------------------------------

def logtesterror(msg):
  """Print error if test fails.
  
  USAGE
    logtesterror(message)
  """
  
  log.critical('ERROR (caching.test): %s' % msg)
  log.critical('Please send this code example and output to ')
  log.critical('Ole.Nielsen@anu.edu.au')
  log.critical()
  log.critical()
  
  raise Exception

#-------------------------------------------------------------
if __name__ == "__main__":
  pass
