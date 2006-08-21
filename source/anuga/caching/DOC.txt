Documentation for module caching.py

Ole.Nielsen@anu.edu.au
-----------------------------------
	
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
    
    Sometimes compression==1 can cause a MemoryException. 
    In that case try the uncompressed version.

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
