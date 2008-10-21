
import unittest
from Numeric import arange, array

from anuga.caching import *
from anuga.caching.dummy_classes_for_testing import Dummy, Dummy_memorytest


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
  
def f_numeric(A, B):
  """Operation on Numeric arrays
  """
  
  return 3.1*A + B + 1
  
  
def f_object(A, B):
  """Operation of objecs of class Dummy
  """
  
  return A.value+B.value, A.another+B.another 
  
  
  
def clear_and_create_cache(Dummy, verbose=False):

  a = cache(Dummy, 'clear', verbose=verbose)
  a = cache(Dummy, args=(9,10),
            verbose=verbose)
      

def retrieve_cache(Dummy, verbose=False):
  if verbose: print 'Check that cache is there'
  
  X = cache(Dummy, args=(9,10), test=1,
            verbose=verbose)      
            
  msg = 'Cached value was not found'
  assert X is not None, msg

      

    

class Test_Caching(unittest.TestCase):
    def setUp(self):
        set_option('verbose', 0)  #Why

        pass

    def tearDown(self):
        pass

    def test_simple(self):
        """Check set_option (and switch stats off)
        """

        set_option('savestat', 0)
        assert options['savestat'] == 0
        set_option('verbose', 0)
        assert options['verbose'] == 0        
        

    def test_basic_caching(self):

        verbose=False
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

        comprange = 2

        for comp in range(comprange):
  
            # Evaluate and store
            #
            T1 = cache(f, (a,b,c,N), {'x':x, 'y':y}, evaluate=1, \
                       compression=comp, verbose=verbose)

            # Retrieve
            #                           
            T2 = cache(f, (a,b,c,N), {'x':x, 'y':y}, compression=comp) 

            # Reference result
            #   
            T3 = f(a,b,c,N,x=x,y=y)  # Compute without caching


            assert T1 == T2, 'Cached result does not match computed result'
            assert T2 == T3, 'Cached result does not match computed result'
            

    def test_caching_of_numeric_arrays(self):
        """test_caching_of_numeric_arrays
        
        Test that Numeric arrays can be recognised by caching even if their id's are different
        """
        
        verbose = False
        
        # Make some test input arguments
        A0 = arange(5)
        B0 = array([1.1, 2.2, 0.0, -5, -5])
        
        A1 = A0.copy()
        B1 = B0.copy()
        
        # Check that their ids are different
        assert id(A0) != id(A1)
        assert id(B0) != id(B1)        
        
        
        # Test caching
        comprange = 2
        for comp in range(comprange):
  
            # Evaluate and store
            T1 = cache(f_numeric, (A0, B0), evaluate=1,
                       compression=comp, verbose=verbose)

            # Retrieve
            T2 = cache(f_numeric, (A1, B1), 
                       compression=comp, test=1, verbose=verbose) 
                       
            # Check for presence of cached result 
            msg = 'Different array objects with same contents were not recognised'            
            assert T2 is not None, msg

            # Reference result
            T3 = f_numeric(A0, B0) # Compute without caching


            assert T1 == T2, 'Cached result does not match computed result'
            assert T2 == T3, 'Cached result does not match computed result'
            


    def test_caching_of_objects(self):
        """test_caching_of_objects
        
        Test that Objecs can be recognised as input variabelse 
        by caching even if their id's are different
        """
    

        verbose = False
        
        # Make some test input arguments
        A0 = Dummy(5, 7)
        B0 = Dummy(2.2, -5)
        
        A1 = A0.copy()
        B1 = B0.copy()
        
        # Check that their ids are different
        assert id(A0) != id(A1)
        assert id(B0) != id(B1)        
        
        
        # Test caching
        comprange = 2
        for comp in range(comprange):
  
            # Evaluate and store
            T1 = cache(f_object, (A0, B0), evaluate=1,
                       compression=comp, verbose=verbose)

            # Retrieve
            T2 = cache(f_object, (A1, B1), 
                       compression=comp, test=1, verbose=verbose) 
                       
            # Check for presence of cached result 
            msg = 'Different objects with same attributes were not recognised'
            assert T2 is not None, msg

            # Reference result
            T3 = f_object(A0, B0) # Compute without caching


            assert T1 == T2, 'Cached result does not match computed result'
            assert T2 == T3, 'Cached result does not match computed result'
            

                        
            
            
    def test_cachefiles(self):
        """Test existence of cachefiles
        """        
        N = 5000  #Make N fairly small here

        a = [1,2]
        b = ('Thou shalt count the number three',4)
        c = {'Five is right out': 6, (7,8): 9}
        x = 3
        y = 'holy hand granate'

        
        FN = cache(f,(a,b,c,N), {'x':x, 'y':y}, verbose=0, \
                  return_filename = 1)


        assert FN[:2] == 'f['

        CD = checkdir(cachedir)
        compression = 1

        (datafile,compressed0) = myopen(CD+FN+'_'+file_types[0],"rb",compression)
        (argsfile,compressed1) = myopen(CD+FN+'_'+file_types[1],"rb",compression)
        (admfile,compressed2) =  myopen(CD+FN+'_'+file_types[2],"rb",compression)

        datafile.close()
        argsfile.close()
        admfile.close()

    def test_test(self):        
        """Test 'test' function when cache is present
        """
        N = 5000  #Make N fairly small here

        a = [1,2]
        b = ('Thou shalt count the number three',4)
        c = {'Five is right out': 6, (7,8): 9}
        x = 3
        y = 'holy hand granate'
        

        T1 = cache(f,(a,b,c,N), {'x':x, 'y':y}, evaluate=1)
        
        T4 = cache(f,(a,b,c,N), {'x':x, 'y':y}, test=1)
        assert T1 == T4, "Option 'test' when cache file present failed"      


    def test_clear(self):        
        """Test that 'clear' works
        """

        N = 5000  #Make N fairly small here

        a = [1,2]
        b = ('Thou shalt count the number three',4)
        c = {'Five is right out': 6, (7,8): 9}
        x = 3
        y = 'holy hand granate'
        

        T1 = cache(f,(a,b,c,N), {'x':x, 'y':y}, evaluate = 1)
        
        
        cache(f, (a,b,c,N), {'x':x, 'y':y}, clear = 1)    

  
        # Test 'test' function when cache is absent
        
        
        T4 = cache(f, (a,b,c,N), {'x':x, 'y':y}, test=1)
        #print 'T4', T4
        assert T4 is None, "Option 'test' when cache absent failed"


    def test_dependencies(self):
        
        # Make a dependency file
        CD = checkdir(cachedir)        

        DepFN = CD + 'testfile.tmp'
        DepFN_wildcard = CD + 'test*.tmp'
        Depfile = open(DepFN,'w')
        Depfile.write('We are the knights who say NI!')
        Depfile.close()

        
        # Test 
        #

        N = 5000  #Make N fairly small here

        a = [1,2]
        b = ('Thou shalt count the number three',4)
        c = {'Five is right out': 6, (7,8): 9}
        x = 3
        y = 'holy hand granate'
        
        T1 = cache(f,(a,b,c,N), {'x':x, 'y':y}, dependencies=DepFN)  
        T2 = cache(f,(a,b,c,N), {'x':x, 'y':y}, dependencies=DepFN)                     
                       
        assert T1 == T2, 'Dependencies do not work'


        # Test basic wildcard dependency
        T3 = cache(f,(a,b,c,N), {'x':x, 'y':y}, dependencies=DepFN_wildcard)                     
    
        assert T1 == T3, 'Dependencies with wildcards do not work'


        # Test that changed timestamp in dependencies triggers recomputation
  
        # Modify dependency file
        Depfile = open(DepFN,'a')
        Depfile.write('You must cut down the mightiest tree in the forest with a Herring')
        Depfile.close()
  
        T3 = cache(f,(a,b,c,N), {'x':x, 'y':y}, dependencies=DepFN, test = 1)
        
        assert T3 is None, 'Changed dependencies not recognised'
  
        # Test recomputation when dependencies have changed
        #
        T3 = cache(f,(a,b,c,N), {'x':x, 'y':y}, dependencies=DepFN)                       
        assert T1 == T3, 'Recomputed value with changed dependencies failed'

    #def test_performance(self):        
    #    """Performance test (with statistics)
    #    Don't really rely on this as it will depend on specific computer. 
    #    """
    #
    #    import time
    #    set_option('savestat', 1)
    #
    #    N = 300000   #Should be large on fast computers...
    #    a = [1,2]
    #    b = ('Thou shalt count the number three',4)
    #    c = {'Five is right out': 6, (7,8): 9}
    #    x = 3
    #    y = 'holy hand granate'
    #     
    #    
    #    tt = time.time()
    #    T1 = cache(f,(a,b,c,N), {'x':x, 'y':y})
    #    t1 = time.time() - tt
    # 
    #    tt = time.time()
    #    T2 = cache(f,(a,b,c,N), {'x':x, 'y':y})
    #    t2 = time.time() - tt
     
    #    assert T1 == T2
    #     assert t1 > t2, 'WARNING: Performance a bit low - this could be specific to current platform. Try to increase N in this test'
    #    #test_OK('Performance test: relative time saved = %s pct' \
    #    #        %str(round((t1-t2)*100/t1,2)))


    def test_statsfile(self):                    
        """Test presence of statistics file
        """
        import os, string
        statsfile  = '.cache_stat'  # Basefilename for cached statistics.
        
        CD = checkdir(cachedir)                
        DIRLIST = os.listdir(CD)
        SF = []
        for FN in DIRLIST:
            if string.find(FN,statsfile) >= 0:
                try:
                    fid = open(CD+FN,'r')
                    fid.close()
                except:
                    raise 'Statistics files cannot be opened'          
  
          

    # def test_network_cachedir(self):

#         #set_option('cachedir', 'H:\\.python_cache\\')
#         set_option('cachedir', 'V:\\2\\cit\\.python_cache\\')
#         set_option('verbose', 1)

        
#         # Make some test input arguments
#         #
#         N = 5000  #Make N fairly small here

#         a = [1,2]
#         b = ('Thou shalt count the number three',4)
#         c = {'Five is right out': 6, (7,8): 9}
#         x = 3
#         y = 'holy hand granate'

#         # Test caching
#         #

#         comprange = 2

#         for comp in range(comprange):
  
#             # Evaluate and store
#             #
#             T1 = cache(f, (a,b,c,N), {'x':x, 'y':y}, evaluate=1, \
#                        compression=comp)

#             # Retrieve
#             #                           
#             T2 = cache(f, (a,b,c,N), {'x':x, 'y':y}, compression=comp) 

#             # Reference result
#             #   
#             T3 = f(a,b,c,N,x=x,y=y)  # Compute without caching


#             assert T1 == T2, 'Cached result does not match computed result'
#             assert T2 == T3, 'Cached result does not match computed result'
              


    def Will_fail_test_objects(self):
      """
      This test shows how instances can't be effectively cached.
      myhash uses hash which uses id which uses the memory address. 
      """
      verbose = True
      #verbose = False

      for i in range(2):
        if verbose: print "clear cache"
        a = cache(Dummy, 'clear')
        
        if verbose: print "cache for first time"
        a = cache(Dummy, args=(9,10), verbose=verbose)
        hash_value = myhash(a)
        
        #print "hash_value",hash_value 
        if verbose: print "cache for second time"
        a = cache(Dummy, args=(9,10), verbose=verbose)
        
        #print "myhash(a)",myhash(a) 
        assert hash_value == myhash(a)


    # This test works in the caching dir and in anuga_core, but not in the
    # anuga_core/source/anuga dir
    # This has to do with pickle (see e.g. http://telin.ugent.be/~slippens/drupal/pickleproblem)
    # The error message is 
    # PicklingError: Can't pickle test_caching.Dummy: it's not the same object as test_caching.Dummy
    #
    # This problem was fixed by moving the class into the separate module
    
    def test_objects_are_created(self):
      """
      This test shows how instances can be created from cache
      as long as input arguments are unchanged.

      Such instances will have different id's and cannot be currently be used as input
      arguments in subsequent caches. However, this is still useful.

      Do it for all combinations of compression

      """

      verbose = False

      for compression_store in [False, True]:
        for compression_retrieve in [False, True]:        
        
          if verbose: print 'clear cache'
          a = cache(Dummy, 'clear')
        
          if verbose: print 'cache for first time'
          a = cache(Dummy, args=(9,10),
                    compression=compression_store,
                    verbose=verbose)
          
          if verbose: print 'Check that cache is there'
          assert cache(Dummy, args=(9,10), test=1,
                       compression=compression_retrieve,
                       verbose=verbose)



    # NOTE (Ole): This test has been commented out because, although the test will pass
    #             inside the caching dir and also at the anuga_core level,
    #             it won't pass at the anuga_core/source/anuga level.
    # It may have to do with the comments above.
    #
    # But this is probably not so important, really
    def Xtest_objects_are_created_memory(self):
      """
      
      This test shows how instances can be created from cache
      as long as input arguments are unchanged - even if the class
      lives in different memory locations.

      This is using cache created in the main program below
      """

      verbose = False

      # Redefine class Dummy_memorytest
      class Dummy_memorytest:
        def __init__(self, value, another):
          self.value = value      

      # Make sure that class has been redefined to another address
      print
      print 'Initial_addr  ', initial_addr
      print 'Redefined addr', `Dummy_memorytest`
      msg = 'Redefined class ended up at same memory location as '
      msg += 'original class making this test irrelevant. Try to run '
      msg += 'it again and see if this error goes away.'
      msg += 'If it persists contact Ole.Nielsen@ga.gov.au'
      assert initial_addr != `Dummy_memorytest`, msg   

      
      retrieve_cache(Dummy_memorytest, verbose=verbose)      
          
# Cache created for use with 'test_objects_are_created_memory'
initial_addr = `Dummy_memorytest`
clear_and_create_cache(Dummy_memorytest, verbose=False)
  
      



     
        


#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Caching,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
