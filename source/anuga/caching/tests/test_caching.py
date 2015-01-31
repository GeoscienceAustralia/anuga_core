
import unittest

from copy import deepcopy
        
from anuga.caching import *
from anuga.caching.dummy_classes_for_testing import Dummy, Dummy_memorytest

import numpy as num


# Define a test function to be cached
#
def f(a, b, c, N, x=0, y='abcdefg'):
    """f(a,b,c,N)
    Do something time consuming and produce a complex result.
    """

    import string
    
    B = []
    for n in range(N):
        s = str(n + 2.0 / (n + 4.0)) + '.a' * 10
    B.append((a, b, c, s, n, x, y))
    return(B)
  
def f_numeric(A, B):
    """Operation on numeric arrays
    """
  
    return 3.1 * A + B + 1
  
  
def f_object(A, B):
    """Operation of objects of class Dummy
    """
  
    return A.value + B.value, A.another + B.another 
  

def f_generic(A):
    return A 
  
def clear_and_create_cache(Dummy, verbose=False):
    
    a = cache(Dummy, 'clear', verbose=verbose)
    a = cache(Dummy, args=(9, 10),
            verbose=verbose)
      

def retrieve_cache(Dummy, verbose=False):
    if verbose: print('Check that cache is there')
  
    X = cache(Dummy, args=(9, 10), test=1,
            verbose=verbose)
    
    msg = 'Cached value was not found'
    assert X is not None, msg
    
class Test_Caching(unittest.TestCase):
    def setUp(self):
        set_option('verbose', 0)  # Why

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

        verbose = False
        # Make some test input arguments
        #
        N = 5000  # Make N fairly small here

        a = [1, 2]
        b = ('Thou shalt count the number three', 4)
        c = {'Five is right out': 6, (7, 8): 9}
        x = 3
        y = 'holy hand granate'

        # Test caching
        #

        comprange = 2

        for comp in range(comprange):
  
            # Evaluate and store
            #
            T1 = cache(f, (a, b, c, N), {'x':x, 'y':y}, evaluate=1, \
                       compression=comp, verbose=verbose)

            # Retrieve
            #                           
            T2 = cache(f, (a, b, c, N), {'x':x, 'y':y}, compression=comp) 

            # Reference result
            #   
            T3 = f(a, b, c, N, x=x, y=y)  # Compute without caching


            assert T1 == T2, 'Cached result does not match computed result'
            assert T2 == T3, 'Cached result does not match computed result'
            

            
    def test_caching_of_numeric_arrays(self):
        """test_caching_of_numeric_arrays
        
        Test that numeric arrays can be recognised by caching even if their id's are different
        """
        
        verbose = False
        
        # Make some test input arguments
        A0 = num.arange(5)
        B0 = num.array([1.1, 2.2, 0.0, -5, -5])
        
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
            T3 = f_numeric(A0, B0)  # Compute without caching


            assert num.alltrue(T1 == T2), 'Cached result does not match computed result'
            assert num.alltrue(T2 == T3), 'Cached result does not match computed result'
            

    def test_hash_collision(self):
        """Test that hash collisons are dealt with correctly"""
        
        verbose = False
        
        # Make test input arguments
        A0 = num.arange(5) * 1.0
        B = ('x', 15)
        
        # Create different A that hashes to the same address (having the same average)
        A1 = num.array([2.0, 2.0, 2.0, 2.0, 2.0])        
        
        assert myhash(A0) == myhash(A1)
            
            
        # Test caching
        comprange = 2
        for comp in range(comprange):
        
            # Clear
            cache(f_numeric, (A0, A0), clear=1,
                  compression=comp, verbose=verbose)        
            cache(f_numeric, (A1, A1), clear=1,
                  compression=comp, verbose=verbose)                          
  
  
            # Evaluate and store
            T1 = cache(f_numeric, (A0, A0), evaluate=1,
                       compression=comp, verbose=verbose)

            
            # Check that A1 doesn't trigger retrieval of the previous result 
            # even though it hashes to the same address
            T2 = cache(f_numeric, (A1, A1),
                       compression=comp, verbose=verbose) 
           
            T1_ref = f_numeric(A0, A0)
            T2_ref = f_numeric(A1, A1)
            T1_ref = f_numeric(A0, A0)
            T2_ref = f_numeric(A1, A1)

            assert num.alltrue(T1 == T1_ref)
            assert num.alltrue(T2 == T2_ref)


    def test_caching_of_dictionaries(self):
        """test_caching_of_dictionaries
        
        Real example from ANUGA that caused some
        hashing problems
        """
    

        verbose = False  # True
        
        D = {'point_attributes': None,
             'use_cache': True,
             'vertex_coordinates': None,
             'verbose': False,
             'max_read_lines': 500,
             'acceptable_overshoot': 1.01,
             'mesh': None,
             'data_origin': None,
             'alpha': 0.02,
             'mesh_origin': None,
             'attribute_name': None,
             'triangles': None}         
        
        DD = deepcopy(D)  # Mangles the dictionary ordering 
        
        assert myhash(DD) == myhash(D)

        # Also test caching now that we are at it
        comprange = 2
        for comp in range(comprange):
  
            # Evaluate and store using D
            T1 = cache(f_generic, D, evaluate=1,
                       compression=comp, verbose=verbose)

            # Retrieve using copy (DD)
            T2 = cache(f_generic, DD,
                       compression=comp, test=1, verbose=verbose) 
                       
            # Check for presence of cached result 
            msg = 'Cached object was not found'            
            assert T2 is not None, msg

            # Reference result
            T3 = f_generic(D)  # Compute without caching

            
            msg = 'Cached result does not match computed result' 
            
            # Compare dictionaries
            for key in T1:
                assert T1[key] == T2[key]
                assert T2[key] == T3[key]                
                
            
           
            

    def test_caching_of_objects(self):
        """test_caching_of_objects
        
        Test that Objecs can be recognised as input variabelse 
        by caching even if their id's are different
        """
    

        verbose = False
        
        # Make some test input arguments
        A0 = Dummy(5, 7)
        B0 = Dummy(2.2, -5)
        
        A0.new_attribute = 'x'  
        B0.new_attribute = 'x'        
        
        A1 = deepcopy(A0)
        B1 = deepcopy(B0)        
        
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
                       compression=comp,
                       test=1,
                       verbose=verbose) 
                       
            # Check for presence of cached result 
            msg = 'Different objects with same attributes were not recognised'
            assert T2 is not None, msg

            # Reference result
            T3 = f_object(A0, B0)  # Compute without caching


            assert T1 == T2, 'Cached result does not match computed result'
            assert T2 == T3, 'Cached result does not match computed result'
            

    def test_caching_of_circular_structures(self):
        """test_caching_of_circular_structures
        
        Test that Caching doesn't recurse infinitely in case of
        circular or self-referencing structures
        """
        
        verbose = False
        
        # Create input argument
        A = Dummy(5, 7)
        B = {'x': 10, 'A': A}
        C = [B, 15]
        A.value = C  # Make it circular
        A.x = [1, 2, C, 5, A]  # More circular and self referential
        
        AA = deepcopy(A)

        # Test caching
        comprange = 2
        for comp in range(comprange):
  
            # Evaluate and store
            T1 = cache(f_generic, A,
                       evaluate=1,
                       compression=comp, verbose=verbose)

            # Retrieve
            T2 = cache(f_generic, AA,
                       compression=comp,
                       test=1, verbose=verbose) 
                       
            # Check for presence of cached result 
            msg = 'Cached object was not found'            
            assert T2 is not None, msg

            # Reference result
            T3 = f_generic(A)  # Compute without caching

            
            msg = 'Cached result does not match computed result' 
            assert str(T1) == str(T2), msg
            assert str(T2) == str(T3), msg
                                    
            
    def test_caching_of_callable_objects(self):
        """test_caching_of_callable_objects(self)
        
        Test that caching will discern between calls of
        two different instances of a callable object with the same input.
        
        This cause problems with add_quantity in ANUGA in changeset:6225
        where a previous instance of Polygon_function was picked up.
        """

        class call:
        
          def __init__(self, a, b):
            self.a = a
            self.b = b
            
          def __call__(self, x):
            return self.a * x + self.b

            
        f1 = call(2, 3)
        f2 = call(5, 7)

        # Check that hash value of callable objects don't change
        # FIXME (Ole): The hash values do appear to change when OS
        # and/or dependencies are upgraded
        if os.name == 'posix' and os.uname()[4] in ['x86_64', 'ia64']:
          # 64 bit hash values
          f1hash = 7079146893884768701
          f2hash = -6995306676314913340

          # Prior to cluster upgrades Feb 2009
          # f1hash = 1914027059797211698 
          # f2hash = 1914027059807087171
        else:
          # 32 bit hash values
          f1hash = -758136387
          f2hash = -11221564     
          
        assert myhash(f1) == f1hash
        assert myhash(f2) == f2hash
        
        bc1 = get_bytecode(f1)
        bc2 = get_bytecode(f2)
        
        msg = 'Byte code should be different'
        assert bc1 != bc2, msg

        
        x = num.arange(10).astype(num.float)
        
        ref1 = f1(x)
        ref2 = f2(x)

        # Clear cache for f1(x) and f2(x) and verify that all is clear
        cache(f1, x, clear=True, verbose=False)
        flag = cache(f1, x, test=True, verbose=False)        
        assert flag is None
        
        cache(f2, x, clear=True, verbose=False)
        flag = cache(f2, x, test=True, verbose=False)        
        assert flag is None        
        
        # Run f1(x) and cache
        res1 = cache(f1, x, verbose=False)
        assert num.allclose(res1, ref1)        
        
        # Test that f1(x) has been cached correctly
        res1 = cache(f1, x, test=True, verbose=False)                
        assert num.allclose(res1, ref1)                
        
        # Test that f2(x) is still clear
        cache(f2, x, clear=True, verbose=False)
        flag = cache(f2, x, test=True, verbose=False)        
        assert flag is None                
        
        # Run f2(x) and test result
        res2 = cache(f2, x, verbose=False)
        msg = 'Wrong result for f2(x)'
        assert num.allclose(res2, ref2), msg
        

        
        
    def test_uniqueness_of_hash_values(self):
        """test_uniqueness_of_hash_values(self):
        
        Test that Caching can handle a realistic 
        complex structure by hashing it consistently and
        uniquely.
        """
        
        verbose = False
        
        # Create input argument
        A = Dummy(5, 7)
        B = {'x': 10, 'A': A}
        C = [B, num.array([1.2, 3, 5, 0.1])]
        A.value = C  # Make it circular

        # Create identical but separate object    
        AA = Dummy(None, None)
        BB = {'A': AA, 'x': 10}
        CC = [BB, num.array([1.200, 3.000, 5.00, 1.0 / 10])]
        AA.value = CC  # Make it circular
        AA.another = 3 + 4        
        
        
        assert myhash(A) == myhash(AA)     
           
           
        
        # Also test caching now that we are at it
        comprange = 2
        for comp in range(comprange):
  
            # Evaluate and store using A
            T1 = cache(f_generic, A, evaluate=1,
                       compression=comp, verbose=verbose)

            # Retrieve using copy (AA)
            T2 = cache(f_generic, AA,
                       compression=comp, test=1, verbose=verbose) 
                       
            # Check for presence of cached result 
            msg = 'Cached object was not found'            
            assert T2 is not None, msg

            # Reference result
            T3 = f_generic(A)  # Compute without caching

            
            msg = 'Cached result does not match computed result' 
            assert str(T1) == str(T2), msg
            assert str(T2) == str(T3), msg
            
           

    def test_caching_of_simple_circular_dictionaries(self):
        """test_caching_of_circular_structures
        
        Test that Caching doesn't recurse infinitely in case of
        circular or self-referencing structures
        """
        
        verbose = False  # True
        
        # Create input argument
        A = {'x': 10, 'B': None}
        B = [A, 15]
        A['B'] = B  # Make it circular
        
        # Test caching
        comprange = 2
        for comp in range(comprange):
  
            # Evaluate and store
            T1 = cache(f_generic, A, evaluate=1,
                       compression=comp, verbose=verbose)
                       

            # Retrieve
            T2 = cache(f_generic, A,
                       compression=comp, test=1, verbose=verbose) 
                       
            # Check for presence of cached result 
            msg = 'Cached object was not found'            
            assert T2 is not None, msg

            # Reference result
            T3 = f_generic(A)  # Compute without caching


            msg = 'Cached result does not match computed result'
            assert str(T1) == str(T2), msg
            assert str(T2) == str(T3), msg

            
            
    def test_cachefiles(self):
        """Test existence of cachefiles
        """        
        N = 5000  # Make N fairly small here

        a = [1, 2]
        b = ('Thou shalt count the number three', 4)
        c = {'Five is right out': 6, (7, 8): 9}
        x = 3
        y = 'holy hand granate'

        
        FN = cache(f, (a, b, c, N), {'x':x, 'y':y}, verbose=0, \
                  return_filename=1)


        assert FN[:2] == 'f['

        CD = checkdir(cachedir)
        compression = 1

        (datafile, compressed0) = myopen(CD + FN + '_' + file_types[0], "rb", compression)
        (argsfile, compressed1) = myopen(CD + FN + '_' + file_types[1], "rb", compression)
        (admfile, compressed2) = myopen(CD + FN + '_' + file_types[2], "rb", compression)

        datafile.close()
        argsfile.close()
        admfile.close()

    def test_test(self):        
        """Test 'test' function when cache is present
        """
        
        verbose = False
        
        N = 5  

        a = [1, 2]
        b = ('Thou shalt count the number three', 4)
        c = {'Five is right out': 6, (7, 8): 9}
        x = 3
        y = 'holy hand granate'
        

        T1 = cache(f, (a, b, c, N), {'x':x, 'y':y},
                   evaluate=1,
                   verbose=verbose)
        
        T2 = cache(f, (a, b, c, N), {'x':x, 'y':y},
                   test=1,
                   verbose=verbose)
                   
                   
        # Check for presence of cached result 
        msg = 'Different objects with same attributes were not recognised'
        assert T2 is not None, msg                   
                   
        assert T1 == T2, "Option 'test' when cache file present failed"      


    def test_clear(self):        
        """Test that 'clear' works
        """

        N = 5000  # Make N fairly small here

        a = [1, 2]
        b = ('Thou shalt count the number three', 4)
        c = {'Five is right out': 6, (7, 8): 9}
        x = 3
        y = 'holy hand granate'
        

        T1 = cache(f, (a, b, c, N), {'x':x, 'y':y}, evaluate=1)
        
        
        cache(f, (a, b, c, N), {'x':x, 'y':y}, clear=1)    

  
        # Test 'test' function when cache is absent
        
        
        T4 = cache(f, (a, b, c, N), {'x':x, 'y':y}, test=1)
        assert T4 is None, "Option 'test' when cache absent failed"


    def test_dependencies(self):
        
        # Make a dependency file
        CD = checkdir(cachedir)        

        DepFN = CD + 'testfile.tmp'
        DepFN_wildcard = CD + 'test*.tmp'
        Depfile = open(DepFN, 'w')
        Depfile.write('We are the knights who say NI!')
        Depfile.close()

        
        # Test 
        #

        N = 5000  # Make N fairly small here

        a = [1, 2]
        b = ('Thou shalt count the number three', 4)
        c = {'Five is right out': 6, (7, 8): 9}
        x = 3
        y = 'holy hand granate'
        
        T1 = cache(f, (a, b, c, N), {'x':x, 'y':y}, dependencies=DepFN)  
        T2 = cache(f, (a, b, c, N), {'x':x, 'y':y}, dependencies=DepFN)                     
                       
        assert T1 == T2, 'Dependencies do not work'


        # Test basic wildcard dependency
        T3 = cache(f, (a, b, c, N), {'x':x, 'y':y}, dependencies=DepFN_wildcard)                     
    
        assert T1 == T3, 'Dependencies with wildcards do not work'


        # Test that changed timestamp in dependencies triggers recomputation
  
        # Modify dependency file
        Depfile = open(DepFN, 'a')
        Depfile.write('You must cut down the mightiest tree in the forest with a Herring')
        Depfile.close()
  
        T3 = cache(f, (a, b, c, N), {'x':x, 'y':y}, dependencies=DepFN, test=1)
        
        assert T3 is None, 'Changed dependencies not recognised'
  
        # Test recomputation when dependencies have changed
        #
        T3 = cache(f, (a, b, c, N), {'x':x, 'y':y}, dependencies=DepFN)                       
        assert T1 == T3, 'Recomputed value with changed dependencies failed'

    # def test_performance(self):        
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
        statsfile = '.cache_stat'  # Basefilename for cached statistics.
        
        CD = checkdir(cachedir)                
        DIRLIST = os.listdir(CD)
        SF = []
        for FN in DIRLIST:
            if string.find(FN, statsfile) >= 0:
                try:
                    fid = open(CD + FN, 'r')
                    fid.close()
                except:
                    raise Exception('Statistics files cannot be opened')
  
          

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
              


    def Will_fail_objects(self):
      """
      This test shows how instances can't be effectively cached.
      myhash uses hash which uses id which uses the memory address. 
      
      This will be a NIL problem once caching can handle instances with different id's and 
      identical content.
      
      The test is disabled.
      """
      

      
      verbose = True
      # verbose = False

      for i in range(2):
        if verbose: print "clear cache"
        a = cache(Dummy, 'clear')
        
        if verbose: print "cache for first time"
        a = cache(Dummy, args=(9, 10), verbose=verbose)
        hash_value = myhash(a)
        
        # print "hash_value",hash_value 
        if verbose: print "cache for second time"
        a = cache(Dummy, args=(9, 10), verbose=verbose)
        
        # print "myhash(a)",myhash(a) 
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
          a_ref = cache(Dummy, args=(9, 10),
                        compression=compression_store,
                        verbose=verbose)
          
          if verbose: print 'Check that cache is there'
          assert cache(Dummy, args=(9, 10), test=1,
                       compression=compression_retrieve,
                       verbose=verbose)
                       
          if verbose: print 'Check cached result'
          a = cache(Dummy, args=(9, 10),
                    compression=compression_store,
                    verbose=verbose)                       
          assert a.__dict__ == a_ref.__dict__



#     # NOTE (Ole): This test has been commented out because, 
#     #             although the test will pass (not anymore!)
#     #             inside the caching dir and also at the anuga_core level,
#     #             it won't pass at the anuga_core/source/anuga level.
#     # It may have to do with the comments above.
#     #
#     # But this is a very nice test to run occasionally within the caching
#     # area
#     def Xtest_objects_are_created_memory(self):
#       """
#       
#       This test shows how instances can be created from cache
#       as long as input arguments are unchanged - even if the class
#       lives in different memory locations.
# 
#       This is using cache created in the main program below
#       """
# 
#       verbose = True  # False
# 
#       # Redefine class Dummy_memorytest
#       class Dummy_memorytest:
#         def __init__(self, value, another):
#           self.value = value      
# 
#       # Make sure that class has been redefined to another address
#       print
#       print 'Initial_addr  ', initial_addr
#       print 'Redefined addr', repr(Dummy_memorytest)
#       msg = 'Redefined class ended up at same memory location as '
#       msg += 'original class making this test irrelevant. Try to run '
#       msg += 'it again and see if this error goes away.'
#       msg += 'If it persists contact Ole.Nielsen@ga.gov.au'
#       assert initial_addr != repr(Dummy_memorytest), msg   
# 
#       
#       retrieve_cache(Dummy_memorytest, verbose=verbose)      
#           
# # Cache created for use with 'test_objects_are_created_memory'
# # initial_addr = `Dummy_memorytest`
# # clear_and_create_cache(Dummy_memorytest, verbose=False)
  


#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Caching, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
