#!/usr/bin/env python
# Test of MPI module 'pypar' for Python
# 
# Run as 
#   python pypartest.py
# or 
#   mpirun -np 2 pypartest.py
# (perhaps try number of processors more than 2)
#
# To verify bandwidth of your architecture please run pytiming (and ctiming) 
#
# OMN, FEB 2002

try:
  import Numeric
except:
  raise 'Module Numeric must be present to run pypar'


print "Importing pypar"
import pypar
methods = dir(pypar)
assert 'Abort' in methods
assert 'Finalize' in methods
assert 'Get_processor_name' in methods
assert 'Wtime' in methods
assert 'rank' in methods
assert 'raw_receive' in methods
assert 'raw_send' in methods
assert 'receive' in methods
assert 'send' in methods
assert 'bcast' in methods
assert 'size' in methods

print "Module pypar imported OK"
#pypar.Barrier()


myid =    pypar.rank()
numproc = pypar.size()
node =    pypar.Get_processor_name()

print "I am processor %d of %d on node %s" %(myid, numproc, node)
pypar.Barrier()


if numproc > 1:
  # Test simple raw communication (arrays, strings and general)
  #
  N = 17 #Number of elements
  
  if myid == 0:
    # Integer arrays
    #
    A = Numeric.array(range(N))
    B = Numeric.zeros(N)    
    pypar.raw_send(A,1)
    pypar.raw_receive(B,numproc-1)
    
    assert Numeric.allclose(A, B)
    print "Raw communication of numeric integer arrays OK"

    # Real arrays
    #
    A = Numeric.array(range(N)).astype('f')
    B = Numeric.zeros(N).astype('f')    
    pypar.raw_send(A,1)
    pypar.raw_receive(B,numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Raw communication of numeric real arrays OK"

    # Strings (< 256 characters)
    #
    A = "and now to something completely different !"
    B = " "*len(A)
    pypar.raw_send(A,1)
    pypar.raw_receive(B,numproc-1)
    
    assert A == B
    print "Raw communication of strings OK"
    
    # A more general structure
    #
    A = ['ABC', (1,2,3.14), {8: 'Monty'}, Numeric.array([13.45, 1.2])]
    B = ['   ', (0,0,0.0), {0: '     '}, Numeric.zeros(2).astype('f')]    
    pypar.raw_send(A,1)
    B = pypar.raw_receive(B,numproc-1)
    
    assert A == B
    print "Raw communication of general structures OK"
    
  else:  
    # Integers
    #
    X = Numeric.zeros(N)
    pypar.raw_receive(X, myid-1)  
    pypar.raw_send(X, (myid+1)%numproc)
  
    # Floats
    #
    X = Numeric.zeros(N).astype('f')
    pypar.raw_receive(X, myid-1)  
    pypar.raw_send(X, (myid+1)%numproc)    

    # Strings
    #
    X = " "*256
    pypar.raw_receive(X, myid-1)  
    pypar.raw_send(X.strip(), (myid+1)%numproc)    

    # General
    #
    X = ['   ', (0,0,0.0), {0: '     '}, Numeric.zeros(2).astype('f')]
    X = pypar.raw_receive(X, myid-1)  
    pypar.raw_send(X, (myid+1)%numproc)    
    

  # Test easy communication  - without buffers (arrays, strings and general)
  #
  N = 17 #Number of elements
  
  if myid == 0:
    # Integer arrays
    #
    A = Numeric.array(range(N))

    pypar.send(A,1)
    B = pypar.receive(numproc-1)
    

    assert Numeric.allclose(A, B)
    print "Simplified communication of numeric integer arrays OK"

    # Real arrays
    #
    A = Numeric.array(range(N)).astype('f')
    pypar.send(A,1)
    B=pypar.receive(numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Simplified communication of numeric real arrays OK"

    # Strings
    #
    A = "and now to something completely different !"
    pypar.send(A,1)
    B=pypar.receive(numproc-1)
    
    assert A == B
    print "Simplified communication of strings OK"
    
    # A more general structure
    #
    A = ['ABC', (1,2,3.14), {8: 'Monty'}, Numeric.array([13.45, 1.2])]
    pypar.send(A,1)
    B = pypar.receive(numproc-1)
    
    assert A == B
    print "Simplified communication of general structures OK"
    
  else:  
    # Integers
    #
    X=pypar.receive(myid-1)  
    pypar.send(X, (myid+1)%numproc)
  
    # Floats
    #
    X=pypar.receive(myid-1)  
    pypar.send(X, (myid+1)%numproc)    

    # Strings
    #
    X=pypar.receive(myid-1)  
    pypar.send(X, (myid+1)%numproc)    

    # General
    #
    X = pypar.receive(myid-1)  
    pypar.send(X, (myid+1)%numproc)    


  # Test broadcast  - with buffers (arrays, strings and general)
  #
  N = 17 #Number of elements
      
  testString = 'test' + str(myid)
  pypar.bcast(testString, 0)
  assert testString == 'test0'
  
  testString = 'test' + str(myid)
  pypar.bcast(testString, numproc-1)
  assert testString == 'test' + str(numproc-1)
  
  if myid == 0:
    print "Broadcast communication of strings OK"
  
  testArray = myid * Numeric.array(range(N))
  pypar.bcast(testArray, 1)
  assert Numeric.allclose(testArray, 1 * testArray)
  
  if myid == 0:    
    print "Broadcast communication of numeric integer array OK"


  testArray = myid * Numeric.array(range(N)).astype('f')
  pypar.bcast(testArray, 1)
  assert Numeric.allclose(testArray, 1 * testArray)
      
  if myid == 0:
    print "Broadcast communication of numeric real array OK"
    
  testGeneral = ['ABC', myid, (1,2,3), {8: 'Monty'}, Numeric.array([13.45, 1.2])]
  
  testGeneral = pypar.bcast(testGeneral, 1)
  
  assert testGeneral ==  ['ABC', 1, (1,2,3), {8: 'Monty'}, Numeric.array([13.45, 1.2])]
  
  if myid == 0:
    print "Broadcast communication of general structures OK"
  





  # Test scatter  - with/without buffers (arrays, strings)
  #
  N = 17 #Number of elements
      
  testString = 'test' + str(myid)
  s_size = 1   
  X = ' '*s_size
  pypar.raw_scatter(testString, s_size, X, 0)
  
  Y = pypar.scatter(testString, s_size, 0)
      
  if myid == 0:
    assert X == 't'
    assert Y == 't'
    print "Scatter communication of strings OK"

  testArray = Numeric.array(range(N))
  s_size = 1   
  X = Numeric.zeros(s_size)
  pypar.raw_scatter(testArray, s_size, X, 0)
  
  Y = pypar.scatter(testArray, s_size, 0)
  

  if myid == 0:
    assert X == [0]
    assert Y == [0]
    print "Scatter communication of numeric integer array OK"


  testArray = Numeric.array(range(N)).astype('f')
  s_size = 1   
  X = Numeric.zeros(s_size).astype('f')
  pypar.raw_scatter(testArray, s_size, X, 0)
    
  Y = pypar.scatter(testArray, s_size, 0)
  
  if myid == 0:
    assert X == [0.0]
    assert Y == [0.0]
    print "Scatter communication of numeric real array OK"



  # Test gather  - with/without buffers (arrays, strings)
  #
  N = 17 #Number of elements
      
  testString = 'AB'
  s_size = 2 # to help test
  X = ' '*(s_size*numproc)
  pypar.raw_gather(testString, s_size, X, 0, 0)

  if myid == 0:
    assert X == 'AB' * numproc

  Y =  pypar.gather(testString, s_size, 0) 
  #print myid, Y
  
  if myid == 0:
    assert X == 'AB' * numproc
    assert Y == 'AB' * numproc
    print "Gather communication of strings OK"
  

  testArray = Numeric.array(range(N))
  s_size = N   
  X = Numeric.zeros(s_size*numproc)
  pypar.raw_gather(testArray, s_size, X, 0, 0)

  Y = pypar.gather(testArray, s_size, 0)
  
  if myid == 0:
    for i in range(numproc):       
      assert Numeric.allclose(testArray, X[(i * s_size): ((i+1)*s_size)])
    print "Gather communication of numeric integer array OK"
    
    
  testArray = Numeric.array(range(N)).astype('f')
  s_size = N   
  X = Numeric.zeros(s_size*numproc).astype('f')
  pypar.raw_gather(testArray, s_size, X, 0, 0)
  
  Y = pypar.gather(testArray, s_size, 0)
    
  if myid == 0:
    for i in range(numproc):       
      assert Numeric.allclose(testArray, X[(i * s_size): ((i+1)*s_size)])
    print "Gather communication of numeric real array OK"
    

  # Test reduce  - with/without buffers (arrays, strings)
  #
  N = 17 #Number of elements
      
  testArray = Numeric.array(range(N))
  s_size = N # to help test
  X = Numeric.zeros(s_size)
  pypar.raw_reduce(testArray, X, s_size, pypar.mpi_SUM, 0, 0)
  if myid == 0:
    print 'SUM', X
  pypar.raw_reduce(testArray, X, s_size, pypar.mpi_MAX, 0, 0)
  if myid == 0:
    print 'MAX', X
  pypar.raw_reduce(testArray, X, s_size, pypar.mpi_MIN, 0, 0)
  if myid == 0:
    print 'MIN', X
  pypar.raw_reduce(testArray, X, s_size, pypar.mpi_PROD, 0, 0)
  if myid == 0:
    print 'PROD', X
  pypar.raw_reduce(testArray, X, s_size, pypar.mpi_LAND, 0, 0)
  if myid == 0:
    print 'LAND', X
  pypar.raw_reduce(testArray, X, s_size, pypar.mpi_BAND, 0, 0)
  if myid == 0:
    print 'BAND', X
  pypar.raw_reduce(testArray, X, s_size, pypar.mpi_LOR, 0, 0)
  if myid == 0:
    print 'LOR', X
  pypar.raw_reduce(testArray, X, s_size, pypar.mpi_BOR, 0, 0)
  if myid == 0:
    print 'BOR', X
  pypar.raw_reduce(testArray, X, s_size, pypar.mpi_LXOR, 0, 0)
  if myid == 0:
    print 'LXOR', X
  pypar.raw_reduce(testArray, X, s_size, pypar.mpi_BXOR, 0, 0)
  if myid == 0:
    print 'BXOR', X
  pypar.raw_reduce(testArray, X, s_size, pypar.mpi_MAXLOC, 0, 0)  
  if myid == 0:
    print 'MAXLOC', X
  pypar.raw_reduce(testArray, X, s_size, pypar.mpi_MINLOC, 0, 0)
  if myid == 0:
    print 'MINLOC', X
  pypar.raw_reduce(testArray, X, s_size, pypar.mpi_REPLACE, 0, 0)
  if myid == 0:
    print 'REPLACE', X
