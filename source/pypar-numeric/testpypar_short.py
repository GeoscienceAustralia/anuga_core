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
# OMN, GPC FEB 2002

import os, sys
print 'PATH:', os.getenv('PATH')
print 'PYTHONPATH:', os.getenv('PYTHONPATH')
print 'Version:', sys.version
try:
  import Numeric
except Exception, e:
  msg = 'Module Numeric must be present to run pypar: %s' %e
  raise msg


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
  print testString
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
  






