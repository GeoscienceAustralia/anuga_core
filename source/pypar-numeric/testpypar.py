#!/usr/bin/env python
# Test of MPI module 'pypar' for Python
# 
# Run as 
#   python testpypar.py
# or 
#   mpirun -np 2 testpypar.py
# (perhaps try number of processors more than 2)
#
# To verify bandwidth of your architecture please run pytiming (and ctiming) 
#
# OMN, GPC FEB 2002

try:
  import Numeric
except:
  raise 'Module Numeric must be present to run pypar'


#print "Importing pypar"
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

#print "Module pypar imported OK"
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

    # Long integer arrays
    #
    A = Numeric.array(range(N)).astype('l')
    B = Numeric.zeros(N)    
    pypar.raw_send(A,1)
    pypar.raw_receive(B,numproc-1)
    
    assert Numeric.allclose(A, B)
    print "Raw communication of numeric long integer arrays OK"

    # Real arrays
    #
    A = Numeric.array(range(N)).astype('f')
    B = Numeric.zeros(N).astype('f')    
    pypar.raw_send(A,1)
    pypar.raw_receive(B,numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Raw communication of numeric real arrays OK"

    # Complex arrays
    #
    A = Numeric.array(range(N)).astype('D')
    A += 3j
    B = Numeric.zeros(N).astype('D')
    pypar.raw_send(A, 1)
    B = pypar.raw_receive(B, numproc-1)

    assert Numeric.allclose(A, B)    
    print "Raw communication of numeric complex arrays OK"

    # Strings (< 256 characters)
    #
    A = "and now to something completely different !"
    B = " "*len(A)
    pypar.raw_send(A,1)
    B, status = pypar.raw_receive(B,numproc-1,return_status=True)

    assert A == B
    print "Raw communication of strings OK"
    
    # A more general structure
    #
    A = ['ABC', (1,2,3.14), {8: 'Monty'}, Numeric.array([13.45, 1.2])]
    B = ['   ', (0,0,0.0), {0: '     '}, Numeric.zeros(2).astype('f')]    
    pypar.raw_send(A,1)
    B, status = pypar.raw_receive(B,numproc-1, return_status=True)

    assert A == B
    print "Raw communication of general structures OK"
    
  else:  
    # Integers
    #
    X = Numeric.zeros(N)
    pypar.raw_receive(X, myid-1)  
    pypar.raw_send(X, (myid+1)%numproc)
  
    # Long integers
    #
    X = Numeric.zeros(N).astype('l')
    pypar.raw_receive(X, myid-1)  
    pypar.raw_send(X, (myid+1)%numproc)    

    # Floats
    #
    X = Numeric.zeros(N).astype('f')
    pypar.raw_receive(X, myid-1)  
    pypar.raw_send(X, (myid+1)%numproc)    

    # Complex
    #
    X = Numeric.zeros(N).astype('D')
    X = pypar.raw_receive(X, myid-1)  
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
    




  #Test (raw communication of) multi dimensional arrays

  M = 13  #Number of elements in dim 1
  N = 17  #Number of elements in higher dims
  
  if myid == 0:
    # 2D real arrays
    #
    A = Numeric.array(range(M*N)).astype('f')
    B = Numeric.zeros(M*N).astype('f')

    A = Numeric.reshape(A, (M,N))
    B = Numeric.reshape(B, (M,N))    
    
    pypar.raw_send(A,1)
    pypar.raw_receive(B,numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Raw communication of 2D real arrays OK"

    # 2D complex arrays
    #
    A = Numeric.array(range(M*N)).astype('D')
    B = Numeric.zeros(M*N).astype('D')

    A = Numeric.reshape(A, (M,N))
    B = Numeric.reshape(B, (M,N))    

    pypar.raw_send(A,1)
    pypar.raw_receive(B,numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Raw communication of 2D complex arrays OK"

    # 3D real arrays
    #
    A = Numeric.array(range(M*N*N)).astype('f')
    B = Numeric.zeros(M*N*N).astype('f')

    A = Numeric.reshape(A, (M,N,N))
    B = Numeric.reshape(B, (M,N,N))    
    
    pypar.raw_send(A,1)
    pypar.raw_receive(B,numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Raw communication of 3D real real arrays OK"

    # 4D real arrays
    #
    A = Numeric.array(range(M*N*N*M)).astype('f')
    B = Numeric.zeros(M*N*N*M).astype('f')

    A = Numeric.reshape(A, (M,N,N,M))
    B = Numeric.reshape(B, (M,N,N,M))    
    
    pypar.raw_send(A,1)
    pypar.raw_receive(B,numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Raw communication of 4D real real arrays OK"

    # 5D real arrays
    #
    A = Numeric.array(range(M*N*2*N*M)).astype('f')
    B = Numeric.zeros(M*N*2*N*M).astype('f')

    A = Numeric.reshape(A, (M,N,2,N,M))
    B = Numeric.reshape(B, (M,N,2,N,M))    
    
    pypar.raw_send(A,1)
    pypar.raw_receive(B,numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Raw communication of 5D real real arrays OK"

    # 5D double arrays
    #
    A = Numeric.array(range(M*N*2*N*M)).astype('d')
    B = Numeric.zeros(M*N*2*N*M).astype('d')

    A = Numeric.reshape(A, (M,N,2,N,M))
    B = Numeric.reshape(B, (M,N,2,N,M))    
    
    pypar.raw_send(A,1)
    pypar.raw_receive(B,numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Raw communication of 5D double arrays OK"

    # 5D complex arrays
    #
    A = Numeric.array(range(M*N*2*N*M)).astype('D')
    B = Numeric.zeros(M*N*2*N*M).astype('D')

    A = Numeric.reshape(A, (M,N,2,N,M))
    B = Numeric.reshape(B, (M,N,2,N,M))    
    
    pypar.raw_send(A,1)
    pypar.raw_receive(B,numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Raw communication of 5D complex arrays OK"
  else:  
    # 2D real arrays
    #
    X = Numeric.zeros(M*N).astype('f')
    X = Numeric.reshape(X, (M,N))
    
    pypar.raw_receive(X, myid-1)  
    pypar.raw_send(X, (myid+1)%numproc)
  
    # 2D complex arrays
    #
    X = Numeric.zeros(M*N).astype('D')
    X = Numeric.reshape(X, (M,N))

    X = pypar.raw_receive(X, myid-1)
    pypar.raw_send(X, (myid+1)%numproc)
  
    # 3D real arrays
    #
    X = Numeric.zeros(M*N*N).astype('f')
    X = Numeric.reshape(X, (M,N,N))
    
    pypar.raw_receive(X, myid-1)  
    pypar.raw_send(X, (myid+1)%numproc)

    # 4D real arrays
    #
    X = Numeric.zeros(M*N*N*M).astype('f')
    X = Numeric.reshape(X, (M,N,N,M))
    
    pypar.raw_receive(X, myid-1)  
    pypar.raw_send(X, (myid+1)%numproc)

    # 5D real arrays
    #
    X = Numeric.zeros(M*N*2*N*M).astype('f')
    X = Numeric.reshape(X, (M,N,2,N,M))
    
    pypar.raw_receive(X, myid-1)  
    pypar.raw_send(X, (myid+1)%numproc)

    # 5D double arrays
    #
    X = Numeric.zeros(M*N*2*N*M).astype('d')
    X = Numeric.reshape(X, (M,N,2,N,M))
    
    pypar.raw_receive(X, myid-1)  
    pypar.raw_send(X, (myid+1)%numproc)

    # 5D complex arrays
    #
    X = Numeric.zeros(M*N*2*N*M).astype('D')
    X = Numeric.reshape(X, (M,N,2,N,M))
    
    pypar.raw_receive(X, myid-1)  
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

    # Long integer arrays
    #
    A = Numeric.array(range(N)).astype('l')
    pypar.send(A,1)
    B=pypar.receive(numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Simplified communication of long integer real arrays OK"

    # Real arrays
    #
    A = Numeric.array(range(N)).astype('f')
    pypar.send(A,1)
    B=pypar.receive(numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Simplified communication of numeric real arrays OK"

    # Complex arrays
    #
    A = Numeric.array(range(N)).astype('D')
    A += 3j
    pypar.send(A,1)
    B=pypar.receive(numproc-1)

    assert Numeric.allclose(A, B)    
    print "Simplified communication of numeric complex arrays OK"


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
  
    # Long integers
    #
    X=pypar.receive(myid-1)  
    pypar.send(X, (myid+1)%numproc)    

    # Floats
    #
    X=pypar.receive(myid-1)  
    pypar.send(X, (myid+1)%numproc)    

    # Complex
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




    
  #Test (easy communication of) multi dimensional arrays

  M = 13  #Number of elements in dim 1
  N = 17  #Number of elements in higher dims
  
  if myid == 0:
    # 2D real arrays
    #
    A = Numeric.array(range(M*N)).astype('f')

    A = Numeric.reshape(A, (M,N))
    
    pypar.send(A,1)
    B = pypar.receive(numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Simplified communication of 2D real arrays OK"

    # 3D real arrays
    #
    A = Numeric.array(range(M*N*N)).astype('f')
    A = Numeric.reshape(A, (M,N,N))
    
    pypar.send(A,1)
    B=pypar.receive(numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Simplified communication of 3D real arrays OK"

    # 4D real arrays
    #
    A = Numeric.array(range(M*N*N*M)).astype('f')
    A = Numeric.reshape(A, (M,N,N,M))
    
    pypar.send(A,1)
    B=pypar.receive(numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Simplified communication of 4D real arrays OK"

    # 5D real arrays
    #
    A = Numeric.array(range(M*N*2*N*M)).astype('f')
    A = Numeric.reshape(A, (M,N,2,N,M))
    
    pypar.send(A,1)
    B=pypar.receive(numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Simplified communication of 5D real arrays OK"

    # 5D double arrays
    #
    A = Numeric.array(range(M*N*2*N*M)).astype('d')
    A = Numeric.reshape(A, (M,N,2,N,M))
    
    pypar.send(A,1)
    B=pypar.receive(numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Simplified communication of 5D double arrays OK"

    # 5D complex arrays
    #
    A = Numeric.array(range(M*N*2*N*M)).astype('D')
    A = Numeric.reshape(A, (M,N,2,N,M))
    
    pypar.send(A,1)
    B=pypar.receive(numproc-1)
    
    assert Numeric.allclose(A, B)    
    print "Simplified communication of 5D complex real arrays OK"

  else:  
    # 2D real arrays
    #
    
    X = pypar.receive(myid-1)  
    pypar.send(X, (myid+1)%numproc)
  
    # 3D real arrays
    #
    X=pypar.receive(myid-1)  
    pypar.send(X, (myid+1)%numproc)

    # 4D real arrays
    #
    X = pypar.receive(myid-1)  
    pypar.send(X, (myid+1)%numproc)

    # 5D real arrays
    #
    X=pypar.receive(myid-1)  
    pypar.send(X, (myid+1)%numproc)

    # 5D double arrays
    #
    X=pypar.receive(myid-1)  
    pypar.send(X, (myid+1)%numproc)

    # 5D complex arrays
    #
    X=pypar.receive(myid-1)  
    pypar.send(X, (myid+1)%numproc)

  # Test broadcast  - with buffers (arrays, strings and general)
  #
      
  testString = ('test' + str(myid)).ljust(10)  #Buffers must have the same length on all procs!
  pypar.bcast(testString, 0)
  assert testString.strip() == 'test0'
  
  testString = ('test' + str(myid)).ljust(10)  #Buffers must have the same length on all procs!
  pypar.bcast(testString, numproc-1)
  assert testString.strip() == 'test' + str(numproc-1)
  
  if myid == 0:
    print "Broadcast communication of strings OK"

  ####################################################  
  N = 17 #Number of elements
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


  M = 13
  testArray = myid * Numeric.array(range(M*N)).astype('f')
  testArray = Numeric.reshape(testArray, (M,N))  
  pypar.bcast(testArray, 1)
  assert Numeric.allclose(testArray, 1 * testArray)
  if myid == 0:
    print "Broadcast communication of 2D numeric real array OK"

  testArray = myid * Numeric.array(range(M*2*N)).astype('f')
  testArray = Numeric.reshape(testArray, (M,2,N))  
  pypar.bcast(testArray, 1)
  assert Numeric.allclose(testArray, 1 * testArray)
  if myid == 0:
    print "Broadcast communication of 3D numeric real array OK"

  testArray = myid * Numeric.array(range(M*2*N)).astype('D')
  testArray = Numeric.reshape(testArray, (M,2,N))  
  pypar.bcast(testArray, 1)
  assert Numeric.allclose(testArray, 1 * testArray)
  if myid == 0:
    print "Broadcast communication of 3D numeric complex array OK"

    
  testGeneral = ['ABC', myid, (1,2,3), {8: 'Monty'}, Numeric.array([13.45, 1.2])]
  
  testGeneral = pypar.bcast(testGeneral, 1)
  
  assert testGeneral ==  ['ABC', 1, (1,2,3), {8: 'Monty'}, Numeric.array([13.45, 1.2])], 'Should have been %s' %testGeneral
  
  if myid == 0:
    print "Broadcast communication of general structures OK"
  



  # Test scatter  - with/without buffers (arrays, strings)
  #
  N = 16 #Number of elements

  NP = N/numproc
  
  testString = 'ABCDEFGHIJKLMNOP'  #Length = 16
  X = ' '*NP

  #print 'P%d: s=%s, r=%s' %(myid, testString, X)
  #pypar.raw_scatter(testString, X, 2)
  #Y = pypar.scatter(testString, 2)


  #assert X==Y, 'X=%s, Y=%s' %(X,Y)
  #assert Y == testString[myid*NP:(myid+1)*NP]
  #assert X == testString[myid*NP:(myid+1)*NP]

  #if myid == 0:    
  #  print "Scatter communication of strings OK"
    

  #Scatter Arrays
  testArray = Numeric.array(range(N))
  X = Numeric.zeros(NP)
  pypar.raw_scatter(testArray, X, 0)
  Y = pypar.scatter(testArray, 0)

  assert Numeric.allclose(X, Y)  
  assert Numeric.allclose(X, testArray[myid*NP:(myid+1)*NP])
  assert Numeric.allclose(Y, testArray[myid*NP:(myid+1)*NP])   

  if myid == 0:
    print "Scatter communication of numeric integer array OK"


  testArray = Numeric.array(range(N)).astype('f')
  X = Numeric.zeros(NP).astype('f')
  pypar.raw_scatter(testArray, X, 0)
    
  Y = pypar.scatter(testArray, 0)

  assert Numeric.allclose(X, Y)  
  assert Numeric.allclose(X, testArray[myid*NP:(myid+1)*NP])
  assert Numeric.allclose(Y, testArray[myid*NP:(myid+1)*NP])   

  if myid == 0:
    print "Scatter communication of numeric real arrays OK"
  #else:
  #  print X, testArray, Y
  #  assert Numeric.allclose(X, Y)

  testArray = Numeric.array(range(N)).astype('D')
  X = Numeric.zeros(NP).astype('D')
  pypar.raw_scatter(testArray, X, 0)
    
  Y = pypar.scatter(testArray, 0)

  assert Numeric.allclose(X, Y)  
  assert Numeric.allclose(X, testArray[myid*NP:(myid+1)*NP])
  assert Numeric.allclose(Y, testArray[myid*NP:(myid+1)*NP])   

  if myid == 0:
    print "Scatter communication of numeric complex array OK"


###################################################
#FIXME: 2D scatter doesn't work yet    
#   M = 16
#   N = 13
#   MP = M/numproc
  
#   testArray = Numeric.array(range(M*N)).astype('D')
#   testArray = Numeric.reshape(testArray, (M,N))    
#   X = Numeric.zeros(MP*N).astype('D')
#   X = Numeric.reshape(X, (MP,N))
  
#   pypar.raw_scatter(testArray, X, 0)
#   Y = pypar.scatter(testArray, 0)
  
#   assert Numeric.allclose(X, Y)  
#   assert Numeric.allclose(X, testArray[myid*MP:(myid+1)*MP,:])
#   assert Numeric.allclose(Y, testArray[myid*MP:(myid+1)*MP,:])   

#   if myid == 0:
#     print "Scatter communication of 2D numeric complex array OK"



  # Test gather  - with/without buffers (arrays, strings)
  #
  N = 17 #Number of elements
      
  testString = 'AB'
  X = '_'*(len(testString)*numproc)    #Blanks caused errors when numproc >= 6 
  pypar.raw_gather(testString, X, 0)

  Y =  pypar.gather(testString, 0) 
  
  if myid == 0:
    assert X == 'AB' * numproc
    assert Y == 'AB' * numproc
    print "Gather communication of strings OK"
  

  testArray = Numeric.array(range(N))
  X = Numeric.zeros(N*numproc)
  pypar.raw_gather(testArray, X, 0)

  Y = pypar.gather(testArray, 0)
  
  if myid == 0:
    for i in range(numproc):       
      assert Numeric.allclose(testArray, X[(i * N): ((i+1)*N)])

    assert Numeric.allclose(X, Y)
    print "Gather communication of numeric integer array OK"
    
    
  testArray = Numeric.array(range(N)).astype('f')
  X = Numeric.zeros(N*numproc).astype('f')
  pypar.raw_gather(testArray, X, 0)
  
  Y = pypar.gather(testArray, 0)
  if myid == 0:
    for i in range(numproc):       
      assert Numeric.allclose(testArray, X[(i * N): ((i+1)*N)])
    assert Numeric.allclose(X, Y)      
    print "Gather communication of numeric real array OK"
    
  
  testArray = Numeric.array(range(N)).astype('D')
  X = Numeric.zeros(N*numproc).astype('D')
  pypar.raw_gather(testArray, X, 0)
  
  Y = pypar.gather(testArray, 0)
  if myid == 0:
    for i in range(numproc):       
      assert Numeric.allclose(testArray, X[(i * N): ((i+1)*N)])
    assert Numeric.allclose(X, Y)      
    print "Gather communication of numeric complex arrays OK"

  M = 13  
  testArray = Numeric.array(range(M*N)).astype('D')
  testArray = Numeric.reshape(testArray, (M,N))
  X = Numeric.zeros(M*N*numproc).astype('D')
  X = Numeric.reshape(X, (M*numproc,N))
  
  pypar.raw_gather(testArray, X, 0)
  
  Y = pypar.gather(testArray, 0)
  if myid == 0:
    for i in range(numproc):       
      assert Numeric.allclose(testArray, X[(i * M): ((i+1)*M), :])
    assert Numeric.allclose(X, Y)      
    print "Gather communication of 2D numeric complex arrays OK"
    
  
  ########################################################
  # Test reduce
  #######################################################
  N = 17 #Number of elements
  
  # Create one (different) array on each processor
  #    
  testArray = Numeric.array(range(N)) * (myid+1)
  #print testArray
  X = Numeric.zeros(N) # Buffer for results

  pypar.raw_reduce(testArray, X, pypar.SUM, 0)
  if myid == 0:
    Y = Numeric.zeros(N)
    for i in range(numproc):
      Y = Y+Numeric.array(range(N))*(i+1)    
    #print X
    #print Y  
    assert Numeric.allclose(X, Y)
    print "Raw reduce using pypar.SUM OK"
        
  pypar.raw_reduce(testArray, X, pypar.MAX, 0, 0)
  if myid == 0:
    Y = Numeric.array(range(N))*numproc
    assert Numeric.allclose(X, Y)
    print "Raw reduce using pypar.MAX OK"

  pypar.raw_reduce(testArray, X, pypar.MIN, 0, 0)
  if myid == 0:
    Y = Numeric.array(range(N))
    assert Numeric.allclose(X, Y)
    print "Raw reduce using pypar.MIN OK"
    
  if numproc <= 20:
    testArray_float = testArray.astype('f')  #Make room for big results
    X_float = X.astype('f')
    pypar.raw_reduce(testArray_float, X_float, pypar.PROD, 0, 0)
    if myid == 0:
      Y = Numeric.ones(N).astype('f')    
      for i in range(numproc):
        Y = Y*Numeric.array(range(N))*(i+1)    
      #print X_float
      #print Y  
      assert Numeric.allclose(X_float, Y)
      print "Raw reduce using pypar.PROD OK"
  else:
    if myid == 0:
      print "Skipping product-reduce - try again with numproc < 20"    

  pypar.raw_reduce(testArray, X, pypar.LAND, 0, 0)
  if myid == 0:  
    Y = Numeric.ones(N)    
    for i in range(numproc):
      Y = Numeric.logical_and(Y, Numeric.array(range(N))*(i+1))  
    assert Numeric.allclose(X, Y)
    print "Raw reduce using pypar.LAND OK"    
    
  pypar.raw_reduce(testArray, X, pypar.BAND, 0, 0)
  if myid == 0:
    Y = Numeric.ones(N)*255  #Neutral element for &   
    for i in range(numproc):
      Y = Numeric.bitwise_and(Y, Numeric.array(range(N))*(i+1))
    assert Numeric.allclose(X, Y)
    print "Raw reduce using pypar.BAND OK"    

  pypar.raw_reduce(testArray, X, pypar.LOR, 0, 0)
  if myid == 0:  
    Y = Numeric.zeros(N)    
    for i in range(numproc):
      Y = Numeric.logical_or(Y, Numeric.array(range(N))*(i+1))  
    assert Numeric.allclose(X, Y)
    print "Raw reduce using pypar.LOR OK"    
  
  pypar.raw_reduce(testArray, X, pypar.BOR, 0, 0)
  if myid == 0:
    Y = Numeric.zeros(N)   #Neutral element for |   
    for i in range(numproc):
      Y = Numeric.bitwise_or(Y, Numeric.array(range(N))*(i+1))
    assert Numeric.allclose(X, Y)
    print "Raw reduce using pypar.BOR OK"    

  pypar.raw_reduce(testArray, X, pypar.LXOR, 0, 0)
  if myid == 0:  
    Y = Numeric.zeros(N)    
    for i in range(numproc):
      Y = Numeric.logical_xor(Y, Numeric.array(range(N))*(i+1))  
    assert Numeric.allclose(X, Y)
    print "Raw reduce using pypar.LXOR OK"    

  pypar.raw_reduce(testArray, X, pypar.BXOR, 0, 0)
  if myid == 0:
    Y = Numeric.zeros(N)   #Neutral element for xor ?   
    for i in range(numproc):
      Y = Numeric.bitwise_xor(Y, Numeric.array(range(N))*(i+1))
    assert Numeric.allclose(X, Y)
    print "Raw reduce using pypar.BXOR OK"    

  # NOT YET SUPPORTED
  #  
  #pypar.raw_reduce(testArray, X, N, pypar.MAXLOC, 0, 0)  
  #if myid == 0:
  #  print 'MAXLOC', X
  #pypar.raw_reduce(testArray, X, N, pypar.MINLOC, 0, 0)
  #if myid == 0:
  #  print 'MINLOC', X
  
  #
  #  FIXME
  # Don't know how to test this (not available on all MPI systems)
  #
  #pypar.raw_reduce(testArray, X, N, pypar.REPLACE, 0, 0)
  #if myid == 0:
  #  print 'REPLACE', X



  # Test status block (simple communication)
  N = 17  
  if myid == 0:
    # Integer arrays
    #
    A = Numeric.array(range(N))

    pypar.send(A,1)
    B, status = pypar.receive(numproc-1, return_status = True)

    #print status

    sz = A.itemsize()
    assert Numeric.allclose(A, B)
    assert len(B) == status.length, 'Reported length == %d should be %d'\
           %(status.length, len(B))
    assert status.size == sz, 'Reported size == %d should be %d'\
           %(status.size, sz)
    assert status.tag == pypar.default_tag, 'Reported tag == %d should be %d'\
           %(status.tag, pypar.default_tag)
    assert status.error == 0
    assert status.source == numproc-1, 'Reported source == %d should be %d'\
           %(status.source, numproc-1)

    print "Status object (numeric integer arrays) OK"
           
    # Real arrays
    #
    A = Numeric.array(range(N)).astype('f')
    pypar.send(A,1)
    B, status = pypar.receive(numproc-1, return_status = True)    

    sz = A.itemsize()    
    assert Numeric.allclose(A, B)
    assert len(B) == status.length, 'Reported length == %d should be %d'\
           %(status.length, len(B))
    assert status.size == sz, 'Reported size == %d should be %d'\
           %(status.size, sz)
    assert status.tag == pypar.default_tag, 'Reported tag == %d should be %d'\
           %(status.tag, pypar.default_tag)
    assert status.error == 0
    assert status.source == numproc-1, 'Reported source == %d should be %d'\
           %(status.source, numproc-1)
    
    print "Status object (numeric real arrays) OK"

    # Strings
    #
    A = "and now to something completely different !"
    pypar.send(A,1)
    B, status = pypar.receive(numproc-1, return_status = True)        

    sz = 1 #Characters are one byte long
    assert A == B
    assert len(B) == status.length, 'Reported length == %d should be %d'\
           %(status.length, len(B))
    assert status.size == sz, 'Reported size == %d should be %d'\
           %(status.size, sz)
    assert status.tag == pypar.default_tag, 'Reported tag == %d should be %d'\
           %(status.tag, pypar.default_tag)
    assert status.error == 0
    assert status.source == numproc-1, 'Reported source == %d should be %d'\
           %(status.source, numproc-1)

    print "Status object (strings) OK"
    
    # A more general structure
    #
    A = ['ABC', (1,2,3.14), {8: 'Monty'}, Numeric.array([13.45, 1.2])]
    pypar.send(A,1)
    B, status = pypar.receive(numproc-1, return_status = True)            
    
    assert A == B

    #Length is the number of characters needed to encode the structure
    #Can't think of a test.

    sz = 1
    assert status.size == sz, 'Reported size == %d should be %d'\
           %(status.size, sz)
    assert status.tag == pypar.default_tag, 'Reported tag == %d should be %d'\
           %(status.tag, pypar.default_tag)
    assert status.error == 0
    assert status.source == numproc-1, 'Reported source == %d should be %d'\
           %(status.source, numproc-1)

    print "Status object (general structures) OK"
    
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








  # Test status block (raw communication)
  N = 17 #Number of elements
  if myid == 0:
    # Integer arrays
    #
    A = Numeric.array(range(N))
    B = Numeric.zeros(N)    
    pypar.raw_send(A,1)
    B, status = pypar.raw_receive(B,numproc-1,return_status=True)
    
    assert Numeric.allclose(A, B)

    sz = A.itemsize()
    assert Numeric.allclose(A, B)
    assert len(B) == status.length, 'Reported length == %d should be %d'\
           %(status.length, len(B))
    assert status.size == sz, 'Reported size == %d should be %d'\
           %(status.size, sz)
    assert status.tag == pypar.default_tag, 'Reported tag == %d should be %d'\
           %(status.tag, pypar.default_tag)
    assert status.error == 0
    assert status.source == numproc-1, 'Reported source == %d should be %d'\
           %(status.source, numproc-1)

    print "Status object (raw numeric integer arrays) OK"


    # Real arrays
    #
    A = Numeric.array(range(N)).astype('f')
    B = Numeric.zeros(N).astype('f')    
    pypar.raw_send(A,1)
    B, status = pypar.raw_receive(B,numproc-1,return_status=True)    
    
    assert Numeric.allclose(A, B)
    sz = A.itemsize()
    assert Numeric.allclose(A, B)
    assert len(B) == status.length, 'Reported length == %d should be %d'\
           %(status.length, len(B))
    assert status.size == sz, 'Reported size == %d should be %d'\
           %(status.size, sz)
    assert status.tag == pypar.default_tag, 'Reported tag == %d should be %d'\
           %(status.tag, pypar.default_tag)
    assert status.error == 0
    assert status.source == numproc-1, 'Reported source == %d should be %d'\
           %(status.source, numproc-1)

    print "Status object (raw numeric real arrays) OK"

    # Strings (< 256 characters)
    #
    A = "and now to something completely different !"
    B = " "*len(A)
    pypar.raw_send(A,1)
    B, status = pypar.raw_receive(B,numproc-1,return_status=True)


    sz = 1 #Characters are one byte long
    assert A == B
    assert len(B) == status.length, 'Reported length == %d should be %d'\
           %(status.length, len(B))
    assert status.size == sz, 'Reported size == %d should be %d'\
           %(status.size, sz)
    assert status.tag == pypar.default_tag, 'Reported tag == %d should be %d'\
           %(status.tag, pypar.default_tag)
    assert status.error == 0
    assert status.source == numproc-1, 'Reported source == %d should be %d'\
           %(status.source, numproc-1)

    print "Status object (raw strings) OK"
    

    
    # A more general structure
    #
    A = ['ABC', (1,2,3.14), {8: 'Monty'}, Numeric.array([13.45, 1.2])]
    B = ['   ', (0,0,0.0), {0: '     '}, Numeric.zeros(2).astype('f')]    
    pypar.raw_send(A,1)
    B, status = pypar.raw_receive(B,numproc-1, return_status=True)


    assert A == B
    sz = 1
    assert status.size == sz, 'Reported size == %d should be %d'\
           %(status.size, sz)
    assert status.tag == pypar.default_tag, 'Reported tag == %d should be %d'\
           %(status.tag, pypar.default_tag)
    assert status.error == 0
    assert status.source == numproc-1, 'Reported source == %d should be %d'\
           %(status.source, numproc-1)

    
    print "Status object (raw general structures) OK"
    
   
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
    X = pypar.raw_receive(X, myid-1)  
    pypar.raw_send(X.strip(), (myid+1)%numproc)    

    # General
    #
    X = ['   ', (0,0,0.0), {0: '     '}, Numeric.zeros(2).astype('f')]
    X = pypar.raw_receive(X, myid-1)  
    pypar.raw_send(X, (myid+1)%numproc)    
    


pypar.Finalize()




