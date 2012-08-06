# Timing of MPI module for Python and estimation of latency and bandwidth
# 
# Send numerical array in a ring from processor 0 to 1 etc back to 0 
# Perform timings and compare different sending strategies 
#
# OMN, OCT 2001

import time, sys, pypar, numpy

# The send/recv routines
#
import pypar

#from mpiext import send_array, receive_array  #most direct methods
#Not available in final install - use bypass form for now

#--------------------------------------------------------------

method = 2   # Use 
             # 0: automatically allocated buffers
             # 1: user-supplied buffers
             # 2: Use bypass - let pypar use direct mpi_ext call
             # 3 (not there): direct call to mpiext (with buffers),
             #    only fractionally better than bypass
                  
vanilla = False  # Force use of vanilla format (1)
consistency_check = False # Check correctness


#--------------------------------------------------------------
# linfit
#
def linfit(x, y):
  """Fit a and b to the model y = ax + b. Return a,b,variance
  """
  
  Sx = Sy = SSoN = SxoN = norm = varest = 0.0
  N = len(x)
  assert len(y) == N, "x and y must have same length"
  
  for i in range(N):
    #print("x,y = %f, %f\n",x[i],y[i])
    Sx  = Sx + x[i]
    Sy  = Sy + y[i]
  
  SxoN = Sx/N
  
  a = 0.0 
  for i in range(N):
    t    = x[i] - SxoN
    SSoN = SSoN + t*t
    a    = a + t*y[i]


  a = a/SSoN            # a = (N Sxy - SxSy)/(NSxx - Sx^2) */
  b = (Sy - Sx*a)/N
  
  # Quality - variance estimate \sum_i r_i^2 /(m-n) 
  for i in range(N):  
    norm = norm + float(x[i])*x[i]
    res = y[i] - a*x[i] - b
    varest = varest + res*res

  varest = varest/norm/(N-2)
  return a, b, varest




#--------------------------------------------------------------
# Main program
#
MAXI  = 10         # Number of blocks 
MAXM  = 500000     # Largest block 
BLOCK = MAXM/MAXI  # Block size 

repeats = 10
msgid = 0
vanilla = 0 #Select vanilla mode (slower but general)

numprocs = pypar.size()
myid = pypar.rank()
processor_name = pypar.get_processor_name()

if myid == 0:
  # Main process - Create message, pass on, verify correctness and log timing
  #
  print "MAXM = %d, number of processors = %d" %(MAXM, numprocs)
  print "Measurements are repeated %d times for reliability" %repeats

if numprocs < 2:
  print "Program needs at least two processors - aborting\n"
  pypar.abort()
   
pypar.barrier() #Synchronize all before timing   
print "I am process %d on %s" %(myid,processor_name)


#Initialise data and timings
#

try:
  from numpy.random import uniform, seed
  seed(17)
  A = uniform(0.0,100.0,MAXM)
except:
  print 'problem with RandomArray'
  from numpy import ones, Float
  A = ones(MAXM).astype('f')
  
elsize = A.itemsize
#print elsize

noelem  = [0]*MAXI
bytes   = [0]*MAXI         
avgtime = [0.0]*MAXI         
mintime = [ 1000000.0]*MAXI      
maxtime = [-1000000.0]*MAXI            




if myid == 0:   
  # Determine timer overhead 
  cpuOH = 1.0;
  for k in range(repeats):   # Repeat to get reliable timings 
    t1 = pypar.time()
    t2 = pypar.time()
    if t2-t1 < cpuOH: cpuOH = t2-t1
    
  print "Timing overhead is %f seconds.\n" %cpuOH         

     
# Pass msg circularly   
for k in range(repeats):
  if myid == 0:
    print "Run %d of %d" %(k+1,repeats)
    
  for i in range(MAXI):
    m=BLOCK*i+1       
   
    noelem[i] = m
   
    pypar.barrier() # Synchronize 
   
    if myid == 0:
      #
      # Main process
      #
      t1 = pypar.time()

      if method==0:
        pypar.send(A[:m], destination=1, tag=msgid, vanilla=vanilla)
        C = pypar.receive(numprocs-1, tag=msgid, vanilla=vanilla)
      elif method == 1:  
        pypar.send(A[:m], use_buffer=True, destination=1,
                   tag=msgid, vanilla=vanilla)
        C = pypar.receive(numprocs-1, buffer=A[:m], tag=msgid, vanilla=vanilla)
      elif method==2:
        pypar.send(A[:m], use_buffer=True, destination=1,
                   tag=msgid, vanilla=vanilla, bypass=True)
        C = pypar.receive(numprocs-1, buffer=A[:m], tag=msgid,
                          vanilla=vanilla, bypass=True)
      else:
        raise 'Unknown method'
        #send_array(A[:m], 1, msgid)    
        #stat = receive_array(A[:m], numprocs-1, msgid)
        #C = A[:m]
        
      t2 = pypar.time() - t1 - cpuOH
      t2 = t2/numprocs
      avgtime[i] = avgtime[i] + t2
      if t2 < mintime[i]: mintime[i] = t2
      if t2 > maxtime[i]: maxtime[i] = t2

      # Uncomment to verify integrity of data
      # However, this may affect accuracy of timings for some reason.
      #
      if consistency_check:
        assert numpy.alltrue(C == A[:m])
    else:
      #
      # Parallel process - get msg and pass it on
      #

      if method==0:
        C = pypar.receive(myid-1, tag=msgid, vanilla=vanilla)
        pypar.send(A[:m], destination=(myid+1)%numprocs,
                   tag=msgid, vanilla=vanilla)
      elif method==1:  
        C = pypar.receive(myid-1, buffer=A[:m], tag=msgid, vanilla=vanilla)
        pypar.send(A[:m], use_buffer=True, destination=(myid+1)%numprocs,
                   tag=msgid, vanilla=vanilla)                        
      elif method==2:
        # Use pypar bypass
        C = pypar.receive(myid-1, buffer=A[:m], tag=msgid,
                          vanilla=vanilla, bypass=True)
        pypar.send(A[:m], use_buffer=True, destination=(myid+1)%numprocs,
                   tag=msgid, vanilla=vanilla, bypass=True)
      else:
        raise 'Unknown method'
        # Use direct mpiext call
        #stat = receive_array(A[:m], myid-1, msgid)                
        #send_array(A[:m], (myid+1)%numprocs, msgid)


# Output stats
#
if myid == 0:
  print "Bytes transferred   time (micro seconds)"
  print "                    min        avg        max "     
  print "----------------------------------------------"
     
  for i in range(MAXI):
    avgtime[i] = avgtime[i]/repeats*1.0e6 #Average micro seconds
    mintime[i] = mintime[i]*1.0e6         #Min micro seconds       
    maxtime[i] = maxtime[i]*1.0e6         #Min micro seconds              
          
    m = noelem[i]
    bytes[i] = elsize*noelem[i]       
      
    print "%10d    %10d %10d %10d" %(bytes[i], mintime[i], avgtime[i], maxtime[i]) 


  Tbw, Tlat, varest = linfit(bytes, mintime)
  print "\nLinear regression on best timings (t = t_l + t_b * bytes):\n",
  print "  t_b = %f\n  t_l = %f" %(Tbw,Tlat)
  print "  Estimated relative variance = %.9f\n" %varest
       
  print "Estimated bandwith (1/t_b):  %.3f Mb/s" %(1.0/Tbw)   
  print "Estimated latency:           %d micro s" %int(mintime[0]-bytes[0]*Tbw)  

   
pypar.finalize()




  
  
