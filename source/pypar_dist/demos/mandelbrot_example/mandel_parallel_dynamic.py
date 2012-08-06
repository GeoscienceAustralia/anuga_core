"""Parallel program computing the Mandelbrot set using dynamic load balancing.
   Simplified version for use with exercise 15.

   To keep master process busy, run with P+1 processes
   where P is the number of physical processors.
   
   Ole Nielsen, SUT 2003
"""

from mandelbrot import calculate_region, balance
from mandelplot import plot
import pypar, numpy


# User definable parameters
kmax = 2**15  # Maximal number of iterations (=number of colors)
M = N = 700   # width = height = N
B = 24        # Number of blocks (first dim)


# Region in complex plane [-2:2]
real_min = -2.0
real_max =  1.0
imag_min = -1.5
imag_max =  1.5

# MPI controls
work_tag = 0
result_tag = 1

#Initialise
t = pypar.time()
P = pypar.size()
p = pypar.rank()
processor_name = pypar.get_processor_name()

print 'Processor %d initialised on node %s' %(p,processor_name)

assert P > 1, 'Must have at least one slave'
assert B > P-1, 'Must have more work packets than slaves'


A = numpy.zeros((M,N), dtype='i')
if p == 0:
    # Create work pool (B blocks)
    # using balanced work partitioning
    workpool = []
    for i in range(B):
        Mlo, Mhi = balance(M, B, i)    
        workpool.append( (Mlo, Mhi) )


    # Distribute initial work to slaves
    w = 0 
    for d in range(1, P):
        pypar.send(workpool[w], destination=d, tag=work_tag)
        w += 1

    #Receive computed work and distribute more
    terminated = 0
    while(terminated < P-1):
        R, status = pypar.receive(pypar.any_source, tag=result_tag,
                                  return_status=True)
        A += R            #Aggregate data        
        d = status.source #Id of slave that just finished        

        if w < len(workpool):
            #Send new work to slave d
            pypar.send(workpool[w], destination=d, tag=work_tag)
            w += 1
        else:
            #Tell slave d to terminate
            pypar.send(None, destination=d, tag=work_tag) 
            terminated += 1
        
    print 'Computed region in %.2f seconds' %(pypar.time()-t)
    try:
        plot(A, kmax)
    except:
        pass    
else:
    while(True):
        #Receive work (or None)
        W = pypar.receive(source=0, tag=work_tag)
        
        if W is None:
            print 'Slave p%d finished: time = %.2f' %(p, pypar.time() - t)
            break

        #Compute allocated work
        A = calculate_region(real_min, real_max, imag_min, imag_max, kmax,
                             M, N, Mlo = W[0], Mhi = W[1])

        #Return result
        pypar.send(A, destination=0, tag=result_tag)

pypar.finalize()







