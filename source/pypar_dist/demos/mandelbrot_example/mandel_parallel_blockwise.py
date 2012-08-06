"""Parallel program computing the Mandelbrot set using static partitioned 
load balancing.

   Due to the uneven work load of different regions,
   this strategy will *not* generally lead to linear speedup

   See mandel_parallel_dynamic.py for a better approach.

   Ole Nielsen, SUT 2003
"""

from mandelbrot import calculate_region, balance
from mandelplot import plot
import pypar

# User definable parameters
kmax = 2**15   # Maximal number of iterations (=number of colors)
M = N = 700    # width = height = N (200 or 700)

# Region in complex plane
real_min = -2.0
real_max =  1.0
imag_min = -1.5
imag_max =  1.5

#Initialise
t = pypar.time()
P = pypar.size()
p = pypar.rank()
processor_name = pypar.get_processor_name()

print 'Processor %d initialised on node %s' %(p,processor_name)


# Balanced work partitioning (row wise)
Mlo, Mhi = pypar.balance(M, P, p)
print 'p%d: [%d, %d], Interval length=%d' %(p, Mlo, Mhi, Mhi-Mlo)

# Parallel computation 
A = calculate_region(real_min, real_max, imag_min, imag_max, kmax,
                     M, N, Mlo = Mlo, Mhi = Mhi)

print 'Processor %d: time = %.2f' %(p, pypar.time() - t)


# Communication phase
if p == 0:
    for d in range(1, P):
        A += pypar.receive(source=d)

    print 'Computed region in %.2f seconds' %(pypar.time()-t)
    try:
        plot(A, kmax)
    except:
        pass    
else:
    pypar.send(A, destination=0)

pypar.finalize()                





