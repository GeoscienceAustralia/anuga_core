"""Parallel program computing the Mandelbrot set using static,
   cyclic load balancing.

   This is probably the best approach for this problem.

   Ole Nielsen, SUT 2003
"""

from mandelbrot import calculate_region_cyclic
from mandelplot import plot
import pypar

# User definable parameters
kmax = 2**15  # Maximal number of iterations (=number of colors)
M = N = 700   # width = height = N (200 or 700)

# Region in complex plane [-2:2]
real_min = -2.0
real_max =  1.0
imag_min = -1.5
imag_max =  1.5

#Initialise
t = pypar.time()
P = pypar.size()
p = pypar.rank()
processor_name = pypar.get_processor_name()

print 'Processor %d initialised on node %s' %(p, processor_name)


# Parallel computation
A = calculate_region_cyclic(real_min, real_max, imag_min,
                            imag_max, kmax,
                            M, N, p, P)

print 'Processor %d: time = %.2f' %(p, pypar.time() - t)


# Communication phase
if p == 0:
    for d in range(1, P):
        print 'Proc 0 receiving from %i' % d
        A += pypar.receive(source=d)
        print 'Received from %i' % d

    print 'Computed region in %.2f seconds' %(pypar.time() - t)
    try:
        plot(A, kmax)
    except:
        pass

else:
    print 'Sending to 0'
    pypar.send(A, destination=0)

pypar.finalize()





