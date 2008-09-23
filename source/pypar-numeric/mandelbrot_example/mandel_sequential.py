#!/usr/bin/env python

"""Sequential program computing the Mandelbrot set.

Ole Nielsen, SUT 2003
"""

from mandelbrot import calculate_region
from mandelplot import plot
import time

# User definable parameters
kmax = 2**15  # Maximal number of iterations (=number of colors)
M = N = 700   # width = height = N (200, 400, 600, 700 are good)

# Region in complex plane
real_min = -2.0
real_max =  1.0
imag_min = -1.5
imag_max =  1.5

# Compute Mandelbrot set
t0 = time.time()
A = calculate_region(real_min, real_max, imag_min, imag_max, kmax, M, N)
print 'Computed region in %.2f seconds' %(time.time()-t0)

# Plot result
plot(A, kmax)




