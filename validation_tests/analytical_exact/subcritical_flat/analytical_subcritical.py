"""
Subcritical flow over flat surface

Steve Roberts, ANU 2014
"""
from numpy import zeros, ones



qA = 4.42  # This is the imposed momentum
hx = 2.0   # This is the water height downstream

def analytic_sol(x):    

    h = ones(len(x))*hx
    z = zeros(len(x))

    return h,z


