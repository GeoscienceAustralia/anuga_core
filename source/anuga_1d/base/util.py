"""
This module contains various auxiliary function used by pyvolution.
"""


import numpy

def mean(x):
    return numpy.sum(x)/len(x)


def gradient(x0, x1, q0, q1):

    if q1-q0 != 0:
        a = (q1-q0)/(x1-x0)
    else:
        a = 0
        
    return a


