#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$15/06/2010 3:08:41 PM$"


import numpy

def uniform_mesh(n, x_0=0.0, x_1=1.0):
    """Create points, and boundary for a uniform mesh with n sub-interval
    ranging from x_0 to x_1
    """

    assert n>0


    points  = x_0 + (x_1 - x_0)*numpy.arange(n+1,dtype=numpy.float)/n
    boundary = {(0, 0): 'left', (n-1, 1): 'right'}

    return points, boundary
