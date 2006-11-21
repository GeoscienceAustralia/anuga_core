"""Example using the module for solving the advection equation


Initial condition is zero, boundary conditions are transmissive everywhere except
for three segments where three different time dependent conditions are applied.
"""
import sys
from os import sep
sys.path.append('..'+sep+'abstract_2d_finite_volumes')

from mesh_factory import rectangular
from advection import Domain, Transmissive_boundary, Dirichlet_boundary, Time_boundary

#Create basic mesh
points, vertices, boundary = rectangular(20, 4)

#Create advection domain
# - Initial condition is zero by default
domain = Domain(points, vertices, boundary, velocity=[1.0, -0.8])
#domain.smooth = False

#Boundaries
T = Transmissive_boundary(domain)
D = Dirichlet_boundary([0.1])

from math import sin, pi
T1 = Time_boundary(domain, f=lambda x: [0.2*(sin(2*x*pi)+1)/2])
T2 = Time_boundary(domain, f=lambda x: [0.2*(sin(4*x*pi)+1)/3])
T3 = Time_boundary(domain, f=lambda x: [0.2*(sin(6*x*pi)+1)/4])

#Modify boundary tags
domain.boundary[(7, 1)] = 't1'
domain.boundary[(5, 2)] = 't2'
domain.boundary[(3, 2)] = 't3'

domain.set_boundary( {'left': T, 'right': T, 'bottom': T, 'top': T,
                      't1': T1, 't2': T2, 't3': T3} )
domain.check_integrity()


###################
# Evolution
for t in domain.evolve(yieldstep = 0.02, finaltime = None):
    domain.write_time()




