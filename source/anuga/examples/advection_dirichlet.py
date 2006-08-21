"""Simple example using the module for solving the advection equation

Initial condition is zero, boudary conditions are dirichlet at the left edge
and transmissive everywhere else
"""
import sys
from os import sep
sys.path.append('..'+sep+'pyvolution')

from mesh_factory import rectangular
from advection import Domain, Transmissive_boundary, Dirichlet_boundary

from Numeric import array

#Create basic mesh
points, vertices, boundary = rectangular(8, 8)

#Create advection domain with direction (1,-1)
# - Initial condition is zero by default
domain = Domain(points, vertices, boundary, velocity=[1.0, -1.0])
domain.smooth = False
domain.visualise = True

#Boundaries
T = Transmissive_boundary(domain)
D = Dirichlet_boundary(array([0.1]))

domain.set_boundary( {'left': D, 'right': T, 'bottom': T, 'top': T} )
domain.check_integrity()

###################
# Evolution
for t in domain.evolve(yieldstep = 0.02, finaltime = 2):
    domain.write_time()



