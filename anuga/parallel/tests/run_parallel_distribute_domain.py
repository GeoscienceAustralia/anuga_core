#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

from builtins import range
from builtins import object
import unittest
import os
import sys

from anuga.utilities.system_tools import get_pathname_from_package

import numpy as num

from anuga.utilities import parallel_abstraction as pypar

#------------------------------------------
# anuga imports
#------------------------------------------
import anuga

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.utilities.util_ext        import double_precision
from anuga.utilities.norms           import l1_norm, l2_norm, linf_norm

from anuga import Domain
from anuga import Reflective_boundary
from anuga import Dirichlet_boundary
from anuga import Time_boundary
from anuga import Transmissive_boundary

from anuga import rectangular_cross
from anuga import create_domain_from_file


from anuga.parallel import distribute, myid, numprocs, finalize

#----------------------------------
# set up MPI to abort on error
#----------------------------------
from anuga.utilities.parallel_abstraction import global_except_hook
import sys
sys.excepthook = global_except_hook

#--------------------------------------------------------------------------
# Setup parameters
#--------------------------------------------------------------------------

mod_path = get_pathname_from_package('anuga.parallel')
mesh_filename = os.path.join(mod_path,'data','merimbula_10785_1.tsh')

verbose = False

#--------------------------------------------------------------------------
# Setup procedures
#--------------------------------------------------------------------------
class Set_Stage(object):
    """Set an initial condition with constant water height, for x < x0
    """

    def __init__(self, x0=0.25, x1=0.5, h=1.0):
        self.x0 = x0
        self.x1 = x1
        self.h = h

    def __call__(self, x, y):
        return self.h * ((x > self.x0) & (x < self.x1))



domain = create_domain_from_file(mesh_filename)
domain.set_quantity('stage', Set_Stage(756000.0, 756500.0, 2.0))

#--------------------------------------------------------------------------
# Create parallel domain if requested
#--------------------------------------------------------------------------

if numprocs > 1:
    if myid == 0 and verbose: print('DISTRIBUTING PARALLEL DOMAIN')
    domain = distribute(domain)

#------------------------------------------------------------------------------
# Setup boundary conditions
# This must currently happen *after* domain has been distributed
#------------------------------------------------------------------------------
domain.store = False
Br = Reflective_boundary(domain)      # Solid reflective wall
domain.set_boundary({'exterior': Br, 'open': Br})

#------------------------------------------------------------------------------
# Setup diagnostic arrays
#------------------------------------------------------------------------------
l1list = []
l2list = []
linflist = []
l1norm = num.zeros(3, float)
l2norm = num.zeros(3, float)
linfnorm = num.zeros(3, float)
recv_norm = num.zeros(3, float)

#------------------------------------------------------------------------------
# Evolution
#------------------------------------------------------------------------------
if numprocs > 1:
    if myid == 0 and verbose: print('PARALLEL EVOLVE')
else:
    if verbose: print('SEQUENTIAL EVOLVE')
    
for t in domain.evolve(yieldstep=1, finaltime=20):
    edges = domain.quantities['stage'].edge_values.take(num.flatnonzero(domain.tri_full_flag),axis=0)
    l1norm[0] = l1_norm(edges[:,0])
    l1norm[1] = l1_norm(edges[:,1])
    l1norm[2] = l1_norm(edges[:,2])
    l2norm[0] = l2_norm(edges[:,0])
    l2norm[1] = l2_norm(edges[:,1])
    l2norm[2] = l2_norm(edges[:,2])
    linfnorm[0] = linf_norm(edges[:,0])
    linfnorm[1] = linf_norm(edges[:,1])
    linfnorm[2] = linf_norm(edges[:,2])
    if numprocs > 1:
        l2norm[0] = pow(l2norm[0], 2)
        l2norm[1] = pow(l2norm[1], 2)
        l2norm[2] = pow(l2norm[2], 2)
        if myid == 0:
            #domain.write_time()

            #print edges[:,1]            
            for p in range(1, numprocs):
                recv_norm = pypar.receive(p)
                l1norm += recv_norm
                recv_norm = pypar.receive(p)
                l2norm += recv_norm
                recv_norm = pypar.receive(p)
                linfnorm[0] = max(linfnorm[0], recv_norm[0])
                linfnorm[1] = max(linfnorm[1], recv_norm[1])
                linfnorm[2] = max(linfnorm[2], recv_norm[2])

            l2norm[0] = pow(l2norm[0], 0.5)
            l2norm[1] = pow(l2norm[1], 0.5)
            l2norm[2] = pow(l2norm[2], 0.5)

            l1list.append(l1norm)                
            l2list.append(l2norm)
            linflist.append(linfnorm)                
        else:
            pypar.send(l1norm, 0)
            pypar.send(l2norm, 0)
            pypar.send(linfnorm, 0)
    else:
        #domain.write_time()
        l1list.append(l1norm)                
        l2list.append(l2norm)
        linflist.append(linfnorm)

# Store results in the appropriate file
if numprocs > 1:
    fid = open('distribute_domain_parallel.txt', 'w')
else: 
    fid = open('distribute_domain_sequential.txt', 'w')

for i in range(len(l1list)):
    fid.write('%f %f %f %f %f %f %f %f %f\n' % (l1list[i][0], l1list[i][1], l1list[i][2],
                                                l2list[i][0], l2list[i][1], l2list[i][2],
                                                linflist[i][0], linflist[i][1], linflist[i][2]))
    
fid.close()
