#########################################################
#
#  Main file for parallel mesh testing.
#
#  This is a modification of the run_parallel_advection.py
# file.
#
#
#  Authors: Linda Stals, Steve Roberts and Matthew Hardy,
# June 2005
#
#
#
#########################################################

import time
import sys


#----------------------------
# Sequential interface
#---------------------------
from anuga import Transmissive_boundary, Reflective_boundary
from anuga import rectangular_cross_domain

#----------------------------
# Parallel interface
#---------------------------
from anuga_parallel.sequential_distribute import sequential_distribute


t0 = time.time()

verbose = True

#--------------------------------------------------------------------------
# Setup Domain only on processor 0
#--------------------------------------------------------------------------
myid = 0
numprocs = 100
if myid == 0:
    length = 2.0
    width = 2.0
    dx = dy = 0.005  # 640,000
    dx = dy = 0.00125
    #dx = dy  = 0.5
    domain = rectangular_cross_domain(int(length/dx), int(width/dy),
                                              len1=length, len2=width, verbose=verbose)


    print domain.number_of_global_triangles
    domain.set_store(True)
    domain.set_quantity('elevation', lambda x,y : -1.0-x )
    domain.set_quantity('stage', 1.0)
    domain.set_flow_algorithm('tsunami')
    domain.set_name('sw_rectangle')
    #domain.print_statistics()

else:
    domain = None

t1 = time.time()

if myid == 0 :
    print 'Create sequential domain ',t1-t0

if myid == 0 and verbose: 
    print 'DISTRIBUTING DOMAIN'
    sys.stdout.flush()
    
#barrier()

#-------------------------------------------------------------------------
# Distribute domain
#----------------------------------------------------------------------
#domain = distribute(domain,verbose=verbose)

sequential_distribute(domain,numprocs, verbose = True)


t2 = time.time()

if myid == 0 :
    print 'Distribute domain ',t2-t1
    
if myid == 0 : print 'after parallel domain'


