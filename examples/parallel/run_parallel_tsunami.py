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
import numpy
import anuga
import math


#----------------------------
# Parallel interface
#---------------------------
from anuga import distribute, myid, numprocs, finalize, barrier
from anuga import rectangular_cross_domain


#------------------------------
# Simulation parameters
#------------------------------
length = 2.0
width = 2.0
refinement_factor = 50
yieldstep = 0.2
finaltime = 2.5

verbose = True

#--------------------------------------------------------------------------
# Setup functions for topograpy etc
#--------------------------------------------------------------------------
scale_me=1.0

def topography(x,y):
	return (-x/2.0 +0.05*numpy.sin((x+y)*200.0))*scale_me

def stagefun(x,y):
    stge=-0.2*scale_me #+0.01*(x>0.9)
    #topo=topography(x,y)
    return stge#*(stge>topo) + (topo)*(stge<=topo)

#--------------------------------------------------------------------------
# Create domains
#--------------------------------------------------------------------------
t0 = time.time()

#--------------------------------------------------------------------------
# Setup Domain only on processor 0
#--------------------------------------------------------------------------
if myid == 0:

    cbrtN = int((numprocs)**(1.0/3.0)*refinement_factor)
    domain = rectangular_cross_domain(cbrtN, cbrtN,
                                      len1=length, len2=width, 
                                      origin=(-length/2, -width/2), 
                                      verbose=verbose)

    domain.set_minimum_allowed_height(0.01)

    domain.set_store(True)
    domain.set_quantity('elevation',topography)     # Use function for elevation
    domain.set_quantity('friction',0.03)            # Constant friction
    domain.set_quantity('stage', stagefun)          # Constant negative initial stage

    domain.set_name('rectangular_tsunami')

    domain.print_statistics()
else:
    domain = None

t1 = time.time()

if myid == 0 :
    print ('Create sequential domain ',t1-t0)

if myid == 0 and verbose: 
    print ('DISTRIBUTING DOMAIN')
    sys.stdout.flush()
    
barrier()

#-------------------------------------------------------------------------
# Distribute domain
#-------------------------------------------------------------------------
domain = distribute(domain,verbose=verbose)


t2 = time.time()

if myid == 0 :
    print ('Distribute domain ',t2-t1)
    
if myid == 0 : print ('after parallel domain')

#Boundaries
T = anuga.Transmissive_boundary(domain)
R = anuga.Reflective_boundary(domain)
D = anuga.Dirichlet_boundary([-0.1*scale_me,0.,0.])

domain.set_boundary( {'left': R, 'right': D, 'bottom': R, 'top': R} )

if myid == 0 : print ('after set_boundary')

barrier()

t0 = time.time()

#===========================================================================
# Main Evolve Loop
#===========================================================================
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    if myid == 0:
        domain.write_time()
        sys.stdout.flush()
						



for p in range(numprocs):
    barrier()
    if myid == p:
        print (50*'=')
        print ('P%g' %(myid))
        print ('That took %.2f seconds' %(time.time()-t0))
        print ('Communication time %.2f seconds'%domain.communication_time)
        print ('Reduction Communication time %.2f seconds'%domain.communication_reduce_time)
        print ('Broadcast time %.2f seconds'%domain.communication_broadcast_time)
        sys.stdout.flush()


if domain.number_of_global_triangles < 50000:
    if myid == 0 :
        print ('Create dump of triangulation for %g triangles' % domain.number_of_global_triangles)
    domain.dump_triangulation(filename="rectangular_cross_%g.png"% numprocs)


if myid==0:
    print (50*'=')
    print('Number of triangles:', domain.number_of_global_triangles)
finalize()
