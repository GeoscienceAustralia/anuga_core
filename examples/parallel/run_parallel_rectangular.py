#########################################################
#
#  Example of running a simple parallel model
#
#  Need mpi setup for your machine 
#
#  To run in parallel on 4 processes, use the following
#
#  mpiexec -np 4 python -u run_parallel_sw_rectangular_cross.py
#
#
#  Note the use of "if myid == 0" to restrict some calculations 
#  to just one processor, in particular the creation of a 
#  full domain on processor 0 which is then distributed to the
#  processors. 
#
#  Authors: 
#  Linda Stals, Steve Roberts and Matthew Hardy - June 2005
#  Steve Roberts - 2018
#
#
#
#########################################################

import time
import sys
import math


#----------------------------
# Sequential interface
#---------------------------
from anuga import Transmissive_boundary, Reflective_boundary
from anuga import rectangular_cross_domain
from anuga import Set_stage

#----------------------------
# Parallel interface
#---------------------------
from anuga import distribute, myid, numprocs, finalize, barrier


t0 = time.time()

#----------------------------
# simulation parameters
#----------------------------
length = 2.0
width = 2.0
refinement_factor = 100
yieldstep = 0.005
finaltime = 0.05

verbose = True

#--------------------------------------------------------------------------
# Setup Domain only on processor 0
#--------------------------------------------------------------------------
if myid == 0:

    cbrtN = int((numprocs)**(1.0/3.0)*refinement_factor)
    domain = rectangular_cross_domain(cbrtN, cbrtN,
                                      len1=length, len2=width, 
                                      origin=(-length/2, -width/2), 
                                      verbose=verbose)


    domain.set_store(True)
    domain.set_quantity('elevation', lambda x,y : -1.0-x )
    domain.set_quantity('stage', 1.0)
    domain.set_flow_algorithm('DE0')
    domain.set_name('sw_rectangle')
    domain.print_statistics()
else:
    domain = None

t1 = time.time()

if myid == 0 :
    print ('Create sequential domain: Time',t1-t0)

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
    print ('Distribute domain: Time ',t2-t1)
    
if myid == 0 : print ('After parallel domain')

#Boundaries
T = Transmissive_boundary(domain)
R = Reflective_boundary(domain)


domain.set_boundary( {'left': R, 'right': R, 'bottom': R, 'top': R} )


if myid == 0 : print ('After set_boundary')

# Let's use a setter to set stage
setter = Set_stage(domain,center=(0.0,0.0), radius=0.5, stage = 2.0)

# evaluate setter
setter()

if myid == 0 : print ('After set quantity')

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



if domain.number_of_global_triangles < 5000:
    if myid == 0 :
        print ('Create dump of triangulation for %g triangles' % domain.number_of_global_triangles)
    domain.dump_triangulation(filename="rectangular_cross_%g.png"% numprocs)

domain.sww_merge(delete_old=True)



if myid==0:
    print(50*"=")
    print('Number of Triangles:', domain.number_of_global_triangles)

finalize()
