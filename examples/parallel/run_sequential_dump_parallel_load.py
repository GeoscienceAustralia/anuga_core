#########################################################
#
#  Example of running a simple parallel model where the 
#  sequential domain is partitioned and dumped as files
#  via sequential_dump and read in via parallel sequential_load
#
#  Note: You can separate the sequential domain creation and dump 
#  in one script and the sequential_load and parallel domain creation
#  script to run the parallel evolution. THe intial sequential domain creation
#  can be run on a computer with a large memory. The partition files are 
#  stored in the directory Partitions
#
#  Run in parallel as follows (on 4 processors)
#
#  mpiexec -np 4 python run_sequential_dump_parallel_load.py
#
#  Note the use of "if myid == 0" to restrict some calculations 
#  to just one processor, in particular the creation and partition of a 
#  full domain on processor 0 which is then the partitions are
#  loaded on all the processors. 
#
#  Authors: 
#  Linda Stals, Steve Roberts and Matthew Hardy - June 2005
#  Steve Roberts - 2022. Adapted to use sequential dump and load
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
from anuga import sequential_distribute_dump, sequential_distribute_load



#---------------------------
# Domain paramters
#---------------------------
#sqrtN = int(math.sqrt(numprocs)*4)
sqrtN = 500



t0 = time.time()

verbose = False

domain_name = 'sw_dump_load'
partition_dir = 'Partitions'

#--------------------------------------------------------------------------
# Setup Domain only on processor 0
#--------------------------------------------------------------------------
if myid == 0:
    length = 2.0
    width = 2.0
    #dx = dy = 0.005
    #dx = dy = 0.00125
    domain = rectangular_cross_domain(sqrtN, sqrtN,
                                      len1=length, len2=width, 
                                      origin=(-length/2, -width/2), 
                                      verbose=verbose)


    domain.set_store(True)
    domain.set_quantity('elevation', lambda x,y : -1.0-x )
    domain.set_quantity('stage', 1.0)
    domain.set_flow_algorithm('DE0')
    domain.set_name(domain_name)
    if verbose: domain.print_statistics()
else:
    domain = None

t1 = time.time()

if myid == 0 :
    creation_time = t1-t0
    print ('Creation of sequential domain: Time =',t1-t0)
    print ('Creation of sequential domain: Number of Triangles =',domain.number_of_global_triangles)

if myid == 0: 
    print ('Dumping partition')
    sys.stdout.flush()
    
barrier()

#-------------------------------------------------------------------------
# Distribute domain
#-------------------------------------------------------------------------
if myid == 0:
    sequential_distribute_dump(domain,numprocs=numprocs, verbose=verbose, partition_dir=partition_dir)


if myid == 0: 
    print ('Loading partitions')
    sys.stdout.flush()

barrier()

domain = sequential_distribute_load(filename=domain_name, partition_dir=partition_dir)


t2 = time.time()

if myid == 0 :
    distribute_time = t2-t1
    print ('Dump and Load Domain: Time ',distribute_time)


#Boundaries
T = Transmissive_boundary(domain)
R = Reflective_boundary(domain)


domain.set_boundary( {'left': R, 'right': R, 'bottom': R, 'top': R, 'ghost': None} )


if myid == 0 : print ('After set_boundary')

# Let's use a setter to set stage
setter = Set_stage(domain,center=(0.0,0.0), radius=0.5, stage = 2.0)

# evaluate setter
setter()

if myid == 0 : print ('After set quantity')

yieldstep = 0.005
finaltime = 0.015

barrier()

t0 = time.time()

#===========================================================================
# Main Evolve Loop
#===========================================================================
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    if myid == 0:
        domain.write_time()
        sys.stdout.flush()



evolve_time = time.time()-t0

if myid == 0 :
    print ('Evolve: Time',evolve_time)

if verbose:
    for p in range(numprocs):
        barrier()
        if myid == p:
            print (50*'=')
            print ('P%g' %(myid))
            print ('That took %.2f seconds' %(evolve_time))
            print ('Communication time %.2f seconds'%domain.communication_time)
            print ('Reduction Communication time %.2f seconds'%domain.communication_reduce_time)
            print ('Broadcast time %.2f seconds'%domain.communication_broadcast_time)
            sys.stdout.flush()



if domain.number_of_global_triangles < 10:
    if myid == 0 :
        print ('Plot triangulation for %g triangles' % domain.number_of_global_triangles)
    domain.dump_triangulation(filename="rectangular_cross_%g.png"% numprocs)

domain.sww_merge(delete_old=True)

if myid == 0:
    print(50*'=')
    print('numprocs, no triangles, creation_time, distribute_time, evolve_time')
    
    msg = "%d,%d,%f,%f,%f"% (numprocs, domain.number_of_global_triangles, creation_time, distribute_time, evolve_time)
    
    print(msg)

finalize()
