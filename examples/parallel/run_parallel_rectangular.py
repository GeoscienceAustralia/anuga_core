#########################################################
#
#  Example of running a simple parallel model
#
#  Need mpi setup for your machine 
#
#  To run in parallel on 4 processes, use the following
#
#  mpiexec -np 4 python -u run_parallel_rectangular.py
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
from xml import dom
import anuga
import nvtx


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
refinement_factor = 100
sqrtN = int((numprocs)**(1.0/2.0)*refinement_factor)

sqrtN = 100
length = 2.0
width = 2.0

yieldstep = 0.005
finaltime = 0.015

fixed_flux_timestep = 0.0

import argparse
parser = argparse.ArgumentParser(description='Rectangular')

parser.add_argument('-ft', '--finaltime', type=float, default=finaltime,
                    help='finaltime')
parser.add_argument('-ys', '--yieldstep', type=float, default=yieldstep,
                    help='yieldstep')
parser.add_argument('-sn', '--sqrtN', type=int, default=sqrtN,
                    help='Size of grid: 500 -> 1_000_000 triangles')
parser.add_argument('-gl', '--ghost_layer', type=int, default=2,
                    help='Size of ghost layer')

parser.add_argument('-fdt', '--fixed_dt', type=float, default=fixed_flux_timestep,
                    help='Set a fixed flux timestep')
parser.add_argument('-ta', '--test_allreduce', action='store_true',
                    help='run fixed timestep with dummy allreduce')
parser.add_argument('-mp', '--multi_processor_mode', type=int, default=0,
                    help='set multiprocessor mode in [0,1,2,3,4]')

parser.add_argument('-v', '--verbose', action='store_true', help='turn on verbosity')

parser.add_argument('-ve', '--evolve_verbose', action='store_true', help='turn on evolve verbosity')

args = parser.parse_args()

if myid == 0: print(args)

multi_processor_mode = args.multi_processor_mode
sqrtN = args.sqrtN
yieldstep = args.yieldstep
finaltime = args.finaltime
verbose = args.verbose
evolve_verbose = args.evolve_verbose
fixed_flux_timestep = args.fixed_dt
test_allreduce = args.test_allreduce

dist_params = {}
dist_params['ghost_layer_width'] = args.ghost_layer

if fixed_flux_timestep == 0.0:
    fixed_flux_timestep = None

#print('fixed_flux_timestep ',fixed_flux_timestep)




#--------------------------------------------------------------------------
# Setup Domain only on processor 0
#--------------------------------------------------------------------------
if myid == 0:

    #nvtx marker
    rng = nvtx.start_range(message="rect_example_creat_time", color="blue")


    domain = rectangular_cross_domain(sqrtN, sqrtN,
                                      len1=length, len2=width, 
                                      origin=(-length/2, -width/2), 
                                      verbose=verbose)


    domain.set_store(True)
    domain.set_quantity('elevation', lambda x,y : -1.0-x )
    domain.set_quantity('stage', 1.0)
    domain.set_flow_algorithm('DE0')
    domain.set_name('sw_rectangle')

    domain.set_multiprocessor_mode(multi_processor_mode)
 
    if verbose: domain.print_statistics()
    # nvtx marker
    nvtx.end_range(rng)

else:
    domain = None

t1 = time.time()

creation_time = t1-t0

if myid == 0 :
    print ('Creation of sequential domain: Time =',t1-t0)
    print ('Creation of sequential domain: Number of Triangles =',domain.number_of_global_triangles)

if myid == 0: 
    print ('DISTRIBUTING DOMAIN')
    sys.stdout.flush()
    
barrier()

#-------------------------------------------------------------------------
# Distribute domain
#-------------------------------------------------------------------------
# nvtx marker
rng = nvtx.start_range(message="rectangular_exam_domain_distr", color="blue")

domain = distribute(domain,verbose=verbose,parameters=dist_params)
# nvtx marker
nvtx.end_range(rng)


# FIXME: THis should be able to be set in the sequential domain
domain.set_fixed_flux_timestep(fixed_flux_timestep)
domain.set_CFL(1.0)
if myid == 0: 
    print('CFL ',domain.CFL)
    print('fixed_flux_timestep ',domain.fixed_flux_timestep)
domain.test_allreduce = test_allreduce

t2 = time.time()

distribute_time = t2-t1

if myid == 0 :
    print ('Distribute domain: Time ',distribute_time)
    
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

# nvtx marker
rng = nvtx.start_range(message="rect_exam_evolve_time", color="blue")

#===========================================================================
# Main Evolve Loop
#===========================================================================
for t in domain.evolve(yieldstep = yieldstep, finaltime = finaltime):
    if myid == 0:
        domain.write_time()
        sys.stdout.flush()

# nvtx marker
nvtx.end_range(rng)
        
        
evolve_time = time.time()-t0

if myid == 0 :
    print ('Evolve: Time',evolve_time)

if evolve_verbose:
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
        print ('Create dump of triangulation for %g triangles' % domain.number_of_global_triangles)
    domain.dump_triangulation(filename="rectangular_cross_%g.png"% numprocs)

# to save time avoid merge
#domain.sww_merge(delete_old=True)


if myid == 0:
    print(80*'=')
    print('np,ntri,ctime,dtime,etime')
    msg = "%d,%d,%f,%f,%f"% (numprocs, domain.number_of_global_triangles, creation_time, distribute_time, evolve_time)
    print(msg)

finalize()
