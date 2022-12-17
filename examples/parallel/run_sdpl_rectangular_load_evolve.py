#########################################################
#
#  run_sdpl_rectangular_load_evolve
#
#  Here sdpl stands for sequential dump, parallel load
#
#  Example of running a parallel simulation where the 
#  initial creation and partitioning of the domain in
#  done on a single processor, and the resulting sub-domains
#  are dumped to files via:
#  
#  python run_sdpl_rectangular_create_partition_dump.py -np N kwargs
#
#  (where kwargs are the extra command line args defining the run)
#  and then the evolve part of the simulation is run in 
#  parallel via:
#
#  mpiexec -np N  python -u run_sdpl_rectangular_load_evolve.py kwargs
#
#
#  Note the use of "if myid == 0" to restrict some calculations 
#  to just one processor, in particular the creation of a 
#  full domain on processor 0 which is then distributed to the
#  processors. 
#
#  Authors: 
#  Linda Stals, Steve Roberts and Matthew Hardy - June 2005
#  Steve Roberts - 2018 - 2022
#
#########################################################

import time
import sys
import math
from xml import dom
import anuga


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

parser.add_argument('-v', '--verbose', action='store_true', help='turn on verbosity')

parser.add_argument('-ve', '--evolve_verbose', action='store_true', help='turn on evolve verbosity')

parser.add_argument('-sww', '--store_sww', action='store_true', help='turn on storing sww file')

args = parser.parse_args()

if myid == 0: print(args)

sqrtN = args.sqrtN
yieldstep = args.yieldstep
finaltime = args.finaltime
verbose = args.verbose
evolve_verbose = args.evolve_verbose
fixed_flux_timestep = args.fixed_dt
test_allreduce = args.test_allreduce
ghost_layer = args.ghost_layer
store_sww = args.store_sww

ncpus = anuga.numprocs

dist_params = {}
dist_params['ghost_layer_width'] = ghost_layer

if fixed_flux_timestep == 0.0:
    fixed_flux_timestep = None

domain_name = f'rect_gl_{ghost_layer}_sqrtn_{sqrtN}_ncpus_{ncpus}'
partition_dir = 'Partitions'

#print('fixed_flux_timestep ',fixed_flux_timestep)

creation_time = 0.0

if myid == 0: 
    print ('Loading partitions')
    sys.stdout.flush()

barrier()

t1 = time.time()

domain = anuga.sequential_distribute_load(filename=domain_name, partition_dir=partition_dir)

t2 = time.time()

barrier()

if myid == 0 :
    distribute_time = t2-t1
    print ('Load Domain: Time ',distribute_time)



# FIXME: THis should be able to be set in the sequential domain
domain.set_fixed_flux_timestep(fixed_flux_timestep)
domain.set_CFL(1.0)
if myid == 0: 
    print('CFL ',domain.CFL)
    print('fixed_flux_timestep ',domain.fixed_flux_timestep)
domain.test_allreduce = test_allreduce


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
    msg = "%d,%d,%f,%f,%f"% (numprocs, domain.number_of_global_triangles, domain.creation_time, distribute_time, evolve_time)
    print(msg)

finalize()
