#########################################################
#
#  Parallel rectangular model for WEAK SCALING / constant-timestep studies.
#
#  Based on run_parallel_rectangular.py, but instead of holding the domain
#  extent fixed (which shrinks the triangles, and hence the timestep, as the
#  grid is refined) this version GROWS the rectangle with the grid so that the
#  triangle size -- and therefore the explicit CFL timestep -- stays constant.
#
#  Why two things must be held constant for dt to be constant
#  ----------------------------------------------------------
#  The explicit shallow-water timestep is dt ~ CFL * dx / (|u| + sqrt(g*h)).
#  To keep dt fixed as the triangle count grows we hold constant:
#    1. the triangle size dx -- by setting length = width = cell_size * sqrtN,
#       so the number of triangles grows (4 * sqrtN^2) while each triangle keeps
#       the same edge length cell_size; and
#    2. the water depth h (which sets the gravity-wave speed sqrt(g*h)) -- by
#       normalising the sloped bed by the extent, elevation = -1 - 2*x/length,
#       so the bed/depth range is identical at every problem size.
#
#  With the default cell_size = 0.002, sqrtN = 1000 gives length = 2.0 and so
#  reproduces run_parallel_rectangular.py exactly at that size; larger sqrtN then
#  keeps dx (and dt) fixed instead of shrinking them.
#
#  Triangle counts (tris = 4 * sqrtN^2):
#    -sn 1000 ->  4,000,000 triangles, extent 2.0 x 2.0
#    -sn 2000 -> 16,000,000 triangles, extent 4.0 x 4.0  (same dx, same dt)
#
#  To run in parallel on 4 processes:
#
#    mpiexec -np 4 python -u run_parallel_rectangular_weak_scaling.py -sn 2000
#
#  Authors:
#  Linda Stals, Steve Roberts and Matthew Hardy - June 2005
#  Steve Roberts - 2018
#
#########################################################

import time
import sys
import math
from xml import dom
import anuga
import numpy as np


#----------------------------
# Sequential interface
#---------------------------
from anuga import rectangular_cross
from anuga import Domain, Mesh
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

# Fixed triangle edge length (resolution). The extent is derived from this and
# sqrtN, so this -- not the extent -- is what stays constant across problem
# sizes, keeping the timestep constant. 0.002 makes -sn 1000 give extent 2.0.
cell_size = 0.002

yieldstep = 0.005
finaltime = 0.015

fixed_flux_timestep = 0.0

import argparse
parser = argparse.ArgumentParser(description='Rectangular (weak scaling / constant timestep)')

parser.add_argument('-ft', '--finaltime', type=float, default=finaltime,
                    help='finaltime')
parser.add_argument('-ys', '--yieldstep', type=float, default=yieldstep,
                    help='yieldstep')
parser.add_argument('-sn', '--sqrtN', type=int, default=sqrtN,
                    help='Size of grid: 1000 -> 4,000,000 triangles (tris = 4*sn^2)')
parser.add_argument('-cs', '--cell_size', type=float, default=cell_size,
                    help='Fixed triangle edge length; extent = cell_size*sqrtN '
                         '(holds dx, and hence the timestep, constant). '
                         'Default 0.002 reproduces run_parallel_rectangular.py at -sn 1000')
parser.add_argument('-gl', '--ghost_layer', type=int, default=2,
                    help='Size of ghost layer')

parser.add_argument('-fdt', '--fixed_dt', type=float, default=fixed_flux_timestep,
                    help='Set a fixed flux timestep')
parser.add_argument('-ta', '--test_allreduce', action='store_true',
                    help='run fixed timestep with dummy allreduce')
parser.add_argument('-mp', '--multi_processor_mode', type=int, default=1,
                    help='set multiprocessor mode in [1,2]')

parser.add_argument('-sww', '--store_sww', action='store_true', help='store sww files')

parser.add_argument('-v', '--verbose', action='store_true', help='turn on verbosity')

parser.add_argument('-ve', '--evolve_verbose', action='store_true', help='turn on evolve verbosity')

parser.add_argument('-ps', '--partition_scheme', type=str, default='metis',
                    help='set partition scheme in [metis, morton, hilbert]')

parser.add_argument('-ro', '--reorder', type=str, default='none',
                    choices=['none', 'hilbert', 'morton', 'rcm', 'metis', 'metis_hilbert', 'metis_rcm'],
                    help='reorder triangles for cache locality before evolving '
                         '(none, hilbert, morton, rcm, metis, metis_hilbert, metis_rcm)')

parser.add_argument('-rn', '--reorder_nprocs', type=int, default=None,
                    help='number of Metis partitions for metis/metis_hilbert/metis_rcm reordering '
                         '(defaults to OMP_NUM_THREADS if not set)')

args = parser.parse_args()

if myid == 0: print(args)

multi_processor_mode = args.multi_processor_mode
sqrtN = args.sqrtN
cell_size = args.cell_size
yieldstep = args.yieldstep
finaltime = args.finaltime
verbose = args.verbose
evolve_verbose = args.evolve_verbose
fixed_flux_timestep = args.fixed_dt
test_allreduce = args.test_allreduce
store_sww = args.store_sww
reorder_method = args.reorder

#----------------------------
# Grow the extent with the grid so the triangle size (and thus the timestep)
# is independent of sqrtN. length/sqrtN == cell_size is the invariant.
#----------------------------
length = cell_size * sqrtN
width  = cell_size * sqrtN

dist_params = {}
dist_params['ghost_layer_width'] = args.ghost_layer
dist_params['partition_scheme'] = args.partition_scheme


if fixed_flux_timestep == 0.0:
    fixed_flux_timestep = None

#print('fixed_flux_timestep ',fixed_flux_timestep)




#--------------------------------------------------------------------------
# Setup Domain only on processor 0
#--------------------------------------------------------------------------
if myid == 0:

    domain = rectangular_cross_domain(sqrtN, sqrtN,
                                      len1=length, len2=width,
                                      origin=(-length/2, -width/2), verbose=verbose)


    domain.set_store(store_sww)
    # Normalise the sloped bed by the extent so the bed (and hence the water
    # depth, and hence the wave speed sqrt(g*h)) range is identical at every
    # problem size. At length == 2 this is exactly the original -1 - x.
    domain.set_quantity('elevation', lambda x,y : -1.0 - 2.0*x/length )
    domain.set_quantity('stage', 1.0)
    domain.set_flow_algorithm('DE0')
    domain.set_name('sw_rectangle_weak_scaling')

    domain.set_multiprocessor_mode(multi_processor_mode)

    if verbose: domain.print_statistics()

else:
    domain = None

t1 = time.time()

creation_time = t1-t0

if myid == 0 :
    print ('Creation of sequential domain: Time =',t1-t0)
    print ('Creation of sequential domain: Number of Triangles =',domain.number_of_global_triangles)
    print ('Extent = %g x %g, cell_size (dx) = %g' % (length, width, cell_size))

if myid == 0:
    print ('DISTRIBUTING DOMAIN')
    sys.stdout.flush()

barrier()

#-------------------------------------------------------------------------
# Distribute domain
#-------------------------------------------------------------------------
domain = distribute(domain,verbose=verbose, parameters=dist_params)


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

#-------------------------------------------------------------------------
# Reorder triangles within each rank for cache locality (optional)
#-------------------------------------------------------------------------
if reorder_method != 'none':
    if myid == 0:
        print('REORDERING DOMAIN using %s' % reorder_method)
    anuga.reorder_domain(domain, method=reorder_method,
                         n_procs=args.reorder_nprocs, verbose=verbose)

#Boundaries
T = Transmissive_boundary(domain)
R = Reflective_boundary(domain)


domain.set_boundary( {'left': R, 'right': R, 'bottom': R, 'top': R} )


if myid == 0 : print ('After set_boundary')

# A fixed-size central disturbance: physically identical at every problem size
# (it stays well resolved -- 0.5/cell_size cells across -- and bounded, so it
# does not change the timestep as the domain grows).
setter = Set_stage(domain,center=(0.0,0.0), radius=0.5, stage = 2.0)

# evaluate setter
setter()

if myid == 0 : print ('After set quantity')

barrier()

t0 = time.time()

if myid == 0 :
    anuga.print_domain_memory_stats(domain)

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
    # Report the realised flux timestep so constant-dt can be confirmed across
    # problem sizes (it should not change with sqrtN for this script).
    dt_min = getattr(domain, 'recorded_min_timestep', None)
    dt_max = getattr(domain, 'recorded_max_timestep', None)
    if dt_min is not None and dt_max is not None:
        print ('Recorded flux timestep: min = %g, max = %g' % (dt_min, dt_max))

if myid == 0 :
    anuga.print_domain_memory_stats(domain)

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
domain.sww_merge(delete_old=True)


if myid == 0:
    print(80*'=')
    print('np,ntri,cell_size,extent,ctime,dtime,etime')
    msg = "%d,%d,%g,%g,%f,%f,%f"% (numprocs, domain.number_of_global_triangles,
                                   cell_size, length, creation_time, distribute_time, evolve_time)
    print(msg)

finalize()
