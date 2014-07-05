""" 
Towradgi Creek 17 August 1998 Storm Event Calibration
Ubuntu Linux 12.04 LTS 64 bit
By Petar Milevski, some revisions by Gareth Davies
"""

#------------------------------------------------------------------------------
# IMPORT NECESSARY MODULES
#------------------------------------------------------------------------------
import time
import numpy

from os.path import join
from anuga import file_function

import anuga

from anuga import myid, numprocs, finalize, barrier
from anuga import sequential_distribute_load

from project import *

#===============================================================================
# Start Simulation
#===============================================================================


args = anuga.get_args()
alg = args.alg
verbose = args.verbose

if myid == 0 and verbose: print 'STARTING PARALLEL SIMULATION'

if myid == 0:
    try:
        os.mkdir(model_output_dir)
    except:
        pass


    
#===============================================================================
# Create sequential domain and partition
#===============================================================================
if myid == 0 and verbose: print 'CREATING PARTITIONED DOMAIN'

if myid == 0:
    from setup_domain_and_partition import setup_domain, setup_partition
    
     
    setup_domain(verbose=verbose)
    
    setup_partition(np=numprocs, verbose=verbose)

barrier()

#===============================================================================
# Setup parallel domains
#===============================================================================
if myid == 0 and verbose: print 'LOADING PARTITIONED DOMAIN'

domain = sequential_distribute_load(filename=join('Partitions',outname), verbose = verbose)
print domain.get_name()
 
 


#===============================================================================
# Create structures such as culverts and bridges
#===============================================================================
if myid == 0 and verbose: print 'CREATING STRUCTURES'
execfile('setup_structures.py')

    
#===============================================================================
# Create rainfall functions associated with polygons
#===============================================================================
if myid == 0 and verbose: print 'CREATING RAINFALL FUNCTIONS'
from setup_rainfall import setup_rainfall

setup_rainfall(domain)

                      
 
#===========================================================================
# SetupBoundary conditions
#===========================================================================
if myid == 0 and verbose: print 'Setting up Boundaries'

func = file_function(join('Forcing','Tide','Pioneer.tms'), quantities='rainfall')
Bd = anuga.Dirichlet_boundary([0,0,0])
Bw = anuga.Time_boundary(domain=domain, function=lambda t: [func(t)[0], 0.0, 0.0])

domain.set_boundary({'west': Bd, 'south': Bd, 'north': Bd, 'east': Bw})

#------------------------------------------------------------------------------
# EVOLVE SYSTEM THROUGH TIME
#------------------------------------------------------------------------------
barrier()

if myid == 0 and verbose: print 'EVOLVE'
    
t0 = time.time()
    
for t in domain.evolve(yieldstep = 300., finaltime = 83700.):
    #if t == 37800.0: #time when bridge deck starts to get submerged, increase n to act as bridge deck, handrail and blockage effects
    ## Try to block all culverts / bridges, as described in the flood study
    #if t == 44100.0: #time when water level drops below bridge deck, bring n back down to existing conditions
    #    print 'Reset friction...'
    if myid == 0:
        domain.write_time()



barrier()
if myid == 0:
    print 'Number of processors %g ' %numprocs
    print 'That took %.2f seconds' %(time.time()-t0)
    print 'Communication time %.2f seconds'%domain.communication_time
    print 'Reduction Communication time %.2f seconds'%domain.communication_reduce_time
    print 'Broadcast time %.2f seconds'%domain.communication_broadcast_time



finalize()
