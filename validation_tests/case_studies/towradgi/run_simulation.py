
from simulation import Simulation
from anuga import myid, finalize, barrier, numprocs
import anuga
import time

args = anuga.get_args()
alg = args.alg
verbose = args.verbose

from project import *

towradgi = Simulation(outname=outname, verbose=verbose)

towradgi.setup_original_domain()
towradgi.setup_structures()
towradgi.setup_rainfall()
towradgi.setup_boundaries()

domain = towradgi.domain

if myid == 0 and verbose: print 'EVOLVE'
    
t0 = time.time()
    
for t in domain.evolve(yieldstep = 1., finaltime = 300):#= 83700.):

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
