#!/usr/bin/env python
#################################################################
# Master/Slave Parallel decomposition sample 
# 
# Run as 
#   python demo3.py
# or 
#   mpirun -np 2 demo3.py
# (perhaps try number of processors more than 2)
#################################################################
#
# To verify bandwidth of your architexture please 
# run pytiming (and ctiming) 
#
# OMN, GPC FEB 2002
#
#

import sys

try:
  import Numeric
except:
  raise 'Module Numeric must be present to run pypar'
  
try:
  import pypar
except:
  raise 'Module pypar must be present to run parallel'

sys.stderr.write("Modules Numeric, pypar imported OK\n")

WORKTAG = 1
DIETAG =  2


def master():
    numCompleted = 0
    
    sys.stderr.write("[MASTER]: I am processor %d of %d on node %s\n" %(MPI_myid, MPI_numproc, MPI_node))
    
    # start slaves distributing the first work slot
    for i in range(1, min(MPI_numproc, numWorks)): 
        work = workList[i]
        pypar.raw_send(work, i, WORKTAG) 
        sys.stderr.write("[MASTER]: sent work '%s' to node '%d'\n" %(work, i))

    # dispach the remaining work slots on dynamic load-balancing policy
    # the quicker to do the job, the more jobs it takes
    for work in workList[MPI_numproc:]:
        result = '  '
        err, status = pypar.raw_receive(result, pypar.any_source, pypar.any_tag, return_status=True) 
        #sys.stderr.write( "[MASTER]: received result '%s' from node '%d'\n" %(result, err[1][0]))
        sys.stderr.write("[MASTER]: received result '%s' from node '%d'\n" %(result, status.source))
        numCompleted += 1
        pypar.raw_send(work, status.source, WORKTAG)
        sys.stderr.write("[MASTER]: sent work '%s' to node '%d'\n" %(work, status.source))
    
    # all works have been dispatched out
    sys.stderr.write("[MASTER]: toDo : %d\n" %numWorks)
    sys.stderr.write("[MASTER]: done : %d\n" %numCompleted)
    
    # I've still to take into the remaining completions   
    while(numCompleted < numWorks): 
        result = '  '
        err, status = pypar.raw_receive(result, pypar.any_source, pypar.any_tag, return_status=True) 
        sys.stderr.write("[MASTER]: received (final) result '%s' from node '%d'\n" %(result, status.source))
        numCompleted += 1
        sys.stderr.write("[MASTER]: %d completed\n" %numCompleted)
        
    sys.stderr.write( "[MASTER]: about to terminate slaves\n")

    # say slaves to stop working
    for i in range(1, MPI_numproc): 
        pypar.raw_send('#', i, DIETAG) 
        sys.stderr.write("[MASTER]: sent (final) work '%s' to node '%d'\n" %(0, i))
        
    return
    
def slave():

    sys.stderr.write( "[SLAVE %d]: I am processor %d of %d on node %s\n" %(MPI_myid, MPI_myid, MPI_numproc, MPI_node))

    while 1:
        result = ' '
        err, status = pypar.raw_receive(result, pypar.any_source, pypar.any_tag, return_status=True) 
        sys.stderr.write("[SLAVE %d]: received work '%s' with tag '%d' from node '%d'\n"\
	      %(MPI_myid, result, status.tag, status.source))
       
        if (status.tag == DIETAG):
            sys.stderr.write("[SLAVE %d]: received termination from node '%d'\n" %(MPI_myid, 0))
            return
        else:
            result = 'X'+result
            pypar.raw_send(result, 0)
            sys.stderr.write("[SLAVE %d]: sent result '%s' to node '%d'\n" %(MPI_myid, result, 0))
            
       

if __name__ == '__main__':
    MPI_myid =    pypar.rank()
    MPI_numproc = pypar.size()
    MPI_node =    pypar.Get_processor_name()

    _workList = ('_dummy_', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j')
    workList = ('_dummy_', 'a', 'b', 'c')
    numWorks = len(workList) - 1
    
    
    #FIXME, better control here
    if MPI_numproc > numWorks or MPI_numproc < 2:
        pypar.Finalize()
	if MPI_myid == 0:
          sys.stderr.write("ERROR: Number of processors must be in the interval [2,%d].\n" %numWorks)
	  
        sys.exit(-1)

    if MPI_myid == 0:
        master()
    else:
        slave()

    pypar.Finalize()
    sys.stderr.write("MPI environment finalized.\n")
        	
