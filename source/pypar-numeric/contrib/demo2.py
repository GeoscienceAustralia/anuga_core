#!/usr/bin/env python
# Test of MPI module 'pypar' for Python
# Demonstrates a 'master-slave' parallel program
# 
# Run as 
#   python demo2.py
# or 
#   mpirun -np 4 demo2.py
#
# GPC, FEB 2002

try:
  import Numeric
except:
  raise 'Module Numeric must be present to run pypar'
  
try:
  import pypar
except:
  raise 'Module pypar must be present to run parallel'

print "Modules Numeric, pypar imported OK"

WORKTAG = 1
DIETAG =  2


def master():
    numCompleted = 0
    
    print "[MASTER]: I am processor %d of %d on node %s" %(MPI_myid, MPI_numproc, MPI_node)
    
    # start slaves distributing the first work slot
    for i in range(1, MPI_numproc): 
        work = workList[i]
        pypar.raw_send(work, i, WORKTAG) 
        print "[MASTER]: sent work '%s' to node '%d'" %(work, i)

    # dispach the remaining work slots on dynamic load-balancing policy
    # the quicker to do the job, the more jobs it takes
    for work in workList[MPI_numproc:]:
        result = '  '
        err = pypar.raw_receive(result, pypar.MPI_ANY_SOURCE, pypar.MPI_ANY_TAG) 
        print "[MASTER]: received result '%s' from node '%d'" %(result, err[1][0])
        numCompleted += 1
        pypar.raw_send(work, err[1][0], WORKTAG)
        print "[MASTER]: sent work '%s' to node '%d'" %(work, err[1][0])
    
    # all works have been dispatched out
    print  "[MASTER]: toDo : %d" %numWorks
    print  "[MASTER]: done : %d" %numCompleted
    
    # I've still to take into the remaining completions   
    while(numCompleted < numWorks): 
        result = '  '
        err = pypar.raw_receive(result, pypar.MPI_ANY_SOURCE, pypar.MPI_ANY_TAG) 
        print "[MASTER]: received (final) result '%s' from node '%d'" %(result, err[1][0])
        numCompleted += 1
        print "[MASTER]: %d completed" %numCompleted
        
    print "[MASTER]: about to terminate slaves"

    # say slaves to stop working
    for i in range(1, MPI_numproc): 
        pypar.raw_send('#', i, DIETAG) 
        print "[MASTER]: sent (final) work '%s' to node '%d'" %(0, i)
        
    return
    
def slave():

    print "[SLAVE %d]: I am processor %d of %d on node %s" %(MPI_myid, MPI_myid, MPI_numproc, MPI_node)

    while 1:
        result = ' '
        err = pypar.raw_receive(result, pypar.MPI_ANY_SOURCE, pypar.MPI_ANY_TAG) 
        print "[SLAVE %d]: received work '%s' with tag '%d' from node '%d'" %(MPI_myid, result, err[1][1], err[1][0])
       
        if (err[1][1] == DIETAG):
            print "[SLAVE %d]: received termination from node '%d'" %(MPI_myid, 0)
            return
        else:
            result = 'X'+result
            pypar.raw_send(result, 0)
            print "[SLAVE %d]: sent result '%s' to node '%d'" %(MPI_myid, result, 0)
            
       

if __name__ == '__main__':
    MPI_myid =    pypar.rank()
    MPI_numproc = pypar.size()
    MPI_node =    pypar.Get_processor_name()

    workList = ('_dummy_', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j')
    _workList = ('_dummy_', 'a', 'b', 'c', 'd', 'e', 'f')
    numWorks = len(workList) - 1

    if MPI_myid == 0:
        master()
    else:
        slave()

    pypar.Finalize()
    print "MPI environment finalized."
        	
