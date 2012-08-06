"""
Master/Slave Parallel decomposition sample

Run as
   python demo3.py
or
   mpirun -np 2 demo3.py
  (perhaps try number of processors more than 2)
  

OMN, GPC FEB 2002
"""

import sys

try:
    import numpy
except:
    raise Exception, 'Module numpy must be present to run pypar'
  
try:
    import pypar
except:
    raise Exception, 'Module pypar must be present to run parallel'


print 'Modules numpy, pypar imported OK'

WORKTAG = 1
DIETAG =  2


def master():
    numCompleted = 0

    print '[MASTER]: I am processor %d of %d on node %s'\
          %(MPI_myid, MPI_numproc, MPI_node)

    # start slaves distributing the first work slot
    for i in range(1, min(MPI_numproc, numWorks)):
        work = workList[i]
        pypar.send(work, destination=i, tag=WORKTAG)
        print '[MASTER]: sent work "%s" to node %d' %(work, i)

    # dispatch the remaining work slots on dynamic load-balancing policy
    # the quicker to do the job, the more jobs it takes
    for work in workList[MPI_numproc:]:
        result, status = pypar.receive(source=pypar.any_source, tag=WORKTAG,
                                       return_status=True)

        print '[MASTER]: received result "%s" from node %d'\
              %(result, status.source)

        numCompleted += 1
        pypar.send(work, destination=status.source, tag=WORKTAG)
        print '[MASTER]: sent work "%s" to node %d' %(work, status.source)

    # all works have been dispatched out
    print '[MASTER]: toDo : %d' %numWorks
    print '[MASTER]: done : %d' %numCompleted

    # I've still to take into the remaining completions
    while (numCompleted < numWorks):
        result, status = pypar.receive(source=pypar.any_source,
                                       tag=WORKTAG,
                                       return_status=True)
        print '[MASTER]: received (final) result "%s" from node %d'\
                  %(result, status.source)
        numCompleted += 1
        print '[MASTER]: %d completed' %numCompleted

    print '[MASTER]: about to terminate slaves'

    # Tell slaves to stop working
    for i in range(1, MPI_numproc):
        pypar.send('#', destination=i, tag=DIETAG)
        print '[MASTER]: sent termination signal to node %d' %(i, )

    return

def slave():

    print '[SLAVE %d]: I am processor %d of %d on node %s'\
                     %(MPI_myid, MPI_myid, MPI_numproc, MPI_node)

    while True:
        result, status = pypar.receive(source=0,
                                       tag=pypar.any_tag,
                                       return_status=True)
        print '[SLAVE %d]: received work "%s" with tag %d from node %d'\
                  %(MPI_myid, result, status.tag, status.source)

        if (status.tag == DIETAG):
            print '[SLAVE %d]: received termination from node %d'\
                             %(MPI_myid, 0)
            return
        else:
            result = 'X'+result
            pypar.send(result, destination=0, tag=WORKTAG)
            print '[SLAVE %d]: sent result "%s" to node %d'\
            %(MPI_myid, result, 0)



if __name__ == '__main__':
    MPI_myid =    pypar.rank()
    MPI_numproc = pypar.size()
    MPI_node =    pypar.get_processor_name()

    workList = ('_dummy_', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j')
    numWorks = len(workList) - 1

    #FIXME, better control here
    if MPI_numproc > numWorks or MPI_numproc < 2:
        pypar.finalize()
        if MPI_myid == 0:
              print 'ERROR: Number of processors must be in the interval [2,%d].'%numWorks
              sys.exit(-1)

    if MPI_myid == 0:
        master()
    else:
        slave()

    pypar.finalize()
    print 'MPI environment finalized.'
    
