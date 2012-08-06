# coding: UTF-8
"""
Simple load balancing with pypar
(based on demo3.py from pypar demo package)
Felix Richter <felix.richter2@uni-rostock.de>
"""
import sys
import time

import numpy
import pypar

PYPAR_WORKTAG = 1
PYPAR_DIETAG = 2

def mprint(txt):
    """
    Print message txt
    with indentation following the node's rank
    """
    import pypar
    
    pre = " " * 8 * pypar.rank()
    if type(txt) != type('dummy'):
        txt = txt.__str__()
    pat = "-%d-"
    print pre + (pat % pypar.rank()) + txt

class PyparWork(object):
    """Abstract base class for ant work to be balanced"""

    def __init__(self):
        pass
    
    def uplink(self, balancer, myid, numprocs, node):
        self.balancer = balancer
        self.pypar_id = myid
        self.pypar_numprocs = numprocs
        self.pypar_node = node
    
    def getNumWorkItems(self):
        pass
    
    def handleWorkResult(self, result, status):
        pass
    
    def calcWorkResult(self, worknum):
        pass
    
    def masterBeforeWork(self):
        """Master node calls this before sending out the work"""
        pass

    def slaveBeforeWork(self):
        """Slave nodes call this before receiving work"""
        pass

    def masterAfterWork(self):
        """Master node calls this after receiving the last work result"""
        pass

    def slaveAfterWork(self):
        """Slave nodes call this after sending the last work result"""
        pass
    
    def msgprint(self, txt):
        pre = " " * 8 * self.pypar_id
        if type(txt) != type('dummy'):
            txt = txt.__str__()
        pat = "-%d-"
        print pre + (pat % self.pypar_id) + txt
        

class PyparBalancer(object):
    """The Load Balancer Class
    Initialize it with a PyparWork-derived class instance
    which describes the actual work to do.
    
    debug == True - more status messages
    """
    
    def __init__(self, work, debug = False):
        self.numprocs = pypar.size()           # Number of processes as specified by mpirun
        self.myid = pypar.rank()               # Id of of this process (myid in [0, numproc-1]) 
        self.node = pypar.get_processor_name() # Host name on which current process is running
        self.debug= debug
        self.work = work

        # Added by Ole Nielsen, 15 May 2008
        if self.numprocs < 2:
            msg = 'PyparBalancer must run on at least 2 processes'
            msg += ' for the Master Slave paradigm to make sense.'
            raise Exception, msg

        
        self.work.uplink(self, self.myid, self.numprocs, self.node)
        
        self.numworks = self.work.getNumWorkItems()
        print "PyparBalancer initialised on proc %d of %d on node %s" %(self.myid, self.numprocs, self.node)

    def master(self):
        numcompleted = 0
        #--- start slaves distributing the first work slot
        for i in range(0, min(self.numprocs-1, self.numworks)): 
            work = i
            slave= i+1
            pypar.send(work, destination=slave, tag=PYPAR_WORKTAG) 
            print '[MASTER ]: sent first work "%s" to node %d' %(work, slave)
    
        # dispatch the remaining work slots on dynamic load-balancing policy
        # the quicker to do the job, the more jobs it takes
        for work in range(self.numprocs-1, self.numworks):
            result, status = pypar.receive(source=pypar.any_source, tag=PYPAR_WORKTAG, return_status=True) 
            print '[MASTER ]: received result from node %d' %(status.source, )
            #print result
            numcompleted += 1
            pypar.send(work, destination=status.source, tag=PYPAR_WORKTAG)
            if self.debug: print '[MASTER ]: sent work "%s" to node %d' %(work, status.source)
            
            self.work.handleWorkResult(result, status)
        
        # all works have been dispatched out
        print '[MASTER ]: ToDo : %d' %self.numworks
        print '[MASTER ]: Done : %d' %numcompleted
        
        # I've still to take into the remaining completions   
        while (numcompleted < self.numworks): 
            result, status = pypar.receive(source=pypar.any_source, tag=PYPAR_WORKTAG, return_status=True) 
            print '[MASTER ]: received (final) result from node %d' % (status.source, )
            print result
            numcompleted += 1
            print '[MASTER ]: %d completed' %numcompleted
            
            self.work.handleWorkResult(result, status)
            
        print '[MASTER ]: about to terminate slaves'
    
        # Tell slaves to stop working
        for i in range(1, self.numprocs): 
            pypar.send('#', destination=i, tag=PYPAR_DIETAG) 
            if self.debug: print '[MASTER ]: sent DIETAG to node %d' %(i,)
            
    
    def slave(self):
        if self.debug: print '[SLAVE %d]: I am processor %d of %d on node %s' % (self.myid, self.myid, self.numprocs, self.node)
        if self.debug: print '[SLAVE %d]: Entering work loop' % (self.myid,)
        while True:
            result, status = pypar.receive(source=0, tag=pypar.any_tag, return_status=True) 
            print '[SLAVE %d]: received work with tag %d from node %d'\
                      %(self.myid, status.tag, status.source)
           
            if (status.tag == PYPAR_DIETAG):
                print '[SLAVE %d]: received termination from node %d' % (self.myid, 0)
                return
            else:
                worknum = result
                if self.debug: print '[SLAVE %d]: work number is %s' % (self.myid, worknum)
                myresult = self.work.calcWorkResult(worknum)
                pypar.send(myresult, destination=0)
                if self.debug: print '[SLAVE %d]: sent result to node %d' % (self.myid, 0)

    def run(self):
        if self.myid == 0:
            self.work.masterBeforeWork()
            self.master()
            self.work.masterAfterWork()
        else:
            self.work.slaveBeforeWork()
            self.slave()
            self.work.slaveAfterWork()
        
        pypar.finalize()
        if self.myid != 0:
            sys.exit()
        # und schluss.

class PyparDemoWork(PyparWork):
    """Example PyparWork implementation"""
    def __init__(self):
        import numpy
        self.worklist = numpy.arange(0.0,20.0)
        self.resultlist = numpy.zeros_like(self.worklist)
        
    def getNumWorkItems(self):
        return len(self.worklist)
    
    def calcWorkResult(self, worknum):
        return [worknum, self.worklist[worknum] + 1]

    def handleWorkResult(self, result, status):
        self.resultlist[result[0]] = result[1]
        
    def masterBeforeWork(self):
        print self.worklist

    def slaveBeforeWork(self):
        pass

    def masterAfterWork(self):
        print self.resultlist

    def slaveAfterWork(self):
        pass
        
if __name__ == "__main__":
    print "-----------------------"
    print "::: PyParBalancer TEST "
    print "-----------------------"
    
    # create instance of work class
    pyparwork = PyparDemoWork()

    # create instance of balancer class,
    # initialize with work class
    balancer = PyparBalancer(pyparwork, True)
    
    # run it
    balancer.run()
