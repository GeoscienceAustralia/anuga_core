# To change this template, choose Tools | Templates
# and open the template in the editor.

import anuga
import numpy
import math
import anuga.structures.inlet

from anuga.utilities.system_tools import log_to_file
from anuga.structures.inlet_operator import Inlet_operator
import parallel_inlet
import pypar

class Parallel_Inlet_operator(Inlet_operator):
    """Inlet Operator - add water to an inlet.
    Sets up the geometry of problem

    Inherit from this class (and overwrite
    discharge_routine method for specific subclasses)

    Input: domain, Two points
    """

    # PETE: This only counts the number of inlets in processor?
    counter = 0

    def __init__(self,
                 domain,
                 line,
                 Q = 0.0,
                 description = None,
                 label = None,
                 logging = False,
                 master_proc = 0,
                 procs = None,
                 verbose = False):

        # TODO: Include statement that excludes non-participating process
        # PETE: Only set if domain actually contains the line, EXIT otherwise
        self.domain = domain
        self.domain.set_fractional_step_operator(self)
        self.line = numpy.array(line, dtype='d')
        self.master_proc = master_proc #PETE: id of the master processor that gathers global data

        if procs is None:
            self.procs = [self.master_proc]
        else:
            self.procs = procs

        self.myid = pypar.rank()

        # should set this up to be a function of time and or space)
        self.Q = Q

        # PETE: Have description mentioning the name of the processor
        if description == None:
            self.description = ' '
        else:
            self.description = description


        if label == None:
            self.label = "inlet_%g" % Inlet_operator.counter + "_P" + str(self.myid)
        else:
            self.label = label + '_%g' % Inlet_operator.counter + "_P" + str(self.myid)


        self.verbose = verbose

        # Keep count of inlet operator

        # Only master proc can update the static counter
        if self.myid == master_proc:
            Inlet_operator.counter += 1

        # PETE: Should the line be global or local? What if enquiry point is elsewhere
        # TODO: Must determine the location of the enquiry point
        self.enquiry_point = 0.5*(self.line[0] + self.line[1])
        # TODO: Check whether the current processor contains enquiry point, tell the other processors
        # who owns it

        self.outward_vector = self.line
        self.inlet = parallel_inlet.Parallel_Inlet(self.domain, self.line, master_proc = master_proc,
                                                    procs = procs, verbose= verbose)

        #TODO: Should the master processor do this?
        self.set_logging(logging)

    def __call__(self):

        volume = 0

        # PETE: The master proc calculates the volume
        if self.myid == self.master_proc:
            timestep = self.domain.get_timestep()

            t = self.domain.get_time()
            Q1 = self.update_Q(t)
            Q2 = self.update_Q(t + timestep)

            volume = 0.5*(Q1+Q2)*timestep

            assert 0.5*(Q1+Q2) >= 0.0, 'Q < 0: Water to be removed from an inlet!'

            #print "Volume to be removed from Inlet = " + str(volume)

        # PETE: this is ok, assume that the master proc for inlet operator is the same as that
        # for the the inlet itself, thus the other processes need not know the volume
        self.inlet.set_stages_evenly(volume)

        # Distribute volume evenly over all cells
        #self.inlet.set_depths_evenly(volume)

    def update_Q(self, t):
        """Virtual method allowing local modifications by writing an
        overriding version in descendant
        """
        # Only one processor should call this unless Q is parallelizable
        if callable(self.Q):
            Q = self.Q(t)[0]
        else:
            Q = self.Q

        return Q

    def statistics(self):

        message = ''

        inlet_stats = self.inlet.statistics()
        

        if self.myid == self.master_proc:
            message  = '=======================================\n'
            message += 'Parallel Inlet Operator: %s\n' % self.label
            message += '=======================================\n'
            
            message += 'Description\n'
            message += '%s' % self.description
            message += '\n'

            message += inlet_stats

            message += '=====================================\n'

        return message


    def print_statistics(self):

        print self.statistics()


    def print_timestepping_statistics(self):

        if self.myid == self.master_proc:
            message = '---------------------------------------------\n'
            message += 'Parallel Inlet report for %s:\n' % self.label
            message += '--------------------------------------------\n'
            message += 'Q [m^3/s]: %.2f\n' % self.Q
        
            print message


    def set_logging(self, flag=True):

        stats = self.statistics()
        self.logging = flag
        
        if self.myid == self.master_proc:
            # If flag is true open file with mode = "w" to form a clean file for logging
            # PETE: Have to open separate file for each processor
            if self.logging:
                self.log_filename = self.label + '.log'
                log_to_file(self.log_filename, stats, mode='w')
                log_to_file(self.log_filename, 'time,Q')

            #log_to_file(self.log_filename, self.culvert_type)


    def timestepping_statistics(self):

        message  = '%.5f, ' % self.domain.get_time()
        message += '%.5f, ' % self.Q

        return message

    def log_timestepping_statistics(self):
        
        if self.myid == self.master_proc:
            if self.logging:
                log_to_file(self.log_filename, self.timestepping_statistics())



    def set_Q(self, Q):

        self.Q = Q

    def get_Q(self):

        return self.Q


    def get_inlet(self):

        return self.inlet

    def get_line(self):

        return self.line

    def get_master_proc(self):
        return self.master_proc

    def __parallel_safe(self):
        return True
