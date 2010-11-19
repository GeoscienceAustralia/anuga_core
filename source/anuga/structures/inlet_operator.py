import anuga
import numpy
import math
import inlet

from anuga.utilities.system_tools import log_to_file


class Inlet_operator:
    """Inlet Operator - add water to an inlet.
    Sets up the geometry of problem
    
    Inherit from this class (and overwrite
    discharge_routine method for specific subclasses)
    
    Input: domain, Two points
    """ 

    counter = 0

    def __init__(self,
                 domain,
                 line,
                 Q = 0.0,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):
        
        self.domain = domain
        self.domain.set_fractional_step_operator(self)
        self.line = numpy.array(line, dtype='d')

        # should set this up to be a function of time and or space)
        self.Q = Q


        if description == None:
            self.description = ' '
        else:
            self.description = description
        

        if label == None:
            self.label = "inlet_%g" % Inlet_operator.counter
        else:
            self.label = label + '_%g' % Inlet_operator.counter


        self.verbose = verbose

        
        # Keep count of inlet operator
        Inlet_operator.counter += 1

        import pdb
        #pdb.set_trace()
        
        self.enquiry_point = 0.5*(self.line[0] + self.line[1])
        self.outward_vector = self.line
        self.inlet = inlet.Inlet(self.domain, self.line, verbose= verbose)

        self.set_logging(logging)

    def __call__(self):

        timestep = self.domain.get_timestep()
        
        Q = self.Q

        assert Q >= 0.0, 'removing water from an inlet!'


        # FIXME (SR): Might be nice to spread the over the inlet so that it is flat


        # Spread Q*timestep amount of water evenly over the inlet region

        areas = self.inlet.get_areas()
        stages = self.inlet.get_stages()

        stages_order = stages.argsort()

        summed_areas = numpy.zeros_like(areas)
        summed_Q = numpy.zeros_like(areas)

        for i,a in enumerate(areas[stages_order]):
            print i,a, stages[stages_order[i]]
            if i == 0:
                summed_areas[i] = a
                summed_Q[i] = stages[stages_order[i]] - stages[stages_order[i]]
            else:
                summed_areas[i] = summed_areas[i-1] + a
                summed_Q[i] = summed_Q[i-1]

        print summed_areas


            



        print stages_order
        print stages
        print areas
        print stages[stages_order]


        new_inlet_depth = self.inlet.get_average_depth() + (Q*timestep/self.inlet.get_area())
        self.inlet.set_depths(new_inlet_depth)

            

    def statistics(self):


        message  = '=====================================\n'
        message += 'Inlet Operator: %s\n' % self.label
        message += '=====================================\n'

        message += 'Description\n'
        message += '%s' % self.description
        message += '\n'
        
        inlet = self.inlet

        message += '-------------------------------------\n'
        message +=  'Inlet %i\n' % i
        message += '-------------------------------------\n'

        message += 'inlet triangle indices and centres\n'
        message += '%s' % inlet.triangle_indices
        message += '\n'
            
        message += '%s' % self.domain.get_centroid_coordinates()[inlet.triangle_indices]
        message += '\n'

        message += 'polyline\n'
        message += '%s' % inlet.polyline
        message += '\n'

        message += '=====================================\n'

        return message


    def print_statistics(self):

        print self.statistics()


    def print_timestepping_statistics(self):

        message = '---------------------------\n'
        message += 'Inlet report for %s:\n' % self.label
        message += '--------------------------\n'
        message += 'Q [m^3/s]: %.2f\n' % self.Q
        

        print message


    def set_logging(self, flag=True):

        self.logging = flag

        # If flag is true open file with mode = "w" to form a clean file for logging
        if self.logging:
            self.log_filename = self.label + '.log'
            log_to_file(self.log_filename, self.statistics(), mode='w')
            log_to_file(self.log_filename, 'time,Q')

            #log_to_file(self.log_filename, self.culvert_type)


    def timestepping_statistics(self):

        message  = '%.5f, ' % self.domain.get_time()
        message += '%.5f, ' % self.Q

        return message

    def log_timestepping_statistics(self):

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

