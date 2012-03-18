import anuga
import numpy
import inlet


class Inlet_operator(anuga.Operator):
    """Inlet Operator - add water to an inlet.
    Sets up the geometry of problem
    
    Inherit from this class (and overwrite
    discharge_routine method for specific subclasses)
    
    Input: domain, Two points
    """ 


    def __init__(self,
                 domain,
                 line,
                 Q = 0.0,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):


        anuga.Operator.__init__(self, domain, description, label, logging, verbose)

        
        self.line = numpy.array(line, dtype='d')

        # should set this up to be a function of time and or space)
        self.Q = Q


        self.enquiry_point = 0.5*(self.line[0] + self.line[1])
        self.outward_vector = self.line
        self.inlet = inlet.Inlet(self.domain, self.line, verbose= verbose)

        self.applied_Q = 0.0

    def __call__(self):

        timestep = self.domain.get_timestep()

        t = self.domain.get_time()

        Q1 = self.update_Q(t)
        Q2 = self.update_Q(t + timestep)

        Q = 0.5*(Q1+Q2)
        volume = Q*timestep
        
        assert volume >= 0.0, 'Q < 0: Water to be removed from an inlet!'

        # store last discharge
        self.applied_Q = Q

        # Distribute volume so as to obtain flat surface
        self.inlet.set_stages_evenly(volume)
        
        # Distribute volume evenly over all cells
        #self.inlet.set_depths_evenly(volume)
        
    def update_Q(self, t):
        """Allowing local modifications of Q
        """
        
        if callable(self.Q):
            Q = self.Q(t)[0]
        else:
            Q = self.Q

        return Q    
  
    def statistics(self):


        message  = '=====================================\n'
        message += 'Inlet Operator: %s\n' % self.label
        message += '=====================================\n'

        message += 'Description\n'
        message += '%s' % self.description
        message += '\n'
        
        inlet = self.inlet

        message += '-------------------------------------\n'
        message +=  'Inlet\n' 
        message += '-------------------------------------\n'

        message += 'inlet triangle indices and centres\n'
        message += '%s' % inlet.triangle_indices
        message += '\n'
            
        message += '%s' % self.domain.get_centroid_coordinates()[inlet.triangle_indices]
        message += '\n'

        message += 'polyline\n'
        message += '%s' % inlet.line
        message += '\n'

        message += '=====================================\n'

        return message


    def timestepping_statistics(self):

        message = '---------------------------\n'
        message += 'Inlet report for %s:\n' % self.label
        message += '--------------------------\n'
        message += 'Q [m^3/s]: %.2f\n' % self.applied_Q

        return message


    def set_Q(self, Q):

        self.Q = Q

    def get_Q(self):

        return self.Q


    def get_inlet(self):

        return self.inlet

    def get_line(self):

        return self.line

