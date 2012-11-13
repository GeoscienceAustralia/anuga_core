import anuga
import numpy
import inlet

from warnings import warn

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
                 default = None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):


        anuga.Operator.__init__(self, domain, description, label, logging, verbose)

        
        self.line = numpy.array(line, dtype='d')

        # should set this up to be a function of time and or space)
        self.Q = Q


        #print self.Q

        self.enquiry_point = 0.5*(self.line[0] + self.line[1])
        self.outward_vector = self.line
        self.inlet = inlet.Inlet(self.domain, self.line, verbose= verbose)

        self.applied_Q = 0.0

        self.set_default(default)


    def __call__(self):

        timestep = self.domain.get_timestep()

        t = self.domain.get_time()



        # Need to run global command on all processors
        current_volume = self.inlet.get_total_water_volume()

        Q1 = self.update_Q(t)
        Q2 = self.update_Q(t + timestep)


        #print Q1,Q2
        Q = 0.5*(Q1+Q2)
        volume = Q*timestep

        #print volume
        
        #print Q, volume

        # store last discharge
        self.applied_Q = Q

        msg =  'Requesting too much water to be removed from an inlet! \n'
        msg += 'current_water_volume = %5.2e Increment volume = %5.2e' % (current_volume, volume)
        import warnings
        if current_volume + volume < 0.0:
            #warnings.warn(msg)
            volume = -current_volume
            self.applied_Q = volume/timestep


        #print 'applied_Q', self.applied_Q
        
        # Distribute volume so as to obtain flat surface
        self.inlet.set_stages_evenly(volume)
        
        # Distribute volume evenly over all cells
        #self.inlet.set_depths_evenly(volume)



    def update_Q(self, t):
        """Allowing local modifications of Q
        """
        from anuga.fit_interpolate.interpolate import Modeltime_too_early, Modeltime_too_late
        
        if callable(self.Q):
            try:
                Q = self.Q(t)
            except Modeltime_too_early, e:
                raise Modeltime_too_early(e)
            except Modeltime_too_late, e:
                Q = self.get_default(t,err_msg=e)
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


    def set_default(self, default=None):
        
        """ Either leave default as None or change it into a function"""

        if default is not None:
            # If it is a constant, make it a function
            if not callable(default):
                tmp = default
                default = lambda t: tmp

            # Check that default_rate is a function of one argument
            try:
                default(0.0)
            except:
                raise Exception(msg)

        self.default = default
        self.default_invoked = False



    def get_default(self,t, err_msg=' '):
        """ Call get_default only if exception
        Modeltime_too_late(msg) has been raised
        """

        from anuga.fit_interpolate.interpolate import Modeltime_too_early, Modeltime_too_late


        if self.default is None:
            msg = '%s: ANUGA is trying to run longer than specified data.\n' %str(err_msg)
            msg += 'You can specify keyword argument default in the '
            msg += 'operator to tell it what to do in the absence of time data.'
            raise Modeltime_too_late(msg)
        else:
            # Pass control to default rate function
            value = self.default(t)

            if self.default_invoked is False:
                # Issue warning the first time
                msg = ('%s\n'
                       'Instead I will use the default rate: %s\n'
                       'Note: Further warnings will be supressed'
                       % (str(err_msg), str(self.default(t))))
                warn(msg)

                # FIXME (Ole): Replace this crude flag with
                # Python's ability to print warnings only once.
                # See http://docs.python.org/lib/warning-filter.html
                self.default_invoked = True

        return value


    def set_Q(self, Q):

        self.Q = Q

    def get_Q(self):

        return self.Q


    def get_inlet(self):

        return self.inlet

    def get_line(self):

        return self.line

