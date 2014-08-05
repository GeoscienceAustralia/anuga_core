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
                 region,
                 Q = 0.0,
                 velocity = None,
                 default = None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):


        anuga.Operator.__init__(self, domain, description, label, logging, verbose)

        self.inlet = inlet.Inlet(self.domain, region, verbose= verbose)

        # should set this up to be a function of time and or space)
        self.Q = Q

        if velocity is not None:
            assert len(velocity)==2

        self.velocity = velocity

        self.applied_Q = 0.0

        self.set_default(default)

        self.activate_logging()


    def __call__(self):

        timestep = self.domain.get_timestep()

        t = self.domain.get_time()


        # Need to run global command on all processors
        current_volume = self.inlet.get_total_water_volume()
        total_area = self.inlet.get_area()

        assert current_volume >= 0.0

        Q1 = self.update_Q(t)
        Q2 = self.update_Q(t + timestep)


        #print Q1,Q2
        Q = 0.5*(Q1+Q2)
        volume = Q*timestep

        #print volume
        
        #print Q, volume

        # store last discharge
        self.applied_Q = Q
 

        
        # Distribute positive volume so as to obtain flat surface otherwise
        # just pull water off to have a uniform depth.
        if volume >= 0.0 :
            self.inlet.set_stages_evenly(volume)
            self.domain.fractional_step_volume_integral+=volume
            if self.velocity is not None:
                depths = self.inlet.get_depths()
                self.inlet.set_xmoms(self.inlet.get_xmoms()+depths*self.velocity[0])
                self.inlet.set_ymoms(self.inlet.get_ymoms()+depths*self.velocity[1])
        elif current_volume + volume >= 0.0 :
            depth = (current_volume + volume)/total_area
            self.inlet.set_depths(depth)
            self.domain.fractional_step_volume_integral+=volume
        else: #extracting too much water!
            self.inlet.set_depths(0.0)
            self.domain.fractional_step_volume_integral-=current_volume
            self.applied_Q = current_volume/timestep

            #msg =  'Requesting too much water to be removed from an inlet! \n'
            #msg += 'current_water_volume = %5.2e Increment volume = %5.2e' % (current_volume, volume)



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

        message += 'region\n'
        message += '%s' % inlet
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


