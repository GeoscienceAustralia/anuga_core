"""
Rate operators (such as rain)

Constraints: See GPL license in the user guide
Version: 1.0 ($Revision: 7731 $)
"""

__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"



from anuga.config import indent
import numpy as num
import anuga.utilities.log as log
from anuga.utilities.function_utils import evaluate_temporal_function


from anuga import Quantity
from anuga.operators.base_operator import Operator
from anuga import Region

class Rate_operator(Operator,Region):
    """
    Add water at certain rate (ms^{-1} = vol/Area/sec) over a
    triangles specified by

    indices: None == all triangles, Empty list [] no triangles

    rate can be a function of time.

    Other units can be used by using the factor argument.

    """

    def __init__(self,
                 domain,
                 rate=0.0,
                 factor=1.0,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 default_rate=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):


        Operator.__init__(self, domain, description, label, logging, verbose)


        Region.__init__(self, domain,
                        indices=indices,
                        polygon=polygon,
                        center=center,
                        radius=radius,
                        verbose=verbose)


        #------------------------------------------
        # Local variables
        #------------------------------------------
        self.factor = factor

        self.rate_callable = False
        self.rate_spatial = False
        
        self.set_rate(rate)
        self.set_default_rate(default_rate)


        self.default_rate_invoked = False    # Flag

        self.set_areas()
        self.set_full_indices()
        
        # Mass tracking
        self.local_influx=0.

    def __call__(self):
        """
        Apply rate to those triangles defined in indices

        indices == [], then don't apply anywhere
        indices == None, then apply everywhere
        otherwise apply for the specific indices
        """

        if self.indices is []:
            return

        t = self.domain.get_time()
        timestep = self.domain.get_timestep()
        factor = self.factor
        indices = self.indices


        if self.rate_spatial:
            if indices is None:
                x = self.coord_c[:,0]
                y = self.coord_c[:,1]
            else:
                x = self.coord_c[indices,0]
                y = self.coord_c[indices,1]

            rate = self.get_spatial_rate(x,y,t)
        elif self.rate_type == 'quantity':
            if indices is None:
                rate  = self.rate.centroid_values
            else:
                rate = self.rate.centroid_values[indices]
        else:
            rate = self.get_non_spatial_rate(t)

        if self.verbose is True:
            log.critical('Rate of %s at time = %.2f = %f'
                         % (self.quantity_name, self.domain.get_time(), rate))

        fid = self.full_indices
        if num.all(rate >= 0.0):
            # Record the local flux for mass conservation tracking
            if indices is None:
                local_rates=factor*timestep*rate
                self.local_influx = (local_rates*self.areas)[fid].sum()
                self.stage_c[:] = self.stage_c[:] + local_rates
            else:
                local_rates = factor*timestep*rate
                self.local_influx=(local_rates*self.areas[fid]).sum()
                self.stage_c[indices] = self.stage_c[indices] \
                       + local_rates
        else: # Be more careful if rate < 0
            if indices is None:
                #self.local_influx=(num.minimum(factor*timestep*rate, self.stage_c[:]-self.elev_c[:])*self.areas)[fid].sum()
                #self.stage_c[:] = num.maximum(self.stage_c  \
                #       + factor*rate*timestep, self.elev_c )
                local_rates = num.maximum(factor*timestep*rate, self.elev_c[:]-self.stage_c[:])
                self.local_influx = (local_rates*self.areas)[fid].sum()
                self.stage_c[:] = self.stage_c + local_rates 
            else:
                #self.local_influx=(num.minimum(factor*timestep*rate, self.stage_c[indices]-self.elev_c[indices])*self.areas)[fid].sum()
                #self.stage_c[indices] = num.maximum(self.stage_c[indices] \
                #       + factor*rate*timestep, self.elev_c[indices])
                
                local_rates = num.maximum(factor*timestep*rate, self.elev_c[indices]-self.stage_c[indices])
                self.local_influx=(local_rates*self.areas)[fid].sum()
                self.stage_c[indices] = self.stage_c[indices] + local_rates 
        # Update mass inflows from fractional steps
        self.domain.fractional_step_volume_integral+=self.local_influx

        return

    def get_non_spatial_rate(self, t=None):
        """Provide a rate to calculate added volume
        """

        if t is None:
            t = self.get_time()


        assert not self.rate_spatial


        rate = evaluate_temporal_function(self.rate, t, default_right_value=self.default_rate)

        if rate is None:
            msg = ('Attribute rate must be specified in '+self.__name__+
                   ' before attempting to call it')
            raise Exception(msg)

        return rate

    def get_spatial_rate(self, x=None, y=None, t=None):
        """Provide a rate to calculate added volume
        only call if self.rate_spatial = True
        """

        assert self.rate_spatial

        if t is None:
            t = self.get_time()

        if x is None:
            assert y is None
            if self.indices is None:
                x = self.coord_c[:,0]
                y = self.coord_c[:,1]
            else:
                x = self.coord_c[self.indices,0]
                y = self.coord_c[self.indices,1]

        assert x is not None
        assert y is not None

        assert isinstance(t, (int, float))
        assert len(x) == len(y)

        #print xy
        #print t

        #print self.rate_type, self.rate_type == 'x,y,t'
        if self.rate_type == 'x,y,t':
            rate = self.rate(x,y,t)
        else:
            rate = self.rate(x,y)

        return rate


    def set_rate(self, rate):
        """Set rate 
        Can change rate while running
        Can be a scalar, or a function of t or x,y or x,y,t or a quantity
        """

        # Test if rate is a quantity
        if isinstance(rate, Quantity):
            self.rate_type = 'quantity'
        else:
            # Possible types are 'scalar', 't', 'x,y' and 'x,y,t'
            from anuga.utilities.function_utils import determine_function_type
            self.rate_type = determine_function_type(rate)


        self.rate = rate
        

        if self.rate_type == 'scalar':
            self.rate_callable = False
            self.rate_spatial = False
        elif self.rate_type == 'quantity':
            self.rate_callable = False
            self.rate_spatial = False
        elif self.rate_type == 't':
            self.rate_callable = True
            self.rate_spatial = False
        else:
            self.rate_callable = True
            self.rate_spatial = True






    def set_areas(self):

        if self.indices is None:
            self.areas = self.domain.areas
            return

        if self.indices is []:
            self.areas = []
            return

        self.areas = self.domain.areas[self.indices]

    def set_full_indices(self):

        if self.indices is None:
            self.full_indices = num.where(self.domain.tri_full_flag ==1)[0]
            return

        if self.indices is []:
            self.full_indices = []
            return

        self.full_indices = num.where(self.domain.tri_full_flag[self.indices] == 1)[0]

    def get_Q(self, full_only=True):
        """ Calculate current overall discharge
        """

        if full_only:
            if self.rate_spatial:
                rate = self.get_spatial_rate() # rate is an array
                fid = self.full_indices
                return num.sum(self.areas[fid]*rate[fid])*self.factor
            elif self.rate_type == 'quantity':
                rate = self.rate.centroid_values # rate is a quantity
                fid = self.full_indices
                return num.sum(self.areas[fid]*rate[fid])*self.factor
            else:
                rate = self.get_non_spatial_rate() # rate is a scalar
                fid = self.full_indices
                return num.sum(self.areas[fid]*rate)*self.factor
        else:
            if self.rate_spatial:
                rate = self.get_spatial_rate() # rate is an array
                return num.sum(self.areas*rate)*self.factor
            elif self.rate_type == 'quantity':
                rate = self.rate.centroid_values # rate is a quantity
                return num.sum(self.areas*rate)*self.factor                
            else:
                rate = self.get_non_spatial_rate() # rate is a scalar
                return num.sum(self.areas*rate)*self.factor

    def set_default_rate(self, default_rate):
        """
        Check and store default_rate
        """
        msg = ('Default_rate must be either None '
               'a scalar, or a function of time.\nI got %s.' % str(default_rate))
        assert (default_rate is None or
                isinstance(default_rate, (int, float)) or
                callable(default_rate)), msg


        #------------------------------------------
        # Allow longer than data
        #------------------------------------------
        if default_rate is not None:
            # If it is a constant, make it a function
            if not callable(default_rate):
                tmp = default_rate
                default_rate = lambda t: tmp

            # Check that default_rate is a function of one argument
            try:
                default_rate(0.0)
            except:
                raise Exception(msg)

        self.default_rate = default_rate


    def parallel_safe(self):
        """Operator is applied independently on each cell and
        so is parallel safe.
        """
        return True

    def statistics(self):

        message = 'You need to implement operator statistics for your operator'
        return message


    def timestepping_statistics(self):

        if self.rate_spatial:
            rate = self.get_spatial_rate()
            min_rate = num.min(rate)
            max_rate = num.max(rate)
            Q = self.get_Q()
            message  = indent + self.label + ': Min rate = %g m/s, Max rate = %g m/s, Total Q = %g m^3/s'% (min_rate,max_rate, Q)
        else:
            rate = self.get_non_spatial_rate()
            Q = self.get_Q()
            message  = indent + self.label + ': Rate = %g m/s, Total Q = %g m^3/s' % (rate, Q)


        return message






#===============================================================================
# Specific Rate Operators for circular region.
#===============================================================================
class Circular_rate_operator(Rate_operator):
    """
    Add water at certain rate (ms^{-1} = vol/Area/sec) over a
    circular region

    rate can be a function of time.

    Other units can be used by using the factor argument.

    """

    def __init__(self, domain,
                 rate=0.0,
                 factor=1.0,
                 center=None,
                 radius=None,
                 default_rate=None,
                 verbose=False):


        Rate_operator.__init__(self,
                               domain,
                               rate=rate,
                               factor=factor,
                               center=center,
                               radius=radius,
                               default_rate=default_rate,
                               verbose=verbose)



#===============================================================================
# Specific Rate Operators for polygonal region.
#===============================================================================
class Polygonal_rate_operator(Rate_operator):
    """
    Add water at certain rate (ms^{-1} = vol/Area/sec) over a
    polygonal region

    rate can be a function of time.

    Other units can be used by using the factor argument.

    """

    def __init__(self, domain,
                 rate=0.0,
                 factor=1.0,
                 polygon=None,
                 default_rate=None,
                 verbose=False):


        Rate_operator.__init__(self,
                               domain,
                               rate=rate,
                               factor=factor,
                               polygon=polygon,
                               default_rate=default_rate,
                               verbose=verbose)
                               


                




