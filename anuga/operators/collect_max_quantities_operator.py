"""
Collect stage quantity stage info


"""

__author__="gareth"
__date__ ="$01/08/2014 4:46:39 PM$"



import numpy as num
from anuga.geometry.polygon import inside_polygon
from anuga.operators.base_operator import Operator
from anuga import Quantity
from anuga.parallel import myid
from anuga.config import velocity_protection
from anuga.config import max_float


class Collect_max_quantities_operator(Operator):
    """
    Simple operator to collect the max stage, depth, speed, and speed*depth during a run.

    Maxima are updated every update_frequency timesteps [any integer >=1 is
    ok], after t exceeds collection_start_time.
   
    In theory this might save time (??), since computing e.g. velocity/momentum etc in python might be expensive
    
    Optionally velocities can be zeroed below velocity_zero_height (defaults to minimum_allowed_height if not required) 
    """

    def __init__(self,
                 domain,
                 update_frequency=1,
                 collection_start_time=0.,
                 velocity_zero_height=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):


        Operator.__init__(self, domain, description, label, logging, verbose)

        self.domain=domain

        #------------------------------------------
        # Setup a quantity to store max_stage
        #------------------------------------------
        self.max_stage = num.zeros(len(domain.centroid_coordinates[:,0])) - max_float
        self.max_depth = num.zeros(len(domain.centroid_coordinates[:,0]))
        self.max_speed = num.zeros(len(domain.centroid_coordinates[:,0]))
        self.max_speedDepth = num.zeros(len(domain.centroid_coordinates[:,0]))

        #------------------------------------------
        # Aliases for stage quantity
        #------------------------------------------
        self.xy=domain.centroid_coordinates
        self.stage  = domain.quantities['stage']
        self.elev   = domain.quantities['elevation']
        self.xmom   = domain.quantities['xmomentum']
        self.ymom   = domain.quantities['ymomentum']

        #------------------------------------------
        # Counter so we don't have to update every timestep
        #------------------------------------------        
        self.counter=0
        assert update_frequency>0, 'Update frequency must be >=1'
        self.update_frequency=update_frequency
        self.collection_start_time = collection_start_time

        #------------------------------------------
        # Can (rarely) get high velocities being recorded
        # (e.g. in 3 cells out of 170000 for a 48 hour simulation,
        #  which are not persistent -- probably just spikes from nearly dry cells)
        # Try to remove this by zeroing velocity in very shallow cells
        #------------------------------------------
        if velocity_zero_height is not None:
            self.velocity_zero_height=velocity_zero_height
        else:
            self.velocity_zero_height=domain.minimum_allowed_height


    def __call__(self):
        """
        Calculate max_quantities at every 'update_frequency' timesteps once time > collection_start_time
        """

        t = self.domain.get_time()

        if(t > self.collection_start_time):
            self.counter+=1

            if(self.counter==self.update_frequency):
            
                self.max_stage=num.maximum(self.max_stage, self.stage.centroid_values)

                momNorm = (self.xmom.centroid_values**2 + self.ymom.centroid_values**2)**0.5
                self.max_speedDepth=num.maximum(self.max_speedDepth, momNorm)
                
                localDepth=num.maximum(self.stage.centroid_values-self.elev.centroid_values, self.domain.minimum_allowed_height) 
                self.max_depth = num.maximum(self.max_depth,localDepth)

                #velMax=(momNorm/(localDepth+velocity_protection/localDepth))*(localDepth>self.domain.minimum_allowed_height)
                velMax=(momNorm/localDepth)*(localDepth>self.velocity_zero_height)
                self.max_speed = num.maximum(self.max_speed, velMax)

                self.counter=0

    def parallel_safe(self):
        """Operator is applied independently on each cell and
        so is parallel safe.
        """
        return True

    def statistics(self):

        message = self.label + ': Collect_max_quantity operator'
        return message


    def timestepping_statistics(self):
        from anuga import indent

        message  = indent + self.label + ': Collecting_max_quantity'
        return message


    def export_max_quantities_to_csv(self, filename_start='Max_Quantities_'):
        """

        Export max-quantities to a csv 

        """
        fullInds=self.domain.tri_full_flag.nonzero()[0]
        outArray=num.vstack([self.xy[fullInds,0]+self.domain.geo_reference.xllcorner, 
                             self.xy[fullInds,1]+self.domain.geo_reference.yllcorner, 
                             self.max_stage[fullInds], self.max_depth[fullInds], 
                             self.max_speed[fullInds], self.max_speedDepth[fullInds]]).transpose()

        outname=filename_start+'P'+str(myid)+'_X_Y_Stage_Depth_Speed_UH_MAX.csv'
        num.savetxt(outname, outArray,delimiter=',')
        return

