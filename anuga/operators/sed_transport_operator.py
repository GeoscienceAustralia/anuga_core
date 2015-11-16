"""
Erosion operators


"""

__author__="steve"
__date__ ="$09/03/2012 4:46:39 PM$"



import numpy as num


from anuga import Domain
from anuga import Quantity
from anuga.operators.base_operator import Operator
from anuga import Region

from erosion_operators import Erosion_operator
from math import sqrt


#===============================================================================
# Specific Erosion operator trying to implement bed shear
#===============================================================================
class Sed_transport_operator(Erosion_operator):
    """
    Sed transport operator based on Bed shear erosion operator


    """

    def __init__(self, domain,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):
                 

        Erosion_operator.__init__(self,
                 domain,
                 indices=indices,
                 polygon=polygon,
                 center=center,
                 radius=radius,
                 description=description,
                 label=label,
                 logging=logging,
                 verbose=verbose)



    def update_quantities(self):
        """Update the vertex values of the quantities to model erosion
        """
        import numpy as num

        
        t = self.get_time()
        dt = self.get_timestep()
        
#         self.conc = 0.01
        self.conc = self.domain.quantities['concentration'].centroid_values


        Elev = self.domain.quantities['elevation']

        Elev.compute_local_gradients()

        self.elev_gx = Elev.x_gradient
        self.elev_gy = Elev.y_gradient


        Stage = self.domain.quantities['stage']

        Stage.compute_local_gradients()

        self.stage_gx = Stage.x_gradient
        self.stage_gy = Stage.y_gradient
        
        self.depth = self.stage_c - self.elev_c
        self.stage_slope = num.sqrt(self.stage_gx**2 + self.stage_gy**2)

        updated = True

        edot = self.erosion()
        ddot = self.deposition()
        porosity = 0.3

        dz = (ddot - edot) / (1 - porosity) * dt


        if self.domain.get_using_discontinuous_elevation():
            self.elev_c = self.elev_c - dz
            self.stage_c = self.elev_c + self.depth

        else:
            # v needs to be stacked to get the right shape (len(ids),3)
            dz = num.vstack((dz,dz,dz)).T
            

            # Ensure we don't erode below self.base level
            self.elev_v = self.elev_v - dz


        return updated

    def erosion(self):
    
        Ke_star = 2.
        D50 = 0.0001
        R = 1.65
        g = 9.81
        
        criticalshear_star = 0.3

        Ke = Ke_star * D50 * sqrt(R * g * D50)
        
        shear_stress_star = self.depth * self.stage_slope / (R * D50)
        
        
        edot = Ke * (shear_stress_star - criticalshear_star)
        edot[edot<0.0] = 0.0
        
        return edot

    def deposition(self):
    
        c1 = 18.
        c2 = 0.4
        mu = 1.0e-6
        rousecoeff = 1.
        R = 1.65
        g = 9.81
        D50 = 0.0001
        
    
        settlingvelocity = ((R * g * D50**2) /
                            ((c1 * mu) + (0.75 * c2 * (R * g * D50**3)**0.5)))

        ddot = rousecoeff * self.conc * settlingvelocity
        
        ddot[ddot<0.0] = 0.0

    
        return ddot

#------------------------------------------------
# Auxilary functions
#------------------------------------------------



def lineno():
    """Returns the current line number in our program."""
    import inspect
    return inspect.currentframe().f_back.f_back.f_lineno



def stage_elev_info(self):
    print 80*"="

    print 'In Evolve: line number ', lineno()
    import inspect
    print inspect.getfile(lineno)

    print 80*"="
    ind = num.array([ 976,  977,  978,  979,  980,  981,  982,  983, 1016, 1017, 1018,
             1019, 1020, 1021, 1022, 1023])
    elev_v = self.get_quantity('elevation').vertex_values
    stage_v = self.get_quantity('stage').vertex_values
    elev_c = self.get_quantity('elevation').centroid_values
    stage_c = self.get_quantity('stage').centroid_values

    from pprint import pprint
    print 'elev_v, elev_c, elev_avg \n'
    pprint( num.concatenate( (elev_v[ind], (elev_c[ind]).reshape(16,1),
                               num.mean(elev_v[ind],axis=1).reshape(16,1)), axis = 1))
    print 'stage_v, stage_c, stage_avg \n'
    pprint( num.concatenate( (stage_v[ind], (stage_c[ind]).reshape(16,1),
                               num.mean(stage_v[ind],axis=1).reshape(16,1)), axis = 1))

    print 80*"="
