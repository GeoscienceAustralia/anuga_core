""" -----------------------------------------------sanddune_erosion_operator V2

This script provides an Anuga operator that simulates the removal of sand
associated with over-topping of a sand dune by a tsunami. It is provided to 
permit the exploration of a scenario where a dune line provides protection for 
land(or water) behind the dune line that could be lost if part of the dune was  
to be over topped by the first or successive waves of the tsunami.

The erosion mechanism assumes a clear water scour (no significant sediment
entrained in the approaching wave that could impact detachment rates and an 
environment where the detached sediment is well within the transport capacity
of the water column in the eroding area(viz no settlement only erosion 
within the target area). 

The erosional relationships used are based on the work by Dr David Froelich. 
Froelich (2002) IMPACT Project Field Tests 1 and 2 "Blind" Simulation by DaveF
Mo-i-Rana Norway, 2002

Erosional parameters of relevance extracted from the above are;

Bed shear  Tau_bed = Wd*G*(n**2)* (m**2) / (d**2.333) Pa (m is absolute momentum)
Where 
Wd = water mass density,               (1000 kg/m3)
Sd = sediment mass density             ( 1800 kg/m3)
G = Accel due to gravity,              (9.8 m/sec/sec)
n = mannings n,                        (sand n = 0.025)
m = abs momentum                       ((mx**2+my**2)**0.5)
d = water depth                        (Anuga stage-elevation m)

Dune erosion will occur when bed shear stress > critical bed shear stress
from Froelich Table 1 
Critical (detachment) bed shear stress for sand Tau_crit is 2.1 Pa
Detachment rate constant Kd is 0.0250 Kg/sec/m2/Pa and
detachment rate Sk = Kd*(tau_bed-Tau_Crit)/Sd m3/sec/m2 (Froelich Eq4)

Note That while Froelich provides erosional properties for other soils, this
operator is specific to the erosion of a sand dune.

The second process simulated by this operator is the collapse, fluidisation and 
removal of sand from the dune system as a consequence of the above (vertically) 
aligned erosion process creating face slopes that would be greater than the
angle of repose.  This is applied as a separate check after the erosion computations 
for a particular timestep have been completed and a new eroded surface computed. 
In this second process, each triangle within the specified erosion zone is checked 
to see if the CG (elevation) of any neighbouring triangles lie below the angle of 
repose line from the current erosion triangles CG . If so then the current triangles 
CG (elevation) that is above the repose angle is lowered to lie at the angle of repose 
relative to the lowest neighbours CG. A check is performed to prevent any reduction 
in elevation that would lower the current triangles CG below the specified base level. 

This is an approximate simplified lower bound approach to the full solution of surface 
level adjustments based on true dip and strike.  While this algorithm allows 
each triangles height in the erosion zone to be reviewed and adjusted for instability at
each timestep, it can not guarantee  all triangles will be fully adjusted back to a 
stable face. Given the small size of each time step it is considered that it can 
however reasonably simulate the collapse of steep sand faces.

It is noted that a third removal process that involves the regressive fluidisation and 
collapse of the downstream face from a rising phreatic surface within the dune from the 
rising sea level. This has not been included in this model as it is considered a 
significantly longer duration process than the other two processes. 

NOTE - NOTE - NOTE - NOTE - NOTE - NOTE - NOTE - NOTE
Anuga is a 2D model and as such can not simulate vertical components of velocity where 
they exist. Flow approaching and in particular down the downstream face of a dune will 
have a significant vertical velocity component that will alter the rate of detachment 
of sand. This code is therefore only an approximation of the realworld process and its 
results should not be used where full simulation of the real world process is the intent.



Script Author E Rigby, ted.rigby@rienco.com.au 
Version       1.00 October 2015
              2.00 June 2016 code modified  to lower triangle CG if lowest neighbour was 
			                 below angle of repose rather than alter neighbour levels if
							 highest was above angle of repose.
"""

from anuga.operators.base_operator import Operator
from anuga import Region
import pdb



class sanddune_erosion_operator(Operator, Region)  :
    
    def __init__(self,
                 domain,
                 threshold= 0.0,
                 base=0.0,
                 indices=None,
                 polygon=None,
                 center=None,
                 radius=None,
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
        self.base = base

       
    def __call__(self):
        """
        Applies  removal of sand by erosion and  collapse due to slope instability to
        those triangles defined in indices by the specified polygon

        indices == [], then don't apply anywhere
        otherwise apply for the specific indices within the erosion area polygon
        """
        
        import math
        import numpy as num

        updated = True             # flag called OK
        
        if self.indices is not []:     # empty list no polygon - return
        
            #-------------------------------------------
            # set some sand and water properties
            #-------------------------------------------
            Wd       = 1000       # water mass density kg/m3
            Sd       = 1800       # sediment mass density kg/m3
            G        = 9.8        # acceleration due to gravity m/sec/sec
            n        = 0.030      # sand mannings n - mostly bare with undulations
            Tau_crit = 2.1        # critical (detachment) bed shear stress Pa
            Kd       = 0.025      # detachment factor Froelich Table 2 Kg/sec/m2/Pa
            Ra       = 34.0       # Repose angle in degrees for dry sand (note wet varies 30-45)
            Rs       = math.tan(Ra*math.pi/180)    # Repose slope mV/mH
    
            #-------------------------------------------
            # get some useful model parameters
            #-------------------------------------------        
            
            dt = self.get_timestep()
    
            #-----------------------------------------------------------------------------------------
            # Compute erosion depths during the timestep and update centroid elevations accordingly 
            # Don't allow erosion below the specified base level.
    		# Note this operator is called for each seperate erosion polygon.
            #-----------------------------------------------------------------------------------------
            
                                                                              
            ind = self.indices                                                   # indices of triangles in polygon
                        
            # store previous timestep centroid heights(depths) for later use
            h = self.stage_c[ind] - self.elev_c[ind]    
            
            m = num.sqrt(self.xmom_c[ind]**2 + self.ymom_c[ind]**2)              # abs Momentum m2/sec vector
            d =(self.stage_c[ind]-self.elev_c[ind])                              # Depth meters vector
                    
            # compute bed shear stress Pa
            Tau_bed  = Wd*G*(n**2)*(m**2)/((d**2.333)+ 0.000001)                 # refer Froelich 2002 vector
            
            # compute de (m) elevation change increment due to scour during timestep
            de = (Kd*(Tau_bed-Tau_crit)/Sd) *dt                                  # refer Froelich 2002 vector
            
            # no scour though if Tau_bed < Tau_crit ie if de is <= 0
            de = num.where(de > 0.0, de, 0.0)                                    # de=de whenever de>0 vector
                            
            # Also ensure we don't erode below the base surface level  
            self.elev_c[ind] = num.maximum(self.elev_c[ind] - de, num.minimum(self.elev_c[ind],self.base[ind]) )   
            
            #-------------------------------------------------------------------------------------------
            # Reduce triangles elevations in any area where erosion has created unstable face slopes, so  
            # that the final face slope lies at or about the angle of repose. This only approximates the 
    		# collapse process as slopes are CG to CG and not specifically oriented to the dip angle.
            #--------------------------------------------------------------------------------------------
            neighbours = self.domain.surrogate_neighbours
            
            neighbourindices = neighbours[ind]           # get the neighbour Indices for each triangle in the erosion zone(s)
                 
            n0 = neighbourindices[:,0]                   # separate into three lists
            n1 = neighbourindices[:,1]
            n2 = neighbourindices[:,2] 
                
            k = n0.shape[0]                              # Erosion poly lEN - num of triangles in poly  
            
            e = num.zeros((k,3))                         # create elev array k triangles, 3 neighbour elev
                
            e[:,0] = self.elev_c[n0]                     # get the elev of each neighbours CG
            e[:,1] = self.elev_c[n1]
            e[:,2] = self.elev_c[n2]
                
            minid0 = num.argmin(e, axis=1)               # get indices of min elev for each triangles neighbours
            
            ident  = num.arange(k)                       # ident is array 0.. k-1, step 1
                
            n0ind = neighbourindices[ident,minid0]       # store the neighbour indices in order lowest elev first (n0)
            
        
            # compute the plan distance lxy from each triangles CG to the lowest neighbours CG 
            # Note when partitioned the neighbour indices can return the traget triangles index if on a ghost boundary so
            # need to deal with lxy = 0 situations         
            lxy = num.sqrt((self.coord_c[ind][:,0]-self.coord_c[n0ind][:,0])**2 + (self.coord_c[ind][:,1]-self.coord_c[n0ind][:,1])**2 )
            
            # compute the slope between each triangle and its lowest neighbour (CG to CG) from neighbour to current triangle 
            # Note when partitioned the neighbour indices can return the current triangles index if on a ghost boundary so
            # need to deal with lxy = 0 situations
            with num.errstate(divide='ignore', invalid='ignore'):            
                s = num.where(lxy>0.0, (self.elev_c[ind]-self.elev_c[n0ind])/lxy, 0.0)   # positive if current is above lowest neighbour
            # Lower the current triangle to be at the angle of repose from the lowest neighbour if s > repose but not below base
            self.elev_c[ind] = num.where(s > Rs, num.maximum(self.elev_c[n0ind]+(Rs*lxy), self.base[ind]), self.elev_c[ind])       
            
    		
                    
          
            # once all triangles processed for erosion and collapse add back height to maintain volumes correctly
            self.stage_c[ind] = self.elev_c[ind] + h
            
        self.domain.update_ghosts(['elevation', 'stage'])
              
        return (updated)
                
    
    def parallel_safe(self):
        """If Operator is applied independently on each cell and
        so is parallel safe.
        """
        
        if self.domain.get_using_discontinuous_elevation():
            return True
        else:
            return False
        
                        
