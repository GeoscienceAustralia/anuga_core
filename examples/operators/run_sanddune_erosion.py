"""sand dune erosion test                          run_sanddune_erosion.py

This tests the sanddune_erosion operator confirming that; 
    1. flow creates erosion when bed shear > critical and that erosion rates 
	are higher in higher bed shear zones (typically higher velocity areas).
	2. that erosion is augmented by collapse of the sand face whenever
	erosion creates slopes > the angle of repose.  
	this process leads to widening of the notch (laterally)and head
    like recession of 
	the main scour zone  as scour acts to steepen the longitudinal 
	bed, triggering collapse of the steep bed and reshaping back 
	to the angle of repose).
	3. that the operator can handle multiple erosion polygons with
    different base levels
	 
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
import anuga
from anuga import  myid, distribute, barrier, finalize
from anuga import Region, Geo_reference
import numpy as num

x0 = 314036.58727982
y0 = 6224951.2960092
geo = Geo_reference(56, x0, y0)

#------------------------------------------------------------------------------
# Function to describe topography as function of X and y
#------------------------------------------------------------------------------
def topography(x,y):
    """Complex topography defined by a function of vectors x and y."""
    print (' Creating topography....')
    
    z = 0.0*(x)                             # horizontal plane 

    N = len(x)
   
    for i in range(N):
        # First notched sand dune across Channel
        if 1010.6 < x[i] <= 1012.0:
            z[i] +=  1.0*(x[i] -1010.6)      # Sloping U/S Face 1:1
        if 1012.0 < x[i] < 1013.0 :
            z[i] +=  1.4                   # Crest of Embankment at +1.4
        if 1013.0 <= x[i] < 1014.4:
            z[i] +=  1.4-1.0*(x[i] -1013.0)  # Sloping D/S Face at 1:1

        # add notch in crest 1m wide by nom 300 deep
        # note sides are near vertical so will collapse back to repose even without erosion

        if 1011.7 <= x[i] <= 1013.3 and 1002.0 <= y[i] <= 1003.0:
            z[i] =  1.1                   # add 300 Notch in Embankment crest  
        # second lower plain sand dune across Channel
        if 1023.0 < x[i] <= 1024.0:
            z[i] +=  1.0*(x[i] -1023.0)      # Sloping U/S Face 1:1        
        if 1024.0 < x[i] < 1025.0 :
            z[i] +=  1.0                   # Crest of Embankment at +1.0
        if 1025.0 <= x[i] < 1026.0:
            z[i] +=  1.0-1.0*(x[i] -1025.0)  # Sloping D/S Face at 1:1      
    return z

#------------------------------------------------------------------------------
# build the check points and erosion polygons now so available to all processors
#------------------------------------------------------------------------------

polygon1    = num.array([ [1010.6, 1000.0], [1014.4, 1000.0], [1014.4, 1005.0], [1010.6, 1005.0] ])
polygon2    = num.array([ [1023.0, 1000.0], [1026.0, 1000.0], [1026.0, 1005.0], [1023.0, 1005.0] ])
poly1nsbase = 0.5  # set poly no scour base level in m
poly2nsbase = 0.3

polygon1 += num.array([x0,y0])
polygon2 += num.array([x0,y0])

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------

if myid == 0:
    print ('>>>>> DUNE EROSION TEST SCRIPT V2')
    print ('>>>>> Setting  up Domain on processor 0...')
    length = 36.
    width = 5.
    dx = dy = 0.1            # Resolution: Length of subdivisions on both axes
    print ('>>>>> Domain has L = %f, W = %f, dx=dy= %f' %(length,width,dx))
    points, vertices, boundary = anuga.rectangular_cross(int(length/dx), int(width/dy), len1=length, len2=width)
    points = points + 1000.0
    
 
    domain = anuga.Domain(points, vertices, boundary, geo_reference = geo)
    domain.set_flow_algorithm('DE0')
    domain.set_name('sanddune_testV2SR') # Output name
    
    domain.set_store_vertices_uniquely(False)
    domain.set_quantities_to_be_stored({'elevation': 2,'stage': 2,'xmomentum': 2,'ymomentum': 2})

    domain.set_quantity('elevation', topography)           # elevation is a function
    domain.set_quantity('friction', 0.01)                  # Constant friction
    domain.set_quantity('stage', expression='elevation')   # Dry initial condition

    print (domain.statistics())
    
    # get the indices of triangles in each erosion poly so can setup nsbase_ in domain
    poly1ind = (Region(domain, polygon=polygon1)).indices
    poly2ind = (Region(domain, polygon=polygon2)).indices
    print ('>>>>> poly1ind is of length ', len(poly1ind), ' and contains triangles', poly1ind[0], ' to ' , poly1ind[-1])
    print ('>>>>> poly2ind is of length ', len(poly2ind), ' and contains triangles', poly2ind[0], ' to ' , poly2ind[-1])
    
    # get the initial model surface elevation
    nsbase_elev_c = domain.get_quantity('elevation').get_values(location='centroids')
    print ('>>>>> nsbase_elev_c is of length ',len(nsbase_elev_c)) 
        
        
    # build the no scour base surface  by combining initial elev where < nsbase and nsbase in each scour poly
    nsbase_elev_c[poly1ind] = num.minimum(nsbase_elev_c[poly1ind], poly1nsbase)
    nsbase_elev_c[poly2ind] = num.minimum(nsbase_elev_c[poly2ind], poly2nsbase)
    
    
    # create new Anuga quantity and assign nsbase_elev_c values to domain so can distribute later
    anuga.Quantity(domain, name='nsbase_elevation', register=True)
    domain.set_quantity('nsbase_elevation',nsbase_elev_c, location='centroids')
    
else:
    
    domain = None

#------------------------------------------------------------------------------
# Print out polygon points and ids as check
#------------------------------------------------------------------------------

if myid == 0:
    print ('>>>>> erosion polygon1 contains ', polygon1)
    print ('>>>>> erosion polygon2 contains ', polygon2)

#------------------------------------------------------------------------------
# Distribute the domain onto the n partitions
#------------------------------------------------------------------------------
domain = distribute(domain)


#------------------------------------------------------------------------------
# Setup boundary conditions on the distributed domain
#------------------------------------------------------------------------------
Bi = anuga.Dirichlet_boundary([1.2, 0, 0])          # Inflow at depth
Br = anuga.Reflective_boundary(domain)              # Solid reflective side walls
Bo = anuga.Dirichlet_boundary([-5, 0, 0])           # uncontrolled outflow 

domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Setup sanddune erosion operator 
#------------------------------------------------------------------------------

if myid == 0:
    print ('>>>>> Setting up Erosion Area(s) to test...')

# power up the erosion operator
from anuga import Sanddune_erosion_operator

# assign ns base elevations to partitioned domains
nsbase_elev_c = domain.get_quantity('nsbase_elevation').get_values(location='centroids')

poly1ind = (Region(domain, polygon=polygon1)).indices
poly2ind = (Region(domain, polygon=polygon2)).indices
indices_union = list(set(poly1ind) | set(poly2ind))

# setup and create operator within polys setting the scour base elevations for each poly
op0 = Sanddune_erosion_operator(domain, base=nsbase_elev_c, indices=indices_union, Ra=45)   # both dunes
#op1 = sanddune_erosion_operator(domain, base=nsbase_elev_c, polygon=polygon1)   # first notched dune
#op2 = sanddune_erosion_operator(domain, base=nsbase_elev_c, polygon=polygon2)   # second plain dune


#------------------------------------------------------------------------------
# Evolve simulation through time  
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=1, duration=6.0):
    if myid == 0: domain.print_timestepping_statistics()
    

#  run completed - tidy up 
barrier()  
if myid == 0: 
    print (' >>>>> Simulation completed successfully ')
    print (' >>>>> Merging the individual cpu sww files and deleting the individual swws once merged')

# Merge the individual parallel swws created by the n processors
barrier()                         # wait foir all processors to complete
domain.sww_merge(delete_old=True)

if myid == 0: 
    print (' >>>>> Finalising the run -- all done')

# Finaise the parallel code and this model run
barrier()                         # wait for all processors to complete
finalize()



   
