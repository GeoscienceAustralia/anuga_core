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
from anuga import Region
import numpy as num


#------------------------------------------------------------------------------
# Function to describe topography as function of X and y
#------------------------------------------------------------------------------
def topography(x,y):
    """Complex topography defined by a function of vectors x and y."""
    print ' Creating topography....'
    
    z = 0.0*(x)                             # horizontal plane 
    	
    N = len(x)
   
    for i in range(N):
        # First notched sand dune across Channel
        if 10.6 < x[i] <= 12.0:
            z[i] +=  1.0*(x[i] -10.6)      # Sloping U/S Face 1:1
        if 12.0 < x[i] < 13.0 :
            z[i] +=  1.4                   # Crest of Embankment at +1.4
        if 13.0 <= x[i] < 14.4:
            z[i] +=  1.4-1.0*(x[i] -13.0)  # Sloping D/S Face at 1:1
			# add notch in crest 1m wide by nom 300 deep
			# note sides are near vertical so will collapse back to repose even without erosion
        if 11.7 <= x[i] <= 13.3 and 2.0 <= y[i] <= 3.0:
            z[i] =  1.1                   # add 300 Notch in Embankment crest  
        # second lower plain sand dune across Channel
        if 23.0 < x[i] <= 24.0:
            z[i] +=  1.0*(x[i] -23.0)      # Sloping U/S Face 1:1        
        if 24.0 < x[i] < 25.0 :
            z[i] +=  1.0                   # Crest of Embankment at +1.0
        if 25.0 <= x[i] < 26.0:
            z[i] +=  1.0-1.0*(x[i] -25.0)  # Sloping D/S Face at 1:1			
    return z

#------------------------------------------------------------------------------
# Setup computational domain
#------------------------------------------------------------------------------
print ' Setting  up Domain ...'
length = 36.
width = 5.
dx = dy = 0.1            # Resolution: Length of subdivisions on both axes
# set points over crest to check longitudinal progress of erosion as simulation evolves
check_dune1_x = [[11.0,2.5],[12.1,2.5],[12.5,2.5],[12.9,2.5],[14.0,2.5]] # US to DS of crest notch
# set points along crest to check lateral progress of erosion as simulation evolves
check_dune1_y = [[12.5,0.5],[12.5,1.5],[12.5,2.5],[12.5,3.5],[12.5,4.5]] # Left to right of middle viewed US

if anuga.myid == 0:
    points, vertices, boundary = anuga.rectangular_cross(int(length/dx), int(width/dy),
                                                   len1=length, len2=width)
    domain = anuga.Domain(points, vertices, boundary)
    domain.set_flow_algorithm('DE0')
    domain.set_name() # Output name
    domain.set_store_vertices_uniquely(True)

    print domain.statistics()

    domain.set_quantities_to_be_stored({'elevation': 2,
                                        'stage': 2,
                                        'xmomentum': 2,
                                        'ymomentum': 2})



    domain.set_quantity('elevation', topography)           # elevation is a function
    domain.set_quantity('friction', 0.01)                  # Constant friction
    domain.set_quantity('stage', expression='elevation')   # Dry initial condition
else:
    domain = None


domain = anuga.distribute(domain)

#------------------------------------------------------------------------------
# Setup boundary conditions
#------------------------------------------------------------------------------
Bi = anuga.Dirichlet_boundary([1.2, 0, 0])          # Inflow at depth
Br = anuga.Reflective_boundary(domain)              # Solid reflective side walls
Bo = anuga.Dirichlet_boundary([-5, 0, 0])           # uncontrolled outflow 

domain.set_boundary({'left': Bi, 'right': Bo, 'top': Br, 'bottom': Br})

#------------------------------------------------------------------------------
# Setup sanddune erosion operator 
#------------------------------------------------------------------------------
print 'Set up Erosion Areas to test...'


polygon1 = [ [10.6, 0.0], [14.4, 0.0], [14.4, 5.0], [10.6, 5.0] ]
polygon2 = [ [23.0, 0.0], [26.0, 0.0], [26.0, 5.0], [23.0, 5.0] ]

poly1nsbase = 0.5  # set poly no scour base level in m
poly2nsbase = 0.3

# get the indices of triangles in each erosion poly
poly1ind = (Region(domain, polygon=polygon1)).indices
poly2ind = (Region(domain, polygon=polygon2)).indices
print 'poly1ind is of length and contains ', len(poly1ind), poly1ind
print 'poly2ind is of length and contains ', len(poly2ind), poly2ind

# get the initial model surface elevation
nsbase_elev_c = domain.get_quantity('elevation').get_values(location='centroids')
print 'nsbase_elev_c is of length ',len(nsbase_elev_c) 


# build the no scour base surface  by combining initial elev where < nsbase and nsbase in each scour poly
nsbase_elev_c[poly1ind] = num.minimum(nsbase_elev_c[poly1ind], poly1nsbase)
nsbase_elev_c[poly2ind] = num.minimum(nsbase_elev_c[poly2ind], poly2nsbase)



# setup and create operator within polys setting the scour base elevations for each poly
op1 = anuga.Sanddune_erosion_operator(domain, base=nsbase_elev_c, polygon=polygon1)   # first notched dune
op2 = anuga.Sanddune_erosion_operator(domain, base=nsbase_elev_c, polygon=polygon2)   # second plain dune




#check_dune1_x_ids = [ domain.get_triangle_containing_point(point) for point in check_dune1_x ]
#check_dune1_y_ids = [ domain.get_triangle_containing_point(point) for point in check_dune1_y ]

#print check_dune1_x_ids
#print check_dune1_y_ids

#------------------------------------------------------------------------------
# Evolve sanddune erosion simulation through time  
#------------------------------------------------------------------------------
for t in domain.evolve(yieldstep=1, duration=600.0):

    if anuga.myid==0:
        domain.print_timestepping_statistics()
    
    #dune1_profile_x = domain.get_quantity('elevation').get_values(interpolation_points=check_dune1_x)
    #dune1_profile_x = domain.get_quantity('elevation').centroid_values[check_dune1_x_ids]
    #print '    Sand profile over first dune crest at time ', t, ' is ' 
    #print '    ', dune1_profile_x
    #dune1_profile_y = domain.get_quantity('elevation').get_values(interpolation_points=check_dune1_y)
    #dune1_profile_y = domain.get_quantity('elevation').centroid_values[check_dune1_y_ids]
    #print '    Sand profile along first dune crest at time ', t, ' is '
    #print '    ', dune1_profile_y
