"""Example of shallow water wave equation.

This is called Netherlands because it shows a dam with a gap in it and
stylised housed behind it and below the water surface.

"""

######################
# Module imports 
#
from shallow_water import Domain, Reflective_boundary, Dirichlet_boundary,\
     Transmissive_boundary, Constant_height

from mesh_factory import rectangular
from Numeric import array
    

class Weir:
    """Set a bathymetry for simple weir with a hole.
    x,y are assumed to be in the unit square
    """
    
    def __init__(self, stage):
        self.inflow_stage = stage

    def __call__(self, x, y):    
        from Numeric import zeros, Float 
    
        N = len(x)
        assert N == len(y)    

        z = zeros(N, Float)
        for i in range(N):
            z[i] = -x[i]/20  #General slope

            #Flattish bit to the left
            if x[i] <= 0.3:
                #z[i] = -x[i]/5
                z[i] = -x[i]/20
                
            
            #Weir
            if x[i] > 0.3 and x[i] < 0.4:
                z[i] = -x[i]/20+1.2

            #Dip
            #if x[i] > 0.6 and x[i] < 0.9:
            #    z[i] = -x[i]/20-0.5  #-y[i]/5

            #Hole in weir
            #if x[i] > 0.3 and x[i] < 0.4 and y[i] > 0.2 and y[i] < 0.4:
            if x[i] > 0.3 and x[i] < 0.4 and y[i] > 0.4 and y[i] < 0.6:
                #z[i] = -x[i]/5
                z[i] = -x[i]/20
            
            #Poles
            #if x[i] > 0.65 and x[i] < 0.8 and y[i] > 0.55 and y[i] < 0.65 or\
            #       x[i] > 0.75 and x[i] < 0.9 and y[i] > 0.35 and y[i] < 0.45:
            #    z[i] = -x[i]/20+0.4

            if (x[i] - 0.72)**2 + (y[i] - 0.6)**2 < 0.05**2:# or\
                   #x[i] > 0.75 and x[i] < 0.9 and y[i] > 0.35 and y[i] < 0.45:
                z[i] = -x[i]/20+0.4

            #Wall
            if x[i] > 0.995:
                z[i] = -x[i]/20+0.3                

        return z/2
        
    
    
######################
# Domain
#

N = 250
#N= 8
N = 16
#N = 4
#N = 102
N = 25
N = 16
N = 60
N = 150 #size = 45000
N = 130 #size = 33800
#N = 60
#N = 40
N = 260
#N = 150
N = 264

N = 600 #Size = 720000
N = 20
#N = 150
N = 110
N = 60

N = 40
#N = 140
#N = 15

print 'Creating domain'
#Create basic mesh
points, vertices, boundary = rectangular(N, N)

#Create shallow water domain
domain = Domain(points, vertices, boundary)
    
domain.check_integrity()
domain.default_order = 2
#domain.beta_h=0

#Output params
domain.smooth = True
domain.reduction = min  #Looks a lot better on top of steep slopes

print "Number of triangles = ", len(domain)


if N > 40:
    domain.visualise = False
    domain.checkpoint = False
    domain.store = True    #Store for visualisation purposes
    domain.format = 'sww'   #Native netcdf visualisation format
    import sys, os
    #FIXME: This was os.path.splitext but caused weird filenames based on root
    base = os.path.basename(sys.argv[0])
    domain.filename, _ = os.path.splitext(base)
else:
    domain.visualise = False
    domain.checkpoint = False
    domain.store = False    

    
#Set bed-slope and friction
inflow_stage = 0.08
manning = 0.02
Z = Weir(inflow_stage)

print 'Field values'
domain.set_quantity('elevation', Z)
domain.set_quantity('friction', manning)


######################
# Boundary conditions
#
print 'Boundaries'
Br = Reflective_boundary(domain)
Bt = Transmissive_boundary(domain)

#Constant inflow
Bd = Dirichlet_boundary([2*inflow_stage, 0.0, 0.0])


#Set boundary conditions
domain.set_boundary({'left': Bd, 'right': Br, 'bottom': Br, 'top': Br})
                  

######################
#Initial condition
#
print 'Initial condition'
domain.set_quantity('stage', Constant_height(Z, 0.))


visualize = True
if visualize:
    from anuga.visualiser import RealtimeVisualiser
    vis = RealtimeVisualiser(domain)
    vis.render_quantity_height("elevation", zScale=100, offset = 5.0, dynamic=False)
    vis.render_quantity_height("stage", zScale =100, dynamic=True, opacity = 0.3)
    vis.colour_height_quantity('stage', (lambda q:q['stage'], -0.5, 0.5))
    vis.start()

time.sleep(2.0)



#Evolve
import time
t0 = time.time()

for t in domain.evolve(yieldstep = 0.5, finaltime = 1.0):
    domain.write_time()

    if visualize: vis.update()
    
if visualize: vis.evolveFinished()

    
print 'That took %.2f seconds' %(time.time()-t0)
print 'time', domain.write_time()

print domain.coordinates
print '*****'
print domain.vertex_coordinates
print '*****'
print domain.quantities['xmomentum'].centroid_values
print '*****'
print domain.quantities['xmomentum'].edge_values
print '*****'
print domain.quantities['stage'].vertex_values
print '*****'
print domain.quantities['stage'].explicit_update

from shallow_water import *

def compute_fluxes_python(domain):
    """Compute all fluxes and the timestep suitable for all volumes
    in domain.

    Compute total flux for each conserved quantity using "flux_function"

    Fluxes across each edge are scaled by edgelengths and summed up
    Resulting flux is then scaled by area and stored in
    explicit_update for each of the three conserved quantities 
    stage, xmomentum and ymomentum    

    The maximal allowable speed computed by the flux_function for each volume
    is converted to a timestep that must not be exceeded. The minimum of
    those is computed as the next overall timestep.

    Post conditions:
      domain.explicit_update is reset to computed flux values
      domain.timestep is set to the largest step satisfying all volumes. 
    """

    import sys
    from Numeric import zeros, Float

    N = domain.number_of_elements
    
    #Shortcuts
    Stage = domain.quantities['stage']
    Xmom = domain.quantities['xmomentum']
    Ymom = domain.quantities['ymomentum']
    Bed = domain.quantities['elevation']        

    #Arrays
    stage = Stage.edge_values
    xmom =  Xmom.edge_values
    ymom =  Ymom.edge_values
    bed =   Bed.edge_values    

    stage_bdry = Stage.boundary_values
    xmom_bdry =  Xmom.boundary_values
    ymom_bdry =  Ymom.boundary_values
    
    flux = zeros((N,3), Float) #Work array for summing up fluxes

    #Loop
    timestep = float(sys.maxint)    
    for k in range(N):

        for i in range(3):
            #Quantities inside volume facing neighbour i
            ql = [stage[k, i], xmom[k, i], ymom[k, i]]
            zl = bed[k, i]

            #Quantities at neighbour on nearest face
            n = domain.neighbours[k,i] 
            if n < 0:
                m = -n-1 #Convert negative flag to index
                qr = [stage_bdry[m], xmom_bdry[m], ymom_bdry[m]]
                zr = zl #Extend bed elevation to boundary
            else:    
                m = domain.neighbour_edges[k,i]
                qr = [stage[n, m], xmom[n, m], ymom[n, m]]
                zr = bed[n, m]                

                
            #Outward pointing normal vector   
            normal = domain.normals[k, 2*i:2*i+2]

            #Flux computation using provided function
            edgeflux, max_speed = flux_function(normal, ql, qr, zl, zr)
        
            flux[k,:] = edgeflux

    return flux
    
flux = compute_fluxes_python(domain)
print 'flux'
print flux

   
# THis was pulled out of 
