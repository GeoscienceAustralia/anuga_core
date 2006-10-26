"""Example of shallow water wave equation.

This is called Netherlands because it shows a dam with a gap in it and
stylised housed behind it and below the water surface.

"""

######################
# Module imports
#
#import rpdb
#rpdb.set_active()

from anuga.shallow_water import Domain, Reflective_boundary, Dirichlet_boundary,\
     Transmissive_boundary, Constant_height, Constant_stage

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
from Numeric import array
from anuga.visualiser.vtk_realtime_visualiser import Visualiser

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


N = 150 #size = 45000
N = 130 #size = 33800
N = 600 #Size = 720000
N = 40



#N = 15

print 'Creating domain'
#Create basic mesh
points, elements, boundary = rectangular_cross(N, N)

#Create shallow water domain
domain = Domain(points, elements, boundary, use_inscribed_circle=True)

domain.check_integrity()

#Setup order and all the beta's for the limiters (these should become defaults
domain.default_order = 2
domain.beta_w      = 1.0
domain.beta_w_dry  = 0.2
domain.beta_uh     = 1.0
domain.beta_uh_dry = 0.2
domain.beta_vh     = 1.0
domain.beta_vh_dry = 0.2

domain.alpha_balance = 10.0

#Output params
domain.smooth = False
domain.reduction = min  #Looks a lot better on top of steep slopes

print "Number of triangles = ", len(domain)
print "Extent = ", domain.get_extent()

#Set bed-slope and friction
inflow_stage = 0.5
manning = 0.03
Z = Weir(inflow_stage)

if N > 150:
    domain.visualise = False
    domain.checkpoint = False
    domain.store = True    #Store for visualisation purposes
    domain.format = 'sww'   #Native netcdf visualisation format
    import sys, os
    #FIXME: This was os.path.splitext but caused weird filenames based on root
    base = os.path.basename(sys.argv[0])
    basename, _ = os.path.splitext(base)
    domain.set_name(basename)
else:
    #domain.initialise_visualiser(rect=[0.0,0.0,1.0,1.0])
    #domain.initialise_visualiser()
    #domain.visualiser.coloring['stage'] = False
    domain.visualise_timer = True
    domain.checkpoint = False
    domain.store = False
    vis = Visualiser(domain,title="netherlands")
    vis.setup['elevation'] = True
    vis.updating['stage'] = True
    vis.qcolor['stage'] = (0.0,0.0,0.8)
    vis.coloring['stage']= False



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
Bd = Dirichlet_boundary([inflow_stage, 0.0, 0.0])


#Set boundary conditions
domain.set_boundary({'left': Bd, 'right': Br, 'bottom': Br, 'top': Br})


######################
#Initial condition
#
print 'Initial condition'
domain.set_quantity('stage', expression='elevation + 0.0')

#Evolve
import time
t0 = time.time()


#from realtime_visualisation_new import Visualiser
#V = Visualiser(domain,title='netherlands')
#V.update_quantity('stage')
#V.update_quantity('elevation')




for t in domain.evolve(yieldstep = 0.005, finaltime = None):
    domain.write_time()
    #domain.write_boundary()
    vis.update()
    print domain.quantities['stage'].get_values(location='centroids',
                                                indices=[0])
    #domain.visualiser.update_quantity('elevation')
    #time.sleep(0.1)
    #raw_input('pause>')
    #V.update_quantity('stage')
    #rpdb.set_active()
    #domain.visualiser.scale_z = 1.0
    #domain.visualiser.update_quantity_color('xmomentum',scale_z = 4.0)
    #integral_label.text='Integral=%10.5e'%domain.quantities['stage'].get_integral()


print 'That took %.2f seconds' %(time.time()-t0)
vis.shutdown()
