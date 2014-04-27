"""Class Domain -
1D interval domains for finite-volume computations of
the shallow water wave equation.

This module contains a specialisation of class Generic_domain from module
generic_domain.py
consisting of methods specific to Channel flow using the Shallow Water Wave Equation

This particular modification of the Domain class implements the ability to
vary the width of the 1D channel that the water flows in. As a result the
conserved variables are different than previous implementations and so are the
equations.

U_t + E_x = S

where

U = [A, Q] = [b*h, u*b*h]
E = [Q, Q^2/A + g*b*h^2/2]
S represents source terms forcing the system
(e.g. gravity, boundary_stree, friction, wind stress, ...)
gravity = -g*b*h*z_x
boundary_stress = 1/2*g*b_x*h^2

and _t, _x, _y denote the derivative with respect to t, x and y respectiely.

The quantities are

symbol    variable name    explanation
A         area             Wetted area = b*h
Q         discharge        flux of water = u*b*h
x         x                horizontal distance from origin [m]
z         elevation        elevation of bed on which flow is modelled [m]
h         height           water height above z [m]
w         stage            absolute water level, w = z+h [m]
u                          speed in the x direction [m/s]
uh        xmomentum        momentum in the x direction [m^2/s]
b         width            width of channel
eta                        mannings friction coefficient [to appear]
nu                         wind stress coefficient [to appear]

The conserved quantities are A, Q
--------------------------------------------------------------------------
For details see e.g.
Christopher Zoppou and Stephen Roberts,
Catastrophic Collapse of Water Supply Reservoirs in Urban Areas,
Journal of Hydraulic Engineering, vol. 127, No. 7 July 1999


John Jakeman, Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou,
Padarn Wilson, Geoscience Australia, 2008
"""


from anuga_1d.base.generic_domain import *
import numpy


#Shallow water domain
class Domain(Generic_domain):

    def __init__(self, coordinates, boundary = None, tagged_elements = None):

        conserved_quantities = ['area', 'discharge']
        evolved_quantities = ['area', 'discharge', 'elevation', 'height', 'velocity','width','stage']
        other_quantities = ['friction']
        Generic_domain.__init__(self,
                                coordinates = coordinates,
                                boundary = boundary,
                                conserved_quantities = conserved_quantities,
                                evolved_quantities = evolved_quantities,
                                other_quantities = other_quantities,
                                tagged_elements = tagged_elements)
        
        from anuga_1d.config import minimum_allowed_height, g, h0
        self.minimum_allowed_height = minimum_allowed_height
        self.g = g
        self.h0 = h0
        self.setstageflag = False
        self.discontinousb = False


        # forcing terms gravity and boundary stress are included in the flux calculation
        #self.forcing_terms.append(gravity)
        #self.forcing_terms.append(boundary_stress)
        #self.forcing_terms.append(manning_friction)
       

        
        #Stored output
        self.store = True
        self.format = 'sww'
        self.smooth = True

        
        #Reduction operation for get_vertex_values
        from anuga_1d.base.util import mean
        self.reduction = mean
        #self.reduction = min  #Looks better near steep slopes

        self.set_quantities_to_be_stored(['area','discharge'])

        self.__doc__ = 'channel_domain'

        self.check_integrity()


    def check_integrity(self):

        #Check that we are solving the shallow water wave equation

        msg = 'First conserved quantity must be "area"'
        assert self.conserved_quantities[0] == 'area', msg
        msg = 'Second conserved quantity must be "discharge"'
        assert self.conserved_quantities[1] == 'discharge', msg

        msg = 'First evolved quantity must be "area"'
        assert self.evolved_quantities[0] == 'area', msg
        msg = 'Second evolved quantity must be "discharge"'
        assert self.evolved_quantities[1] == 'discharge', msg
        msg = 'Third evolved quantity must be "elevation"'
        assert self.evolved_quantities[2] == 'elevation', msg
        msg = 'Fourth evolved quantity must be "height"'
        assert self.evolved_quantities[3] == 'height', msg
        msg = 'Fifth evolved quantity must be "velocity"'
        assert self.evolved_quantities[4] == 'velocity', msg
        msg = 'Sixth evolved quantity must be "width"'
        assert self.evolved_quantities[5] == 'width', msg
        msg = 'Seventh evolved quantity must be "stage"'
        assert self.evolved_quantities[6] == 'stage', msg

        Generic_domain.check_integrity(self)

    def compute_fluxes(self):
        #Call correct module function
        #(either from this module or C-extension)
        compute_fluxes_channel(self)

    def distribute_to_vertices_and_edges(self):
        #Call correct module function
        #(either from this module or C-extension)
        #distribute_to_vertices_and_edges_limit_s_v_h(self)
        distribute_to_vertices_and_edges_limit_s_v(self)

    def update_derived_quantites(self):

        pass



    def initialize_plotting(self,
                            stage_lim = [-1.0, 40.0],
                            width_lim = [-1.0, 10.0],
                            velocity_lim = [-5.0, 5.0]):

        import pylab
        pylab.ion()


        x = self.get_vertices().flatten()

        z = self.quantities['elevation'].vertex_values.flatten()
        w = self.quantities['stage'].vertex_values.flatten()
        h = self.quantities['height'].vertex_values.flatten()
        v = self.quantities['velocity'].vertex_values.flatten()
        b = self.quantities['width'].vertex_values.flatten()

        print x.shape
        print z.shape

        #-------------------------------
        # Top plot
        #-------------------------------
        self.plot1 = pylab.subplot(311)

        self.zplot, = pylab.plot(x, z, 'k')
        self.wplot, = pylab.plot(x, w, 'k')

        self.plot1.set_ylim(stage_lim)
        #pylab.xlabel('Position')
        pylab.ylabel('Bed/Stage')

        #-------------------------------
        # Middle Plot
        #-------------------------------
        self.plot2 = pylab.subplot(312)

        self.bplot, = pylab.plot(x, b, 'k')

        self.plot2.set_ylim(width_lim)
        #pylab.xlabel('Position')
        pylab.ylabel('Width')

        #-------------------------------
        # Bottom Plot
        #-------------------------------
        self.plot3 = pylab.subplot(313)

        self.vplot, = pylab.plot(x, v, 'k')

        self.plot3.set_ylim(velocity_lim)

        pylab.xlabel('Position')
        pylab.ylabel('Velocity')


    def update_plotting(self):

        import pylab

        #x = self.get_vertices().flatten()
        z = self.quantities['elevation'].vertex_values.flatten()
        w = self.quantities['stage'].vertex_values.flatten()
        h = self.quantities['height'].vertex_values.flatten()
        v = self.quantities['velocity'].vertex_values.flatten()
        b = self.quantities['width'].vertex_values.flatten()


        self.zplot.set_ydata(z)
        self.wplot.set_ydata(w)
        self.bplot.set_ydata(b)
        self.vplot.set_ydata(v)

        pylab.draw()


    def hold_plotting(self,save=None):

        self.update_plotting()
        import pylab
        
        pylab.ioff()

        if save != None:
            file = save+".pdf"
            pylab.savefig(file)

        pylab.show()



    def finalize_plotting(self):

        pass



#=============== End of Channel Domain ===============================

#-----------------------------------
# Compute flux definition with channel
#-----------------------------------
def compute_fluxes_channel(domain):
    import sys
    timestep = float(sys.maxint)

    area       = domain.quantities['area']
    discharge  = domain.quantities['discharge']
    bed        = domain.quantities['elevation']
    height     = domain.quantities['height']
    velocity   = domain.quantities['velocity']
    width      = domain.quantities['width']


    from anuga_1d.channel.channel_domain_ext import compute_fluxes_channel_ext
    domain.flux_timestep = compute_fluxes_channel_ext(timestep,domain,area,discharge,bed,height,velocity,width)

#-----------------------------------------------------------------------
# Distribute to verticies with stage, velocity and channel geometry
# reconstructed and then extrapolated.
#-----------------------------------------------------------------------
def distribute_to_vertices_and_edges_limit_s_v(domain):
    import sys
    from anuga_1d.config import epsilon, h0

    N = domain.number_of_elements

    #Shortcuts
    area      = domain.quantities['area']
    discharge = domain.quantities['discharge']
    bed       = domain.quantities['elevation']
    height    = domain.quantities['height']
    velocity  = domain.quantities['velocity']
    width     = domain.quantities['width']
    stage     = domain.quantities['stage']

    #Arrays   
    a_C   = area.centroid_values
    d_C   = discharge.centroid_values
    z_C   = bed.centroid_values
    h_C   = height.centroid_values
    u_C   = velocity.centroid_values
    b_C   = width.centroid_values
    w_C   = stage.centroid_values

    # Calculate height, velocity and stage.
    # Here we assume the conserved quantities and the channel geometry
    # (i.e. bed and width) have been accurately computed in the previous
    # timestep. 
    h_C[:] = numpy.where(a_C > 0.0, a_C/b_C, 0.0)
    u_C[:] = numpy.where(a_C > 0.0, d_C/a_C, 0.0)

    w_C[:] = h_C + z_C

    #print w_C

    # Extrapolate velocity and stage as well as channel geometry.
    for name in ['velocity', 'stage', 'elevation', 'width']:
        Q = domain.quantities[name]
        if domain.order == 1:
            Q.extrapolate_first_order()
        elif domain.order == 2:
            Q.extrapolate_second_order()
        else:
            raise 'Unknown order'

    # Stage, bed, width and velocity have been extrapolated
    w_V  = stage.vertex_values
    u_V  = velocity.vertex_values
    z_V  = bed.vertex_values
    b_V  = width.vertex_values

    # Need to update these vertex_values
    a_V  = area.vertex_values
    h_V  = height.vertex_values
    d_V  = discharge.vertex_values

    # Calculate height and fix up negatives. The main idea here is
    # fix up the wet/dry interface.
    h_V[:,:] = w_V - z_V

    h_0 = numpy.where(h_V[:,0] < 0.0, 0.0, h_V[:,0])
    h_1 = numpy.where(h_V[:,0] < 0.0, h_V[:,1]+h_V[:,0], h_V[:,1])

    h_V[:,0] = h_0
    h_V[:,1] = h_1


    h_0 = numpy.where(h_V[:,1] < 0.0, h_V[:,1]+h_V[:,0], h_V[:,0])
    h_1 = numpy.where(h_V[:,1] < 0.0, 0.0, h_V[:,1])

    h_V[:,0] = h_0
    h_V[:,1] = h_1


    # Protect against negative and small heights. If we set h to zero
    # we better do the same with velocity (i.e. no water, no velocity).
    h_V[:,:] = numpy.where (h_V <= h0, 0.0, h_V)
    u_V[:,:] = numpy.where (h_V <= h0, 0.0, u_V)


    # Clean up conserved quantities
    w_V[:] = z_V + h_V
    a_V[:] = b_V * h_V
    d_V[:] = u_V * a_V

    
    return

#-----------------------------------------------------------------------
# Distribute to verticies with stage, height and velocity reconstructed 
# and then extrapolated.
# In this method, we extrapolate the stage and height out to the vertices.
# The bed, although given as initial data to the problem, is reconstructed
# from the stage and height. This ensures consistency of the reconstructed
# quantities (i.e. w = z + h) as well as protecting against negative
# heights.
#-----------------------------------------------------------------------
def distribute_to_vertices_and_edges_limit_s_v_h(domain):
    import sys
    from anuga_1d.config import epsilon, h0

    N = domain.number_of_elements

    #Shortcuts
    area      = domain.quantities['area']
    discharge = domain.quantities['discharge']
    bed       = domain.quantities['elevation']
    height    = domain.quantities['height']
    velocity  = domain.quantities['velocity']
    width     = domain.quantities['width']
    stage     = domain.quantities['stage']

    #Arrays   
    a_C   = area.centroid_values
    d_C   = discharge.centroid_values
    z_C   = bed.centroid_values
    h_C   = height.centroid_values
    u_C   = velocity.centroid_values
    b_C   = width.centroid_values
    w_C   = stage.centroid_values

    # Construct h,u,w from the conserved quantities after protecting
    # conserved quantities from becoming too small.
    a_C[:] = numpy.where( (a_C>h0), a_C, 0.0 )
    d_C[:] = numpy.where( (a_C>h0), d_C, 0.0 )
    h_C[:] = numpy.where( (b_C>h0), a_C/(b_C + h0/b_C), 0.0 )
    u_C[:] = numpy.where( (a_C>h0), d_C/(a_C + h0/a_C), 0.0 )    

    # Set the stage
    w_C[:] = h_C + z_C		

    # Extrapolate "fundamental" quantities.
    # All other quantities will be reconstructed from these.
    for name in ['velocity', 'stage',  'height', 'width']:
        Q = domain.quantities[name]
        if domain.order == 1:
            Q.extrapolate_first_order()
        elif domain.order == 2:
            Q.extrapolate_second_order()
        else:
            raise 'Unknown order'

    # These quantities have been extrapolated.
    u_V  = velocity.vertex_values
    w_V  = stage.vertex_values
    h_V  = height.vertex_values
    b_V  = width.vertex_values

    # These need to be reconstructed
    a_V  = area.vertex_values
    d_V  = discharge.vertex_values
    z_V  = bed.vertex_values

    # Reconstruct bed from stage and height.
    z_V[:] = w_V-h_V

    # Now reconstruct our conserved quantities from the above
    # reconstructed quantities.
    a_V[:] = b_V*h_V
    d_V[:] = u_V*a_V

    return


#--------------------------------------------------------
#Boundaries - specific to the channel_domain
#--------------------------------------------------------
class Reflective_boundary(Boundary):
    """Reflective boundary returns same conserved quantities as
    those present in its neighbour volume but reflected.

    This class is specific to the shallow water equation as it
    works with the momentum quantities assumed to be the second
    and third conserved quantities.
    """

    def __init__(self, domain = None):
        Boundary.__init__(self)

        if domain is None:
            msg = 'Domain must be specified for reflective boundary'
            raise msg

        #Handy shorthands
        self.normals  = domain.normals
        self.area     = domain.quantities['area'].vertex_values
        self.discharge     = domain.quantities['discharge'].vertex_values
        self.bed      = domain.quantities['elevation'].vertex_values
        self.height   = domain.quantities['height'].vertex_values
        self.velocity = domain.quantities['velocity'].vertex_values
        self.width    = domain.quantities['width'].vertex_values
        self.stage    = domain.quantities['stage'].vertex_values

        self.evolved_quantities = numpy.zeros(7, numpy.float)

    def __repr__(self):
        return 'Reflective_boundary'


    def evaluate(self, vol_id, edge_id):
        """Reflective boundaries reverses the outward momentum
        of the volume they serve.
        """

        q = self.evolved_quantities
        q[0] =  self.area[vol_id, edge_id]
        q[1] = -self.discharge[vol_id, edge_id]
        q[2] =  self.bed[vol_id, edge_id]
        q[3] =  self.height[vol_id, edge_id]
        q[4] = -self.velocity[vol_id, edge_id]
        q[5] =  self.width[vol_id,edge_id]
        q[6] =  self.stage[vol_id,edge_id]

        return q

class Dirichlet_boundary(Boundary):
    """Dirichlet boundary returns constant values for the
    conserved quantities
if k>5 and k<15:
            print discharge_ud[k],-g*zx*avg_h*avg_b
        discharge_ud[k] +=-g*zx*avg_h*avg_b    """


    def __init__(self, evolved_quantities=None):
        Boundary.__init__(self)

        if evolved_quantities is None:
            msg = 'Must specify one value for each evolved quantity'
            raise msg

        assert len(evolved_quantities) == 7

        self.evolved_quantities=numpy.array(evolved_quantities,numpy.float)

    def __repr__(self):
        return 'Dirichlet boundary (%s)' %self.evolved_quantities

    def evaluate(self, vol_id=None, edge_id=None):
        return self.evolved_quantities


#----------------------------
#Standard forcing terms:
#---------------------------
def gravity(domain):
    """Apply gravitational pull in the presence of bed slope
    """

    from util import gradient
    from Numeric import zeros, Float, array, sum



    Area      = domain.quantities['area']
    Discharge  = domain.quantities['discharge']
    Elevation = domain.quantities['elevation']
    Height    = domain.quantities['height']
    Width     = domain.quantities['width']

    discharge_ud  = Discharge.explicit_update


       
    h = Height.vertex_values
    b = Width.vertex_values
    a = Area.vertex_values
    z = Elevation.vertex_values

    x = domain.get_vertex_coordinates()
    g = domain.g
    for k in range(domain.number_of_elements):
        avg_h = 0.5*(h[k,0] + h[k,1])
        avg_b = 0.5*(b[k,0] + b[k,1])
        
        #Compute bed slope
        x0, x1 = x[k,:]
        z0, z1 = z[k,:]
        zx = gradient(x0, x1, z0, z1)
       
        #Update momentum (explicit update is reset to source values)
        discharge_ud[k]+= -g*zx*avg_h*avg_b
      

def boundary_stress(domain):

    from util import gradient
    from Numeric import zeros, Float, array, sum

    Area     = domain.quantities['area']
    Discharge  = domain.quantities['discharge']
    Elevation = domain.quantities['elevation']
    Height    = domain.quantities['height']
    Width     = domain.quantities['width']

    discharge_ud  = Discharge.explicit_update
 
    h = Height.vertex_values
    b = Width.vertex_values
    a = Area.vertex_values
    z = Elevation.vertex_values

    x = domain.get_vertex_coordinates()
    g = domain.g

    for k in range(domain.number_of_elements):
        avg_h = 0.5*(h[k,0] + h[k,1])
        

        #Compute bed slope
        x0, x1 = x[k,:]
        b0, b1 = b[k,:]
        bx = gradient(x0, x1, b0, b1)
        
        #Update momentum (explicit update is reset to source values)
        discharge_ud[k] += 0.5*g*bx*avg_h*avg_h
        #stage_ud[k] = 0.0
 
 
def manning_friction(domain):
    """Apply (Manning) friction to water momentum
    """

    from math import sqrt

    w = domain.quantities['stage'].centroid_values
    z = domain.quantities['elevation'].centroid_values
    h = w-z

    uh = domain.quantities['xmomentum'].centroid_values
    eta = domain.quantities['friction'].centroid_values

    xmom_update = domain.quantities['xmomentum'].semi_implicit_update

    N = domain.number_of_elements
    eps = domain.minimum_allowed_height
    g = domain.g

    for k in range(N):
        if eta[k] >= eps:
            if h[k] >= eps:
            	#S = -g * eta[k]**2 * sqrt((uh[k]**2 + vh[k]**2))
                S = -g * eta[k]**2 * uh[k]
            	S /= h[k]**(7.0/3)

            	#Update momentum
            	xmom_update[k] += S*uh[k]
            	#ymom_update[k] += S*vh[k]

def linear_friction(domain):
    """Apply linear friction to water momentum

    Assumes quantity: 'linear_friction' to be present
    """

    from math import sqrt

    w = domain.quantities['stage'].centroid_values
    z = domain.quantities['elevation'].centroid_values
    h = w-z

    uh = domain.quantities['xmomentum'].centroid_values
    tau = domain.quantities['linear_friction'].centroid_values

    xmom_update = domain.quantities['xmomentum'].semi_implicit_update

    N = domain.number_of_elements
    eps = domain.minimum_allowed_height

    for k in range(N):
        if tau[k] >= eps:
            if h[k] >= eps:
            	S = -tau[k]/h[k]

            	#Update momentum
            	xmom_update[k] += S*uh[k]




def linearb(domain):

    bC = domain.quantities['width'].vertex_values
    
    for i in range(len(bC)-1):
        temp= 0.5*(bC[i,1]+bC[i+1,0])
        bC[i,1]=temp
        bC[i+1,0]=temp



