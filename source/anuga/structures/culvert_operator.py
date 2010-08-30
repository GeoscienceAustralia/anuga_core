from anuga.geometry.polygon import inside_polygon, polygon_area
from anuga.config import g
import anuga.utilities.log as log

from boyd_box_culvert import Boyd_box_culvert

class Culvert_operator:
    """Culvert flow - transfer water from one rectangular box to another.
    Sets up the geometry of problem
    
    This is the base class for culverts. Inherit from this class (and overwrite
    compute_discharge method for specific subclasses)
    
    Input: Two points, pipe_size (either diameter or width, height), 
    mannings_rougness,
    """	

    def __init__(self,
                 domain,
                 end_point0, 
                 end_point1,
                 width,
                 height=None,
                 verbose=False):
        
        self.domain = domain
        self.domain.set_fractional_step_operator(self)
        end_points = [end_point0, end_point1]
        
        if height is None:
            height = width

        self.width = width
        self.height = height
       
        self.culvert = Boyd_box_culvert(self.domain, end_points, self.width, self.height)
        print self.culvert
        self.routine = self.culvert.routine
        print self.routine
        self.inlets = self.culvert.get_inlets()
   
        self.print_stats()


    def __call__(self):

        timestep = self.domain.get_timestep()
        
        Q, barrel_speed, culvert_outlet_depth = self.routine()

        inflow  = self.routine.get_inflow()
        outflow = self.routine.get_outflow()

        outflow_direction = - outflow.outward_culvert_vector

        outflow_momentum_flux = barrel_speed**2*culvert_outlet_depth*outflow_direction


        print Q, barrel_speed, culvert_outlet_depth, outflow_momentum_flux

        #FIXME (SR) Check whether we need to mult/divide by inlet area
        inflow_transfer =  Q*timestep/inflow.get_area()

        outflow_transfer = Q*timestep/outflow.get_area()



        inflow.set_heights(inflow.get_average_height() - inflow_transfer)

        inflow.set_xmoms(0.0)
        inflow.set_ymoms(0.0)

        #u = outflow.get_xvelocities()
        #v = outflow.get_yvelocities()

        outflow.set_heights(outflow.get_average_height() + outflow_transfer)
        #outflow.set_xmoms(outflow.get_xmoms() + timestep*outflow_momentum_flux[0] )
        #outflow.set_ymoms(outflow.get_ymoms() + timestep*outflow_momentum_flux[1] )

    def print_stats(self):

        print '====================================='
        print 'Generic Culvert Operator'
        print '====================================='
        
        for i, inlet in enumerate(self.inlets):
            print '-------------------------------------'
            print 'Inlet %i' % i
            print '-------------------------------------'

            print 'inlet triangle indices and centres'
            print inlet.triangle_indices
            print self.domain.get_centroid_coordinates()[inlet.triangle_indices]
        
            print 'polygon'
            print inlet.polygon

        print '====================================='

                        
# FIXME(Ole): Write in C and reuse this function by similar code
# in interpolate.py
def interpolate_linearly(x, xvec, yvec):

    msg = 'Input to function interpolate_linearly could not be converted '
    msg += 'to numerical scalar: x = %s' % str(x)
    try:
        x = float(x)
    except:
        raise Exception, msg


    # Check bounds
    if x < xvec[0]:
        msg = 'Value provided = %.2f, interpolation minimum = %.2f.'\
            % (x, xvec[0])
        raise Below_interval, msg

    if x > xvec[-1]:
        msg = 'Value provided = %.2f, interpolation maximum = %.2f.'\
            %(x, xvec[-1])
        raise Above_interval, msg


    # Find appropriate slot within bounds
    i = 0
    while x > xvec[i]: i += 1


    x0 = xvec[i-1]
    x1 = xvec[i]
    alpha = (x - x0)/(x1 - x0)

    y0 = yvec[i-1]
    y1 = yvec[i]
    y = alpha*y1 + (1-alpha)*y0

    return y



def read_culvert_description(culvert_description_filename):

    # Read description file
    fid = open(culvert_description_filename)

    read_rating_curve_data = False
    rating_curve = []
    for i, line in enumerate(fid.readlines()):

        if read_rating_curve_data is True:
            fields = line.split(',')
            head_difference = float(fields[0].strip())
            flow_rate = float(fields[1].strip())
            barrel_velocity = float(fields[2].strip())

            rating_curve.append([head_difference, flow_rate, barrel_velocity])

        if i == 0:
            # Header
            continue
        if i == 1:
            # Metadata
            fields = line.split(',')
            label=fields[0].strip()
            type=fields[1].strip().lower()
            assert type in ['box', 'pipe']

            width=float(fields[2].strip())
            height=float(fields[3].strip())
            length=float(fields[4].strip())
            number_of_barrels=int(fields[5].strip())
            #fields[6] refers to losses
            description=fields[7].strip()

        if line.strip() == '': continue # Skip blanks

        if line.startswith('Rating'):
            read_rating_curve_data = True
            # Flow data follows

    fid.close()

    return label, type, width, height, length, number_of_barrels, description, rating_curve

