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
                 apron=None,
                 manning=0.013,
                 enquiry_gap=0.2,
                 use_momentum_jet=True,
                 use_velocity_head=True,
                 verbose=False):
        
        self.domain = domain
        self.domain.set_fractional_step_operator(self)
        self.end_points = [end_point0, end_point1]
        
        if height is None:
            height = width

        if apron is None:
            apron = width

        self.width  = width
        self.height = height
        self.apron  = apron
        self.manning = manning
        self.enquiry_gap = enquiry_gap
        self.verbose = verbose

        self.use_momentum_jet = use_momentum_jet
        self.use_velocity_head= use_velocity_head
       
        self.culvert = Boyd_box_culvert(self.domain,
                                        self.end_points,
                                        self.width,
                                        self.height,
                                        self.apron,
                                        self.manning,
                                        self.enquiry_gap,
                                        self.use_velocity_head,
                                        self.verbose)
        
        self.routine = self.culvert.routine

        self.inlets  = self.culvert.get_inlets()

        if self.verbose:
            self.print_stats()


    def __call__(self):

        timestep = self.domain.get_timestep()
        
        Q, barrel_speed, outlet_depth = self.routine()


        inflow  = self.routine.get_inflow()
        outflow = self.routine.get_outflow()


        old_inflow_height = inflow.get_average_height()
        old_inflow_xmom = inflow.get_average_xmom()
        old_inflow_ymom = inflow.get_average_ymom()
            
        if old_inflow_height > 0.0 :
                Qstar = Q/old_inflow_height
        else:
                Qstar = 0.0

        factor = 1.0/(1.0 + Qstar*timestep/inflow.get_area())

            
            
        new_inflow_height = old_inflow_height*factor
        new_inflow_xmom = old_inflow_xmom*factor
        new_inflow_ymom = old_inflow_ymom*factor
            

        inflow.set_heights(new_inflow_height)

        #inflow.set_xmoms(Q/inflow.get_area())
        #inflow.set_ymoms(0.0)


        inflow.set_xmoms(new_inflow_xmom)
        inflow.set_ymoms(new_inflow_ymom)


        loss = (old_inflow_height - new_inflow_height)*inflow.get_area()

            
        # set outflow
        if old_inflow_height > 0.0 :
                timestep_star = timestep*new_inflow_height/old_inflow_height
        else:
            timestep_star = 0.0

            
        outflow_extra_height = Q*timestep_star/outflow.get_area()
        outflow_direction = - outflow.outward_culvert_vector
        outflow_extra_momentum = outflow_extra_height*barrel_speed*outflow_direction
            

        gain = outflow_extra_height*outflow.get_area()
            
        #print Q, Q*timestep, barrel_speed, outlet_depth, Qstar, factor, timestep_star
        #print '  ', loss, gain



        new_outflow_height = outflow.get_average_height() + outflow_extra_height


        if self.use_momentum_jet :
            # FIXME (SR) Review momentum to account for possible hydraulic jumps at outlet
            #new_outflow_xmom = outflow.get_average_xmom() + outflow_extra_momentum[0]
            #new_outflow_ymom = outflow.get_average_ymom() + outflow_extra_momentum[1]

            new_outflow_xmom = barrel_speed*new_outflow_height*outflow_direction[0]
            new_outflow_ymom = barrel_speed*new_outflow_height*outflow_direction[1]

        else:
            #new_outflow_xmom = outflow.get_average_xmom()
            #new_outflow_ymom = outflow.get_average_ymom()

            new_outflow_xmom = 0.0
            new_outflow_ymom = 0.0


        outflow.set_heights(new_outflow_height)
        outflow.set_xmoms(new_outflow_xmom)
        outflow.set_ymoms(new_outflow_ymom)


            
        #print '   outflow volume ',outflow.get_total_water_volume()

    def print_stats(self):

        print '====================================='
        print 'Generic Culvert Operator'
        print '====================================='

        print 'Culvert'
        print self.culvert

        print 'Culvert Routine'
        print self.routine
        
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

