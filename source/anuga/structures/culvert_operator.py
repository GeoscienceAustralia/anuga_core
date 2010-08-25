import sys

from anuga.shallow_water.forcing import Inflow, General_forcing
from anuga.utilities.system_tools import log_to_file
from anuga.geometry.polygon import inside_polygon, is_inside_polygon
from anuga.geometry.polygon import plot_polygons, polygon_area

from anuga.utilities.numerical_tools import mean
from anuga.utilities.numerical_tools import ensure_numeric, sign
        
from anuga.config import g, epsilon
from anuga.config import minimum_allowed_height, velocity_protection        
import anuga.utilities.log as log

import inlet

import numpy as num
import math

class Below_interval(Exception): pass 
class Above_interval(Exception): pass


class Generic_box_culvert:
    """Culvert flow - transfer water from one rectangular box to another.
    Sets up the geometry of problem
    
    This is the base class for culverts. Inherit from this class (and overwrite
    compute_discharge method for specific subclasses)
    
    Input: Two points, pipe_size (either diameter or width, height), 
    mannings_rougness,
    """	

    def __init__(self,
                 domain,
                 end_point0=None, 
                 end_point1=None,
                 enquiry_gap_factor=0.2,
                 width=None,
                 height=None,
                 verbose=False):
        
        # Input check
        
        self.domain = domain

        self.domain.set_fractional_step_operator(self)

        self.end_points = [end_point0, end_point1]
        self.enquiry_gap_factor = enquiry_gap_factor
        
        if height is None:
            height = width

        self.width = width
        self.height = height
        
        self.verbose=verbose
        self.filename = None
       
        # Create the fundamental culvert polygons and create inlet objects
        self.create_culvert_polygons()

        #FIXME (SR) Put this into a foe loop to deal with more inlets
        self.inlets = []
        polygon0 = self.inlet_polygons[0]
        enquiry_pt0 = self.enquiry_points[0]
        inlet0_vector = self.culvert_vector
        self.inlets.append(inlet.Inlet(self.domain, polygon0, enquiry_pt0, inlet0_vector))

        polygon1 = self.inlet_polygons[1]
        enquiry_pt1 = self.enquiry_points[1]
        inlet1_vector = - self.culvert_vector
        self.inlets.append(inlet.Inlet(self.domain, polygon1, enquiry_pt1, inlet1_vector))
 

   
        self.print_stats()


    def __call__(self):

        # Time stuff
        time     = self.domain.get_time()
        timestep = self.domain.get_timestep()



        inflow  = self.inlets[0]
        outflow = self.inlets[1]

        # Determine flow direction based on total energy difference
        delta_total_energy = inflow.get_average_total_energy() - outflow.get_average_total_energy()

        if delta_total_energy < 0:
            inflow  = self.inlets[1]
            outflow = self.inlets[0]
            delta_total_energy = -delta_total_energy

        delta_z = inflow.get_average_elevation() - outflow.get_average_elevation()
        culvert_slope = delta_z/self.culvert_length

        # Determine controlling energy (driving head) for culvert
        if inflow.get_average_specific_energy() > delta_total_energy:
            # Outlet control
            driving_head = delta_total_energy
        else:
            # Inlet control
            driving_head = inflow.get_average_specific_energy()
            

        # Transfer
        from culvert_routines import boyd_generalised_culvert_model
        Q, barrel_velocity, culvert_outlet_depth =\
                    boyd_generalised_culvert_model(inflow.get_average_height(),
                                         outflow.get_average_height(),
                                         inflow.get_average_speed(),
                                         outflow.get_average_speed(),
                                         inflow.get_average_specific_energy(),
                                         delta_total_energy,
                                         g,
                                         culvert_length=self.culvert_length,
                                         culvert_width=self.width,
                                         culvert_height=self.height,
                                         culvert_type='box',
                                         manning=0.01)

        transfer_water = Q*timestep


        inflow.set_heights(inflow.get_average_height() - transfer_water)
        inflow.set_xmoms(0.0)
        inflow.set_ymoms(0.0)


        outflow.set_heights(outflow.get_average_height() + transfer_water)
        outflow.set_xmoms(0.0)
        outflow.set_ymoms(0.0)



    def print_stats(self):

        print '====================================='
        print 'Generic Culvert Operator'
        print '====================================='
        print "enquiry_gap_factor"
        print self.enquiry_gap_factor
        
        for i, inlet in enumerate(self.inlets):
            print '-------------------------------------'
            print 'Inlet %i' % i
            print '-------------------------------------'

            print 'inlet triangle indices and centres'
            print inlet.triangle_indices[i]
            print self.domain.get_centroid_coordinates()[inlet.triangle_indices[i]]
        
            print 'polygon'
            print inlet.polygon

            print 'enquiry_point'
            print inlet.enquiry_point

        print '====================================='





    def create_culvert_polygons(self):

        """Create polygons at the end of a culvert inlet and outlet.
        At either end two polygons will be created; one for the actual flow to pass through and one a little further away
        for enquiring the total energy at both ends of the culvert and transferring flow.
        """

        # Calculate geometry
        x0, y0 = self.end_points[0]
        x1, y1 = self.end_points[1]

        dx = x1 - x0
        dy = y1 - y0

        self.culvert_vector = num.array([dx, dy])
        self.culvert_length = math.sqrt(num.sum(self.culvert_vector**2))
        assert self.culvert_length > 0.0, 'The length of culvert is less than 0'

        # Unit direction vector and normal
        self.culvert_vector /= self.culvert_length                      # Unit vector in culvert direction
        self.culvert_normal = num.array([-dy, dx])/self.culvert_length  # Normal vector

        # Short hands
        w = 0.5*self.width*self.culvert_normal # Perpendicular vector of 1/2 width
        h = self.height*self.culvert_vector    # Vector of length=height in the
                             # direction of the culvert
        gap = (1 + self.enquiry_gap_factor)*h

        self.inlet_polygons = []
        self.enquiry_points = []

        # Build exchange polygon and enquiry points 0 and 1
        for i in [0, 1]:
            i0 = (2*i-1)
            p0 = self.end_points[i] + w
            p1 = self.end_points[i] - w
            p2 = p1 + i0*h
            p3 = p0 + i0*h
            self.inlet_polygons.append(num.array([p0, p1, p2, p3]))
            self.enquiry_points.append(self.end_points[i] + i0*gap)

        # Check that enquiry points are outside inlet polygons
        for i in [0,1]:
            polygon = self.inlet_polygons[i]
            # FIXME (SR) Probably should calculate the area of all the triangles
            # associated with this polygon, as there is likely to be some
            # inconsistency between triangles and ploygon
            area = polygon_area(polygon)
            

            msg = 'Polygon %s ' %(polygon)
            msg += ' has area = %f' % area
            assert area > 0.0, msg

            for j in [0,1]:
                point = self.enquiry_points[j]
                msg = 'Enquiry point falls inside a culvert polygon.'

                assert not inside_polygon(point, polygon), msg

    

                        
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

