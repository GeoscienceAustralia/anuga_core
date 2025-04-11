
from anuga.geometry.polygon import inside_polygon, is_inside_polygon, line_intersect
from anuga.config import velocity_protection, g
import math

import numpy as num

from . import inlet

class Inlet_enquiry(inlet.Inlet):
    """Contains information associated with each inlet plus an enquiry point
    """

    def __init__(self,
                 domain,
                 region,
                 enquiry_pt,
                 invert_elevation = None,
                 outward_culvert_vector=None,
                 verbose=False):
        """
        region can be a line or a polygon or a region object
        """

        #print region
        #print enquiry_pt
        
        inlet.Inlet.__init__(self, domain, region, verbose)


        
        self.enquiry_pt = enquiry_pt
        self.invert_elevation = invert_elevation
        self.outward_culvert_vector = outward_culvert_vector

        self.compute_enquiry_index()


    def compute_enquiry_index(self):

        # Get boundary (in absolute coordinates)
        bounding_polygon = self.domain_bounding_polygon
        #domain_centroids = self.domain.get_centroid_coordinates(absolute=True)
        #vertex_coordinates = self.domain.get_vertex_coordinates(absolute=True)

                
        point = self.enquiry_pt
        msg = 'Enquiry Point %s ' %  str(point)
        msg += ' did not fall within the domain boundary.'
        assert is_inside_polygon(point, bounding_polygon), msg

        try:
            self.enquiry_index = self.domain.get_triangle_containing_point(self.enquiry_pt)
        except:
            msg = "Enquiry point %s doesn't intersect mesh, maybe inside a building, try reducing enquiry_gap" % str(self.enquiry_pt)
            raise Exception(msg)

        
        if self.enquiry_index in self.triangle_indices:
            msg = 'Enquiry point %s' % (self.enquiry_pt)
            msg += ' is in an inlet triangle'
            import warnings
            warnings.warn(msg)
            

    def get_enquiry_position(self):

        return self.domain.get_centroid_coordinates(absolute=True)[self.enquiry_index]

    def get_enquiry_stage(self):

        return self.domain.quantities['stage'].centroid_values[self.enquiry_index]


    def get_enquiry_xmom(self):

        return self.domain.quantities['xmomentum'].centroid_values[self.enquiry_index]

    def get_enquiry_ymom(self):

        return self.domain.quantities['ymomentum'].centroid_values[self.enquiry_index]


    def get_enquiry_elevation(self):

        return self.domain.quantities['elevation'].centroid_values[self.enquiry_index]

    def get_enquiry_depth(self):

        return max(self.get_enquiry_stage() - self.get_enquiry_invert_elevation(), 0.0)


    def get_enquiry_water_depth(self):

        return self.get_enquiry_stage() - self.get_enquiry_elevation()


    def get_enquiry_invert_elevation(self):

        if  self.invert_elevation is None:
            return self.get_enquiry_elevation()
        else:
            return self.invert_elevation


    def get_enquiry_velocity(self):

            depth = self.get_enquiry_water_depth()
            u = depth*self.get_enquiry_xmom()/(depth**2 + velocity_protection)
            v = depth*self.get_enquiry_ymom()/(depth**2 + velocity_protection)

            return u, v


    def get_enquiry_xvelocity(self):

            depth = self.get_enquiry_water_depth()
            return depth*self.get_enquiry_xmom()/(depth**2 + velocity_protection)

    def get_enquiry_yvelocity(self):

            depth = self.get_enquiry_water_depth()
            return depth*self.get_enquiry_ymom()/(depth**2 + velocity_protection)


    def get_enquiry_speed(self):

            u, v = self.get_enquiry_velocity()

            return math.sqrt(u**2 + v**2)


    def get_enquiry_velocity_head(self):

        # velocity head will be zero if flowing out of inlet


        if self.domain.use_new_velocity_head:
            u, v   = self.get_enquiry_velocity()
            n1, n2 = self.outward_culvert_vector
            normal_speed = min(u*n1 + v*n2, 0.0)

            velocity_head = 0.5*normal_speed**2/g
        else:
            velocity_head = 0.5*self.get_enquiry_speed()**2/g

        return velocity_head

        #return 0.5*self.get_enquiry_speed()**2/g


    def get_enquiry_total_energy(self):

        return self.get_enquiry_velocity_head() + self.get_enquiry_stage()


    def get_enquiry_specific_energy(self):

        return self.get_enquiry_velocity_head() + self.get_enquiry_depth()


