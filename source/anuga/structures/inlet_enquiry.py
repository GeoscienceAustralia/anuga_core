from anuga.geometry.polygon import inside_polygon, is_inside_polygon, line_intersect
from anuga.config import velocity_protection, g
import math

import numpy as num

import inlet

class Inlet_enquiry(inlet.Inlet):
    """Contains information associated with each inlet plus an enquiry point
    """

    def __init__(self, domain, polyline, enquiry_pt,  outward_culvert_vector=None, verbose=False):



        inlet.Inlet.__init__(self, domain, polyline, verbose)


        self.enquiry_pt = enquiry_pt
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
            msg += 'is in an inlet triangle'
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

        return self.get_enquiry_stage() - self.get_enquiry_elevation()


    def get_enquiry_velocity(self):

            depth = self.get_enquiry_depth()
            u = self.get_enquiry_xmom()/(depth + velocity_protection/depth)
            v = self.get_enquiry_ymom()/(depth + velocity_protection/depth)

            return u, v


    def get_enquiry_xvelocity(self):

            depth = self.get_enquiry_depth()
            return self.get_enquiry_xmom()/(depth + velocity_protection/depth)

    def get_enquiry_yvelocity(self):

            depth = self.get_enquiry_depth()
            return self.get_enquiry_ymom()/(depth + velocity_protection/depth)


    def get_enquiry_speed(self):

            u, v = self.get_enquiry_velocity()

            return math.sqrt(u**2 + v**2)


    def get_enquiry_velocity_head(self):

        return 0.5*self.get_enquiry_speed()**2/g


    def get_enquiry_total_energy(self):

        return self.get_enquiry_velocity_head() + self.get_enquiry_stage()


    def get_enquiry_specific_energy(self):

        return self.get_enquiry_velocity_head() + self.get_enquiry_depth()


