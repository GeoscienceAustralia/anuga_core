from anuga.geometry.polygon import inside_polygon, is_inside_polygon
from anuga.config import velocity_protection 

import numpy as num

class Inlet:
    """Contains information associated with each inlet
    """

    def __init__(self, domain, polygon, enquiry_point, inlet_vector):

        self.domain = domain
        self.domain_bounding_polygon = self.domain.get_boundary_polygon()
        self.polygon = polygon
        self.enquiry_point = enquiry_point
        self.inlet_vector = inlet_vector

        # FIXME (SR) Using get_triangle_containing_point which needs to be sped up
        self.enquiry_index = self.domain.get_triangle_containing_point(self.enquiry_point)

        self.compute_triangle_indices()
        self.compute_area()


    def compute_triangle_indices(self):

        # Get boundary (in absolute coordinates)
        bounding_polygon = self.domain_bounding_polygon
        domain_centroids = self.domain.get_centroid_coordinates(absolute=True)

        self.inlet_polygon = self.polygon

        # Check that polygon lies within the mesh.
        for point in self.inlet_polygon:
                msg = 'Point %s in polygon for forcing term' %  str(point)
                msg += ' did not fall within the domain boundary.'
                assert is_inside_polygon(point, bounding_polygon), msg

        self.triangle_indices = inside_polygon(domain_centroids, self.inlet_polygon)

        if len(self.triangle_indices) == 0:
            region = 'Inlet polygon=%s' % (self.inlet_polygon)
            msg = 'No triangles have been identified in '
            msg += 'specified region: %s' % region
            raise Exception, msg


    def compute_area(self):
        
        # Compute inlet area as the sum of areas of triangles identified
        # by polygon. Must be called after compute_inlet_triangle_indices().
        if len(self.triangle_indices) == 0:
            region = 'Inlet polygon=%s' % (self.inlet_polygon)
            msg = 'No triangles have been identified in '
            msg += 'specified region: %s' % region
            raise Exception, msg
        
        self.area = 0.0
        for j in self.triangle_indices:
            self.area += self.domain.areas[j]

        msg = 'Inlet exchange area has area = %f' % self.area
        assert self.area > 0.0


    def get_areas(self):
        
        # Must be called after compute_inlet_triangle_indices().
        return self.domain.areas.take(self.triangle_indices)
    
        
    def get_stages(self):
        
        return self.domain.quantities['stage'].centroid_values.take(self.triangle_indices)
        
    def get_average_stage(self):
        
        return num.sum(self.get_stages())/self.triangle_indices.size
       
        
    def get_elevations(self):    
        
        return self.domain.quantities['elevation'].centroid_values.take(self.triangle_indices)
        
    def get_average_elevation(self):
        
        return num.sum(self.get_elevations())/self.triangle_indices.size
    
    
    def get_xmoms(self):
    
        return self.domain.quantities['xmomentum'].centroid_values.take(self.triangle_indices)
        
    def get_average_xmom(self):
        
        return num.sum(self.get_xmoms())/self.triangle_indices.size
        
    
    def get_ymoms(self):
        
        return self.domain.quantities['ymomentum'].centroid_values.take(self.triangle_indices)
 
 
    def get_average_ymom(self):
        
        return num.sum(self.get_ymoms())/self.triangle_indices.size
 
    
    def get_heights(self):
    
        return self.get_stages() - self.get_elevations()
    
    
    def get_total_water_volume(self):
       
       return num.sum(self.get_heights*self.get_areas())
    
    
    def get_average_height(self):
    
        return self.get_total_water_volume()/self.area
        
        
    def get_velocities(self):
        
            depths = self.get_stages() - self.get_elevations()
            u = self.get_xmoms()/(depths + velocity_protection/depths)
            v = self.get_ymoms()/(depths + velocity_protection/depths)
            
            return u, v
            
            
    def get_average_velocities(self):
 
            depths = self.get_stages() - self.get_elevations()
            u = self.get_xmoms()/(depths + velocity_protection/depths)
            v = self.get_ymoms()/(depths + velocity_protection/depths)
            
            return num.sum(u)/self.triangle_indices.size, num.sum(v)/self.triangle_indices.size
