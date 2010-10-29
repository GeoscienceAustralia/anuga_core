from anuga.geometry.polygon import inside_polygon, is_inside_polygon, line_intersect
from anuga.config import velocity_protection, g
import math

import numpy as num

class Inlet:
    """Contains information associated with each inlet
    """

    def __init__(self, domain, line, verbose=False):

        self.domain = domain
        self.domain_bounding_polygon = self.domain.get_boundary_polygon()
        self.line = line
        self.verbose = verbose

        self.compute_triangle_indices()
        self.compute_area()


    def compute_triangle_indices(self):

        # Get boundary (in absolute coordinates)
        bounding_polygon = self.domain_bounding_polygon
        domain_centroids = self.domain.get_centroid_coordinates(absolute=True)
        vertex_coordinates = self.domain.get_vertex_coordinates(absolute=True)

        # Check that line lies within the mesh.
        for point in self.line:  
                msg = 'Point %s ' %  str(point)
                msg += ' did not fall within the domain boundary.'
                assert is_inside_polygon(point, bounding_polygon), msg
                


        self.triangle_indices = line_intersect(vertex_coordinates, self.line)

        if len(self.triangle_indices) == 0:
            msg = 'Inlet line=%s ' % (self.line)
            msg += 'No triangles intersecting line '
            raise Exception, msg



    def compute_area(self):
        
        # Compute inlet area as the sum of areas of triangles identified
        # by line. Must be called after compute_inlet_triangle_indices().
        if len(self.triangle_indices) == 0:
            region = 'Inlet line=%s' % (self.inlet_line)
            msg = 'No triangles have been identified in region '
            raise Exception, msg
        
        self.area = 0.0
        for j in self.triangle_indices:
            self.area += self.domain.areas[j]

        msg = 'Inlet exchange area has area = %f' % self.area
        assert self.area > 0.0


    def get_area(self):

        return self.area

    
    def get_areas(self):
        
        # Must be called after compute_inlet_triangle_indices().
        return self.domain.areas.take(self.triangle_indices)
    
        
    def get_stages(self):
        
        return self.domain.quantities['stage'].centroid_values.take(self.triangle_indices)
        
        
    def get_average_stage(self):

        return num.sum(self.get_stages()*self.get_areas())/self.area
        
    def get_elevations(self):    
        
        return self.domain.quantities['elevation'].centroid_values.take(self.triangle_indices)
        
    def get_average_elevation(self):

        return num.sum(self.get_elevations()*self.get_areas())/self.area
    
    
    def get_xmoms(self):
    
        return self.domain.quantities['xmomentum'].centroid_values.take(self.triangle_indices)
        
        
    def get_average_xmom(self):

        return num.sum(self.get_xmoms()*self.get_areas())/self.area
        
    
    def get_ymoms(self):
        
        return self.domain.quantities['ymomentum'].centroid_values.take(self.triangle_indices)
 
 
    def get_average_ymom(self):
        
        return num.sum(self.get_ymoms()*self.get_areas())/self.area
    

    def get_heights(self):
    
        return self.get_stages() - self.get_elevations()
    
    
    def get_total_water_volume(self):
       
       return num.sum(self.get_heights()*self.get_areas())
  

    def get_average_height(self):
    
        return self.get_total_water_volume()/self.area
        
        
    def get_velocities(self):
        
            heights = self.get_heights()
            u = self.get_xmoms()/(heights + velocity_protection/heights)
            v = self.get_ymoms()/(heights + velocity_protection/heights)
            
            return u, v


    def get_xvelocities(self):

            heights = self.get_heights()
            return self.get_xmoms()/(heights + velocity_protection/heights)

    def get_yvelocities(self):

            heights = self.get_heights()
            return self.get_ymoms()/(heights + velocity_protection/heights)
            
            
    def get_average_speed(self):
 
            u, v = self.get_velocities()
            
            average_u = num.sum(u*self.get_areas())/self.area
            average_v = num.sum(v*self.get_areas())/self.area
            
            return math.sqrt(average_u**2 + average_v**2)


    def get_average_velocity_head(self):

        return 0.5*self.get_average_speed()**2/g


    def get_average_total_energy(self):
        
        return self.get_average_velocity_head() + self.get_average_stage()
        
    
    def get_average_specific_energy(self):
        
        return self.get_average_velocity_head() + self.get_average_height()



    def set_heights(self,height):

        self.domain.quantities['stage'].centroid_values.put(self.triangle_indices, self.get_elevations() + height)


    def set_stages(self,stage):

        self.domain.quantities['stage'].centroid_values.put(self.triangle_indices, stage)


    def set_xmoms(self,xmom):

        self.domain.quantities['xmomentum'].centroid_values.put(self.triangle_indices, xmom)


    def set_ymoms(self,ymom):

        self.domain.quantities['ymomentum'].centroid_values.put(self.triangle_indices, ymom)


    def set_elevations(self,elevation):

        self.domain.quantities['elevation'].centroid_values.put(self.triangle_indices, elevation)
