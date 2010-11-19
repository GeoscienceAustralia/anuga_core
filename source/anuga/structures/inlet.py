import anuga.geometry.polygon
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
        self.compute_inlet_length()



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


    def compute_inlet_length(self):
        """ Compute the length of the inlet (as
        defined by the input line
        """

        point0 = self.line[0]
        point1 = self.line[1]

        self.inlet_length = anuga.geometry.polygon.line_length(self.line)


    def get_inlet_length(self):

        return self.inlet_length

    def get_line(self):

        return self.line
        
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
    

    def get_depths(self):
    
        return self.get_stages() - self.get_elevations()
    
    
    def get_total_water_volume(self):
       
       return num.sum(self.get_depths()*self.get_areas())
  

    def get_average_depth(self):
    
        return self.get_total_water_volume()/self.area
        
        
    def get_velocities(self):
        
            depths = self.get_depths()
            u = self.get_xmoms()/(depths + velocity_protection/depths)
            v = self.get_ymoms()/(depths + velocity_protection/depths)
            
            return u, v


    def get_xvelocities(self):

            depths = self.get_depths()
            return self.get_xmoms()/(depths + velocity_protection/depths)

    def get_yvelocities(self):

            depths = self.get_depths()
            return self.get_ymoms()/(depths + velocity_protection/depths)
            
            
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
        
        return self.get_average_velocity_head() + self.get_average_depth()



    def set_depths(self,depth):

        self.domain.quantities['stage'].centroid_values.put(self.triangle_indices, self.get_elevations() + depth)


    def set_stages(self,stage):

        self.domain.quantities['stage'].centroid_values.put(self.triangle_indices, stage)


    def set_xmoms(self,xmom):

        self.domain.quantities['xmomentum'].centroid_values.put(self.triangle_indices, xmom)


    def set_ymoms(self,ymom):

        self.domain.quantities['ymomentum'].centroid_values.put(self.triangle_indices, ymom)


    def set_elevations(self,elevation):

        self.domain.quantities['elevation'].centroid_values.put(self.triangle_indices, elevation)

    def set_stages_evenly(self,volume):
        """ Distribute volume of water over
        inlet exchange region so that stage is level
        """

        areas = self.get_areas()
        stages = self.get_stages()

        stages_order = stages.argsort()

        summed_areas = num.zeros_like(areas)
        summed_volume = num.zeros_like(areas)
        #diff_stage = num.zeros_like(areas)

        for i,a in enumerate(areas[stages_order]):
            #print i,a, stages[stages_order[i]]
            if i == 0:
                summed_areas[i] = a
                summed_volume[i] = 0.0
                #diff_stage[i] = 0.0
            else:
                summed_areas[i] = summed_areas[i-1] + a

                summed_volume[i] = summed_volume[i-1] + summed_areas[i-1]*\
                    (stages[stages_order[i]] - stages[stages_order[i-1]])

                #diff_stage[i] = stages[stages_order[i]] - stages[stages_order[i-1]]

        #print 'amount ',amount
        #print summed_areas
        #print summed_amount
        #print diff_stage

            
        #index = len(summed_amount)
        for i,a in enumerate(summed_volume):
            #print 'a ',a
            if volume > a :
                index = i

        #print index

        #print stages_order
        #print stages
        
        depth = (volume - summed_volume[index])/summed_areas[index]
		
        stages[stages_order[0:index+1]] = stages[stages_order[index]]+depth

        self.set_stages(stages) 

    def set_depths_evenly(self,volume):
        """ Distribute volume over all exchange
        cells with equal depth of water
        """
	    
        new_depth = self.get_average_depth() + (volume/self.get_area())
        self.set_depths(new_depth)

