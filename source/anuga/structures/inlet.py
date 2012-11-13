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
        self.line = num.asarray(line, dtype=num.float64)
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

        # accumulate areas of cells ordered by stage
        summed_areas = num.cumsum(areas[stages_order])
        
        # accumulate the volume need to fill cells
        summed_volume = num.zeros_like(areas)       
        summed_volume[1:] = num.cumsum(summed_areas[:-1]*num.diff(stages[stages_order]))


        if volume>= 0.0 :
            # find the number of cells which will be filled
            #print summed_volume
            #print summed_volume<=volume

            #Test = summed_volume<=volume

            #print num.any(Test)
            #print num.nonzero(summed_volume<=volume)
            #print num.nonzero(summed_volume<=volume)[0]

            index = num.nonzero(summed_volume<=volume)[0][-1]

            # calculate stage needed to fill chosen cells with given volume of water
            depth = (volume - summed_volume[index])/summed_areas[index]
            stages[stages_order[0:index+1]] = stages[stages_order[index]]+depth

        else:
            if summed_volume[-1]>= -volume :
                depth = (summed_volume[-1] + volume)/summed_areas[-1]
                stages[:] = depth
            else:
                #assert summed_volume[-1] >= -volume
                import warnings
                warnings.warn('summed_volume < volume to be extracted')
                stages[:] = 0.0

        self.set_stages(stages)




    def set_stages_evenly_new(self,volume):
        """ Distribute volume of water over
        inlet exchange region so that stage is level
        """
        # WARNING: requires synchronization, must be called by all procs associated
        # with this inlet

        centroid_coordinates = self.domain.get_full_centroid_coordinates(absolute=True)
        areas = self.get_areas()
        stages = self.get_stages()

        stages_order = stages.argsort()

        # PETE: send stages and areas, apply merging procedure

        total_stages = len(stages)


        s_areas = areas
        s_stages = stages
        s_stages_order = stages_order


        pos = 0
        summed_volume = 0.
        summed_areas = 0.
        prev_stage = 0.
        num_stages = 0.
        first = True

        for i in self.procs:
            pos[i] = 0

        while num_stages < total_stages:
            # Determine current minimum stage of all the processors in s_stages
            num_stages = num_stages + 1
            current_stage = num.finfo(num.float32).max
            index = -1


            if pos >= len(s_stages):
                continue

            if s_stages[s_stages_order[pos]] < current_stage:
                current_stage = s_stages[i][s_stages_order[i][pos[i]]]
                index = i

            # If first iteration, then only update summed_areas, position, and prev|current stage

            if first:
                first = False
                summed_areas = s_areas[index][s_stages_order[index][pos[index]]]
                pos[index] = pos[index] + 1
                prev_stage = current_stage
                continue

            assert index >= 0, "Index out of bounds"

            # Update summed volume and summed areas
            tmp_volume = summed_volume + (summed_areas * (current_stage - prev_stage))

            # Terminate if volume exceeded
            if tmp_volume >= volume:
                break

            summed_areas = summed_areas + s_areas[index][s_stages_order[index][pos[index]]]
            pos[index] = pos[index] + 1
            summed_volume = tmp_volume

        # Update position of index processor and current stage
        prev_stage = current_stage

        # Calculate new stage
        new_stage = prev_stage + (volume - summed_volume) / summed_areas


        # Update own depth
        stages[stages_order[0:pos]] = new_stage



        self.set_stages(stages)

        stages = self.get_stages()
        stages_order = stages.argsort()


    def set_depths_evenly(self,volume):
        """ Distribute volume over all exchange
        cells with equal depth of water
        """
	    
        new_depth = self.get_average_depth() + (volume/self.get_area())
        self.set_depths(new_depth)

