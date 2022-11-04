

import anuga.geometry.polygon
from anuga.geometry.polygon import inside_polygon, is_inside_polygon, line_intersect
from anuga.config import velocity_protection, g
import math

import numpy as num
from anuga.structures.inlet import Inlet
import warnings
from anuga import Region

class Parallel_Inlet(Inlet):
    """Contains information associated with each inlet
    """

    """
    Parallel inlet:

    master_proc - coordinates all processors associated with this inlet
    usually the processors with domains which contains parts of this inlet.

    procs - is the list of all processors associated with this inlet.

    (We assume that the above arguments are determined correctly by the parallel_operator_factory)
    """

    def __init__(self, domain, poly, master_proc = 0, procs = None, verbose=False):

        self.domain = domain
        self.verbose = verbose

        # poly can be either a line, polygon or a regions
        if isinstance(poly,Region):
            self.region = poly
        else:
            self.region = Region(domain,poly=poly,expand_polygon=True)

        if self.region.get_type() == 'line':
            self.line = True
        else:
            self.line = False

        self.master_proc = master_proc

        if procs is None:
            self.procs = [self.master_proc]
        else:
            self.procs = procs

        from anuga.utilities import parallel_abstraction as pypar
        self.myid = pypar.rank()

        self.triangle_indices = self.region.get_indices(full_only=True)
        self.compute_area()
        #self.compute_inlet_length()

    def compute_area(self):

        # Compute inlet area as the sum of areas of triangles identified
        # by line. Must be called after compute_inlet_triangle_indices().
        if len(self.triangle_indices) == 0:
            region = 'Inlet line=%s' % (self.line)
            msg = 'No triangles have been identified in region '
            print("WARNING: " + msg)

        self.area = 0.0
        for j in self.triangle_indices:
            self.area += self.domain.areas[j]

        msg = 'Inlet exchange area has area = %f' % self.area
        assert self.area >= 0.0

    def compute_inlet_length(self):
        """ Compute the length of the inlet within this domain (as
        defined by the input line
        """

        point0 = self.line[0]
        point1 = self.line[1]

        self.inlet_length = anuga.geometry.polygon.line_length(self.line)


    def get_inlet_length(self):
        # LOCAL
        msg = "Warning: compute_inlet_length not implemented"
        warnings.warn(msg)
        return self.inlet_length

    def get_line(self):
        return self.line

    def get_area(self):
        # LOCAL
        return self.area

    def get_global_area(self):
        # GLOBAL: Master processor gathers area from all child processors, and returns value

        # WARNING: requires synchronization, must be called by all procs associated
        # with this inlet

        from anuga.utilities import parallel_abstraction as pypar
        local_area = self.area
        area = local_area

        if self.myid == self.master_proc:

            for i in self.procs:
                if i == self.master_proc: continue

                val = pypar.receive(i)
                area = area + val
        else:
            pypar.send(area, self.master_proc)

        return area


    def get_areas(self):
        # Must be called after compute_inlet_triangle_indices().
        # LOCAL

        return self.domain.areas.take(self.triangle_indices)


    def get_stages(self):
        # LOCAL

        return self.domain.quantities['stage'].centroid_values.take(self.triangle_indices)


    def get_average_stage(self):
        # LOCAL

        return num.sum(self.get_stages()*self.get_areas())/self.area

    def get_global_average_stage(self):
        # GLOBAL: Master processor gathers stages from all child processors, and returns average

        # WARNING: requires synchronization, must be called by all procs associated
        # with this inlet

        from anuga.utilities import parallel_abstraction as pypar
        local_stage = num.sum(self.get_stages()*self.get_areas())
        global_area = self.get_global_area()



        global_stage = local_stage

        if self.myid == self.master_proc:
            for i in self.procs:
                if i == self.master_proc: continue

                val = pypar.receive(i)
                global_stage = global_stage + val
        else:
            pypar.send(local_stage, self.master_proc)


        if global_area > 0.0:
            return global_stage/global_area
        else:
            return 0.0

    def get_elevations(self):
        # LOCAL
        return self.domain.quantities['elevation'].centroid_values.take(self.triangle_indices)

    def get_average_elevation(self):
        # LOCAL

        if self.area > 0:
            return num.sum(self.get_elevations()*self.get_areas())/self.area
        else:
            return 0.0

    def get_global_average_elevation(self):
        # GLOBAL: Master processor gathers elevations from all child processors, and returns average

        # WARNING: requires synchronization, must be called by all procs associated
        # with this inlet

        from anuga.utilities import parallel_abstraction as pypar
        local_elevation = num.sum(self.get_elevations()*self.get_areas())
        global_area = self.get_global_area()



        global_elevation = local_elevation

        if self.myid == self.master_proc:
            for i in self.procs:
                if i == self.master_proc: continue

                val = pypar.receive(i)
                global_elevation = global_elevation + val
        else:
            pypar.send(local_elevation, self.master_proc)


        if global_area > 0.0:
            return global_elevation/global_area
        else:
            return 0.0

    def get_xmoms(self):
        # LOCAL
        return self.domain.quantities['xmomentum'].centroid_values.take(self.triangle_indices)


    def get_average_xmom(self):
        # LOCAL

        if self.area > 0:
            return num.sum(self.get_xmoms()*self.get_areas())/self.area
        else:
            return 0.0

    def get_global_average_xmom(self):
        # GLOBAL: master proc gathers all xmom values and returns average
        # WARNING: requires synchronization, must be called by all procs associated
        # with this inlet

        from anuga.utilities import parallel_abstraction as pypar
        global_area = self.get_global_area()
        local_xmoms = num.sum(self.get_xmoms()*self.get_areas())
        global_xmoms = local_xmoms

        if self.myid == self.master_proc:
            for i in self.procs:
                if i == self.master_proc: continue

                val = pypar.receive(i)
                global_xmoms = global_xmoms + val
        else:
            pypar.send(local_xmoms, self.master_proc)


        if global_area > 0.0:
            return global_xmoms/global_area
        else:
            return 0.0


    def get_ymoms(self):
        # LOCAL
        return self.domain.quantities['ymomentum'].centroid_values.take(self.triangle_indices)


    def get_average_ymom(self):
        # LOCAL
        return num.sum(self.get_ymoms()*self.get_areas())/self.area

    def get_global_average_ymom(self):
        # GLOBAL: master proc gathers all ymom values and returns average
        # WARNING: requires synchronization, must be called by all procs associated
        # with this inlet

        from anuga.utilities import parallel_abstraction as pypar
        global_area = self.get_global_area()
        local_ymoms = num.sum(self.get_ymoms()*self.get_areas())
        global_ymoms = local_ymoms

        if self.myid == self.master_proc:
            for i in self.procs:
                if i == self.master_proc: continue

                val = pypar.receive(i)
                global_ymoms = global_ymoms + val
        else:
            pypar.send(local_ymoms, self.master_proc)


        if global_area > 0.0:
            return global_ymoms/global_area
        else:
            return 0.0

    def get_depths(self):
        # LOCAL
        return self.get_stages() - self.get_elevations()


    def get_total_water_volume(self):
        # LOCAL
       return num.sum(self.get_depths()*self.get_areas())

    def get_global_total_water_volume(self):
        # GLOBAL: master proc gathers total water volumes from each proc and returns average
        # WARNING: requires synchronization, must be called by all procs associated
        # with this inlet

        from anuga.utilities import parallel_abstraction as pypar
        local_volume = num.sum(self.get_depths()*self.get_areas())
        volume = local_volume

        if self.myid == self.master_proc:

            for i in self.procs:
                if i == self.master_proc: continue

                val = pypar.receive(i)
                volume = volume + val
        else:
            pypar.send(volume, self.master_proc)

        return volume

    def get_average_depth(self):
        # LOCAL

        if self.area > 0.0:
            return self.get_total_water_volume()/self.area
        else:
            return 0.0

    def get_global_average_depth(self):
        # GLOBAL: master proc gathers all depth values and returns average
        # WARNING: requires synchronization, must be called by all procs associated
        # with this inlet

        area = self.get_global_area()
        total_water_volume = self.get_global_total_water_volume()


        if area > 0.0:
            return total_water_volume/ area
        else:
            return 0.0


    def get_velocities(self):
        #LOCAL
        depths = self.get_depths()
        u = depths*self.get_xmoms()/(depths**2 + velocity_protection)
        v = depths*self.get_ymoms()/(depths**2 + velocity_protection)

        return u, v


    def get_xvelocities(self):
        #LOCAL
        depths = self.get_depths()
        return depths*self.get_xmoms()/(depths**2 + velocity_protection)

    def get_yvelocities(self):
        #LOCAL
        depths = self.get_depths()
        return depths*self.get_ymoms()/(depths**2 + velocity_protection)


    def get_average_speed(self):
        #LOCAL
        u, v = self.get_velocities()

        average_u = num.sum(u*self.get_areas())/self.area
        average_v = num.sum(v*self.get_areas())/self.area

        return math.sqrt(average_u**2 + average_v**2)


    def get_average_velocity_head(self):
        #LOCAL
        return 0.5*self.get_average_speed()**2/g


    def get_average_total_energy(self):
        #LOCAL
        return self.get_average_velocity_head() + self.get_average_stage()


    def get_average_specific_energy(self):
        #LOCAL
        return self.get_average_velocity_head() + self.get_average_depth()


# Set routines (ALL LOCAL)

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
        # WARNING: requires synchronization, must be called by all procs associated
        # with this inlet

        from anuga.utilities import parallel_abstraction as pypar
        centroid_coordinates = self.domain.get_full_centroid_coordinates(absolute=True)
        areas = self.get_areas()
        stages = self.get_stages()

        stages_order = stages.argsort()

        # PETE: send stages and areas, apply merging procedure

        s_areas = {}
        s_stages = {}
        s_stages_order = {}
        total_stages = len(stages)

        if self.myid == self.master_proc:
            s_areas[self.myid] = areas
            s_stages[self.myid] = stages
            s_stages_order[self.myid] = stages_order

            # Recieve areas, stages, and stages order
            for i in self.procs:
                if i != self.master_proc:
                    s_areas[i] = pypar.receive(i)
                    s_stages[i] = pypar.receive(i)
                    s_stages_order[i] = pypar.receive(i)
                    total_stages = total_stages + len(s_stages[i])

        else:
            # Send areas, stages, and stages order to master proc of inlet
            pypar.send(areas, self.master_proc)
            pypar.send(stages, self.master_proc)
            pypar.send(stages_order, self.master_proc)

        # merge sorted stage order
        if self.myid == self.master_proc:
            pos = {}
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

                for i in self.procs:
                    if pos[i] >= len(s_stages[i]):
                        continue

                    if s_stages[i][s_stages_order[i][pos[i]]] < current_stage:
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
            new_stage = prev_stage + (volume - summed_volume)/ summed_areas

            # Send postion and new stage to all processors
            for i in self.procs:
                if i != self.master_proc:
                    pypar.send(pos[i], i)
                    pypar.send(new_stage, i)

            # Update own depth
            stages[stages_order[0:pos[self.myid]]] = new_stage
        else:
            pos = pypar.receive(self.master_proc)
            new_stage = pypar.receive(self.master_proc)
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

    def get_master_proc(self):
        return self.master_proc

    def parallel_safe(self):

        return True

    def statistics(self):
        # WARNING: requires synchronization, must be called by all procs associated
        # with this inlet

        from anuga.utilities import parallel_abstraction as pypar

        message = ''

        tri_indices = {}

        if self.myid == self.master_proc:
            tri_indices[self.myid] = self.triangle_indices

            for proc in self.procs:
                if proc == self.master_proc: continue

                tri_indices[proc] = pypar.receive(proc)

        else:
            pypar.send(self.triangle_indices, self.master_proc)


        if self.myid == self.master_proc:
            message += '=====================================\n'
            message +=  'Inlet\n'
            message += '=====================================\n'

            for proc in self.procs:
                message += '======> inlet triangle indices and centres and elevation at P%d\n' %(proc)
                message += '%s' % tri_indices[proc]
                message += '\n'

                message += '%s' % self.domain.get_centroid_coordinates()[tri_indices[proc]]
                message += '\n'

                elev = self.domain.quantities['elevation'].centroid_values[tri_indices[proc]]
                message += '%s' % elev
                message += '\n'

                try:
                    elevation_difference = elev.max() - elev.min()
                except ValueError:
                    elevation_difference = 0.0

                if not num.allclose(elevation_difference, 0.):
                    message += 'Elevation range of ' + str(elevation_difference)
                    message += 'Warning: Non-constant inlet elevation can lead to well-balancing problems'

                try:
                    # If the inlet does not have an enquiry point this will
                    # fail gracefully
                    message += '\n'
                    message += 'Enquiry point:'
                    message += '%s' % self.domain.get_centroid_coordinates()[self.enquiry_index]
                    message += '\n'
                    message += 'Enquiry Index:'
                    message += '%s' % self.enquiry_index
                    message += '\n'
                except:
                    pass


            message += 'line\n'
            message += '%s' % self.line
            message += '\n'

        return message

__author__="pete"
__date__ ="$16/08/2011 6:49:42 PM$"

if __name__ == "__main__":
    print("Hello World")
