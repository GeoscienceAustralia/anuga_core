# To change this template, choose Tools | Templates
# and open the template in the editor.

import anuga.geometry.polygon
from anuga.geometry.polygon import inside_polygon, is_inside_polygon, line_intersect
from anuga.config import velocity_protection, g
import math

import numpy as num
from anuga.structures.inlet import Inlet
import pypar

class Parallel_Inlet(Inlet):
    """Contains information associated with each inlet
    """

    def __init__(self, domain, line, master_proc = 0, procs = None, verbose=False):

        self.domain = domain
        self.line = line
        self.verbose = verbose
        self.master_proc = master_proc #Master processor where global data is gathered

        if procs is None:
            self.procs = [self.master_proc]
        else:
            self.procs = procs

        self.myid = pypar.rank()

        self.compute_triangle_indices()
        self.compute_area()
        self.compute_inlet_length()



    def compute_triangle_indices(self):

        # PETE: Note that all of these are fine but keep in mind that they are local to a domain.

        # Get boundary (in absolute coordinates)
        vertex_coordinates = self.domain.get_full_vertex_coordinates(absolute=True)

        # PETE: Eliminate ghost triangle indices, we can assume that it is in the other inlet

        self.triangle_indices = line_intersect(vertex_coordinates, self.line)

        #print "P%d has %d inlet triangles" %(self.myid, len(self.triangle_indices))

        #print "Triangle Indices:"
        #
        for i in self.triangle_indices:
            assert self.domain.tri_full_flag[i] == 1

        # This effectively does the checks already
        if len(self.triangle_indices) == 0:
            msg = 'Inlet line=%s ' % (self.line)
            msg += 'No triangles intersecting line (Only an enquiry point?)'
            print "WARNING: " + msg
            #raise Exception, msg



    def compute_area(self):

        # Compute inlet area as the sum of areas of triangles identified
        # by line. Must be called after compute_inlet_triangle_indices().
        if len(self.triangle_indices) == 0:
            region = 'Inlet line=%s' % (self.line)
            msg = 'No triangles have been identified in region '
            print "WARNING: " + msg
            #raise Exception, msg

        self.area = 0.0
        for j in self.triangle_indices:
            self.area += self.domain.areas[j]

        # PETE: Do a reduction operation to tally up the areas? Must be asynchronous
        # Can we assume that this will be called roughly at the same time?
        # At this point this calculates the local area

        msg = 'Inlet exchange area has area = %f' % self.area
        assert self.area >= 0.0
        #assert self.area > 0.0


    def compute_inlet_length(self):
        """ Compute the length of the inlet (as
        defined by the input line
        """

        # PETE: This is ok, I think this is independent of the domain?

        point0 = self.line[0]
        point1 = self.line[1]

        #TODO: Go through each point in the line, only count the lenght as the one within the
        #bounding polygon

        self.inlet_length = anuga.geometry.polygon.line_length(self.line)


    def get_inlet_length(self):
        # LOCAL
        msg = "Warning: compute_inlet_length not implemented"
        warnings.warn(msg)
        return self.inlet_length

    def get_line(self):
        # LOCAL
        return self.line

    def get_area(self):
        # LOCAL
        return self.area

    def get_global_area(self):
        # Master processor gathers area from all child processors, and returns value
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
        # LOCAL
        # Must be called after compute_inlet_triangle_indices().
        return self.domain.areas.take(self.triangle_indices)


    def get_stages(self):
        # LOCAL
        # PETE: Do we provide all the stages, is it ok if we just provide the local stages?
        # Are there any dependencies?
        return self.domain.quantities['stage'].centroid_values.take(self.triangle_indices)


    def get_average_stage(self):
        # LOCAL
        return num.sum(self.get_stages()*self.get_areas())/self.area

    def get_global_average_stage(self):
        # LOCAL
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


        return global_stage/global_area

    def get_elevations(self):
        # LOCAL
        return self.domain.quantities['elevation'].centroid_values.take(self.triangle_indices)

    def get_average_elevation(self):
        # LOCAL
        return num.sum(self.get_elevations()*self.get_areas())/self.area


    def get_xmoms(self):
        # LOCAL
        return self.domain.quantities['xmomentum'].centroid_values.take(self.triangle_indices)


    def get_average_xmom(self):
        # LOCAL
        return num.sum(self.get_xmoms()*self.get_areas())/self.area

    def get_global_average_xmom(self):
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


        return global_xmoms/global_area

    def get_ymoms(self):
        # LOCAL
        return self.domain.quantities['ymomentum'].centroid_values.take(self.triangle_indices)


    def get_average_ymom(self):
        # LOCAL
        return num.sum(self.get_ymoms()*self.get_areas())/self.area

    def get_global_average_ymom(self):
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


        return global_ymoms/global_area

    def get_depths(self):
        # LOCAL
        return self.get_stages() - self.get_elevations()


    def get_total_water_volume(self):
        # LOCAL
       return num.sum(self.get_depths()*self.get_areas())

    def get_global_total_water_volume(self):
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
        return self.get_total_water_volume()/self.area

    def get_global_average_depth(self):
        area = self.get_global_area()
        total_water_volume = self.get_global_total_water_volume()

        return total_water_volume / area


    def get_velocities(self):
        #LOCAL
            depths = self.get_depths()
            u = self.get_xmoms()/(depths + velocity_protection/depths)
            v = self.get_ymoms()/(depths + velocity_protection/depths)

            return u, v


    def get_xvelocities(self):
        #LOCAL
            depths = self.get_depths()
            return self.get_xmoms()/(depths + velocity_protection/depths)

    def get_yvelocities(self):
        #LOCAL
            depths = self.get_depths()
            return self.get_ymoms()/(depths + velocity_protection/depths)


    def get_average_speed(self):
        #LOCAL
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


# Set routines

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

    def debug_set_stages_evenly(self,volume):
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

        # find the number of cells which will be filled
        index = num.nonzero(summed_volume<=volume)[0][-1]

        # calculate stage needed to fill chosen cells with given volume of water
        depth = (volume - summed_volume[index])/summed_areas[index]
        new_stage = stages[stages_order[index]]+depth


        #print "Summed Volume = " + str(summed_volume)
        #print "Summed Area = " + str(summed_areas)
        #print "Ordered Stages = " + str(stages[stages_order[:]])
        #print "New Stage = " + str(new_stage) + " Volume = " + str(volume) +  " Summed Volume = " + str(summed_volume[index]) + " Index = " + str(index+1)
        #print "NS = " + str(new_stage) + " SAr = " + str(summed_areas[index]) + " Vol = " + str(volume) + " SVol = " + str(summed_volume[index]) + " I = " + str(index)


        stages[stages_order[0:index+1]] = new_stage
        #stages[stages_order[0:index+1]] = stages[stages_order[index]]+depth

        

        self.set_stages(stages)

    def set_stages_evenly(self,volume):
        """ Distribute volume of water over
        inlet exchange region so that stage is level
        """
        # PETE: THIS must be global and in parallel - this does not appear to set the stage for the part
        # above the volume
        #

        '''
        if pypar.size() == 1:
            self.debug_set_stages_evenly(volume)
            return
        '''
        centroid_coordinates = self.domain.get_full_centroid_coordinates(absolute=True)
        areas = self.get_areas()
        stages = self.get_stages()

        stages_order = stages.argsort()

        # PETE: send stages and areas, apply merging procedure

        s_areas = {}
        s_stages = {}
        s_stages_order = {}
        total_stages = len(stages)

        #for i in stages_order:
        #    print "[%d, %f, %s]" %(self.myid, stages[i], centroid_coordinates[self.triangle_indices[i]])

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
                    #print "Recieved from P%d" %(i)
                    #print str(s_stages[i])
                    total_stages = total_stages + len(s_stages[i])

        else:
            # Send areas, stages, and stages order to master
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
            #sa_debug = []
            #sv_debug = []
            #s_debug = []

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
                        #s_debug.append(current_stage)
                        index = i

                # If first iteration, then only update summed_areas, and current stage
                #print "(%d, %f, %s)" %(index, current_stage, centroid_coordinates[self.triangle_indices[s_stages_order[index][pos[index]]]])

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
                #sa_debug.append(summed_areas)
                #sv_debug.append(summed_volume)
                # Update position of index processor and current stage
                prev_stage = current_stage

            # Calculate new stage
            new_stage = prev_stage + (volume - summed_volume) / summed_areas

            #print "Ordered Stages = " + str(stages[stages_order[:]])
            #print "NS = " + str(new_stage) + " SAr = " + str(summed_areas) + " Vol = " + str(volume) + " SVol = " + str(summed_volume) + " I = " + str(pos[self.myid])

            # Send postion and new stage to all processors
            for i in self.procs:
                if i != self.master_proc:
                    pypar.send(pos[i], i)
                    pypar.send(new_stage, i)

            # Update own depth
            #print "P%d, pos = %d, new_stage = %f" %(self.myid, pos[self.myid], new_stage)
            stages[stages_order[0:pos[self.myid]]] = new_stage
        else:
            pos = pypar.receive(self.master_proc)
            new_stage = pypar.receive(self.master_proc)
            stages[stages_order[0:pos]] = new_stage
            #print "P%d, pos = %d, new_stage = %f" %(self.myid, pos, new_stage)
            #print str(stages)

        self.set_stages(stages)

        stages = self.get_stages()
        stages_order = stages.argsort()

        #for i in stages_order:
            #print "eeee: [%d, %f, %s]" %(self.myid, stages[i], centroid_coordinates[self.triangle_indices[i]])

    def set_depths_evenly(self,volume):
        """ Distribute volume over all exchange
        cells with equal depth of water
        """
        # Is this correct?
        new_depth = self.get_average_depth() + (volume/self.get_area())
        self.set_depths(new_depth)

    def get_master_proc(self):
        return self.master_proc

    def __parallel_safe(self):

        return True

    def statistics(self):
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
                message += '======> inlet triangle indices and centres at P%d\n' %(proc)
                message += '%s' % tri_indices[proc]
                message += '\n'
                
                message += '%s' % self.domain.get_centroid_coordinates()[tri_indices[proc]]
                message += '\n'
                
            message += 'line\n'
            message += '%s' % self.line
            message += '\n'

        return message

__author__="pete"
__date__ ="$16/08/2011 6:49:42 PM$"

if __name__ == "__main__":
    print "Hello World"
