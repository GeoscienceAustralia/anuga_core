
from anuga.geometry.polygon import inside_polygon, is_inside_polygon, line_intersect
from anuga.config import velocity_protection, g
import math

import numpy as num

from . import parallel_inlet

class Parallel_Inlet_enquiry(parallel_inlet.Parallel_Inlet):
    """Contains information associated with each inlet plus an enquiry point
    """

    """
    master_proc - index of the processor which coordinates all processors 
    associated with this inlet operator.
    procs - list of all processors associated with this inlet operator
    enquiry_proc - processor containing inlet enquiry point
    """

    def __init__(self, domain, polyline, enquiry_pt,
                 invert_elevation = None,
                 outward_culvert_vector=None,
                 master_proc = 0,
                 procs = None,
                 enquiry_proc = -1,
                 verbose=False):

   
        parallel_inlet.Parallel_Inlet.__init__(self, domain, polyline,
                                                master_proc = master_proc, procs = procs, verbose=verbose)

        from anuga.utilities import parallel_abstraction as pypar

        self.enquiry_pt = enquiry_pt
        self.invert_elevation = invert_elevation
        self.outward_culvert_vector = outward_culvert_vector
        self.master_proc = master_proc

        if procs is None:
            self.procs = [self.master_proc]
        else:
            self.procs = procs

        self.enquiry_proc = enquiry_proc # Processor of domain containing enquiry point
        self.compute_enquiry_index()


    def compute_enquiry_index(self):

        # Get boundary (in absolute coordinates)
        vertex_coordinates = self.domain.get_full_vertex_coordinates(absolute=True)

        point = self.enquiry_pt
        has_enq_point = False

        try:
            k = self.domain.get_triangle_containing_point(point)
                
            if self.domain.tri_full_flag[k] == 1:
                has_enq_point = True
            else:
                has_enq_point = False
        except:
            has_enq_point = False

        if has_enq_point:
            self.enquiry_index = self.domain.get_triangle_containing_point(self.enquiry_pt)
   
            if self.enquiry_index in self.triangle_indices:
                msg = 'Enquiry point %s' % (self.enquiry_pt)
                msg += 'is in an inlet triangle'
                import warnings
                warnings.warn(msg)


            if self.enquiry_proc >= 0: 
                assert self.enquiry_proc == self.myid, "Specified enquiry proc does not match actual enquiry proc"
            self.enquiry_proc = self.myid
            assert self.enquiry_index >= 0, "Enquiry point inside polygon, but no triangle index found"
        else:
            self.enquiry_index = -1


    def get_enquiry_position(self):
        # WARNING: Must be called by processor containing inlet enquiry point to have effect

        if self.enquiry_index >= 0:
            return self.domain.get_centroid_coordinates(absolute=True)[self.enquiry_index]
        else:
            return None

    def get_enquiry_stage(self):
        # WARNING: Must be called by processor containing inlet enquiry point to have effect

        if self.enquiry_index >= 0:
            return self.domain.quantities['stage'].centroid_values[self.enquiry_index]
        else:
            return None

        return None

    def get_enquiry_xmom(self):
        # WARNING: Must be called by processor containing inlet enquiry point to have effect

        if self.enquiry_index >= 0:
            return self.domain.quantities['xmomentum'].centroid_values[self.enquiry_index]
        else:
            return None

    def get_enquiry_ymom(self):
        # WARNING: Must be called by processor containing inlet enquiry point to have effect

        if self.enquiry_index >= 0:
            return self.domain.quantities['ymomentum'].centroid_values[self.enquiry_index]
        else:
            return None

    def get_enquiry_elevation(self):
        # WARNING: Must be called by processor containing inlet enquiry point to have effect

        if self.enquiry_index >= 0:
            return self.domain.quantities['elevation'].centroid_values[self.enquiry_index]
        else:
            return None


    def get_enquiry_depth(self):
        # WARNING: Must be called by processor containing inlet enquiry point to have effect

        if self.enquiry_index >= 0:
            return max(self.get_enquiry_stage() - self.get_enquiry_invert_elevation(), 0.0)
        else:
            return None

    def get_enquiry_water_depth(self):
        # WARNING: Must be called by processor containing inlet enquiry point to have effect

        if self.enquiry_index >= 0:
            return self.get_enquiry_stage() - self.get_enquiry_elevation()
        else:
            return None

    def get_enquiry_invert_elevation(self):
        # WARNING: Must be called by processor containing inlet enquiry point to have effect

        if self.enquiry_index >= 0:
            if  self.invert_elevation is None:
                return self.get_enquiry_elevation()
            else:
                return self.invert_elevation
        else:
            return None

    def get_enquiry_velocity(self):
        # WARNING: Must be called by processor containing inlet enquiry point to have effect

        if self.enquiry_index >= 0:
            depth = self.get_enquiry_water_depth()
            u = depth*self.get_enquiry_xmom()/(depth**2 + velocity_protection)
            v = depth*self.get_enquiry_ymom()/(depth**2 + velocity_protection)

            return u, v
        else:
            return None


    def get_enquiry_xvelocity(self):
        # WARNING: Must be called by processor containing inlet enquiry point to have effect

        if self.enquiry_index >= 0:
            depth = self.get_enquiry_water_depth()
            return depth*self.get_enquiry_xmom()/(depth**2 + velocity_protection)
        else:
            return None

    def get_enquiry_yvelocity(self):
        # WARNING: Must be called by processor containing inlet enquiry point to have effect

        if self.enquiry_index >= 0:
            depth = self.get_enquiry_water_depth()
            return depth*self.get_enquiry_ymom()/(depth**2 + velocity_protection)
        else:
            return None

    def get_enquiry_speed(self):
        # WARNING: Must be called by processor containing inlet enquiry point to have effect

        if self.enquiry_index >= 0:
            u, v = self.get_enquiry_velocity()

            return math.sqrt(u**2 + v**2)
        else:
            return None


    def get_enquiry_velocity_head(self):
        # WARNING: Must be called by processor containing inlet enquiry point to have effect

        if self.enquiry_index >= 0:
            if self.domain.use_new_velocity_head:
                u, v   = self.get_enquiry_velocity()
                n1, n2 = self.outward_culvert_vector
                normal_speed = min(u*n1 + v*n2, 0.0)

                velocity_head = 0.5*normal_speed**2/g
            else:
                velocity_head = 0.5*self.get_enquiry_speed()**2/g

            return velocity_head
        else:
            return None


    def get_enquiry_total_energy(self):
        # WARNING: Must be called by processor containing inlet enquiry point to have effect

        if self.enquiry_index >= 0:
            return self.get_enquiry_velocity_head() + self.get_enquiry_stage()
        else:
            return None


    def get_enquiry_specific_energy(self):
        # WARNING: Must be called by processor containing inlet enquiry point to have effect

        if self.enquiry_index >= 0:
            return self.get_enquiry_velocity_head() + self.get_enquiry_depth()
        else:
            return None


    def get_master_proc(self):
        return self.master_proc


    def get_enquiry_proc(self):
        return self.enquiry_proc
