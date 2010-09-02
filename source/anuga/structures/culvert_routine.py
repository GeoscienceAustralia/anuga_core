from anuga.config import velocity_protection
from anuga.utilities.numerical_tools import safe_acos as acos

from math import pi, sqrt, sin, cos
from anuga.config import g

class Culvert_routine:
    """Collection of culvert routines for use with Culvert_operator

    This module holds various routines to determine FLOW through CULVERTS and SIMPLE BRIDGES

    Usage:

    NOTE:
     Inlet control:  self.delta_total_energy > self.inflow.get_average_specific_energy()
     Outlet control: self.delta_total_energy < self.inflow.get_average_specific_energy()
     where total energy is (w + 0.5*v^2/g) and
     specific energy is (h + 0.5*v^2/g)
    """

    def __init__(self, culvert, manning=0.0):
        
        self.inlets = culvert.inlets
        
       
        self.culvert_length = culvert.get_culvert_length()
        self.culvert_width = culvert.get_culvert_width()
        self.culvert_height = culvert.get_culvert_height()
        self.sum_loss = 0.0
        self.max_velocity = 10.0
        self.manning = manning
        self.log_filename = None

        self.use_velocity_head = True
        

        
        self.determine_inflow()

        #delta_z = self.self.inflow.get_average_elevation() - self.self.outflow.get_average_elevation()
        #culvert_slope = delta_z/self.culvert.get_self.culvert_length()

        # Determine controlling energy (driving head) for culvert
        #if self.self.inflow.get_average_specific_energy() > self.self.delta_total_energy:
        #    # Outlet control
        #    driving_head = self.delta_total_energy
        #else:
            # Inlet control
        #    driving_head = self.inflow.get_average_specific_energy()
            


    def determine_inflow(self):
        # Determine flow direction based on total energy difference

        if self.use_velocity_head:
            self.delta_total_energy = self.inlets[0].get_enquiry_total_energy() - self.inlets[1].get_enquiry_total_energy()
        else:
            self.delta_total_energy = self.inlets[0].get_enquiry_stage() - self.inlets[1].get_enquiry_stage()


        self.inflow  = self.inlets[0]
        self.outflow = self.inlets[1]
        

        if self.delta_total_energy < 0:
            self.inflow  = self.inlets[1]
            self.outflow = self.inlets[0]
            self.delta_total_energy = -self.delta_total_energy


    def get_inflow(self):
        
        return self.inflow
        
        
    def get_outflow(self):
    
        return self.outflow
    






