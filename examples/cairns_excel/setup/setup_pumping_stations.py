
#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

Setup pumping stations (using internal boundary operator)  

Gareth Davies, Geoscience Australia 2014+

"""


import numpy
import anuga
import scipy.interpolate
from anuga.utilities import spatialInputUtil as su
from anuga.parallel import myid, numprocs, barrier, finalize
from anuga.parallel.parallel_api import pypar_available
if pypar_available:
    from anuga.parallel.parallel_operator_factory import Internal_boundary_operator 
else:
    from anuga.structures.internal_boundary_operator import Internal_boundary_operator 
from anuga.structures.internal_boundary_functions import pumping_station_function

def setup_pumping_stations(domain, project):
    """
    Extract pumping station data from project class
    and apply the internal boundary operator

    """
    pumping_station_data = project.pumping_station_data

    for i in range(len(pumping_station_data)):

        # Extract pumping station parameters
        ps = pumping_station_data[i]

        label = ps[0]
        pump_capacity = ps[1]
        pump_rate_of_increase = ps[2]
        pump_rate_of_decrease = ps[3]
        hw_to_start_pumping = ps[4]
        hw_to_stop_pumping = ps[5]

        # Pumping station basin polygon + elevation are used elsewhere

        exchange_line_0 = su.read_polygon(ps[8])
        exchange_line_1 = su.read_polygon(ps[9])
        exchange_lines = [exchange_line_0, exchange_line_1]

        smoothing_timescale = ps[10]

        print 'Need to implement elevation data adjustments'

        # Function which computes Q
        pump_behaviour = pumping_station_function(
            domain=domain,
            pump_capacity=pump_capacity,
            hw_to_start_pumping=hw_to_start_pumping,
            hw_to_stop_pumping=hw_to_stop_pumping,
            initial_pump_rate=0.,
            pump_rate_of_increase=pump_rate_of_increase,
            pump_rate_of_decrease=pump_rate_of_decrease
            )
        
        # Add operator as side-effect of this operation
        pumping_station = Internal_boundary_operator(
            domain,
            pump_behaviour,
            exchange_lines=exchange_lines,
            enquiry_gap=0.,
            apron = 0.0,
            smoothing_timescale=smoothing_timescale,
            compute_discharge_implicitly=False,
            logging=True,
            label=label,
            verbose=True)

    return
