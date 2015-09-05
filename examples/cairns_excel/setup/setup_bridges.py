
#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

Setup bridges (using hecras internal boundary rating curves)  

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
from anuga.structures.internal_boundary_functions import hecras_internal_boundary_function

def setup_bridges(domain, project):
    """
    Extract bridge data from project class
    and apply the internal boundary operator

    Note that the bridge deck information
    is applied when setting the elevation
    """
    bridge_data = project.bridge_data

    for i in range(len(bridge_data)):

        # Extract bridge parameters
        bd = bridge_data[i]

        label = bd[0]
        # bd[1] and bd[2] are used elsewhere to set deck elevation
        exchange_line_0 = su.read_polygon(bd[3])
        exchange_line_1 = su.read_polygon(bd[4])
        exchange_lines = [exchange_line_0, exchange_line_1]
        enquiry_gap = bd[5]
        #apron = bd[6]
        #assert apron==0.0, 'Apron must be zero until parallel apron issues fixed'
        internal_boundary_curve_file = bd[6]
        vertical_datum_offset = bd[7]
        smoothing_timescale = bd[8]



        # Function which computes Q
        rating_curve = hecras_internal_boundary_function(
            internal_boundary_curves_file=internal_boundary_curve_file,
            allow_sign_reversal=True,
            vertical_datum_offset=vertical_datum_offset)

        # Add operator as side-effect of this operation
        bridge = Internal_boundary_operator(
            domain,
            rating_curve,
            exchange_lines=exchange_lines,
            enquiry_gap=enquiry_gap,
            apron = 0.0,
            zero_outflow_momentum=False,
            smoothing_timescale=smoothing_timescale,
            logging=True,
            label=label,
            verbose=True)

    return
