#!/usr/bin/python
"""

Set up inflows

Read a discharge input timeseries, and impose it

Gareth Davies, Geoscience Australia 2014+

"""


import numpy
import anuga
from anuga.utilities import spatialInputUtil as su
from anuga.scenario.timeseries import make_timeseries_function


def setup_inlets(domain, project):
    """
    Add inlets to domain
    """

    inlet_data = project.inlet_data

    for i in range(len(inlet_data)):
        name = inlet_data[i][0]
        line_file = inlet_data[i][1]
        timeseries_file = inlet_data[i][2]
        start_time = inlet_data[i][3]

        # Add inlet
        timeseries = numpy.genfromtxt(
            timeseries_file, delimiter=',', skip_header=1)

        # Adjust start time
        timeseries[:, 0] = timeseries[:, 0] - start_time

        # Make discharge function
        # (held at the first value before the series starts, so a negative
        # starttime gives a burn-in at the initial inflow)
        qfun = make_timeseries_function(
            timeseries[:, 0], timeseries[:, 1], kind='linear')

        # Make cross-section line
        line = su.read_polygon(line_file)

        anuga.Inlet_operator(domain, line, qfun, label='Inlet: ' + str(name))
    return
