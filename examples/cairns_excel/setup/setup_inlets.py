#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

Set up inflows

Read a discharge input timeseries, and impose it

Gareth Davies, Geoscience Australia 2014+

"""


import numpy
import anuga
from anuga.utilities import spatialInputUtil as su
import scipy
import scipy.interpolate


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
        qfun = scipy.interpolate.interp1d(
            timeseries[:, 0], timeseries[:, 1], kind='linear')

        # Make cross-section line
        line = su.read_polygon(line_file)

        anuga.Inlet_operator(domain, line, qfun, label='Inlet: ' + str(name))
    return
