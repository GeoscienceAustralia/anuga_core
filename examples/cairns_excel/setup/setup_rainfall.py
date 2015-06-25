#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

Set rainfall

Gareth Davies, Geoscience Australia 2014+

"""


import numpy
import anuga
import scipy.interpolate
from anuga.utilities import spatialInputUtil as su


def setup_rainfall(domain, project):
    """
    Function to add rainfall operators to the domain
    """
    rain_data = project.rain_data

    for i in range(len(rain_data)):
        timeseries_file = rain_data[i][0]
        start_time = rain_data[i][1]
        interpolation_type = rain_data[i][2]

        # Get polygon defining rainfall extent
        if ((len(rain_data[i]) >= 4) and (rain_data[i][3] != 'All')):
            polygon = su.read_polygon(rain_data[i][3])
        else:
            polygon = None

        if len(rain_data[i]) >=5:
            multiplier = rain_data[i][4]*1.0
        else:
            multiplier = 1.0

        rain_timeseries = scipy.genfromtxt(
            timeseries_file, delimiter=',', skip_header=1)

        # Adjust starttime
        rain_timeseries[:, 0] = rain_timeseries[:, 0] - start_time

        # Convert units to m/s (from mm/hr)
        rain_timeseries[:, 1] = rain_timeseries[:, 1] / (3600. * 1000.) * multiplier

        # Sanity check
        assert rain_timeseries[:, 1].min() >= 0., 'Negative rainfall input'

        # Make interpolation function and add to ANUGA as operator
        if rain_timeseries[:, 1].max() >= 0.:
            myrain = scipy.interpolate.interp1d(
                rain_timeseries[:, 0], rain_timeseries[:, 1],
                kind=interpolation_type)
            anuga.operators.rate_operators.Rate_operator(
                domain, rate=myrain, polygon=polygon, label=timeseries_file)

    return
