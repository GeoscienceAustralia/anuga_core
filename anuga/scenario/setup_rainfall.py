#!/usr/bin/python
"""

Set rainfall

Gareth Davies, Geoscience Australia 2014+

"""


import numpy
import anuga
from anuga.utilities import spatialInputUtil as su
from anuga.scenario.timeseries import make_timeseries_function


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

        # ANUGA timeseries files (.tms) hold a 'rate' quantity that is read with
        # anuga.file_function rather than parsed as CSV. In this case the
        # multiplier is applied as the Rate_operator unit factor (e.g. 1e-3 for
        # mm-based .tms), matching the classic gauge-rainfall workflow.
        if str(timeseries_file).lower().endswith('.tms'):
            rate_fun = anuga.file_function(timeseries_file, quantities='rate')
            anuga.operators.rate_operators.Rate_operator(
                domain, rate=rate_fun, factor=multiplier, polygon=polygon,
                default_rate=0.0, label=timeseries_file)
            continue

        rain_timeseries = numpy.genfromtxt(
            timeseries_file, delimiter=',', skip_header=1)

        # Adjust starttime
        rain_timeseries[:, 0] = rain_timeseries[:, 0] - start_time

        # Convert units to m/s (from mm/hr)
        rain_timeseries[:, 1] = rain_timeseries[:, 1] / (3600. * 1000.) * multiplier

        # Sanity check
        assert rain_timeseries[:, 1].min() >= 0., 'Negative rainfall input'

        # Make interpolation function and add to ANUGA as operator
        if rain_timeseries[:, 1].max() >= 0.:
            # No rain before the series starts (as for .tms files)
            myrain = make_timeseries_function(
                rain_timeseries[:, 0], rain_timeseries[:, 1],
                kind=interpolation_type, before=0.0)
            anuga.operators.rate_operators.Rate_operator(
                domain, rate=myrain, polygon=polygon, label=timeseries_file)

    return
