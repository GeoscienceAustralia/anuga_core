#!/usr/bin/python
"""

Set the boundary conditions

Gareth Davies, Geoscience Australia 2014+

"""


import sys
import traceback
import numpy
import anuga
import anuga.utilities.spatialInputUtil as su
from anuga.scenario.timeseries import make_timeseries_function
import anuga.shallow_water.boundaries as asb


def setup_boundary_conditions(domain, project):
    """
    Edit this if new types of boundary conditions need to be supported
    """

    try:
        _setup_boundary_conditions_inner(domain, project)
    except Exception:
        sys.__stdout__.write('\nERROR in setup_boundary_conditions:\n')
        traceback.print_exc(file=sys.__stdout__)
        sys.__stdout__.flush()
        raise


def _setup_boundary_conditions_inner(domain, project):
    # Dict to hold { boundary_tag: boundary_condition_function }
    boundary_tags_and_conditions = {}

    for i in range(len(project.boundary_data)):
        bd = project.boundary_data[i]
        boundary_tag = bd[0]
        boundary_condition_type = bd[1]

        if boundary_condition_type == 'Reflective':
            # Simple reflective boundary
            boundary_tags_and_conditions[boundary_tag] = \
                anuga.Reflective_boundary(domain)
        elif ((boundary_condition_type == 'Stage') or
              (boundary_condition_type == 'Flather_Stage')):
            # Here we read a timeseries, offset the starttime, and then pass
            # that function as a set stage transmissive or flather boundary
            boundary_data_file = bd[2]
            start_time = bd[3]
            # Read data
            boundary_data = numpy.genfromtxt(
                boundary_data_file,
                delimiter=',',
                skip_header=1)
            # Make start time = 0
            boundary_data[:, 0] -= start_time
            # Make interpolation function (held at the first stage before
            # the series starts, for a negative starttime)
            stage_time_fun = make_timeseries_function(
                boundary_data[:, 0],
                boundary_data[:, 1])

            # Make boundary condition
            if boundary_condition_type == 'Stage':
                boundary_function = \
                    asb.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(
                        domain,
                        stage_time_fun)

            elif boundary_condition_type == 'Flather_Stage':
                boundary_function = \
                    asb.Flather_external_stage_zero_velocity_boundary(
                        domain,
                        stage_time_fun)

            boundary_tags_and_conditions[boundary_tag] = \
                boundary_function

        else:
            msg = 'Boundary condition type ' + boundary_condition_type +\
                  ' for tag ' + boundary_tag + ' is not implemented'
            raise NotImplementedError(msg)

    # Check that all tags in boundary_info have been set

    for tag in project.boundary_tags:
        if tag not in boundary_tags_and_conditions:
            msg = 'Need to set boundary_tags_and_conditions for tag = ' + tag
            raise Exception(msg)

    # Set the boundary

    domain.set_boundary(boundary_tags_and_conditions)

    return
