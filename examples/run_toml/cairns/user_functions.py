#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Home for various 'user-defined' functions.

The user can define any new functions they need in this script, and make
small modifications to the main scripts to use it.

This is usually nicer than making detailed hacky changes to the main scripts

If you have a nice function, consider pushing it to the ANUGA source


Gareth Davies, Geoscience Australia 2014+
"""

import os
import sys
import anuga
from anuga.parallel import myid, numprocs, finalize, barrier


def print_velocity_statistics(domain, max_quantities):
    """
    Utility for printing velocity stats in the evolve loop
    """

    # Print velocity statistics

    for i in range(numprocs):
        if myid == i:

            # Peak speed info

            xx = domain.quantities['xmomentum'].centroid_values
            yy = domain.quantities['ymomentum'].centroid_values
            dd = domain.quantities['stage'].centroid_values \
                - domain.quantities['elevation'].centroid_values
            dd = dd * (dd > 1.0e-03) + 1.0e-03 * (dd <= 1.0e-03)
            vv = 1 / dd * (xx ** 2 + yy ** 2) ** 0.5
            vv = vv * (dd > 1.0e-03)
            print('    Processor ', myid)
            print('    @ Peak velocity is: ', vv.max(), vv.argmax())
            print('     &- MaxSpeedHistory: ', \
                max_quantities.max_speed.max())
            print('     %- FUF: ', domain.flux_update_frequency.mean())
        else:
            pass
        barrier()

    # Make a newline

    if myid == 0:
        print('')

    return


def print_operator_inputs(domain):
    """

        Utility for printing discharge information for some operators

    """

    # Shorthand notation
    operators = domain.fractional_step_operators

    # Rainfall first
    if myid == 0:
        for i in range(len(operators)):
            if hasattr(operators[i], 'rate'):
                print('    Operator ' + operators[i].label + \
                      ' rate = ' + str(operators[i].rate(domain.time)))

    barrier()

    # Inlets
    for i in range(len(operators)):
        if hasattr(operators[i], 'applied_Q'):
            print('    Operator ' + operators[i].label + \
                  ' Q = ' + str(operators[i].applied_Q))
    barrier()

    if myid == 0:
        print(' ')

    return
