#!/usr/bin/env python3
"""Optional evolve-loop hooks for the complex anuga_toml_run example.

anuga_toml_run imports a ``user_functions`` module from the scenario directory
(if present) and calls these hooks each yieldstep when the matching
``report_*`` flag is set in the TOML [project] section:

    report_peak_velocity_statistics -> print_velocity_statistics(domain, max_q)
    report_operator_statistics      -> print_operator_inputs(domain)

Delete this file (or set the flags false) to run without the extra reporting.
"""

from anuga.parallel import myid, numprocs, barrier


def print_velocity_statistics(domain, max_quantities):
    """Print the peak flow speed on each rank, plus the running max so far."""
    for i in range(numprocs):
        if myid == i:
            xx = domain.quantities['xmomentum'].centroid_values
            yy = domain.quantities['ymomentum'].centroid_values
            depth = (domain.quantities['stage'].centroid_values
                     - domain.quantities['elevation'].centroid_values)
            wet = depth > 1.0e-03
            d = depth * wet + 1.0e-03 * (~wet)
            speed = (xx ** 2 + yy ** 2) ** 0.5 / d * wet
            msg = f'    [rank {myid}] peak speed = {speed.max():.3f} m/s'
            running_max = getattr(max_quantities, 'max_speed', None)
            if running_max is not None:
                msg += f'   (running max {running_max.max():.3f} m/s)'
            print(msg)
        barrier()
    if myid == 0:
        print('')


def print_operator_inputs(domain):
    """Print the current rate of each fractional-step operator (rain/inlets)."""
    if myid == 0:
        for op in domain.fractional_step_operators:
            label = getattr(op, 'label', op.__class__.__name__)
            if hasattr(op, 'rate'):
                try:
                    rate = op.rate(domain.get_time())
                except TypeError:
                    rate = op.rate
                print(f'    operator {label!r}: rate = {rate}')
    barrier()
