"""

Time series helpers shared by the scenario setup modules.

"""

import scipy.interpolate


def make_timeseries_function(times, values, kind='linear', before='first'):
    """Return f(t) interpolating *values* at *times*.

    Before the first time f returns *before*: the first value when
    before='first', otherwise the number given. That lets a scenario start
    before its forcing data does (a negative [project] starttime used as a
    burn-in), with the inflow or tide held at its first value, or the rain
    set to zero. After the last time f raises ValueError as before, so a
    series that runs out early is still reported.
    """
    f = scipy.interpolate.interp1d(times, values, kind=kind)
    t0 = times[0]
    left = values[0] if before == 'first' else float(before)

    def fun(t):
        if t < t0:
            return left
        return f(t)

    return fun
