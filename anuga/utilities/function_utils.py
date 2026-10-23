""" Generic function utilities to test type of function
"""


import numpy as num
from anuga.fit_interpolate.interpolate import Modeltime_too_early
from anuga.fit_interpolate.interpolate import Modeltime_too_late

import warnings
#warnings.filterwarnings('default')

def determine_function_type(function):
    """Test type of function, either
    - scalar
    - function of t (scalar)
    - function of two arrays x,y
    - function of two arrays x,y and t (scalar)"""

    #------------------------------------------
    # Check that argument is at least a scalar
    # or csllable
    #------------------------------------------
    msg = "Input argument must be a scalar, or a function or None"

    if function is None:
        return None



    assert (isinstance(function, (int, float, list)) or
            isinstance(function, num.ndarray) or
            callable(function)), msg


    if callable(function):

        # test if temporal
        x = num.array([0.0, 1.0])
        y = num.array([0.0, 2.0])
        t = 0.0

        #function(x,y,t)
        try:
            function(x,y,t)
        except TypeError:
            #print 'Problem calling with three arguments'
            try:
                function(x,y)
            except TypeError:
                #print 'Problem calling with 2 array arguments'
                try:
                    function(t)
                except TypeError:
                    #print 'problem calling with one scalar argument'
                    msg = 'Input argument cannot be called as f(t), f(x,y) or f(x,y,t)'
                    raise Exception(msg)
                except ValueError:
                    #print 'problem calling out of range'
                    return 't'
                except Modeltime_too_early:
                    #print 'test argument out of range'
                    return 't'
                except Modeltime_too_late:
                    #print 'test argument out of range'
                    return 't'

                else:
                    return 't'
            except ValueError:
                return 'x,y'
            else:
                return 'x,y'
        except ValueError:
            return 'x,y,t'
        else:
            return 'x,y,t'

    elif isinstance(function, (int,float)):
        return 'scalar'

    elif isinstance(function, (list, num.ndarray)):
        return 'array'


def evaluate_temporal_function(function, t, default_left_value=None, default_right_value=None):

    if  callable(function):
        try:
            result = function(t)
        except Modeltime_too_early as e:

            if default_left_value is None:
                msg = '%s: Trying to evaluate function earlier than specified in the data set.\n' %str(e)
                raise Modeltime_too_early(msg)
            else:
                # Pass control to default left function

                warnings.warn('Using default_left_value')
                if callable(default_left_value):
                    result = default_left_value(t)
                else:
                    result = default_left_value

        except Modeltime_too_late as e:
            if default_right_value is None:
                msg = '%s: Trying to evaluate function later than specified in the data set.\n' %str(e)
                raise Modeltime_too_late(msg)
            else:
                # Pass control to default right function

                warnings.warn('Using default_right_value')
                if callable(default_right_value):
                    result = default_right_value(t)
                else:
                    result = default_right_value
    else:
        result = function

    return result








def evaluate_file_function_all_points(F, t):
    """Evaluate a ``file_function`` at time *t* for every interpolation point.

    ``file_function`` objects are evaluated one point at a time with
    ``F(t, point_id=i)``.  For a spatial field applied at every mesh centroid
    (or node) each timestep that is prohibitively slow in Python, so this
    helper interpolates the object's precomputed time series directly.

    Parameters
    ----------
    F : Interpolation_function
        The object returned by :func:`anuga.file_function`.
    t : float
        Model time.  Must lie within ``F.time``.

    Returns
    -------
    ndarray
        Shape ``(n_quantities, n_points)`` for a spatial file function, or
        ``(n_quantities,)`` for a time-only one (e.g. a ``.tms`` file).
        Quantities are ordered as in ``F.quantity_names``.
    """
    times = getattr(F, 'time', None)
    values = getattr(F, 'precomputed_values', None)
    names = getattr(F, 'quantity_names', None)
    if times is None or values is None or names is None:
        raise TypeError('evaluate_file_function_all_points expects the object '
                        'returned by anuga.file_function; got %s'
                        % type(F).__name__)

    times = num.asarray(times, dtype=float)
    msg = ('Model time %.16f is not contained in function domain '
           '[%.16f:%.16f]' % (t, times[0], times[-1]))
    if t < times[0]:
        raise Modeltime_too_early(msg)
    if t > times[-1]:
        raise Modeltime_too_late(msg)

    # Bracket t: times[i0] <= t <= times[i1]
    i1 = int(num.searchsorted(times, t, side='left'))
    i1 = min(i1, len(times) - 1)
    i0 = max(i1 - 1, 0)
    if times[i1] == times[i0]:
        ratio = 0.0
    else:
        ratio = (t - times[i0]) / (times[i1] - times[i0])

    out = []
    for name in names:
        Q = num.asarray(values[name], dtype=float)
        q0 = Q[i0]
        q = q0 if ratio == 0.0 else q0 + ratio * (Q[i1] - q0)
        out.append(q)
    return num.array(out)
