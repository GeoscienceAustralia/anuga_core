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
        t =0.0

        try:
            function(x,y,t)
        except:
            #print 'Problem calling with three arguments'
            try:
                function(x,y)
            except:
                #print 'Problem calling with 2 array arguments'
                try:
                    function(t)
                except:
                    #print 'problem calling with one scalar argument'
                    msg = 'Input argument cannot be called as f(t), f(x,y) or f(x,y,t)'
                    raise Exception, msg
                else:
                    return 't'
            else:
                return 'x,y'
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
        except Modeltime_too_early, e:

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

        except Modeltime_too_late, e:
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






