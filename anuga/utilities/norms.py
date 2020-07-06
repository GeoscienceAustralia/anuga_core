"""Definitions of common norm functions, for consistency checks.
"""

from builtins import map
from math import fabs,sqrt
from functools import reduce

def l1_norm(vector):
    """L_1 norm of a vector"""
    return reduce(lambda p,q : p + q, list(map(fabs, vector)), 0.0)

def l2_norm(vector):
    """L_2 norm of a vector"""
    return sqrt(reduce(lambda p,q : p + q, [x ** 2 for x in vector], 0.0))

def linf_norm(vector):
    """L_\infty norm of a vector"""
    return max(list(map(fabs, vector + [0.0])))
