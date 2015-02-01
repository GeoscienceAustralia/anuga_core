"""Definitions of common norm functions, for consistency checks.
"""

from math import fabs,sqrt

def l1_norm(vector):
    """L_1 norm of a vector"""
    return reduce(lambda p,q : p + q, map(fabs, vector), 0.0)

def l2_norm(vector):
    """L_2 norm of a vector"""
    return sqrt(reduce(lambda p,q : p + q, map(lambda x : x ** 2, vector), 0.0))

def linf_norm(vector):
    """L_\infty norm of a vector"""
    return max(map(fabs, vector + [0.0]))
