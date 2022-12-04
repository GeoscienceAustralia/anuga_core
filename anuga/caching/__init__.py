"""
    Generic caching module.
    
    Allows for the disk caching of the results of any function. If a function
    is cached, its return values will be stored on the local hard drive.
    If the function is called with identical parameters in the future, the
    cached result will be returned.
"""


from .caching import *

from numpy._pytesttester import PytestTester
test = PytestTester(__name__)
del PytestTester




