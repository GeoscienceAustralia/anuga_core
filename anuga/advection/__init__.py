"""
    A simple implementation of the shallow wave equation, mainly for test
    purposes.
"""


from .advection import Advection_Domain


from numpy._pytesttester import PytestTester
test = PytestTester(__name__)
del PytestTester

