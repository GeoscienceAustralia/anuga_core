"""Make directory available as a Python package
"""

from anuga.geospatial_data.geospatial_data import *


from numpy._pytesttester import PytestTester
test = PytestTester(__name__)
del PytestTester

