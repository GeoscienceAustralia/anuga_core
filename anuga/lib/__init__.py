"""Make directory available as a Python package
"""

from numpy._pytesttester import PytestTester
test = PytestTester(__name__)
del PytestTester






