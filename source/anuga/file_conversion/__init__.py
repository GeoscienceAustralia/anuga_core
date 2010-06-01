"""Make directory available as a Python package
"""

# Add path of package to PYTHONPATH to allow C-extensions to be loaded
import sys
sys.path += __path__

