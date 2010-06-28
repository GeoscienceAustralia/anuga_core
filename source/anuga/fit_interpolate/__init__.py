"""
    2D mesh fitting and interpolation.
    
    Maps quantity data over a 2D mesh. It calculates a smooth gradation of 
    data over the mesh, and allows data to be sampled at any given point.
"""

#Add path of package to PYTHONPATH to allow C-extensions to be loaded
import sys
sys.path += __path__



