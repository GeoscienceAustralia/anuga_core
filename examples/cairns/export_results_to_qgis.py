import os
import sys
import project

import anuga

scenario = 'fixed_wave'
#scenario = 'slide'
name = 'cairns_' + scenario

print 'output dir:', name

"""
Produce .tif files extracting results of cairns simulation.

Can be viewed with qgis (a project is available cairns_fixed_wave.qgs)
"""
    
from anuga.utilities.plot_utils import Make_Geotif
Make_Geotif(swwFile=name+'.sww', 
             output_quantities=['stage', 'depth', 'velocity', 'elevation'],
             myTimeStep='max', CellSize=1000.0, 
             lower_left=None, upper_right=None,
             EPSG_CODE=32355, 
             proj4string=None,
             velocity_extrapolation=True,
             min_allowed_height=1.0e-05,
             output_dir='.',
             bounding_polygon=project.bounding_polygon,
             verbose=True)

Make_Geotif(swwFile=name+'.sww', 
             output_quantities=['stage'],
             myTimeStep=0,
             CellSize=1000.0, 
             lower_left=None, upper_right=None,
             EPSG_CODE=32355, 
             proj4string=None,
             velocity_extrapolation=True,
             min_allowed_height=1.0e-05,
             output_dir='.',
             bounding_polygon=project.bounding_polygon,
             verbose=True)
