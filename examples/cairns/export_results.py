import os
import sys
import project

import anuga

scenario = 'fixed_wave'
name = 'cairns_' + scenario

print 'output dir:', name
which_var = 0

if which_var == 0:    # Stage
    outname = name + '_stage'
    quantityname = 'stage'

if which_var == 1:    # Absolute Momentum
    outname = name + '_momentum'
    quantityname = '(xmomentum**2 + ymomentum**2)**0.5'    #Absolute momentum

if which_var == 2:    # Depth
    outname = name + '_depth'
    quantityname = 'stage-elevation'  #Depth

if which_var == 3:    # Speed
    outname = name + '_speed'
    quantityname = '(xmomentum**2 + ymomentum**2)**0.5/(stage-elevation+1.e-3)'  #Speed

if which_var == 4:    # Elevation
    outname = name + '_elevation'
    quantityname = 'elevation'  #Elevation

print 'start sww2dem'

anuga.sww2dem(name+'.sww',
        outname+'.asc',
        quantity=quantityname,
        cellsize=100,      
        #easting_min=project.eastingmin,
        #easting_max=project.eastingmax,
        #northing_min=project.northingmin,
        #northing_max=project.northingmax,        
        reduction=max, 
        verbose=True)
