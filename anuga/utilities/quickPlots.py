

from builtins import str

"""""

Quick code to plot ANUGA outputs

"""""

from anuga.utilities import plot_utils as util
import numpy
import matplotlib
matplotlib.use('Agg') # This let's us run on NCI
from matplotlib import pyplot
import argparse
import os
    


def quickPlots(swwFile=None, ocean_land_threshold=None, fig_dir=None, figScale=None):
    """
        Routine to make a set of 'quick and dirty' plots of initial conditions + maxima. 
        Useful for preliminary check on ANUGA outputs
    """

    if(swwFile is None):
        parser.print_help()
        print(' ')
        raise Exception('Must give an sww file')

    # Make directory for figures
    try:
        os.mkdir(fig_dir)
    except:
        'Cannot make directory'
        pass


    # Read in sww file
    p=util.get_output(swwFile)
    p2=util.get_centroids(p,velocity_extrapolation=True)

    xRange=p2.x.max()-p2.x.min()
    yRange=p2.y.max()-p2.y.min()

    figSize=(figScale, figScale*yRange//xRange)

    # Use spatial coordinates
    x=p2.x+p.xllcorner
    y=p2.y+p.yllcorner

    # Plot friction
    try:
        pyplot.figure(figsize=figSize)
        pyplot.scatter(x,y,c=p2.friction,edgecolors='none')
        pyplot.gca().set_aspect('equal')
        pyplot.title('Friction')
        pyplot.colorbar()
        pyplot.savefig(fig_dir+'/Friction.png')
        pyplot.close()
    except:
        print('Cannot plot friction')

    # Plot elevation
    try:
        pyplot.figure(figsize=figSize)
        pyplot.scatter(x,y,c=p2.elev,edgecolors='none')
        pyplot.gca().set_aspect('equal')
        pyplot.title('Elevation')
        pyplot.colorbar()
        pyplot.savefig(fig_dir+'/Elevation.png')
        pyplot.close()
    except:
        print('Cannot plot elevation')

    # Plot Initial Stage (where elevation<ocean_land_threshold)
    pyplot.figure(figsize=figSize)
    pyplot.scatter(x,y,c=p2.stage[0,:]*(p2.elev<ocean_land_threshold),edgecolors='none')
    pyplot.gca().set_aspect('equal')
    pyplot.title('Initial Stage (zero where elevation > '+ str(ocean_land_threshold) +' )')
    pyplot.colorbar()
    pyplot.savefig(fig_dir+'/Initial_stage.png')
    pyplot.close()

    # Plot Initial Depth 
    pyplot.figure(figsize=figSize)
    pyplot.scatter(x,y,c=p2.height[0,:],edgecolors='none')
    pyplot.gca().set_aspect('equal')
    pyplot.title('Initial Depth')
    pyplot.colorbar()
    pyplot.savefig(fig_dir+'/Initial_depth.png')
    pyplot.close()

    # Initial Speed
    pyplot.figure(figsize=figSize)
    pyplot.scatter(x,y,c=p2.vel[0,:],edgecolors='none')
    pyplot.gca().set_aspect('equal')
    pyplot.title('Initial speed ')
    pyplot.colorbar()
    pyplot.savefig(fig_dir+'/Initial_speed.png')
    pyplot.close()

    # Triangle areas
    triA=util.triangle_areas(p)
    tri_len_scale=(triA*2)**0.5
    pyplot.figure(figsize=figSize)
    pyplot.scatter(x,y,c=tri_len_scale,edgecolors='none')
    pyplot.gca().set_aspect('equal')
    pyplot.title('Mesh triangle side length scale (2*area)^0.5')
    pyplot.colorbar()
    pyplot.savefig(fig_dir+'/Triangle_side_length.png')
    pyplot.close()

    # Max stage in wet areas
    maxStage=p2.stage.max(axis=0)
    pyplot.figure(figsize=figSize)
    pyplot.scatter(x,y,c=maxStage*(maxStage>p2.elev), edgecolors='none')
    pyplot.gca().set_aspect('equal')
    pyplot.title('Max stage (zeroed in dry areas)')
    pyplot.colorbar()
    pyplot.savefig(fig_dir+'/Max_stage.png')
    pyplot.close()

    # Max depth
    maxDepth=p2.height.max(axis=0)
    pyplot.figure(figsize=figSize)
    pyplot.scatter(x,y,c=maxDepth, edgecolors='none')
    pyplot.gca().set_aspect('equal')
    pyplot.title('Max depth')
    pyplot.colorbar()
    pyplot.savefig(fig_dir+'/Max_depth.png')
    pyplot.close()

    # Max speed in wet areas
    maxSpeed=p2.vel.max(axis=0)
    pyplot.figure(figsize=figSize)
    pyplot.scatter(x,y,c=maxSpeed, edgecolors='none')
    pyplot.gca().set_aspect('equal')
    pyplot.title('Max speed')
    pyplot.colorbar()
    pyplot.savefig(fig_dir+'/Max_speed.png')
    pyplot.close()

####################################################

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Make quick and dirty plots from sww file')
    parser.add_argument('-sww', type=str, default = None,
                   help='Name of sww file')
    parser.add_argument('-fig_dir', type=str, default="FIGS",
                   help='Directory to save figures in')
    parser.add_argument('-ocean_land_threshold', type=float, default=9.0e+20,
                   help='Elevation dividing "ocean" from "land". "Land" areas have stage=0.0 for the initial condition plot, which can assist in seeing the stage in wet areas. For example, for tsunami studies I might set ocean_land_threshold=0.0 to see the initial water surface perturbation')
    parser.add_argument('-figScale', type=float, default=12,
                   help='Increase this to increase the figure size')
    args = parser.parse_args()

    ## INPUTS ##
    swwFile=args.sww #'mySww.sww'
    ocean_land_threshold=args.ocean_land_threshold #0. # Above this elevation, we have 'land'
    fig_dir=args.fig_dir
    figScale=args.figScale

    try:
        quickPlots(swwFile, ocean_land_threshold, fig_dir, figScale)
    except:
        print('Error in quick_plots')
