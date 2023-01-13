""" 
Towradgi Creek 17 August 1998 Storm Event Calibration
By Petar Milevski, some revisions by Gareth Davies
"""

#------------------------------------------------------------------------------
# IMPORT NECESSARY MODULES
#------------------------------------------------------------------------------


import anuga
import time
import numpy
import os

from os.path import join, isdir

from anuga import Polygon_function


def read_polygon_list(poly_list):
    # Alternative to read_polygon_dir -- allows us to control order of polygons
    from anuga import read_polygon
    
    result = []
    for i in range(len(poly_list)):
        result.append((read_polygon(poly_list[i][0]) , poly_list[i][1]))
    return result



#--------------------------------------------------------------------------
# Setup Domain
#--------------------------------------------------------------------------
def setup_domain(simulation):
    """
    Takes input from a simulation object which 
    provides parameters and arguments, and returns a
    sequential domain
    """
    
    args = simulation.args
    verbose = args.verbose
    alg = args.alg

    if not isdir('DEM_bridges'):
        msg = """
################################################################################
#
# Could not the find data directories
#
# You can download these directories using the data_download.py script.
# This will download over 120 MB of data!
#
################################################################################
"""
        raise Exception(msg)
    
    N = args.N
    S = args.S
    E = args.E
    W = args.W
    
    from catchment_info import create_catchment_list
    from catchment_info import create_manning_list
    
    CatchmentList = create_catchment_list(simulation)
    ManningList = create_manning_list(simulation)
    
    #---------------------------------------------------------------------------
    # CREATING MESH
    #---------------------------------------------------------------------------
    
    bounding_polygon = [[W, S], [E, S], [E, N], [W, N]]

    interior_regions = read_polygon_list(CatchmentList)

    # FIXME: Have these in a shapefile / other file and read them in    
    breaklines=[[[306612.336559443,6193708.75358065],
                 [306604.441364239,6193693.17994946]],
                [[306977.886673843,6193753.44134088],
                 [306978.027867398,6193710.94208076]],
                [[306956.001672788,6193750.89985688],
                 [306956.707640564,6193706.14149989]],
                [[306627.303076293,6193697.45809624],
                 [306620.525785644,6193683.62112783]],
                [[307236.83565407,6193741.01630802],
                 [307231.682089306,6193721.03741996]],
                [[307224.975395434,6193742.71063068],
                 [307220.880782334,6193723.36711362]],
                [[307624.764946969,6193615.98941489],
                 [307617.98765632,6193601.44647871]],
                [[307613.328268998,6193623.19028621],
                 [307607.751123568,6193610.97704368]]]

    #---------------------------------------------------------------------------
    # SETUP COMPUTATIONAL DOMAIN
    #---------------------------------------------------------------------------
    domain = anuga.create_domain_from_regions(bounding_polygon,
                                      boundary_tags={'south': [0], 'east': [1],
                                                     'north': [2], 'west': [3]},
                                      maximum_triangle_area=args.maximum_triangle_area,
                                      interior_regions=interior_regions,
                                      breaklines=breaklines,
                                      mesh_filename=args.meshname,
                                      use_cache=False,
                                      verbose=False)
    

    domain.set_flow_algorithm(alg)

    if(not domain.get_using_discontinuous_elevation()):
        raise Exception('This model run relies on a discontinuous elevation solver (because of how topography is set up)')

    domain.set_datadir(args.model_output_dir)
    domain.set_name(args.outname)
        
    print (domain.statistics())
    
    #------------------------------------------------------------------------------
    # APPLY MANNING'S ROUGHNESSES
    #------------------------------------------------------------------------------
    
    if verbose: print ('Calculating complicated polygon friction function')
    friction_list = read_polygon_list(ManningList)
    domain.set_quantity('friction', Polygon_function(friction_list, default=args.base_friction, geo_reference=domain.geo_reference))
    
    # Set a Initial Water Level over the Domain
    domain.set_quantity('stage', 0)
        
    if verbose: print('READING %s' % args.basename+'.npy')
    elev_xyz = numpy.load(args.basename+'.npy')

    # Use nearest-neighbour interpolation of elevation
    if verbose: print('CREATING nearest neighbour interpolator')
    from anuga.utilities.quantity_setting_functions import make_nearestNeighbour_quantity_function
    elev_fun_wrapper = make_nearestNeighbour_quantity_function(elev_xyz, domain)

    if verbose: print ('Applying elevation interpolation function')    
    domain.set_quantity('elevation', elev_fun_wrapper, location='centroids')

    
    return domain


