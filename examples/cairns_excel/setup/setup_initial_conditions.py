#!/usr/bin/python
# -*- coding: utf-8 -*-

"""

Make initial condition

Gareth Davies, Geoscience Australia 2014+

"""


import numpy
import glob
import anuga
from anuga.utilities import quantity_setting_functions as qs
from anuga.utilities import spatialInputUtil as su

# Local function
from setup.spatially_averaged_function import \
    make_spatially_averaged_function


def setup_initial_conditions(domain, project):
    """
    Set the initial conditions for ANUGA

    INPUT: domain  = the anuga domain
    project = the project module

    """

    # #########################################################################
    # INSTRUCTIONS ON USAGE
    #
    # Here we set the initial quantities, optionally applying different types
    # of interpolation in different polygons.
    #
    # To do this, each quantity needs an ordered list of [polygon, value] pairs
    #  (where we will set xy points in 'polygon' using 'value').
    #
    # Possible values of 'polygon' are:
    #  1) An ANUGA polygon
    #  2) A filename with .shp extension, or a csv file in anuga_polygon format
    #  (with any extension)
    #  3) None -- apply nowhere
    #  4) 'Extent' (only valid if value is a raster filename) -- apply over the
    #  full raster extent
    #  5) 'All' -- apply everywhere -- must be the last polygon-value pair
    #
    # Possible values of 'value' are:
    # 1) A constant
    # 2) A 3-column numpy array with x,y,z (then nearest neighbour
    # interpolation is used)
    # 3) A gdal-compatible raster (pixel-lookup is used)
    # 4) A comma separated variable file with x,y,z data. The file MUST have a
    # .csv or .txt extension
    # 5) Any function f(x,y). Note the function must act on 'ANUGA-transformed'
    # xy coordinates which have
    #    min(x), min(y) = (0.,0.)
    #
    # If polygons overlap, ones earlier in the list take preference
    # e.g.
    #   [ [polygon_1, value_1], [polygon_2, value_2],
    #     [polygon_3, value_3], ['All', value_4] ]
    # will set all points in polygon_1 using value_1, then all points in
    # polygon_2 which were not in polygon_1 using value_2, then all points in
    # polygon_3 which were not previously set using value_3, then everything
    # remaining with value_4
    #
    #
    # For more information, see (in anuga.utilities.quantity_setting_functions)
    #   composite_quantity_setting_function
    #

    # #########################################################################
    #
    # Set the quantities here
    #
    # #########################################################################

    # To allow for 'discontinuous' elevation or other quantities, set
    # location='centroids'

    ##########################################################################
    def quick_set_quantity(quantity_name, quantity_data, domain,
                           quantity_clip_range, quantity_mean,
                           quantity_additions, location='vertices',
                           mean_type='mean'):
        """Convenience function to set the quantities
            and reduce boilerplate code.
            Sets the quantity + takes care of quantity additions
        """

        # Make the function
        quantity_function = \
            qs.composite_quantity_setting_function(
                quantity_data,
                domain,
                quantity_clip_range,
                nan_treatment='fall_through',
                default_k_nearest_neighbours=2)
        # Spatially average the function if required
        if quantity_mean is not None:
            grid_spacing = [quantity_mean, quantity_mean]
            quantity_function = make_spatially_averaged_function(
                quantity_function, domain, approx_grid_spacing=grid_spacing,
                averaging=mean_type)

        print quantity_name, quantity_data, domain, quantity_clip_range, quantity_mean, quantity_additions

        print quantity_function(numpy.array([0.0]),numpy.array([0.0]))
        # Set the quantity
        domain.set_quantity(quantity_name, quantity_function,
                            location=location)

        # Treat additions
        quantity_addition_function = \
            qs.composite_quantity_setting_function(
                quantity_additions, domain, nan_treatment='fall_through')
        domain.add_quantity(quantity_name, quantity_addition_function,
            location=location)
    ##########################################################################

    # Elevation
    quick_set_quantity(quantity_name='elevation',
                       quantity_data=project.elevation_data,
                       domain=domain,
                       quantity_clip_range=project.elevation_clip_range,
                       quantity_mean=project.elevation_mean,
                       quantity_additions=project.elevation_additions,
                       location='centroids',
                       mean_type='mean')

    # Friction -- if averaging is used, harmonic mean is probably a better
    # averaging method, although no method will be perfect.
    # Why? Say there are many points in a cell, with constant Sf, depth d, but
    # variable n. Pointwise,
    #     U = 1/n * Sf^(0.5) * d^(2/3)
    # Say Sf is fixed (by topography for steady uniform flow)
    # Say d^(2/3) is ~ 1.
    # Then mean(U) and mean(Ud) are preserved by taking n = 1/mean(1/n)
    quick_set_quantity(quantity_name='friction',
                       quantity_data=project.friction_data,
                       domain=domain,
                       quantity_clip_range=project.friction_clip_range,
                       quantity_mean=project.friction_mean,
                       quantity_additions=project.friction_additions,
                       location='centroids',
                       mean_type='harmonic_mean')
    # Stage
    quick_set_quantity(quantity_name='stage',
                       quantity_data=project.stage_data,
                       domain=domain,
                       quantity_clip_range=project.stage_clip_range,
                       quantity_mean=project.stage_mean,
                       quantity_additions=project.stage_additions,
                       location='centroids',
                       mean_type='mean')

    # xmomentum
    quick_set_quantity(quantity_name='xmomentum',
                       quantity_data=project.xmomentum_data,
                       domain=domain,
                       quantity_clip_range=project.xmomentum_clip_range,
                       quantity_mean=project.xmomentum_mean,
                       quantity_additions=project.xmomentum_additions,
                       location='centroids',
                       mean_type='mean')

    # ymomentum
    quick_set_quantity(quantity_name='ymomentum',
                       quantity_data=project.ymomentum_data,
                       domain=domain,
                       quantity_clip_range=project.ymomentum_clip_range,
                       quantity_mean=project.ymomentum_mean,
                       quantity_additions=project.ymomentum_additions,
                       location='centroids',
                       mean_type='mean')

    return
