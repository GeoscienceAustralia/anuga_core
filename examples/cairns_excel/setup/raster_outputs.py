#!/usr/bin/python
# -*- coding: utf-8 -*-
############################################################

# Get the modules we need

"""

Code to make spatial outputs (GeoTifs)

Gareth Davies, Geoscience Australia 2014+

"""


import sys
import os
import glob
import numpy
from anuga.utilities import plot_utils as util
from anuga.utilities import spatialInputUtil as su
from anuga.file.netcdf import NetCDFFile


# @@@@@@@@@@@@@@@@@@#
# STOP EDITING HERE
# @@@@@@@@@@@@@@@@@@#
#
#

def make_resampled_elevation(
    elevation_raster,
    raster_extent,
    cell_size,
    clip_polygon,
    clip_polygon_layer,
    proj4string,
    tif_output_dir,
):
    """
    Call gdalwarp to get the resampled elevation data
    """

    import platform
    ostype = platform.system()

    if ostype == 'Linux':
        gdalwarp = 'gdalwarp'
        gdal_rasterize = 'gdal_rasterize'
    elif ostype == 'Windows':
        gdalwarp = 'gdalwarp.exe'
        gdalwarp = 'gdal_rasterize.exe'
    else:
        raise Exception('Unrecognized OS')

    rast_xmin = raster_extent[0]
    rast_xmax = raster_extent[1]
    rast_ymin = raster_extent[2]
    rast_ymax = raster_extent[3]
    xres = cell_size
    yres = cell_size

    # allSRS='EPSG:'+str(EPSG_CODE)

    all_srs = proj4string
    new_dem = tif_output_dir + 'LIDAR_resampled.tif'
    new_mask = tif_output_dir + 'Mask_resampled.tif'

    extent_info = ' -te ' + str(rast_xmin) + ' ' + str(rast_ymin) + ' ' \
        + str(rast_xmax) + ' ' + str(rast_ymax)
    res_info = ' -tr ' + str(xres) + ' ' + str(yres)

    # mkdir -p ElevationToPlot

    gdalwarp_command = gdalwarp + ' -s_srs ' + all_srs + ' -t_srs ' \
        + all_srs + extent_info + ' ' + res_info + ' ' \
        + elevation_raster + ' ' + new_dem + ' -overwrite'

    os.system(gdalwarp_command)

    import shutil
    shutil.copyfile(new_dem, new_mask)

    gdal_rasterize_command = gdal_rasterize + ' -a_srs ' + all_srs \
        + ' -burn 1 -l ' + clip_polygon_layer + ' ' + extent_info + ' ' \
        + res_info + ' ' + clip_polygon + ' ' + new_mask
    print gdal_rasterize_command
    os.system(gdal_rasterize_command)

    return


#################################################################

def gdal_calc_command(
    stage,
    depth,
    elevation,
    mask,
    depth_threshold,
):
    """
    Call gdal_calc to compute depth by draping
    stage over the high res elevation data
    """

    gdal_calc = 'gdal_calc.py '
    gdal_input = '-A ' + stage + ' -B ' + depth + ' -C ' + elevation \
        + ' -D ' + mask
    calc_command = " --calc '(A-C)*(A>C)*(B>" + str(depth_threshold) \
        + ")*D' "

    # rastOut=$swwTifDir$swwBase'_depth_high_res.tif'

    rast_out = os.path.dirname(depth) + '/HighRes_' \
        + os.path.basename(depth)

    gdalcalc_command = gdal_calc + gdal_input + calc_command \
        + '--outfile ' + rast_out
    print gdalcalc_command
    os.system(gdalcalc_command)

    return


############################################################

def make_me_some_tifs(
    sww_file,
    bounding_polygon,
    proj4string,
    my_time_step='collected_max',
    tif_output_subdir='/TIFS/',
    cell_size=5.0,
    k_nearest_neighbours=1,
    make_highres_drape_plot=False,
    elevation_raster=None,
    depth_threshold=None,
    clip_polygon=None,
    clip_polygon_layer=None,
    creation_options=['COMPRESS=DEFLATE'],
):
    """

    ### INPUT DATA
    - swwFile -- Full path name of sww to read outputs from
    - bounding_polygon -- ANUGA's bounding polygon, or another polygon to clip
      the rasters to, [in ANUGA's polygon format].
    - proj4string defining the coordinate system
    - my_time_step -- option from util.Make_Geotif
      use 'max' to plot the maxima
      use [0, 5, 10] to plot the first, sixth and eleventh output time step
      use 'collected_max' to read the csv outputs from
        collect_max_quantities_operator
    - tif_outputdir -- Write outputs to this folder inside the swwFile's
      directory (MUST INCLUDE TRAILING SLASH /)
    - cell_size -- Desired raster cellSize (m)
    - k_nearest_neighbours -- use inverse distance weighted interpolation with this many neighbours
    - make_highres_drape_plot -- True/False, Make a high-res drape plot?
    - elevation_raster -- Filename of elevation raster for 'high-res-drape'
      depth plot [from subtracting stage from high res topography]
    - depth_threshold -- Depth threshold for high-res-drape' plot (m). If
      ANUGA's triangle has depth
      < depth_threshold, then depth=0 in all high-res cells in that triangle
    - clipPolygon -- Polygon to clip 'high-res-drape' plot. Must be
      provided if make_highres_drape_plot==True (can use bounding polygon or
      another choice)
    - clipPolygonLayer -- Layer for above as used by gdal (usually shapefile
      name without .shp)
    - creation_options -- list of gdal tif creation options

    ## OUTPUT
    Nothing is returned, but tifs are made in tif_output_subdir inside the
    swwFile directory

    """

    # Convert utm_zone to proj4string
    # p=Proj(proj='utm', south=(utm_zone<0.),
    #        zone=abs(utm_zone), ellps='WGS84')
    # proj4string = p.srs

    tif_output_dir = os.path.dirname(sww_file) + tif_output_subdir

    # Make the geotifs

    if my_time_step == 'collected_max':

        # Read max quantity output files inputs

        max_qfiles = glob.glob(os.path.dirname(sww_file)
                               + '/*_UH_MAX.csv')
        if len(max_qfiles) == 0:
            raise Exception(
                'Cannot find any files containing collected maxima')
        for i in range(len(max_qfiles)):
            if i == 0:
                max_qs = numpy.genfromtxt(max_qfiles[i], delimiter=',')
            else:
                extra_data = numpy.genfromtxt(max_qfiles[i], delimiter=',')
                max_qs = \
                    numpy.vstack([max_qs, extra_data])

        # Make the geotiff's

        for (i, quant) in enumerate(['stage', 'depth', 'velocity',
                                     'depthIntegratedVelocity']):

            # FIXME: The interpolation function is remade for every quantity,
            # since only a 3 column array can be passed to Make_Geotif Could
            # make it fast (do it only once) by changing Make_Geotif

            tmp_arr = max_qs[:, [0, 1, i + 2]]

            util.Make_Geotif(
                tmp_arr,
                output_quantities=[quant + '_MAX'],
                CellSize=cell_size,
                proj4string=proj4string,
                verbose=True,
                bounding_polygon=bounding_polygon,
                output_dir=tif_output_dir,
                creation_options=creation_options,
                k_nearest_neighbours=k_nearest_neighbours,
            )

        # Also plot elevation + friction
        # Try to reduce memory demands by only extracting first timestep
        fid = NetCDFFile(sww_file)

        # Make xy coordinates (note -- max_quantity_collector outputs might
        # have repeated x,y at parallel ghost cells)
        x_v = fid.variables['x'][:] + fid.xllcorner
        y_v = fid.variables['y'][:] + fid.yllcorner
        vols = fid.variables['volumes'][:]
        xc = (1. / 3.) * (x_v[vols[:, 0]] + x_v[vols[:, 1]] + x_v[vols[:, 2]])
        yc = (1. / 3.) * (y_v[vols[:, 0]] + y_v[vols[:, 1]] + y_v[vols[:, 2]])

        for (i, quant) in enumerate(['elevation_c', 'friction_c']):

            # Get the quantity if it exists
            if fid.variables.has_key(quant):
                quant_values = fid.variables[quant]
                # If multi time-steps, only get first timestep
                if(len(quant_values.shape) > 1):
                    quant_values = quant_values[0, :]
            else:
                # Set quantity to nan if it is not stored
                quant_values = xc * 0. + numpy.nan

            tmp_arr = numpy.vstack([xc, yc, quant_values]).transpose()

            util.Make_Geotif(
                tmp_arr,
                output_quantities=[quant + '_INITIAL'],
                CellSize=cell_size,
                proj4string=proj4string,
                verbose=True,
                bounding_polygon=bounding_polygon,
                output_dir=tif_output_dir,
                creation_options=creation_options,
                k_nearest_neighbours=k_nearest_neighbours,
            )

    else:
        util.Make_Geotif(
            sww_file,
            myTimeStep=my_time_step,
            output_quantities=[
                'depth',
                'stage',
                'elevation',
                'velocity',
                'depthIntegratedVelocity',
                'friction',
            ],
            CellSize=cell_size,
            proj4string=proj4string,
            verbose=True,
            bounding_polygon=bounding_polygon,
            output_dir=tif_output_dir,
            creation_options=creation_options,
            k_nearest_neighbours=k_nearest_neighbours,

        )

    # Early finish

    if not make_highres_drape_plot:
        return

    # Get extent of geotifs

    sample_rast = glob.glob(tif_output_dir + '*.tif')[0]
    raster_extent = su.getRasterExtent(sample_rast)

    # Make the resampled elevation data

    make_resampled_elevation(
        elevation_raster,
        raster_extent,
        cell_size,
        clip_polygon,
        clip_polygon_layer,
        proj4string,
        tif_output_dir,
    )

    elevation = glob.glob(tif_output_dir + 'LIDAR_resampled*')[0]
    mask = glob.glob(tif_output_dir + 'Mask_resampled*')[0]
    if my_time_step == 'collected_max':
        depth = glob.glob(tif_output_dir + 'PointData_depth_*')[0]
        stage = glob.glob(tif_output_dir + 'PointData_stage_*')[0]
    else:
        depth = glob.glob(tif_output_dir + '*_depth_max.tif')[0]
        stage = glob.glob(tif_output_dir + '*_stage_max.tif')[0]

    # Call gdal_calc

    gdal_calc_command(stage, depth, elevation, mask, depth_threshold)

    return
