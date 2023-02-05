__author__ = 'Allen Zhi Li'
__date__= '2020/06/08'

import numpy as np

from osgeo import gdal

from pyproj import Proj, CRS, transform

from affine import Affine

def tif2array_lat_long(filename, variable_name='elevation',
              easting_min=None, easting_max=None,
              northing_min=None, northing_max=None,
             use_cache=False, verbose=False, proj=None, points=None):
    
    import os


    raster= gdal.Open(filename)
    ncols= raster.RasterXSize
    nrows= raster.RasterYSize
    # x_origin, x_res, _, y_origin, _, y_res= raster.GetGeoTransform()
    NODATA_value= raster.GetRasterBand(1).GetNoDataValue()
    Z= raster.ReadAsArray()
    # treat nan with 0 for now
    Z= np.where(Z==NODATA_value, 0, Z)
    maxRows, maxCols= Z.shape
    try:
        src_georeference= CRS(raster.GetProjection())
    except:
        src_georeference= CRS('EPSG:4326')
    # print src_georeference
    UTM= CRS(proj)
    # print points
    utm_to_84_lons, utm_to_84_lats= transform(UTM,src_georeference,points[:,0], points[:,1])
    transformer= Affine.from_gdal(*raster.GetGeoTransform())
    ilocs= np.array(~ transformer * (utm_to_84_lats,utm_to_84_lons))
    # print utm_to_84_lats, utm_to_84_lons
    icols= ilocs[0,:].astype(int); irows= ilocs[1,:].astype(int)
    # print Z.shape, icols.max(), irows.max()

    #safe return
    # tobe_return= np.zeros(len(points)) * np.nan
    if (icols<maxCols).all() and (irows<maxRows).all():
        return Z[irows, icols]
    elif (icols-3<maxCols).all() and (irows<maxRows).all():
        mask= (icols>=maxCols)
        icols[mask]= maxCols-1
        return Z[irows,icols]
    elif (icols<maxCols).all() and (irows-3<maxRows).all():
        mask= (irows>=maxRows)
        irows[mask]= maxRows-1
        return Z[irows,icols]        
    else:
        #if exceed the boundary by a little bit, then get the nearest neighbors
        msg = 'the input file is inside the boundary, please crop with a larger extent'
        raise ValueError(msg)
