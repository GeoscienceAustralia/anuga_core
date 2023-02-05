import numpy as np

from osgeo import gdal

def tif2array(filename, variable_name='elevation',
              easting_min=None, easting_max=None,
              northing_min=None, northing_max=None,
             use_cache=False, verbose=False,):
    
    import os
    raster= gdal.Open(filename)
    ncols= raster.RasterXSize
    nrows= raster.RasterYSize
    x_origin, x_res, _, y_origin, _, y_res= raster.GetGeoTransform()
    NODATA_value= raster.GetRasterBand(1).GetNoDataValue()
    Z= raster.ReadAsArray()
    Z= np.where(Z==NODATA_value, np.nan, Z)
    
    if y_res<0:
        x= np.linspace(x_origin, x_origin+(ncols-1)*x_res, ncols)
        y= np.linspace(y_origin+(nrows-1)*y_res, y_origin, nrows)
        Z= np.flip(Z, axis=0)
        Z= Z.transpose()
    elif y_res>=0:
        x= np.linspace(x_origin, x_origin+(ncols-1)*x_res, ncols)
        y= np.linspace(y_origin, y_origin+(nrows-1)*y_res, nrows)
        Z= Z.transpose()
    
    
    return x, y, Z
    
    