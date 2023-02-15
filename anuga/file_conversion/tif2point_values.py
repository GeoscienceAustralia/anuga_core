__author__ = 'Allen Zhi Li'
__date__= '2020/06/08'

# Adapted by Stephen Roberts 2023

import numpy as np
from osgeo import gdal
from pyproj import Proj, CRS, transform
from affine import Affine

def tif2point_values(filename, zone=None, south=True, points=None, verbose=False):

    raster= gdal.Open(filename)
    ncols= raster.RasterXSize
    nrows= raster.RasterYSize
 
    NODATA_value= raster.GetRasterBand(1).GetNoDataValue()
    Z= raster.ReadAsArray()
    # treat nan with 0 for now
    Z= np.where(Z==NODATA_value, 0, Z)
    maxRows, maxCols= Z.shape

    # Would expect tif file to be lat long projection ie 'EPSG:4326'
    src_georeference= CRS(raster.GetProjection())
    #print(src_georeference)
    
    if zone == -1 :
        raise Exception('Need to specify zone for domain')
        
    UTM = CRS.from_dict({'proj': 'utm', 'zone': zone, 'south': south})

    #points_lons, points_lats= transform(UTM, src_georeference, points[:,0], points[:,1])

    from pyproj import Transformer
    transformer = Transformer.from_crs(UTM, src_georeference)
    points_lons, points_lats = transformer.transform(points[:,0], points[:,1])

    import numpy
    #assert numpy.allclose(points_lons,points_lons_1)
    #assert numpy.allclose(points_lats,points_lats_1)

    transformer= Affine.from_gdal(*raster.GetGeoTransform())
    ilocs= np.array(~ transformer * (points_lats,points_lons))

    icols= ilocs[0,:].astype(int); irows= ilocs[1,:].astype(int)

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
        msg = 'the input file is inside the boundary, please crop tif file with a larger extent'
        raise ValueError(msg)
