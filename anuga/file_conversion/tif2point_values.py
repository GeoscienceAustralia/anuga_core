__author__ = 'Allen Zhi Li'
__date__= '2020/06/08'

# Adapted by Stephen Roberts 2023

from pprint import pprint

def tif2point_values(filename, zone=None, south=True, points=None, verbose=False):

    import numpy as np
    from osgeo import gdal, osr
    from pyproj import Proj, CRS, transform
    from affine import Affine

    raster= gdal.Open(filename)
    ncols= raster.RasterXSize
    nrows= raster.RasterYSize


    tif_proj = osr.SpatialReference(wkt=raster.GetProjection())
    #print(proj.GetAttrValue('AUTHORITY',1))
    tif_epsg = tif_proj.GetAuthorityCode(None)

    #print(tif_epsg)
 
    NODATA_value= raster.GetRasterBand(1).GetNoDataValue()
    Z= raster.ReadAsArray()
    # treat nan with 0 for now
    Z= np.where(Z==NODATA_value, 0, Z)
    maxRows, maxCols= Z.shape

    # CRS for input points assumed UTM defined by zone and whether south or not
    points_utm = CRS.from_dict({'proj': 'utm', 'zone': zone, 'south': south})
    #print(points_utm)
    #points_lons, points_lats= transform(UTM, src_georeference, points[:,0], points[:,1])

    #print(tif_epsg)
    #print(south)

    if tif_epsg == '4326':
        # tif file is lat long projection ie 'EPSG:4326'
        tif_georeference= CRS(raster.GetProjection())
        #print(tif_georeference)

        from pyproj import Transformer
        transformer = Transformer.from_crs(points_utm, tif_georeference)
        points_lat, points_lon = transformer.transform(points[:,0], points[:,1])

        #print(transformer)
        #pprint(points)
        #print('points_lat')
        #pprint(points_lat)
        #print('points_lon')
        #pprint(points_lon)

        import numpy
        #assert numpy.allclose(points_lons,points_lons_1)
        #assert numpy.allclose(points_lats,points_lats_1)

        affine_transform= Affine.from_gdal(*raster.GetGeoTransform())
        #print(affine_transform)
        #print(~affine_transform)
        #print(affine_transform * (points_lon,points_lat))
        #print(~ affine_transform * (points_lon,points_lat))

        ilocs= np.array(~ affine_transform * (points_lon,points_lat))

    elif (tif_epsg == str(32600 + int(zone))) and not south:
        # no need for transformation
        affine_transform= Affine.from_gdal(*raster.GetGeoTransform())
        ilocs= np.array(~ affine_transform * (points[:,0],points[:,1]))

    elif (tif_epsg == str(32700 + int(zone))) and south:
        # no need for transformation
        affine_transform= Affine.from_gdal(*raster.GetGeoTransform())
        ilocs= np.array(~ affine_transform * (points[:,0],points[:,1]))

    elif (tif_epsg == str(7800 + int(zone)))  and south:
        # no need for transformation
        affine_transform= Affine.from_gdal(*raster.GetGeoTransform())
        ilocs= np.array(~ affine_transform * (points[:,0],points[:,1]))
        
    else:
        msg = 'zone and hemisphere of tif not the same as zone and hemisphere of points'
        raise Exception(msg)

    icols= ilocs[0,:].astype(int); irows= ilocs[1,:].astype(int)

    #pprint(icols)
    #pprint(irows)

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
        msg = 'points outside the extent of the source tif file, please crop tif file with a larger extent'
        raise ValueError(msg)
