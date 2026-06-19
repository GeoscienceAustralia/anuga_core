__author__ = 'Allen Zhi Li'
__date__= '2020/06/08'

# Adapted by Stephen Roberts 2023
# Migrated from osgeo/gdal to rasterio/pyproj 2026

from pprint import pprint

def tif2point_values(filename, zone=None, south=True, points=None, verbose=False):

    import numpy as np
    import rasterio
    from pyproj import CRS, Transformer

    with rasterio.open(filename) as raster:
        ncols = raster.width
        nrows = raster.height
        tif_crs = raster.crs
        tif_epsg = str(tif_crs.to_epsg()) if tif_crs is not None else None
        NODATA_value = raster.nodata
        Z = raster.read(1)
        affine_transform = raster.transform

    # treat nodata with 0 for now
    if NODATA_value is not None:
        Z = np.where(Z == NODATA_value, 0, Z)
    maxRows, maxCols = Z.shape

    # CRS for input points assumed UTM defined by zone and whether south or not
    points_utm = CRS.from_dict({'proj': 'utm', 'zone': zone, 'south': south})

    wgs84_utm_north = {zone: 32600 + zone for zone in range(1, 61)}
    wgs84_utm_south = {zone: 32700 + zone for zone in range(1, 61)}

    nad83_utm_north = {
        10: 26910, 11: 26911, 12: 26912, 13: 26913,
        14: 26914, 15: 26915, 16: 26916, 17: 26917,
        18: 26918, 19: 26919, 20: 26920, 21: 26921,
        22: 26922, 23: 26923,
    }

    if tif_epsg == '4326':
        # tif file is lat long projection ie 'EPSG:4326'
        tif_georeference = CRS.from_epsg(4326)

        transformer = Transformer.from_crs(points_utm, tif_georeference)
        # pyproj dispatches to _transform_point (scalar path) when given a
        # 1-element array, which fails with numpy >= 1.25.  Pass plain Python
        # lists so pyproj always uses the array path, then convert back.
        _lat, _lon = transformer.transform(
            points[:, 0].tolist(), points[:, 1].tolist())
        points_lat = np.asarray(_lat)
        points_lon = np.asarray(_lon)

        ilocs = np.array(~affine_transform * (points_lon, points_lat))

    elif not south:
        zone = int(zone)
        same_utm = (
            tif_epsg == str(wgs84_utm_north.get(zone)) or
            tif_epsg == str(nad83_utm_north.get(zone))
        )
        if same_utm:
            ilocs = np.array(~affine_transform * (points[:, 0], points[:, 1]))
        else:
            raise Exception("zone and hemisphere of tif not the same as zone and hemisphere of points")

    elif south:
        zone = int(zone)
        same_utm = tif_epsg == str(wgs84_utm_south.get(zone))
        if same_utm:
            ilocs = np.array(~affine_transform * (points[:, 0], points[:, 1]))
        else:
            raise Exception("zone and hemisphere of tif not the same as zone and hemisphere of points")

    else:
        msg = 'zone and hemisphere of tif not the same as zone and hemisphere of points'
        raise Exception(msg)

    icols = ilocs[0, :].astype(int)
    irows = ilocs[1, :].astype(int)

    if (icols < maxCols).all() and (irows < maxRows).all():
        return Z[irows, icols]
    elif (icols - 3 < maxCols).all() and (irows < maxRows).all():
        mask = (icols >= maxCols)
        icols[mask] = maxCols - 1
        return Z[irows, icols]
    elif (icols < maxCols).all() and (irows - 3 < maxRows).all():
        mask = (irows >= maxRows)
        irows[mask] = maxRows - 1
        return Z[irows, icols]
    else:
        msg = 'points outside the extent of the source tif file, please crop tif file with a larger extent'
        raise ValueError(msg)
