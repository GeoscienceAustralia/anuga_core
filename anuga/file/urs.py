
from struct import pack, unpack
import array as p_array
import numpy as num


from anuga.coordinate_transforms.geo_reference import Geo_reference

from anuga.geospatial_data.geospatial_data import ensure_absolute, \
                                                    Geospatial_data

from anuga.coordinate_transforms.redfearn import redfearn
from anuga.utilities import system_tools 

class Read_urs(object):
    """
    Read the info in URS mux files.

    for the quantities here's a correlation between the file names and
    what they mean;
    z-mux is height above sea level, m
    e-mux is velocity is Eastern direction, m/s
    n-mux is velocity is Northern direction, m/s
    """

    def __init__(self, urs_file):
        self.iterated = False
        columns = 3                         # long, lat , depth
        mux_file = open(urs_file, 'rb')

        # Number of points/stations
        (self.points_num,) = unpack('i', mux_file.read(4))

        # nt, int - Number of time steps
        (self.time_step_count,) = unpack('i', mux_file.read(4))
        #dt, float - time step, seconds
        (self.time_step,) = unpack('f', mux_file.read(4))
        msg = "Bad data in the urs file."
        if self.points_num < 0:
            mux_file.close()
            raise ANUGAError(msg)
        if self.time_step_count < 0:
            mux_file.close()
            raise ANUGAError(msg)
        if self.time_step < 0:
            mux_file.close()
            raise ANUGAError(msg)

        # The depth is in meters, and it is the distance from the ocean
        # to the sea bottom.
        lonlatdep = p_array.array('f')

        if system_tools.major_version == 2:
            # This one works in Python2.7 but not in Python3.8.
            lonlatdep.read(mux_file, columns * self.points_num)
        elif system_tools.major_version == 3:
            # In Python3 we open the file, read it and assign it to the array
            #print(mux_file, columns, self.points_num, columns * self.points_num)
            
            # FIXME (Ole): The number for in probably the item size of a float.
            # This needs to be looked into, but the tests pass now.
            data = mux_file.read(4 * columns * self.points_num)
            #print('data', data, len(data))
            lonlatdep.frombytes(data)
        else:
            raise Exception('Unknown python version: %' % system_tools.version)

        #print('lonlatdep', lonlatdep)
        lonlatdep = num.array(lonlatdep, dtype=float)
        lonlatdep = num.reshape(lonlatdep, (self.points_num, columns))
        self.lonlatdep = lonlatdep

        self.mux_file = mux_file
        # check this array

    def __iter__(self):
        """
        iterate over quantity data which is with respect to time.

        Note: You can only iterate once over an object

        returns quantity infomation for each time slice
        """

        msg =  "You can only interate once over a urs file."
        assert not self.iterated, msg

        self.iter_time_step = 0
        self.iterated = True

        return self

    def __next__(self):
        if self.time_step_count == self.iter_time_step:
            self.close()
            raise StopIteration

        #Read in a time slice from mux file
        hz_p_array = p_array.array('f')

        if system_tools.major_version == 2:
            # This one works in Python2.7 but not in Python3.8.
            hz_p_array.read(self.mux_file, self.points_num)            
        elif system_tools.major_version == 3:
            # In Python3 we open the file, read it and assign it to the array
            #print(mux_file, columns, self.points_num, columns * self.points_num)
            
            # FIXME (Ole): The number for in probably the item size of a float.
            # This needs to be looked into, but the tests pass now.
            data = self.mux_file.read(4 * self.points_num)
            #print('data', data, len(data))
            hz_p_array.frombytes(data) 
        else:
            raise Exception('Unknown python version: %' % system_tools.version)

        
        hz_p = num.array(hz_p_array, dtype=float)
        self.iter_time_step += 1

        return hz_p

    def close(self):
        self.mux_file.close()
   



### PRODUCING THE POINTS NEEDED FILE ###

# Ones used for FESA 2007 results
#LL_LAT = -50.0
#LL_LONG = 80.0
#GRID_SPACING = 1.0/60.0
#LAT_AMOUNT = 4800
#LONG_AMOUNT = 3600


def save_boundary_as_urs(file_name, boundary_polygon, zone,
                              ll_lat, ll_long,
                              grid_spacing,
                              lat_amount, long_amount,
                              isSouthernHemisphere=True,
                              export_csv=False, use_cache=False,
                              verbose=False):
    """
    Given the info to replicate the URS grid and a polygon output
    a file that specifies the cloud of boundary points for URS.

    This creates a .urs file.  This is in the format used by URS;
    1st line is the number of points,
    each line after represents a point,in lats and longs.

    Note: The polygon cannot cross zones or hemispheres.

    A work-a-round for different zones or hemispheres is to run this twice,
    once for each zone, and then combine the output.

    file_name - name of the urs file produced for David.
    boundary_polygon - a list of points that describes a polygon.
                      The last point is assumed ot join the first point.
                      This is in UTM (lat long would be better though)

     This is info about the URS model that needs to be inputted.

    ll_lat - lower left latitude, in decimal degrees
    ll-long - lower left longitude, in decimal degrees
    grid_spacing - in deciamal degrees
    lat_amount - number of latitudes
    long_amount- number of longs

    Don't add the file extension.  It will be added.
    """

    geo = calculate_boundary_points(boundary_polygon, zone, ll_lat, ll_long,
                            grid_spacing,
                            lat_amount, long_amount, isSouthernHemisphere,
                            use_cache, verbose)

    if not file_name[-4:] == ".urs":
        file_name += ".urs"

    geo.export_points_file(file_name, isSouthHemisphere=isSouthernHemisphere)

    if export_csv:
        if file_name[-4:] == ".urs":
            file_name = file_name[:-4] + ".csv"
        geo.export_points_file(file_name)

    return geo


def calculate_boundary_points(boundary_polygon, zone, ll_lat,
                      ll_long, grid_spacing,
                      lat_amount, long_amount, isSouthHemisphere=True,
                      use_cache=False, verbose=False):
    args = (boundary_polygon,
            zone, ll_lat,
            ll_long, grid_spacing,
            lat_amount, long_amount, isSouthHemisphere)
    kwargs = {}

    if use_cache is True:
        try:
            from anuga.caching import cache
        except:
            msg = 'Caching was requested, but caching module' \
                  'could not be imported'
            raise Exception(msg)

        geo = cache(_calculate_boundary_points,
                    args, kwargs,
                    verbose=verbose,
                    compression=False)
    else:
        geo = _calculate_boundary_points(*args, **kwargs)

    return geo


def _calculate_boundary_points(boundary_polygon,
                       zone, ll_lat,
                       ll_long, grid_spacing,
                       lat_amount, long_amount,
                       isSouthHemisphere):
    """
    boundary_polygon - a list of points that describes a polygon.
                      The last point is assumed ot join the first point.
                      This is in UTM (lat long would b better though)

    ll_lat - lower left latitude, in decimal degrees
    ll-long - lower left longitude, in decimal degrees
    grid_spacing - in decimal degrees

    """
    
    msg = "grid_spacing can not be zero"
    assert not grid_spacing == 0, msg

    a = boundary_polygon

    # List of segments.  Each segment is two points.
    segs = [i and [a[i-1], a[i]] or [a[len(a)-1], a[0]] for i in range(len(a))]


    # convert the segs to Lat's and longs.
    # Don't assume the zone of the segments is the same as the lower left
    # corner of the lat long data!!  They can easily be in different zones
    lat_long_set = frozenset()
    for seg in segs:

        points_lat_long = points_needed(seg, ll_lat, ll_long, grid_spacing,
                                        lat_amount, long_amount, zone,
                                        isSouthHemisphere)

        lat_long_set |= frozenset(points_lat_long)


    if lat_long_set == frozenset([]):
        msg = "URS region specified and polygon does not overlap."
        raise ValueError(msg)

    # Warning there is no info in geospatial saying the hemisphere of
    # these points.  There should be.
    geo = Geospatial_data(data_points=list(lat_long_set),
                          points_are_lats_longs=True)

    return geo


def points_needed(seg, ll_lat, ll_long, grid_spacing,
                  lat_amount, long_amount, zone,
                  isSouthHemisphere):
    """
    seg is two points, in UTM
    return a list of the points, in lats and longs that are needed to
    interpolate any point on the segment.
    """

    from math import sqrt

    geo_reference = Geo_reference(zone=zone)
    geo = Geospatial_data(seg, geo_reference=geo_reference)
    seg_lat_long = geo.get_data_points(as_lat_long=True,
                                       isSouthHemisphere=isSouthHemisphere)

    # 1.415 = 2^0.5, rounded up....
    sqrt_2_rounded_up = 1.415
    buffer = sqrt_2_rounded_up * grid_spacing

    max_lat = max(seg_lat_long[0][0], seg_lat_long[1][0]) + buffer
    max_long = max(seg_lat_long[0][1], seg_lat_long[1][1]) + buffer
    min_lat = min(seg_lat_long[0][0], seg_lat_long[1][0]) - buffer
    min_long = min(seg_lat_long[0][1], seg_lat_long[1][1]) - buffer

    first_row = (min_long - ll_long)/ grid_spacing

    # To round up
    first_row_long = int(round(first_row + 0.5))

    last_row = (max_long - ll_long)/ grid_spacing # round down
    last_row_long = int(round(last_row))

    first_row = (min_lat - ll_lat)/ grid_spacing
    # To round up
    first_row_lat = int(round(first_row + 0.5))

    last_row = (max_lat - ll_lat)/ grid_spacing # round down
    last_row_lat = int(round(last_row))

    max_distance = 157147.4112 * grid_spacing
    points_lat_long = []

    # Create a list of the lat long points to include.
    for index_lat in range(first_row_lat, last_row_lat + 1):
        for index_long in range(first_row_long, last_row_long + 1):

            lat = ll_lat + index_lat*grid_spacing
            long = ll_long + index_long*grid_spacing


            #filter here to keep good points
            if keep_point(lat, long, seg, max_distance):
                points_lat_long.append((lat, long)) #must be hashable


    # Now that we have these points, lets throw ones out that are too far away

    return points_lat_long
       


def keep_point(lat, long, seg, max_distance):
    """
    seg is two points, UTM
    """

    from math import sqrt

    _ , x0, y0 = redfearn(lat, long)
    x1 = seg[0][0]
    y1 = seg[0][1]
    x2 = seg[1][0]
    y2 = seg[1][1]
    x2_1 = x2-x1
    y2_1 = y2-y1


    num = (x2_1)*(x2_1)+(y2_1)*(y2_1)
    if sqrt(num) == 0 and abs(num) == 0:
        return True
    else:
        d = abs((x2_1)*(y1-y0)-(x1-x0)*(y2_1))/num
        return d <= max_distance




