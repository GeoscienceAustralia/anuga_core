from builtins import map
from builtins import str
from anuga.file.urs import Read_urs

def urs2txt(basename_in, location_index=None):
    """
    Not finished or tested
    """

    files_in = [basename_in + WAVEHEIGHT_MUX_LABEL,
                basename_in + EAST_VELOCITY_LABEL,
                basename_in + NORTH_VELOCITY_LABEL]
    quantities = ['HA','UA','VA']

    d = ","

    # instantiate urs_points of the three mux files.
    mux = {}
    for quantity, file in zip(quantities, files_in):
        mux[quantity] = Read_urs(file)

    # Could check that the depth is the same. (hashing)

    # handle to a mux file to do depth stuff
    a_mux = mux[quantities[0]]

    # Convert to utm
    latitudes = a_mux.lonlatdep[:,1]
    longitudes = a_mux.lonlatdep[:,0]
    points_utm, zone = \
        convert_from_latlon_to_utm(latitudes=latitudes, longitudes=longitudes)
    depths = a_mux.lonlatdep[:,2]

    # open the output text file, start writing.
    fid = open(basename_in + '.txt', 'w')

    fid.write("zone: " + str(zone) + "\n")

    if location_index is not None:
        #Title
        li = location_index
        fid.write('location_index' + d + 'lat' + d + 'long' + d +
                  'Easting' + d + 'Northing' + '\n')
        fid.write(str(li) + d + str(latitudes[li]) + d +
                  str(longitudes[li]) + d + str(points_utm[li][0]) + d +
                  str(points_utm[li][0o1]) + '\n')

    # the non-time dependent stuff
    #Title
    fid.write('location_index' + d + 'lat' + d + 'long' + d +
              'Easting' + d + 'Northing' + d + 'depth m' + '\n')
    i = 0
    for depth, point_utm, lat, long in zip(depths, points_utm,
                                           latitudes, longitudes):

        fid.write(str(i) + d + str(lat) + d + str(long) + d +
                  str(point_utm[0]) + d + str(point_utm[0o1]) + d +
                  str(depth) + '\n')
        i += 1

    #Time dependent
    if location_index is not None:
        time_step = a_mux.time_step
        i = 0
        #Title
        fid.write('time' + d + 'HA depth m' + d + 'UA momentum East x m/sec' +
                  d + 'VA momentum North y m/sec' + '\n')
        for HA, UA, VA in zip(mux['HA'], mux['UA'], mux['VA']):
            fid.write(str(i*time_step) + d + str(HA[location_index]) + d +
                      str(UA[location_index]) + d +
                      str(VA[location_index]) + '\n')
            i += 1

