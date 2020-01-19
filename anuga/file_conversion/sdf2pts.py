
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.geospatial_data.geospatial_data import Geospatial_data

def _read_hecras_cross_sections(lines):
    """Return block of surface lines for each cross section
    Starts with SURFACE LINE,
    Ends with END CROSS-SECTION
    """

    points = []

    """
    Test comment
    """
    reading_surface = False
    for i, line in enumerate(lines):
        if len(line.strip()) == 0:    #Ignore blanks
            continue

        if lines[i].strip().startswith('SURFACE LINE'):
            reading_surface = True
            continue

        if lines[i].strip().startswith('END') and reading_surface:
            yield points
            reading_surface = False
            points = []

        if reading_surface:
            fields = line.strip().split(',')
            easting = float(fields[0])
            northing = float(fields[1])
            elevation = float(fields[2])
            points.append([easting, northing, elevation])


def sdf2pts(name_in, name_out=None, verbose=False):
    """
    Read HEC-RAS Elevation datal from the following ASCII format (.sdf)

    basename_in Sterm of input filename.
    basename_out Sterm of output filename.
    verbose True if this function is to be verbose.

    Example:

# RAS export file created on Mon 15Aug2005 11:42
# by HEC-RAS Version 3.1.1

BEGIN HEADER:
  UNITS: METRIC
  DTM TYPE: TIN
  DTM: v:\\1\\cit\\perth_topo\\river_tin
  STREAM LAYER: c:\\local\\hecras\\21_02_03\\up_canning_cent3d.shp
  CROSS-SECTION LAYER: c:\\local\\hecras\\21_02_03\\up_can_xs3d.shp
  MAP PROJECTION: UTM
  PROJECTION ZONE: 50
  DATUM: AGD66
  VERTICAL DATUM:
  NUMBER OF REACHES:  19
  NUMBER OF CROSS-SECTIONS:  14206
END HEADER:

Only the SURFACE LINE data of the following form will be utilised
  CROSS-SECTION:
    STREAM ID:Southern-Wungong
    REACH ID:Southern-Wungong
    STATION:19040.*
    CUT LINE:
      405548.671603161 , 6438142.7594925
      405734.536092045 , 6438326.10404912
      405745.130459356 , 6438331.48627354
      405813.89633823 , 6438368.6272789
    SURFACE LINE:
     405548.67,   6438142.76,   35.37
     405552.24,   6438146.28,   35.41
     405554.78,   6438148.78,   35.44
     405555.80,   6438149.79,   35.44
     405559.37,   6438153.31,   35.45
     405560.88,   6438154.81,   35.44
     405562.93,   6438156.83,   35.42
     405566.50,   6438160.35,   35.38
     405566.99,   6438160.83,   35.37
     ...
   END CROSS-SECTION

    Convert to NetCDF pts format which is

    points:  (Nx2) float array
    elevation: N float array
    """

    import os
    from anuga.file.netcdf import NetCDFFile

    if name_in[-4:] != '.sdf':
        raise IOError('Input file %s should be of type .sdf.' % name_in)

    if name_out is None:
        name_out = name_in[:-4] + '.pts'
    elif name_out[-4:] != '.pts':
        raise IOError('Input file %s should be of type .pts.' % name_out)

    # Get ASCII file
    infile = open(name_in, 'r')

    if verbose: log.critical('Reading DEM from %s' % (root + '.sdf'))

    lines = infile.readlines()
    infile.close()

    if verbose: log.critical('Converting to pts format')

    # Scan through the header, picking up stuff we need.
    i = 0
    while lines[i].strip() == '' or lines[i].strip().startswith('#'):
        i += 1

    assert lines[i].strip().upper() == 'BEGIN HEADER:'
    i += 1

    assert lines[i].strip().upper().startswith('UNITS:')
    units = lines[i].strip().split()[1]
    i += 1

    assert lines[i].strip().upper().startswith('DTM TYPE:')
    i += 1

    assert lines[i].strip().upper().startswith('DTM:')
    i += 1

    assert lines[i].strip().upper().startswith('STREAM')
    i += 1

    assert lines[i].strip().upper().startswith('CROSS')
    i += 1

    assert lines[i].strip().upper().startswith('MAP PROJECTION:')
    projection = lines[i].strip().split(':')[1]
    i += 1

    assert lines[i].strip().upper().startswith('PROJECTION ZONE:')
    zone = int(lines[i].strip().split(':')[1])
    i += 1

    assert lines[i].strip().upper().startswith('DATUM:')
    datum = lines[i].strip().split(':')[1]
    i += 1

    assert lines[i].strip().upper().startswith('VERTICAL DATUM:')
    i += 1

    assert lines[i].strip().upper().startswith('NUMBER OF REACHES:')
    i += 1

    assert lines[i].strip().upper().startswith('NUMBER OF CROSS-SECTIONS:')
    number_of_cross_sections = int(lines[i].strip().split(':')[1])
    i += 1

    # Now read all points
    points = []
    elevation = []
    for j, entries in enumerate(_read_hecras_cross_sections(lines[i:])):
        for k, entry in enumerate(entries):
            points.append(entry[:2])
            elevation.append(entry[2])

    msg = 'Actual #number_of_cross_sections == %d, Reported as %d'\
          %(j+1, number_of_cross_sections)
    assert j+1 == number_of_cross_sections, msg

    # Get output file, write PTS data
    if name_out is None:
        ptsname = name_in[:-4] + '.pts'
    else:
        ptsname = name_out

    geo_ref = Geo_reference(zone, 0, 0, datum, projection, units)
    geo = Geospatial_data(points, {"elevation":elevation},
                          verbose=verbose, geo_reference=geo_ref)
    geo.export_points_file(ptsname)
