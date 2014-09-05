

def xya2pts(filename, verbose=False):
    """Convert a file with xy and elevation data to NetCDF pts file.
    """
    
    from anuga.geospatial_data import Geospatial_data
    
    if verbose: print 'Creating', filename
    
    # Read the ascii (.xya) version of this file,
    # make it comma separated and invert the bathymetry
    # (Below mean sea level should be negative)
    infile = open(filename[:-4] + '.xya')
    out_filename = filename[:-4] + '.pts'

    points = []
    attribute = []
    for line in infile.readlines()[1:]: #Skip first line (the header)
        fields = line.strip().split(',')

        x = float(fields[0])
        y = float(fields[1])
        z = float(fields[2]) # Bathymetry is inverted in original file
        
        points.append([x,y])
        attribute.append(z)
    infile.close()

    # Convert to geospatial data and store as NetCDF
    G = Geospatial_data(data_points=points,
                        attributes=attribute)
    G.export_points_file(out_filename)
    
    