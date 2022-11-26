#!/usr/bin/env python

import numpy as num

celltype_map = {'IEEE4ByteReal': num.float32, 'IEEE8ByteReal': num.float64}


def write_ermapper_grid(ofile, data, header = {}):
    """
    write_ermapper_grid(ofile, data, header = {}):

    Function to write a 2D numeric array to an ERMapper grid.  There are a series of conventions adopted within
    this code, specifically:
    1)  The registration coordinate for the data is the SW (or lower-left) corner of the data
    2)  The registration coordinates refer to cell centres
    3)  The data is a 2D numeric array with the NW-most data in element (0,0) and the SE-most data in element (N,M)
        where N is the last line and M is the last column
    4)  There has been no testng of the use of a rotated grid.  Best to keep data in an NS orientation

    Input Parameters:
    ofile:      string - filename for output (note the output will consist of two files
                ofile and ofile.ers.  Either of these can be entered into this function
    data:       array - 2D array containing the data to be output to the grid
    header:     dictionary - contains spatial information about the grid, in particular:
                    header['datum'] datum for the data ('"GDA94"')
                    header['projection'] - either '"GEOGRAPHIC"' or '"PROJECTED"'
                    header['coordinatetype'] - either 'EN' (for eastings/northings) or
                                                      'LL' (for lat/long)
                    header['rotation'] - rotation of grid ('0:0:0.0')
                    header['celltype'] - data type for writing data ('IEEE4ByteReal')
                    header['nullcellvalue'] - value for null cells ('-99999')
                    header['xdimension'] - cell size in x-dir in units dictated by 'coordinatetype' ('100')
                    header['registrationcellx'] == '0'
                    header['ydimension'] - cell size in y-dir in units dictated by 'coordinatetype' ('100')
                    header['longitude'] - co-ordinate of registration cell ('0:0:0')
                    header['latitude'] - co-ordinate of registration line ('0:0:0')
                    header['nrofbands'] - number of bands ('1')
                    header['value'] - name of grid ('"Default_Band"')

                    Some entries are determined automatically from the data
                    header['nroflines'] - number of lines in data
                    header['nrofcellsperline'] - number of columns in data
                    header['registrationcelly'] == last line of data

    Written by Trevor Dhu, Geoscience Australia 2005
    """
    # extract filenames for header and data files from ofile
    ers_index = ofile.find('.ers')
    if ers_index > 0:
        data_file = ofile[0:ers_index]
        header_file = ofile
    else:
        data_file = ofile
        header_file = ofile + '.ers'


    # Check that the data is a 2 dimensional array
    data_size = num.shape(data)
    assert len(data_size) == 2

    header['nroflines'] = str(data_size[0])
    header['nrofcellsperline'] = str(data_size[1])


    header = create_default_header(header)
    write_ermapper_header(header_file, header)
    write_ermapper_data(data, data_file, data_format = header['celltype'])


def read_ermapper_grid(ifile):
    ers_index = ifile.find('.ers')
    if ers_index > 0:
        data_file = ifile[0:ers_index]
        header_file = ifile
    else:
        data_file = ifile
        header_file = ifile + '.ers'

    header = read_ermapper_header(header_file)

    nroflines = int(header['nroflines'])
    nrofcellsperlines = int(header['nrofcellsperline'])
    data = read_ermapper_data(data_file)
    data = num.reshape(data,(nroflines,nrofcellsperlines))
    return data


def write_ermapper_header(ofile, header = {}):

    header = create_default_header(header)
    # Determine if the dataset is in lats/longs or eastings/northings and set header parameters
    # accordingly
    if header['coordinatetype'] == 'LL':
        X_Class = 'Longitude'
        Y_Class = 'Latitude'
    elif header['coordinatetype'] == 'EN':
        X_Class = 'Eastings'
        Y_Class = 'Northings'

    # open the header file for writing to
    fid = open(ofile,'wt')

    # Begin writing the header
    fid.write('DatasetHeader Begin\n')
    fid.write('\tVersion\t\t= "6.4"\n')
    fid.write('\tDatasetType\t= ERStorage\n')
    fid.write('\tDataType\t= Raster\n')
    fid.write('\tByteOrder\t= LSBFirst\n')

    # Write the coordinate space information
    fid.write('\tCoordinateSpace Begin\n')
    fid.write('\t\tDatum\t\t\t = ' + header['datum'] + '\n')
    fid.write('\t\tProjection\t\t = ' + header['projection'] + '\n')
    fid.write('\t\tCoordinateType\t = ' + header['coordinatetype'] + '\n')
    fid.write('\t\tRotation\t\t = ' + header['rotation'] + '\n')
    fid.write('\t\tUnits\t\t = ' + header['units'] + '\n')
    fid.write('\tCoordinateSpace End\n')

    # Write the raster information
    fid.write('\tRasterInfo Begin\n')
    fid.write('\t\tCellType\t\t\t = ' + header['celltype'] + '\n')
    fid.write('\t\tNullCellValue\t\t = ' + header['nullcellvalue'] + '\n')
    fid.write('\t\tRegistrationCellX\t\t = ' + header['registrationcellx'] +'\n')
    fid.write('\t\tRegistrationCellY\t\t = ' + header['registrationcelly'] +'\n')
    # Write the cellsize information
    fid.write('\t\tCellInfo Begin\n')
    fid.write('\t\t\tXDimension\t\t\t = ' + header['xdimension'] + '\n')
    fid.write('\t\t\tYDimension\t\t\t = ' + header['ydimension'] + '\n')
    fid.write('\t\tCellInfo End\n')
    # Continue with wrting the raster information
    fid.write('\t\tNrOfLines\t\t\t = ' + header['nroflines'] + '\n')
    fid.write('\t\tNrOfCellsPerLine\t = ' + header['nrofcellsperline'] + '\n')
    # Write the registration coordinate information
    fid.write('\t\tRegistrationCoord Begin\n')
    fid.write('\t\t\t' + X_Class + '\t\t\t = ' + header[X_Class.lower()] + '\n')
    fid.write('\t\t\t' + Y_Class + '\t\t\t = ' + header[Y_Class.lower()] + '\n')
    fid.write('\t\tRegistrationCoord End\n')
    # Continue with wrting the raster information
    fid.write('\t\tNrOfBands\t\t\t = ' + header['nrofbands'] + '\n')
    fid.write('\t\tBandID Begin\n')
    fid.write('\t\t\tValue\t\t\t\t = ' + header['value'] + '\n')
    fid.write('\t\tBandID End\n')
    fid.write('\tRasterInfo End\n')
    fid.write('DatasetHeader End\n')

    fid.close()

def read_ermapper_header(ifile):
    # function for reading an ERMapper header from file
    header = {}

    fid = open(ifile,'rt')
    header_string = fid.readlines()
    fid.close()

    for line in header_string:
        if line.find('=') > 0:
            tmp_string = line.strip().split('=')
            header[tmp_string[0].strip().lower()]= tmp_string[1].strip()

    return header

def write_ermapper_data(grid, ofile, data_format=num.float32):


    try:
        data_format = celltype_map[data_format]
    except:
        pass


    #if isinstance(data_format, str):
    #    #celltype_map is defined at top of code
    #    if celltype_map.has_key(data_format):
    #        data_format = celltype_map[data_format]
    #    else:
    #        msg = 'Format %s is not yet defined by celltype_map' %data_format
    #        raise Exception(msg)


    # Convert the array to data_format (default format is Float32)
    grid_as_float = grid.astype(data_format)

    # Convert array to a string for writing to output file
    output_string = grid_as_float.tobytes()

    # open output file in a binary format and write the output string
    fid = open(ofile,'wb')
    fid.write(output_string)
    fid.close()


def read_ermapper_data(ifile, data_format = num.float32):
    # open input file in a binary format and read the input string
    fid = open(ifile,'rb')
    input_string = fid.read()
    fid.close()

    # convert input string to required format (Note default format is num.float32)
    grid_as_float = num.frombuffer(input_string,data_format)
    return grid_as_float

def create_default_header(header = {}):
    # fill any blanks in a header dictionary with default values
    # input parameters:
    # header:   a dictionary containing fields that are not meant
    #           to be filled with default values


    if 'datum' not in header:
        header['datum'] = '"GDA94"'
    if 'projection' not in header:
        header['projection'] = '"GEOGRAPHIC"'
    if 'coordinatetype' not in header:
        header['coordinatetype'] = 'LL'
    if 'rotation' not in header:
        header['rotation'] = '0:0:0.0'
    if 'units' not in header:
        header['units'] = '"METERS"'
    if 'celltype' not in header:
        header['celltype'] = 'IEEE4ByteReal'
    if 'nullcellvalue' not in header:
        header['nullcellvalue'] = '-99999'
    if 'xdimension' not in header:
        header['xdimension'] = '100'
    if 'latitude' not in header:
        header['latitude'] = '0:0:0'
    if 'longitude' not in header:
        header['longitude'] = '0:0:0'
    if 'ydimension' not in header:
        header['ydimension'] = '100'
    if 'nroflines' not in header:
        header['nroflines'] = '3'
    if 'nrofcellsperline' not in header:
        header['nrofcellsperline'] = '4'
    if 'registrationcellx' not in header:
        header['registrationcellx'] = '0'
    if 'registrationcelly' not in header:
        header['registrationcelly'] = str(int(header['nroflines'])-1)
    if 'nrofbands' not in header:
        header['nrofbands'] = '1'
    if 'value' not in header:
        header['value'] = '"Default_Band"'


    return header
