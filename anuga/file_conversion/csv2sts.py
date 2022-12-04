#!/usr/bin/env python

"""
    Module to convert a .csv file to an .sts file.
    
    You can use the function from your own scripts, or call this script on the 
    command line.
    
    Use this module to convert an arbitrary csv file to an sts file.
    The csv file must have a header containing the column names. This will
    be converted to a NetCDF .sts file containing the same data, with the
    addition of optional latitude and longitude values.

    For example, the following .csv file:
    
        time stage
        0 4
        1 150.66667
        2 150.83334
        3 151.
        4 151.16667
        5 -34.
        6 -34.16667
        7 -34.33333
        8 -34.5
        9 -1.
        10 -5.
        11 -9.
        12 -13.
    
    Using this command:
        python csv2sts --lat 14 --lon 56 infile.csv foo.sts
    Will be converted to the following .sts file:

        netcdf foo {
        dimensions:
            number_of_timesteps = 13 ;
        variables:
            double stage(number_of_timesteps) ;
            double time(number_of_timesteps) ;

        // global attributes:
                :latitude = 14. ;
                :longitude = 56. ;
        data:

         stage = 4, 150.66667, 150.83334, 151, 151.16667, -34, -34.16667,
                -34.33333, -34.5, -1, -5, -9, -13 ;

         time = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 ;
        }
    
    As of June 2010 this module has a pylint quality rating of 9.30/10.
"""


import sys
import getopt
from anuga.utilities import log
from anuga.file.netcdf import NetCDFFile
from anuga.file.csv_file import load_csv_as_dict
from anuga.config import netcdf_mode_w, netcdf_float


def csv2sts(infile, outfile, latitude = None, longitude = None,
                    verbose = False):
    """
        Take a csv file and convert it to an sts file.
        
        May be used for timeseries, or any other data.
    """
        
    timeseries_data, col_names = load_csv_as_dict(infile, delimiter=' ')
    
    if not col_names:
        raise IOError('csv2sts: file %s is empty or unreadable.' % infile)
    
    if verbose:
        log.critical('csv2sts input data:')
        for col in col_names:
            log.critical('column ' + col + ':')
            log.critical(timeseries_data[col])        

    data_len = len(list(timeseries_data.values())[0])
    if verbose:
        log.critical('   data length = %d.' % data_len)
    
    fid = NetCDFFile(outfile, netcdf_mode_w)

    fid.createDimension('number_of_timesteps', data_len)

    if latitude:
        fid.latitude = latitude
        
    if longitude:
        fid.longitude = longitude
    
    for col in col_names:
        fid.createVariable(col, netcdf_float, ('number_of_timesteps',))
        
        fid.variables[col][:] = timeseries_data[col]

    fid.close()

                 

######
# Script is being run from command line.
#

def usage():
    """ Display usage of this module from the comand line. """
    print('csv2sts - convert a csv file to an sts file.')
    print('Usage: csv2sts [-hv] [--help] [--verbose]', end=' ')
    print('[-x --lat --latitude <degrees>]', end=' ')
    print('[-y --lon --longitude <degrees>] <in.csv> <out.sts>')
    print('eg:')
    print('python csv2sts.py -v --lat 10 --lon 20 infile.csv sts_out.sts')
    print()

def main(argv):                         
    """ Script is being run from the command line. """
    lat = None
    lon = None
    verbose = False
    
    try:                                
        long_parms = ["help", "verbose", \
                        "lat=", "lon=", "latitude=", "longitude="]
        opts, args = getopt.getopt(argv, "hvx:y:", long_parms)
    except getopt.GetoptError:   
        usage()     
        sys.exit(2)
    for opt, arg in opts:    
        if opt in ("-h", "--help"):
            usage()                     
            sys.exit()          
        elif opt in ("-x", "--lat", "--latitude"):
            lat = float(arg) 
        elif opt in ("-y", "--lon", "--longitude"):
            lon = float(arg)
        if opt in ("-v", "--verbose"):
            verbose = True
                            
    if len(args) != 2:
        usage()     
        sys.exit(2)        
    
    infile = args[0]    
    outfile = args[1]
    
    if verbose:
        msg = 'csv2sts: converting %s to %s' % (infile, outfile)
        if lat and lon:
            msg += ' with lat = %d, lon = %d...' % (lat, lon)
        print(msg)
    csv2sts(infile, outfile, lat, lon)
    
    if verbose:
        print('done!')



    
if __name__ == "__main__":
    main(sys.argv[1:])       
