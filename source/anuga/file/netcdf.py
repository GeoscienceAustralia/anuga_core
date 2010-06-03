""" This module is responsible for loading and saving NetCDF NC files
"""

import os, sys
import csv
import string
import shutil
from struct import unpack
import array as p_array
from os import sep, path, remove, mkdir, access, F_OK, W_OK, getcwd

import numpy as num

from Scientific.IO.NetCDF import NetCDFFile
from os.path import exists, basename, join
from os import getcwd

from anuga.coordinate_transforms.redfearn import redfearn, \
     convert_from_latlon_to_utm
from anuga.coordinate_transforms.geo_reference import Geo_reference, \
     write_NetCDF_georeference, ensure_geo_reference
from anuga.geospatial_data.geospatial_data import Geospatial_data,\
     ensure_absolute
from anuga.config import minimum_storable_height as \
     default_minimum_storable_height
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.config import netcdf_float, netcdf_float32, netcdf_int
from anuga.config import max_float
from anuga.utilities.numerical_tools import ensure_numeric,  mean
from anuga.caching.caching import myhash
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.abstract_2d_finite_volumes.pmesh2domain import \
     pmesh_to_domain_instance
from anuga.abstract_2d_finite_volumes.util import get_revision_number, \
     remove_lone_verts, sww2timeseries, get_centroid_values

from anuga.abstract_2d_finite_volumes.neighbour_mesh import segment_midpoints
from anuga.load_mesh.loadASCII import export_mesh_file
from anuga.geometry.polygon import intersection
from anuga.file_conversion.sww2dem import sww2dem

from anuga.utilities.system_tools import get_vars_in_expression
import anuga.utilities.log as log

from anuga.utilities.file_utils import create_filename,\
                        get_all_swwfiles
from anuga.file.csv_file import load_csv_as_dict
from sww_file import Read_sww, Write_sww

from anuga.anuga_exceptions import DataMissingValuesError, \
                DataFileNotOpenError, DataTimeError, DataDomainError, \
                NewQuantity


# Definitions of various NetCDF dimension names, etc.
lon_name = 'LON'
lat_name = 'LAT'
time_name = 'TIME'
precision = netcdf_float # So if we want to change the precision its done here

##
# @brief Clas for a NetCDF data file writer.
class Write_nc:
    """Write an nc file.

    Note, this should be checked to meet cdc netcdf conventions for gridded
    data. http://www.cdc.noaa.gov/cdc/conventions/cdc_netcdf_standard.shtml
    """

    ##
    # @brief Instantiate a Write_nc instance.
    # @param quantity_name 
    # @param file_name 
    # @param time_step_count The number of time steps.
    # @param time_step The time_step size.
    # @param lon 
    # @param lat 
    def __init__(self,
                 quantity_name,
                 file_name,
                 time_step_count,
                 time_step,
                 lon,
                 lat):
        """Instantiate a Write_nc instance (NetCDF file writer).

        time_step_count is the number of time steps.
        time_step is the time step size

        pre-condition: quantity_name must be 'HA', 'UA'or 'VA'.
        """

        self.quantity_name = quantity_name
        quantity_units = {'HA':'CENTIMETERS',
                          'UA':'CENTIMETERS/SECOND',
                          'VA':'CENTIMETERS/SECOND'}

        multiplier_dic = {'HA':100.0,   # To convert from m to cm
                          'UA':100.0,   #             and m/s to cm/sec
                          'VA':-100.0}  # MUX files have positive x in the
                                        # Southern direction.  This corrects
                                        # for it, when writing nc files.

        self.quantity_multiplier =  multiplier_dic[self.quantity_name]

        #self.file_name = file_name
        self.time_step_count = time_step_count
        self.time_step = time_step

        # NetCDF file definition
        self.outfile = NetCDFFile(file_name, netcdf_mode_w)
        outfile = self.outfile

        #Create new file
        nc_lon_lat_header(outfile, lon, lat)

        # TIME
        outfile.createDimension(time_name, None)
        outfile.createVariable(time_name, precision, (time_name,))

        #QUANTITY
        outfile.createVariable(self.quantity_name, precision,
                               (time_name, lat_name, lon_name))
        outfile.variables[self.quantity_name].missing_value = -1.e+034
        outfile.variables[self.quantity_name].units = \
                                 quantity_units[self.quantity_name]
        outfile.variables[lon_name][:]= ensure_numeric(lon)
        outfile.variables[lat_name][:]= ensure_numeric(lat)

        #Assume no one will be wanting to read this, while we are writing
        #outfile.close()

    ##
    # @brief Write a time-step of quantity data.
    # @param quantity_slice The data to be stored for this time-step.
    def store_timestep(self, quantity_slice):
        """Write a time slice of quantity info

        quantity_slice is the data to be stored at this time step
        """

        # Get the variables
        time = self.outfile.variables[time_name]
        quantity = self.outfile.variables[self.quantity_name]

        # get index oflice to write
        i = len(time)

        #Store time
        time[i] = i * self.time_step    #self.domain.time
        quantity[i,:] = quantity_slice * self.quantity_multiplier

    ##
    # @brief Close file underlying the class instance.
    def close(self):
        self.outfile.close()



##
# @brief Write an NC elevation file.
# @param file_out Path to the output file.
# @param lon ??
# @param lat ??
# @param depth_vector The elevation data to write.
def write_elevation_nc(file_out, lon, lat, depth_vector):
    """Write an nc elevation file."""

    # NetCDF file definition
    outfile = NetCDFFile(file_out, netcdf_mode_w)

    #Create new file
    nc_lon_lat_header(outfile, lon, lat)

    # ELEVATION
    zname = 'ELEVATION'
    outfile.createVariable(zname, precision, (lat_name, lon_name))
    outfile.variables[zname].units = 'CENTIMETERS'
    outfile.variables[zname].missing_value = -1.e+034

    outfile.variables[lon_name][:] = ensure_numeric(lon)
    outfile.variables[lat_name][:] = ensure_numeric(lat)

    depth = num.reshape(depth_vector, (len(lat), len(lon)))
    outfile.variables[zname][:] = depth

    outfile.close()


##
# @brief Write lat/lon headers to a NetCDF file.
# @param outfile Handle to open file to write to.
# @param lon An iterable of the longitudes.
# @param lat An iterable of the latitudes.
# @note Defines lat/long dimensions and variables. Sets various attributes:
#          .point_spacing  and  .units
#       and writes lat/lon data.

def nc_lon_lat_header(outfile, lon, lat):
    """Write lat/lon headers to a NetCDF file.

    outfile is the netcdf file handle.
    lon - a list/array of the longitudes
    lat - a list/array of the latitudes
    """

    outfile.institution = 'Geoscience Australia'
    outfile.description = 'Converted from URS binary C'

    # Longitude
    outfile.createDimension(lon_name, len(lon))
    outfile.createVariable(lon_name, precision, (lon_name,))
    outfile.variables[lon_name].point_spacing = 'uneven'
    outfile.variables[lon_name].units = 'degrees_east'
    outfile.variables[lon_name].assignValue(lon)

    # Latitude
    outfile.createDimension(lat_name, len(lat))
    outfile.createVariable(lat_name, precision, (lat_name,))
    outfile.variables[lat_name].point_spacing = 'uneven'
    outfile.variables[lat_name].units = 'degrees_north'
    outfile.variables[lat_name].assignValue(lat)




##
# @brief Filter data file, selecting timesteps first:step:last.
# @param filename1 Data file to filter.
# @param filename2 File to write filtered timesteps to.
# @param first First timestep.
# @param last Last timestep.
# @param step Timestep stride.
def filter_netcdf(filename1, filename2, first=0, last=None, step=1):
    """Filter data file, selecting timesteps first:step:last.
    
    Read netcdf filename1, pick timesteps first:step:last and save to
    nettcdf file filename2
    """

    from Scientific.IO.NetCDF import NetCDFFile

    # Get NetCDF
    infile = NetCDFFile(filename1, netcdf_mode_r)  #Open existing file for read
    outfile = NetCDFFile(filename2, netcdf_mode_w)  #Open new file

    # Copy dimensions
    for d in infile.dimensions:
        outfile.createDimension(d, infile.dimensions[d])

    # Copy variable definitions
    for name in infile.variables:
        var = infile.variables[name]
        outfile.createVariable(name, var.dtype.char, var.dimensions)

    # Copy the static variables
    for name in infile.variables:
        if name == 'time' or name == 'stage':
            pass
        else:
            outfile.variables[name][:] = infile.variables[name][:]

    # Copy selected timesteps
    time = infile.variables['time']
    stage = infile.variables['stage']

    newtime = outfile.variables['time']
    newstage = outfile.variables['stage']

    if last is None:
        last = len(time)

    selection = range(first, last, step)
    for i, j in enumerate(selection):
        log.critical('Copying timestep %d of %d (%f)'
                     % (j, last-first, time[j]))
        newtime[i] = time[j]
        newstage[i,:] = stage[j,:]

    # Close
    infile.close()
    outfile.close()

