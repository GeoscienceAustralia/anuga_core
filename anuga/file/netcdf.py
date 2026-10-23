"""Read and write NetCDF NC files (netCDF4 only)."""
import numpy as num

from anuga.coordinate_transforms.redfearn import \
     convert_from_latlon_to_utm
from anuga.config import minimum_storable_height as \
     default_minimum_storable_height
from anuga.config import netcdf_mode_r, netcdf_mode_w
from anuga.config import netcdf_float, netcdf_int
from anuga.utilities.numerical_tools import ensure_numeric, mean

import anuga.utilities.log as log


# Definitions of various NetCDF dimension names, etc.
lon_name = 'LON'
lat_name = 'LAT'
time_name = 'TIME'
precision = netcdf_float


def NetCDFFile(file_name, netcdf_mode=netcdf_mode_r):
    """Return a netCDF4 Dataset handle for *file_name*.

    Parameters
    ----------
    file_name : str
        Path to the NetCDF file.
    netcdf_mode : str
        Open mode: ``'r'`` (read), ``'w'`` (write), ``'a'`` (append), or the
        legacy alias ``'wl'`` which maps to write mode with NETCDF3_64BIT format.

    Returns
    -------
    netCDF4.Dataset
    """
    from netCDF4 import Dataset
    if netcdf_mode == 'wl':
        return Dataset(file_name, 'w', format='NETCDF3_64BIT')
    return Dataset(file_name, netcdf_mode, format='NETCDF3_64BIT')



class Write_nc:
    """Write an nc file.

    Note, this should be checked to meet cdc netcdf conventions for gridded
    data. http://www.cdc.noaa.gov/cdc/conventions/cdc_netcdf_standard.shtml
    """

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

    def store_timestep(self, quantity_slice):
        """Write a time slice of quantity info

        quantity_slice is the data to be stored at this time step
        """

        # Get the variables
        time = self.outfile.variables[time_name]
        quantity = self.outfile.variables[self.quantity_name]

        # get index of slice to write
        i = len(time)

        #Store time
        time[i] = i * self.time_step
        quantity[i,:] = quantity_slice * self.quantity_multiplier

    def close(self):
        """ Close file underlying the class instance. """
        self.outfile.close()



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

    depth = num.reshape(depth_vector, (len(lat), len(lon)))
    outfile.variables[zname][:] = depth

    outfile.close()


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
    outfile.variables[lon_name][:] = ensure_numeric(lon)

    # Latitude
    outfile.createDimension(lat_name, len(lat))
    outfile.createVariable(lat_name, precision, (lat_name,))
    outfile.variables[lat_name].point_spacing = 'uneven'
    outfile.variables[lat_name].units = 'degrees_north'
    outfile.variables[lat_name][:] = ensure_numeric(lat)
