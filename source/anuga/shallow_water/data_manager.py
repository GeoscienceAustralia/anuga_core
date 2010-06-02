"""datamanager.py - input output for AnuGA


This module takes care of reading and writing datafiles such as topograhies,
model output, etc


Formats used within AnuGA:

.sww: Netcdf format for storing model output f(t,x,y)
.tms: Netcdf format for storing time series f(t)

.csv: ASCII format for storing arbitrary points and associated attributes
.pts: NetCDF format for storing arbitrary points and associated attributes

.asc: ASCII format of regular DEMs as output from ArcView
.prj: Associated ArcView file giving more meta data for asc format
.ers: ERMapper header format of regular DEMs for ArcView

.dem: NetCDF representation of regular DEM data

.tsh: ASCII format for storing meshes and associated boundary and region info
.msh: NetCDF format for storing meshes and associated boundary and region info

.nc: Native ferret NetCDF format
.geo: Houdinis ascii geometry format (?)


A typical dataflow can be described as follows

Manually created files:
ASC, PRJ:     Digital elevation models (gridded)
TSH:          Triangular meshes (e.g. created from anuga.pmesh)
NC            Model outputs for use as boundary conditions (e.g from MOST)


AUTOMATICALLY CREATED FILES:

ASC, PRJ  ->  DEM  ->  PTS: Conversion of DEM's to native pts file

NC -> SWW: Conversion of MOST bundary files to boundary sww

PTS + TSH -> TSH with elevation: Least squares fit

TSH -> SWW:  Conversion of TSH to sww viewable using Swollen

TSH + Boundary SWW -> SWW: Simluation using abstract_2d_finite_volumes

"""

# This file was reverted from changeset:5484 to changeset:5470 on 10th July
# by Ole.

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



def csv2building_polygons(file_name,
                          floor_height=3,
                          clipping_polygons=None):
    """
    Convert CSV files of the form:

    easting,northing,id,floors
    422664.22,870785.46,2,0
    422672.48,870780.14,2,0
    422668.17,870772.62,2,0
    422660.35,870777.17,2,0
    422664.22,870785.46,2,0
    422661.30,871215.06,3,1
    422667.50,871215.70,3,1
    422668.30,871204.86,3,1
    422662.21,871204.33,3,1
    422661.30,871215.06,3,1

    to a dictionary of polygons with id as key.
    The associated number of floors are converted to m above MSL and 
    returned as a separate dictionary also keyed by id.
    
    Optional parameter floor_height is the height of each building story.
    Optional parameter clipping_olygons is a list of polygons selecting
    buildings. Any building not in these polygons will be omitted.
    
    See csv2polygons for more details
    """

    polygons, values = csv2polygons(file_name,
                                    value_name='floors',
                                    clipping_polygons=None)    

    
    heights = {}
    for key in values.keys():
        v = float(values[key])
        heights[key] = v*floor_height
        
    return polygons, heights                
            

##
# @brief Convert CSV file into a dictionary of polygons and associated values.
# @param filename The path to the file to read, value_name name for the 4th column
def csv2polygons(file_name,
                 value_name='value',
                 clipping_polygons=None):
    """
    Convert CSV files of the form:

    easting,northing,id,value
    422664.22,870785.46,2,0
    422672.48,870780.14,2,0
    422668.17,870772.62,2,0
    422660.35,870777.17,2,0
    422664.22,870785.46,2,0
    422661.30,871215.06,3,1
    422667.50,871215.70,3,1
    422668.30,871204.86,3,1
    422662.21,871204.33,3,1
    422661.30,871215.06,3,1

    to a dictionary of polygons with id as key.
    The associated values are returned as a separate dictionary also keyed by id.


    easting: x coordinate relative to zone implied by the model
    northing: y coordinate relative to zone implied by the model    
    id: tag for polygon comprising points with this tag
    value: numeral associated with each polygon. These must be the same for all points in each polygon.
   
    The last header, value, can take on other names such as roughness, floors, etc - or it can be omitted 
    in which case the returned values will be None
    
    Eastings and Northings will be returned as floating point values while
    id and values will be returned as strings.

    Optional argument: clipping_polygons will select only those polygons that are
    fully within one or more of the clipping_polygons. In other words any polygon from
    the csv file which has at least one point not inside one of the clipping polygons
    will be excluded 
    
    See underlying function load_csv_as_dict for more details.
    """

    X, _ = load_csv_as_dict(file_name)

    msg = 'Polygon csv file must have 3 or 4 columns'
    assert len(X.keys()) in [3, 4], msg
    
    msg = 'Did not find expected column header: easting'
    assert 'easting' in X.keys(), msg
    
    msg = 'Did not find expected column header: northing'    
    assert 'northing' in X.keys(), northing
    
    msg = 'Did not find expected column header: northing'        
    assert 'id' in X.keys(), msg
    
    if value_name is not None:
        msg = 'Did not find expected column header: %s' % value_name        
        assert value_name in X.keys(), msg    
    
    polygons = {}
    if len(X.keys()) == 4:
        values = {}
    else:
        values = None

    # Loop through entries and compose polygons
    excluded_polygons={}
    past_ids = {}
    last_id = None
    for i, id in enumerate(X['id']):

        # Check for duplicate polygons
        if id in past_ids:
            msg = 'Polygon %s was duplicated in line %d' % (id, i)
            raise Exception, msg
        
        if id not in polygons:
            # Start new polygon
            polygons[id] = []
            if values is not None:
                values[id] = X[value_name][i]

            # Keep track of previous polygon ids
            if last_id is not None:
                past_ids[last_id] = i
            
        # Append this point to current polygon
        point = [float(X['easting'][i]), float(X['northing'][i])]

        if clipping_polygons is not None:
            exclude=True
            for clipping_polygon in clipping_polygons:
                if inside_polygon(point, clipping_polygon):
                    exclude=False
                    break
                
            if exclude is True:
                excluded_polygons[id]=True

        polygons[id].append(point)    
            
        # Check that value is the same across each polygon
        msg = 'Values must be the same across each polygon.'
        msg += 'I got %s in line %d but it should have been %s' % (X[value_name][i], i, values[id])
        assert values[id] == X[value_name][i], msg

        last_id = id

    # Weed out polygons that were not wholly inside clipping polygons
    for id in excluded_polygons:
        del polygons[id]
        
    return polygons, values


            




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


##
# @brief 
# @param basename_in 
# @param extra_name_out 
# @param quantities 
# @param timestep 
# @param reduction 
# @param cellsize 
# @param number_of_decimal_places 
# @param NODATA_value 
# @param easting_min 
# @param easting_max 
# @param northing_min 
# @param northing_max 
# @param verbose 
# @param origin 
# @param datum 
# @param format 
# @return 
def export_grid(basename_in, extra_name_out=None,
                quantities=None, # defaults to elevation
                reduction=None,
                cellsize=10,
                number_of_decimal_places=None,
                NODATA_value=-9999,
                easting_min=None,
                easting_max=None,
                northing_min=None,
                northing_max=None,
                verbose=False,
                origin=None,
                datum='WGS84',
                format='ers'):
    """Wrapper for sww2dem.
    See sww2dem to find out what most of the parameters do.

    Quantities is a list of quantities.  Each quantity will be
    calculated for each sww file.

    This returns the basenames of the files returned, which is made up
    of the dir and all of the file name, except the extension.

    This function returns the names of the files produced.

    It will also produce as many output files as there are input sww files.
    """

    if quantities is None:
        quantities = ['elevation']

    if type(quantities) is str:
            quantities = [quantities]

    # How many sww files are there?
    dir, base = os.path.split(basename_in)

    iterate_over = get_all_swwfiles(dir, base, verbose)

    if dir == "":
        dir = "." # Unix compatibility

    files_out = []
    for sww_file in iterate_over:
        for quantity in quantities:
            if extra_name_out is None:
                basename_out = sww_file + '_' + quantity
            else:
                basename_out = sww_file + '_' + quantity + '_' + extra_name_out

            file_out = sww2dem(dir+sep+sww_file, dir+sep+basename_out,
                               quantity,
                               reduction,
                               cellsize,
                               number_of_decimal_places,
                               NODATA_value,
                               easting_min,
                               easting_max,
                               northing_min,
                               northing_max,
                               verbose,
                               origin,
                               datum,
                               format)
            files_out.append(file_out)
    return files_out


##
# @brief Read DEM file, decimate data, write new DEM file.
# @param basename_in Stem of the input filename.
# @param stencil 
# @param cellsize_new New cell size to resample on.
# @param basename_out Stem of the output filename.
# @param verbose True if this function is to be verbose.
def decimate_dem(basename_in, stencil, cellsize_new, basename_out=None,
                 verbose=False):
    """Read Digitial Elevation model from the following NetCDF format (.dem)

    Example:

    ncols         3121
    nrows         1800
    xllcorner     722000
    yllcorner     5893000
    cellsize      25
    NODATA_value  -9999
    138.3698 137.4194 136.5062 135.5558 ..........

    Decimate data to cellsize_new using stencil and write to NetCDF dem format.
    """

    import os
    from Scientific.IO.NetCDF import NetCDFFile

    root = basename_in
    inname = root + '.dem'

    #Open existing netcdf file to read
    infile = NetCDFFile(inname, netcdf_mode_r)

    if verbose: log.critical('Reading DEM from %s' % inname)

    # Read metadata (convert from numpy.int32 to int where appropriate)
    ncols = int(infile.ncols[0])
    nrows = int(infile.nrows[0])
    xllcorner = infile.xllcorner[0]
    yllcorner = infile.yllcorner[0]
    cellsize = int(infile.cellsize[0])
    NODATA_value = int(infile.NODATA_value[0])
    zone = int(infile.zone[0])
    false_easting = infile.false_easting[0]
    false_northing = infile.false_northing[0]
    projection = infile.projection
    datum = infile.datum
    units = infile.units

    dem_elevation = infile.variables['elevation']

    #Get output file name
    if basename_out == None:
        outname = root + '_' + repr(cellsize_new) + '.dem'
    else:
        outname = basename_out + '.dem'

    if verbose: log.critical('Write decimated NetCDF file to %s' % outname)

    #Determine some dimensions for decimated grid
    (nrows_stencil, ncols_stencil) = stencil.shape
    x_offset = ncols_stencil / 2
    y_offset = nrows_stencil / 2
    cellsize_ratio = int(cellsize_new / cellsize)
    ncols_new = 1 + (ncols - ncols_stencil) / cellsize_ratio
    nrows_new = 1 + (nrows - nrows_stencil) / cellsize_ratio

    #print type(ncols_new), ncols_new
    
    #Open netcdf file for output
    outfile = NetCDFFile(outname, netcdf_mode_w)

    #Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = 'NetCDF DEM format for compact and portable ' \
                          'storage of spatial point data'

    #Georeferencing
    outfile.zone = zone
    outfile.projection = projection
    outfile.datum = datum
    outfile.units = units

    outfile.cellsize = cellsize_new
    outfile.NODATA_value = NODATA_value
    outfile.false_easting = false_easting
    outfile.false_northing = false_northing

    outfile.xllcorner = xllcorner + (x_offset * cellsize)
    outfile.yllcorner = yllcorner + (y_offset * cellsize)
    outfile.ncols = ncols_new
    outfile.nrows = nrows_new

    # dimension definition
    #print nrows_new, ncols_new, nrows_new*ncols_new
    #print type(nrows_new), type(ncols_new), type(nrows_new*ncols_new)
    outfile.createDimension('number_of_points', nrows_new*ncols_new)

    # variable definition
    outfile.createVariable('elevation', netcdf_float, ('number_of_points',))

    # Get handle to the variable
    elevation = outfile.variables['elevation']

    dem_elevation_r = num.reshape(dem_elevation, (nrows, ncols))

    #Store data
    global_index = 0
    for i in range(nrows_new):
        if verbose: log.critical('Processing row %d of %d' % (i, nrows_new))

        lower_index = global_index
        telev = num.zeros(ncols_new, num.float)
        local_index = 0
        trow = i * cellsize_ratio

        for j in range(ncols_new):
            tcol = j * cellsize_ratio
            tmp = dem_elevation_r[trow:trow+nrows_stencil,
                                  tcol:tcol+ncols_stencil]

            #if dem contains 1 or more NODATA_values set value in
            #decimated dem to NODATA_value, else compute decimated
            #value using stencil
            if num.sum(num.sum(num.equal(tmp, NODATA_value))) > 0:
                telev[local_index] = NODATA_value
            else:
                telev[local_index] = num.sum(num.sum(tmp * stencil))

            global_index += 1
            local_index += 1

        upper_index = global_index

        elevation[lower_index:upper_index] = telev

    assert global_index == nrows_new*ncols_new, \
           'index not equal to number of points'

    infile.close()
    outfile.close()


    ####  URS 2 SWW  ###

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
# @brief 
# @param long_lat_dep 
# @return A tuple (long, lat, quantity).
# @note The latitude is the fastest varying dimension - in mux files.
def lon_lat2grid(long_lat_dep):
    """
    given a list of points that are assumed to be an a grid,
    return the long's and lat's of the grid.
    long_lat_dep is an array where each row is a position.
    The first column is longitudes.
    The second column is latitudes.

    The latitude is the fastest varying dimension - in mux files
    """

    LONG = 0
    LAT = 1
    QUANTITY = 2

    long_lat_dep = ensure_numeric(long_lat_dep, num.float)

    num_points = long_lat_dep.shape[0]
    this_rows_long = long_lat_dep[0,LONG]

    # Count the length of unique latitudes
    i = 0
    while long_lat_dep[i,LONG] == this_rows_long and i < num_points:
        i += 1

    # determine the lats and longsfrom the grid
    lat = long_lat_dep[:i, LAT]
    long = long_lat_dep[::i, LONG]

    lenlong = len(long)
    lenlat = len(lat)

    msg = 'Input data is not gridded'
    assert num_points % lenlat == 0, msg
    assert num_points % lenlong == 0, msg

    # Test that data is gridded
    for i in range(lenlong):
        msg = 'Data is not gridded.  It must be for this operation'
        first = i * lenlat
        last = first + lenlat

        assert num.allclose(long_lat_dep[first:last,LAT], lat), msg
        assert num.allclose(long_lat_dep[first:last,LONG], long[i]), msg

    msg = 'Out of range latitudes/longitudes'
    for l in lat:assert -90 < l < 90 , msg
    for l in long:assert -180 < l < 180 , msg

    # Changing quantity from lat being the fastest varying dimension to
    # long being the fastest varying dimension
    # FIXME - make this faster/do this a better way
    # use numeric transpose, after reshaping the quantity vector
    quantity = num.zeros(num_points, num.float)

    for lat_i, _ in enumerate(lat):
        for long_i, _ in enumerate(long):
            q_index = lat_i*lenlong + long_i
            lld_index = long_i*lenlat + lat_i
            temp = long_lat_dep[lld_index, QUANTITY]
            quantity[q_index] = temp

    return long, lat, quantity

################################################################################
# END URS 2 SWW
################################################################################




##
# @brief Store keyword params into a CSV file.
# @param verbose True if this function is to be verbose.
# @param kwargs Dictionary of keyword args to store.
# @note If kwargs dict contains 'file_name' key, that has the output filename.
#       If not, make up a filename in the output directory.
def store_parameters(verbose=False, **kwargs):
    """
    Store "kwargs" into a temp csv file, if "completed" is in kwargs,
    csv file is kwargs[file_name] else it is kwargs[output_dir]+details_temp.csv

    Must have a file_name keyword arg, this is what is writing to.
    might be a better way to do this using CSV module Writer and writeDict.

    writes file to "output_dir" unless "completed" is in kwargs, then
    it writes to "file_name" kwargs
    """

    import types

    # Check that kwargs is a dictionary
    if type(kwargs) != types.DictType:
        raise TypeError

    # is 'completed' in kwargs?
    completed = kwargs.has_key('completed')

    # get file name and removes from dict and assert that a file_name exists
    if completed:
        try:
            file = str(kwargs['file_name'])
        except:
            raise 'kwargs must have file_name'
    else:
        # write temp file in output directory
        try:
            file = str(kwargs['output_dir']) + 'detail_temp.csv'
        except:
            raise 'kwargs must have output_dir'

    # extracts the header info and the new line info
    line = ''
    header = ''
    count = 0
    keys = kwargs.keys()
    keys.sort()

    # used the sorted keys to create the header and line data
    for k in keys:
        header += str(k)
        line += str(kwargs[k])
        count += 1
        if count < len(kwargs):
            header += ','
            line += ','
    header += '\n'
    line += '\n'

    # checks the header info, if the same, then write, if not create a new file
    # try to open!
    try:
        fid = open(file, 'r')
        file_header = fid.readline()
        fid.close()
        if verbose: log.critical('read file header %s' % file_header)
    except:
        msg = 'try to create new file: %s' % file
        if verbose: log.critical(msg)
        #tries to open file, maybe directory is bad
        try:
            fid = open(file, 'w')
            fid.write(header)
            fid.close()
            file_header=header
        except:
            msg = 'cannot create new file: %s' % file
            raise Exception, msg

    # if header is same or this is a new file
    if file_header == str(header):
        fid = open(file, 'a')
        fid.write(line)
        fid.close()
    else:
        # backup plan,
        # if header is different and has completed will append info to
        # end of details_temp.cvs file in output directory
        file = str(kwargs['output_dir']) + 'detail_temp.csv'
        fid = open(file, 'a')
        fid.write(header)
        fid.write(line)
        fid.close()

        if verbose:
            log.critical('file %s', file_header.strip('\n'))
            log.critical('head %s', header.strip('\n'))
        if file_header.strip('\n') == str(header):
            log.critical('they equal')

        msg = 'WARNING: File header does not match input info, ' \
              'the input variables have changed, suggest you change file name'
        log.critical(msg)



##
# @brief Convert points to a polygon (?)
# @param points_file The points file.
# @param minimum_triangle_angle ??
# @return 
def points2polygon(points_file, minimum_triangle_angle=3.0):
    """
    WARNING: This function is not fully working.

    Function to return a polygon returned from alpha shape, given a points file.

    WARNING: Alpha shape returns multiple polygons, but this function only
             returns one polygon.
    """

    from anuga.pmesh.mesh import Mesh, importMeshFromFile
    from anuga.shallow_water import Domain
    from anuga.pmesh.mesh_interface import create_mesh_from_regions

    mesh = importMeshFromFile(points_file)
    mesh.auto_segment()
    mesh.exportASCIIsegmentoutlinefile("outline.tsh")
    mesh2 = importMeshFromFile("outline.tsh")
    mesh2.generate_mesh(maximum_triangle_area=1000000000,
                        minimum_triangle_angle=minimum_triangle_angle,
                        verbose=False)
    mesh2.export_mesh_file('outline_meshed.tsh')
    domain = Domain("outline_meshed.tsh", use_cache = False)
    polygon =  domain.get_boundary_polygon()
    return polygon


################################################################################

if __name__ == "__main__":
    # setting umask from config to force permissions for all files and
    # directories created to the same. (it was noticed the "mpirun" doesn't
    # honour the umask set in your .bashrc etc file)

    from config import umask
    import os
    os.umask(umask)
