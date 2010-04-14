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
import exceptions
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
from anuga.utilities.anuga_exceptions import ANUGAError
from anuga.shallow_water import Domain
from anuga.abstract_2d_finite_volumes.pmesh2domain import \
     pmesh_to_domain_instance
from anuga.abstract_2d_finite_volumes.util import get_revision_number, \
     remove_lone_verts, sww2timeseries, get_centroid_values

from anuga.abstract_2d_finite_volumes.neighbour_mesh import segment_midpoints
from anuga.load_mesh.loadASCII import export_mesh_file
from anuga.utilities.polygon import intersection

from anuga.utilities.system_tools import get_vars_in_expression
import anuga.utilities.log as log


# Default block size for sww2dem()
DEFAULT_BLOCK_SIZE = 10000

######
# Exception classes
######

class TitleValueError(exceptions.Exception): pass
class DataMissingValuesError(exceptions.Exception): pass
class DataFileNotOpenError(exceptions.Exception): pass
class DataTimeError(exceptions.Exception): pass
class DataDomainError(exceptions.Exception): pass
class NewQuantity(exceptions.Exception): pass


######
# formula mappings
######

quantity_formula = {'momentum':'(xmomentum**2 + ymomentum**2)**0.5',
                    'depth':'stage-elevation',
                    'speed': \
 '(xmomentum**2 + ymomentum**2)**0.5/(stage-elevation+1.e-6/(stage-elevation))'}


##
# @brief Convert a possible filename into a standard form.
# @param s Filename to process.
# @return The new filename string.
def make_filename(s):
    """Transform argument string into a Sexsuitable filename
    """

    s = s.strip()
    s = s.replace(' ', '_')
    s = s.replace('(', '')
    s = s.replace(')', '')
    s = s.replace('__', '_')

    return s


##
# @brief Check that a specified filesystem directory path exists.
# @param path The dirstory path to check.
# @param verbose True if this function is to be verbose.
# @note If directory path doesn't exist, it will be created.
def check_dir(path, verbose=None):
    """Check that specified path exists.
    If path does not exist it will be created if possible

    USAGE:
       checkdir(path, verbose):

    ARGUMENTS:
        path -- Directory
        verbose -- Flag verbose output (default: None)

    RETURN VALUE:
        Verified path including trailing separator
    """

    import os.path

    if sys.platform in ['nt', 'dos', 'win32', 'what else?']:
        unix = 0
    else:
        unix = 1

    # add terminal separator, if it's not already there
    if path[-1] != os.sep:
        path = path + os.sep

    # expand ~ or ~username in path
    path = os.path.expanduser(path)

    # create directory if required
    if not (os.access(path, os.R_OK and os.W_OK) or path == ''):
        try:
            exitcode = os.mkdir(path)

            # Change access rights if possible
            if unix:
                exitcode = os.system('chmod 775 ' + path)
            else:
                pass  # FIXME: What about access rights under Windows?

            if verbose: log.critical('MESSAGE: Directory %s created.' % path)
        except:
            log.critical('WARNING: Directory %s could not be created.' % path)
            if unix:
                path = '/tmp/'
            else:
                path = 'C:' + os.sep

            log.critical("Using directory '%s' instead" % path)

    return path


##
# @brief Delete directory and all sub-directories.
# @param path Path to the directory to delete.
def del_dir(path):
    """Recursively delete directory path and all its contents
    """

    if os.path.isdir(path):
        for file in os.listdir(path):
            X = os.path.join(path, file)

            if os.path.isdir(X) and not os.path.islink(X):
                del_dir(X)
            else:
                try:
                    os.remove(X)
                except:
                    log.critical("Could not remove file %s" % X)

        os.rmdir(path)


##
# @brief ??
# @param path 
# @param __func__ 
# @param verbose True if this function is to be verbose.
# @note ANOTHER OPTION, IF NEED IN THE FUTURE, Nick B 7/2007
def rmgeneric(path, func, verbose=False):
    ERROR_STR= """Error removing %(path)s, %(error)s """

    try:
        func(path)
        if verbose: log.critical('Removed %s' % path)
    except OSError, (errno, strerror):
        log.critical(ERROR_STR % {'path' : path, 'error': strerror })


##
# @brief Remove directory and all sub-directories.
# @param path Filesystem path to directory to remove.
# @param verbose True if this function is to be verbose.
def removeall(path, verbose=False):
    if not os.path.isdir(path):
        return

    for x in os.listdir(path):
        fullpath = os.path.join(path, x)
        if os.path.isfile(fullpath):
            f = os.remove
            rmgeneric(fullpath, f)
        elif os.path.isdir(fullpath):
            removeall(fullpath)
            f = os.rmdir
            rmgeneric(fullpath, f, verbose)


##
# @brief Create a standard filename.
# @param datadir Directory where file is to be created.
# @param filename Filename 'stem'.
# @param format Format of the file, becomes filename extension.
# @param size Size of file, becomes part of filename.
# @param time Time (float), becomes part of filename.
# @return The complete filename path, including directory.
# @note The containing directory is created, if necessary.
def create_filename(datadir, filename, format, size=None, time=None):
    FN = check_dir(datadir) + filename

    if size is not None:
        FN += '_size%d' % size

    if time is not None:
        FN += '_time%.2f' % time

    FN += '.' + format

    return FN


##
# @brief Get all files with a standard name and a given set of attributes.
# @param datadir Directory files must be in.
# @param filename Filename stem.
# @param format Filename extension.
# @param size Filename size.
# @return A list of fielnames (including directory) that match the attributes.
def get_files(datadir, filename, format, size):
    """Get all file (names) with given name, size and format
    """

    import glob

    dir = check_dir(datadir)
    pattern = dir + os.sep + filename + '_size=%d*.%s' % (size, format)

    return glob.glob(pattern)


##
# @brief Generic class for storing output to e.g. visualisation or checkpointing
class Data_format:
    """Generic interface to data formats
    """

    ##
    # @brief Instantiate this instance.
    # @param domain 
    # @param extension 
    # @param mode The mode of the underlying file.
    def __init__(self, domain, extension, mode=netcdf_mode_w):
        assert mode[0] in ['r', 'w', 'a'], \
               "Mode %s must be either:\n" % mode + \
               "   'w' (write)\n" + \
               "   'r' (read)\n" + \
               "   'a' (append)"

        #Create filename
        self.filename = create_filename(domain.get_datadir(),
                                        domain.get_name(), extension)

        self.timestep = 0
        self.domain = domain

        # Exclude ghosts in case this is a parallel domain
        self.number_of_nodes = domain.number_of_full_nodes
        self.number_of_volumes = domain.number_of_full_triangles
        #self.number_of_volumes = len(domain)

        #FIXME: Should we have a general set_precision function?


##
# @brief Class for storing output to e.g. visualisation
class SWW_file(Data_format):
    """Interface to native NetCDF format (.sww) for storing model output

    There are two kinds of data

    1: Constant data: Vertex coordinates and field values. Stored once
    2: Variable data: Conserved quantities. Stored once per timestep.

    All data is assumed to reside at vertex locations.
    """

    ##
    # @brief Instantiate this instance.
    # @param domain ??
    # @param mode Mode of the underlying data file.
    # @param max_size ??
    # @param recursion ??
    # @note Prepare the underlying data file if mode starts with 'w'.
    def __init__(self, domain, 
                 mode=netcdf_mode_w, max_size=2000000000, recursion=False):
        from Scientific.IO.NetCDF import NetCDFFile

        self.precision = netcdf_float32 # Use single precision for quantities
        self.recursion = recursion
        self.mode = mode
        if hasattr(domain, 'max_size'):
            self.max_size = domain.max_size # File size max is 2Gig
        else:
            self.max_size = max_size
        if hasattr(domain, 'minimum_storable_height'):
            self.minimum_storable_height = domain.minimum_storable_height
        else:
            self.minimum_storable_height = default_minimum_storable_height

        # Call parent constructor
        Data_format.__init__(self, domain, 'sww', mode)

        # Get static and dynamic quantities from domain
        static_quantities = []
        dynamic_quantities = []
        
        for q in domain.quantities_to_be_stored:
            flag = domain.quantities_to_be_stored[q]
        
            msg = 'Quantity %s is requested to be stored ' % q
            msg += 'but it does not exist in domain.quantities'
            assert q in domain.quantities, msg
        
            assert flag in [1,2]
            if flag == 1: static_quantities.append(q)
            if flag == 2: dynamic_quantities.append(q)                
                       
        
        # NetCDF file definition
        fid = NetCDFFile(self.filename, mode)
        if mode[0] == 'w':
            description = 'Output from anuga.abstract_2d_finite_volumes ' \
                          'suitable for plotting'
                          
            self.writer = Write_sww(static_quantities, dynamic_quantities)
            self.writer.store_header(fid,
                                     domain.starttime,
                                     self.number_of_volumes,
                                     self.domain.number_of_full_nodes,
                                     description=description,
                                     smoothing=domain.smooth,
                                     order=domain.default_order,
                                     sww_precision=self.precision)

            # Extra optional information
            if hasattr(domain, 'texture'):
                fid.texture = domain.texture

            if domain.quantities_to_be_monitored is not None:
                fid.createDimension('singleton', 1)
                fid.createDimension('two', 2)

                poly = domain.monitor_polygon
                if poly is not None:
                    N = len(poly)
                    fid.createDimension('polygon_length', N)
                    fid.createVariable('extrema.polygon',
                                       self.precision,
                                       ('polygon_length', 'two'))
                    fid.variables['extrema.polygon'][:] = poly

                interval = domain.monitor_time_interval
                if interval is not None:
                    fid.createVariable('extrema.time_interval',
                                       self.precision,
                                       ('two',))
                    fid.variables['extrema.time_interval'][:] = interval

                for q in domain.quantities_to_be_monitored:
                    fid.createVariable(q + '.extrema', self.precision,
                                       ('numbers_in_range',))
                    fid.createVariable(q + '.min_location', self.precision,
                                       ('numbers_in_range',))
                    fid.createVariable(q + '.max_location', self.precision,
                                       ('numbers_in_range',))
                    fid.createVariable(q + '.min_time', self.precision,
                                       ('singleton',))
                    fid.createVariable(q + '.max_time', self.precision,
                                       ('singleton',))

        fid.close()

    ##
    # @brief Store connectivity data into the underlying data file.
    def store_connectivity(self):
        """Store information about nodes, triangles and static quantities

        Writes x,y coordinates of triangles and their connectivity.
        
        Store also any quantity that has been identified as static.
        """

        # FIXME: Change name to reflect the fact thta this function 
        # stores both connectivity (triangulation) and static quantities
        
        from Scientific.IO.NetCDF import NetCDFFile

        domain = self.domain

        # append to the NetCDF file
        fid = NetCDFFile(self.filename, netcdf_mode_a)

        # Get X, Y from one (any) of the quantities
        Q = domain.quantities.values()[0]
        X,Y,_,V = Q.get_vertex_values(xy=True, precision=self.precision)

        # store the connectivity data
        points = num.concatenate((X[:,num.newaxis],Y[:,num.newaxis]), axis=1)
        self.writer.store_triangulation(fid,
                                        points,
                                        V.astype(num.float32),
                                        points_georeference=\
                                            domain.geo_reference)


        # Get names of static quantities
        static_quantities = {}
        for name in self.writer.static_quantities:
            Q = domain.quantities[name]
            A, _ = Q.get_vertex_values(xy=False, 
                                       precision=self.precision)
            static_quantities[name] = A
        
        # Store static quantities        
        self.writer.store_static_quantities(fid, **static_quantities)
                                            
        fid.close()

    ##
    # @brief Store time and time dependent quantities 
    # to the underlying data file.
    def store_timestep(self):
        """Store time and time dependent quantities
        """

        from Scientific.IO.NetCDF import NetCDFFile
        import types
        from time import sleep
        from os import stat

        # Get NetCDF
        retries = 0
        file_open = False
        while not file_open and retries < 10:
            try:
                # Open existing file
                fid = NetCDFFile(self.filename, netcdf_mode_a)
            except IOError:
                # This could happen if someone was reading the file.
                # In that case, wait a while and try again
                msg = 'Warning (store_timestep): File %s could not be opened' \
                      % self.filename
                msg += ' - trying step %s again' % self.domain.time
                log.critical(msg)
                retries += 1
                sleep(1)
            else:
                file_open = True

        if not file_open:
            msg = 'File %s could not be opened for append' % self.filename
            raise DataFileNotOpenError, msg

        # Check to see if the file is already too big:
        time = fid.variables['time']
        i = len(time) + 1
        file_size = stat(self.filename)[6]
        file_size_increase = file_size / i
        if file_size + file_size_increase > self.max_size * 2**self.recursion:
            # In order to get the file name and start time correct,
            # I change the domain.filename and domain.starttime.
            # This is the only way to do this without changing
            # other modules (I think).

            # Write a filename addon that won't break the anuga viewers
            # (10.sww is bad)
            filename_ext = '_time_%s' % self.domain.time
            filename_ext = filename_ext.replace('.', '_')

            # Remember the old filename, then give domain a
            # name with the extension
            old_domain_filename = self.domain.get_name()
            if not self.recursion:
                self.domain.set_name(old_domain_filename + filename_ext)

            # Temporarily change the domain starttime to the current time
            old_domain_starttime = self.domain.starttime
            self.domain.starttime = self.domain.get_time()

            # Build a new data_structure.
            next_data_structure = SWW_file(self.domain, mode=self.mode,
                                           max_size=self.max_size,
                                           recursion=self.recursion+1)
            if not self.recursion:
                log.critical('    file_size = %s' % file_size)
                log.critical('    saving file to %s'
                             % next_data_structure.filename) 

            # Set up the new data_structure
            self.domain.writer = next_data_structure

            # Store connectivity and first timestep
            next_data_structure.store_connectivity()
            next_data_structure.store_timestep()
            fid.sync()
            fid.close()

            # Restore the old starttime and filename
            self.domain.starttime = old_domain_starttime
            self.domain.set_name(old_domain_filename)
        else:
            self.recursion = False
            domain = self.domain

            # Get the variables
            time = fid.variables['time']
            i = len(time)
             
            if 'stage' in self.writer.dynamic_quantities:            
                # Select only those values for stage, 
                # xmomentum and ymomentum (if stored) where 
                # depth exceeds minimum_storable_height
                #
                # In this branch it is assumed that elevation
                # is also available as a quantity            
            
                Q = domain.quantities['stage']
                w, _ = Q.get_vertex_values(xy=False)
                
                Q = domain.quantities['elevation']
                z, _ = Q.get_vertex_values(xy=False)                
                
                storable_indices = (w-z >= self.minimum_storable_height)
            else:
                # Very unlikely branch
                storable_indices = None # This means take all
            
            
            # Now store dynamic quantities
            dynamic_quantities = {}
            for name in self.writer.dynamic_quantities:
                netcdf_array = fid.variables[name]
                
                Q = domain.quantities[name]
                A, _ = Q.get_vertex_values(xy=False,
                                           precision=self.precision)

                if storable_indices is not None:
                    if name == 'stage':
                        A = num.choose(storable_indices, (z, A))

                    if name in ['xmomentum', 'ymomentum']:
                        # Get xmomentum where depth exceeds 
                        # minimum_storable_height
                        
                        # Define a zero vector of same size and type as A
                        # for use with momenta
                        null = num.zeros(num.size(A), A.dtype.char)
                        A = num.choose(storable_indices, (null, A))
                
                dynamic_quantities[name] = A
                
                                        
            # Store dynamic quantities
            self.writer.store_quantities(fid,
                                         time=self.domain.time,
                                         sww_precision=self.precision,
                                         **dynamic_quantities)


            # Update extrema if requested
            domain = self.domain
            if domain.quantities_to_be_monitored is not None:
                for q, info in domain.quantities_to_be_monitored.items():
                    if info['min'] is not None:
                        fid.variables[q + '.extrema'][0] = info['min']
                        fid.variables[q + '.min_location'][:] = \
                                        info['min_location']
                        fid.variables[q + '.min_time'][0] = info['min_time']

                    if info['max'] is not None:
                        fid.variables[q + '.extrema'][1] = info['max']
                        fid.variables[q + '.max_location'][:] = \
                                        info['max_location']
                        fid.variables[q + '.max_time'][0] = info['max_time']

            # Flush and close
            fid.sync()
            fid.close()


##
# @brief Class to open an sww file so that domain can be populated with quantity values 
class Read_sww:

    def __init__(self, source):
        """The source parameter is assumed to be a NetCDF sww file.
        """

        self.source = source

        self.frame_number = 0

        fin = NetCDFFile(self.source, 'r')

        self.time = num.array(fin.variables['time'], num.float)
        self.last_frame_number = self.time.shape[0] - 1

        self.frames = num.arange(self.last_frame_number+1)

        fin.close()
        
        self.read_mesh()

        self.quantities = {}

        self.read_quantities()


    def read_mesh(self):
        fin = NetCDFFile(self.source, 'r')

        self.vertices = num.array(fin.variables['volumes'], num.int)
        
        self.x = x = num.array(fin.variables['x'], num.float)
        self.y = y = num.array(fin.variables['y'], num.float)

        assert len(self.x) == len(self.y)
        
        self.xmin = num.min(x)
        self.xmax = num.max(x)
        self.ymin = num.min(y)
        self.ymax = num.max(y)



        fin.close()
        
    def read_quantities(self, frame_number=0):

        assert frame_number >= 0 and frame_number <= self.last_frame_number

        self.frame_number = frame_number

        M = len(self.x)/3
        
        fin = NetCDFFile(self.source, 'r')
        
        for q in filter(lambda n:n != 'x' and n != 'y' and n != 'time' and n != 'volumes' and \
                        '_range' not in n, \
                        fin.variables.keys()):
            if len(fin.variables[q].shape) == 1: # Not a time-varying quantity
                self.quantities[q] = num.ravel(num.array(fin.variables[q], num.float)).reshape(M,3)
            else: # Time-varying, get the current timestep data
                self.quantities[q] = num.array(fin.variables[q][self.frame_number], num.float).reshape(M,3)
        fin.close()
        return self.quantities

    def get_bounds(self):
        return [self.xmin, self.xmax, self.ymin, self.ymax]

    def get_last_frame_number(self):
        return self.last_frame_number

    def get_time(self):
        return self.time[self.frame_number]


##
# @brief Class for handling checkpoints data
# @note This is not operational at the moment
class CPT_file(Data_format):
    """Interface to native NetCDF format (.cpt) to be 
    used for checkpointing (one day)
    """

    ##
    # @brief Initialize this instantiation.
    # @param domain ??
    # @param mode Mode of underlying data file (default WRITE).
    def __init__(self, domain, mode=netcdf_mode_w):
        from Scientific.IO.NetCDF import NetCDFFile

        self.precision = netcdf_float #Use full precision

        Data_format.__init__(self, domain, 'sww', mode)

        # NetCDF file definition
        fid = NetCDFFile(self.filename, mode)
        if mode[0] == 'w':
            # Create new file
            fid.institution = 'Geoscience Australia'
            fid.description = 'Checkpoint data'
            #fid.smooth = domain.smooth
            fid.order = domain.default_order

            # Dimension definitions
            fid.createDimension('number_of_volumes', self.number_of_volumes)
            fid.createDimension('number_of_vertices', 3)

            # Store info at all vertices (no smoothing)
            fid.createDimension('number_of_points', 3*self.number_of_volumes)
            fid.createDimension('number_of_timesteps', None) #extensible

            # Variable definitions

            # Mesh
            fid.createVariable('x', self.precision, ('number_of_points',))
            fid.createVariable('y', self.precision, ('number_of_points',))


            fid.createVariable('volumes', netcdf_int, ('number_of_volumes',
                                                       'number_of_vertices'))

            fid.createVariable('time', self.precision, ('number_of_timesteps',))

            #Allocate space for all quantities
            for name in domain.quantities.keys():
                fid.createVariable(name, self.precision,
                                   ('number_of_timesteps',
                                    'number_of_points'))

        fid.close()

    ##
    # @brief Store connectivity data to underlying data file.
    def store_checkpoint(self):
        """Write x,y coordinates of triangles.
        Write connectivity (
        constituting
        the bed elevation.
        """

        from Scientific.IO.NetCDF import NetCDFFile

        domain = self.domain

        #Get NetCDF
        fid = NetCDFFile(self.filename, netcdf_mode_a)

        # Get the variables
        x = fid.variables['x']
        y = fid.variables['y']

        volumes = fid.variables['volumes']

        # Get X, Y and bed elevation Z
        Q = domain.quantities['elevation']
        X,Y,Z,V = Q.get_vertex_values(xy=True, precision=self.precision)

        x[:] = X.astype(self.precision)
        y[:] = Y.astype(self.precision)
        z[:] = Z.astype(self.precision)

        volumes[:] = V

        fid.close()

    ##
    # @brief Store time and named quantities to underlying data file.
    # @param name 
    def store_timestep(self, name):
        """Store time and named quantity to file
        """

        from Scientific.IO.NetCDF import NetCDFFile
        from time import sleep

        #Get NetCDF
        retries = 0
        file_open = False
        while not file_open and retries < 10:
            try:
                fid = NetCDFFile(self.filename, netcdf_mode_a)
            except IOError:
                #This could happen if someone was reading the file.
                #In that case, wait a while and try again
                msg = 'Warning (store_timestep): File %s could not be opened' \
                      ' - trying again' % self.filename
                log.critical(msg)
                retries += 1
                sleep(1)
            else:
                file_open = True

        if not file_open:
            msg = 'File %s could not be opened for append' % self.filename
            raise DataFileNotOpenError, msg

        domain = self.domain

        # Get the variables
        time = fid.variables['time']
        stage = fid.variables['stage']
        i = len(time)

        #Store stage
        time[i] = self.domain.time

        # Get quantity
        Q = domain.quantities[name]
        A,V = Q.get_vertex_values(xy=False, precision=self.precision)

        stage[i,:] = A.astype(self.precision)

        #Flush and close
        fid.sync()
        fid.close()


##
# @brief Class for National Exposure Database storage (NEXIS).

LAT_TITLE = 'LATITUDE'
LONG_TITLE = 'LONGITUDE'
X_TITLE = 'x'
Y_TITLE = 'y'

class Exposure_csv:

    ##
    # @brief Instantiate this instance.
    # @param file_name Name of underlying data file.
    # @param latitude_title ??
    # @param longitude_title ??
    # @param is_x_y_locations ??
    # @param x_title ??
    # @param y_title ??
    # @param refine_polygon ??
    # @param title_check_list ??
    def __init__(self,file_name, latitude_title=LAT_TITLE,
                 longitude_title=LONG_TITLE, is_x_y_locations=None,
                 x_title=X_TITLE, y_title=Y_TITLE,
                 refine_polygon=None, title_check_list=None):
        """
        This class is for handling the exposure csv file.
        It reads the file in and converts the lats and longs to a geospatial
        data object.
        Use the methods to read and write columns.

        The format of the csv files it reads is;
           The first row is a title row.
           comma's are the delimiters
           each column is a 'set' of data

        Feel free to use/expand it to read other csv files.

        It is not for adding and deleting rows

        Can geospatial handle string attributes? It's not made for them.
        Currently it can't load and save string att's.

        So just use geospatial to hold the x, y and georef? Bad, since
        different att's are in diferent structures.  Not so bad, the info
        to write if the .csv file is saved is in attribute_dic

        The location info is in the geospatial attribute.
        """

        self._file_name = file_name
        self._geospatial = None #

        # self._attribute_dic is a dictionary.
        #The keys are the column titles.
        #The values are lists of column data

        # self._title_index_dic is a dictionary.
        #The keys are the column titles.
        #The values are the index positions of file columns.
        self._attribute_dic, self._title_index_dic = \
            csv2dict(self._file_name, title_check_list=title_check_list)
        try:
            #Have code here that handles caps or lower
            lats = self._attribute_dic[latitude_title]
            longs = self._attribute_dic[longitude_title]
        except KeyError:
            # maybe a warning..
            #Let's see if this works..
            if False != is_x_y_locations:
                is_x_y_locations = True
            pass
        else:
            self._geospatial = Geospatial_data(latitudes=lats,
                                               longitudes=longs)

        if is_x_y_locations is True:
            if self._geospatial is not None:
                pass #fixme throw an error
            try:
                xs = self._attribute_dic[x_title]
                ys = self._attribute_dic[y_title]
                points = [[float(i),float(j)] for i,j in map(None,xs,ys)]
            except KeyError:
                # maybe a warning..
                msg = "Could not find location information."
                raise TitleValueError, msg
            else:
                self._geospatial = Geospatial_data(data_points=points)

        # create a list of points that are in the refining_polygon
        # described by a list of indexes representing the points

    ##
    # @brief Create a comparison method.
    # @param self This object.
    # @param other The other object.
    # @return True if objects are 'same'.
    def __cmp__(self, other):
        #check that 'other' is an instance of this class
        if isinstance(self, type(other)):
            result = cmp(self._attribute_dic, other._attribute_dic)
            if result <> 0:
                return result

            # The order of the columns is important. Therefore..
            result = cmp(self._title_index_dic, other._title_index_dic)
            if result <> 0:
                return result
            for self_ls, other_ls in map(None, self._attribute_dic,
                                         other._attribute_dic):
                result = cmp(self._attribute_dic[self_ls],
                             other._attribute_dic[other_ls])
                if result <> 0:
                    return result
            return 0
        else:
            return 1

    ##
    # @brief Get a list of column values given a column name.
    # @param column_name The name of the column to get values from.
    # @param use_refind_polygon Unused??
    def get_column(self, column_name, use_refind_polygon=False):
        """
        Given a column name return a list of the column values

        Note, the type of the values will be String!
        do this to change a list of strings to a list of floats
        time = [float(x) for x in time]

        Not implemented:
        if use_refind_polygon is True, only return values in the
        refined polygon
        """

        if not self._attribute_dic.has_key(column_name):
            msg = 'There is no column called %s!' % column_name
            raise TitleValueError, msg

        return self._attribute_dic[column_name]

    ##
    # @brief ??
    # @param value_column_name ??
    # @param known_column_name ??
    # @param known_values ??
    # @param use_refind_polygon ??
    def get_value(self, value_column_name, known_column_name,
                  known_values, use_refind_polygon=False):
        """
        Do linear interpolation on the known_colum, using the known_value,
        to return a value of the column_value_name.
        """

        pass

    ##
    # @brief Get a geospatial object that describes the locations.
    # @param use_refind_polygon Unused??
    def get_location(self, use_refind_polygon=False):
        """
        Return a geospatial object which describes the
        locations of the location file.

        Note, if there is not location info, this returns None.

        Not implemented:
        if use_refind_polygon is True, only return values in the
        refined polygon
        """

        return self._geospatial

    ##
    # @brief Add column to 'end' of CSV data.
    # @param column_name The new column name.
    # @param column_values The new column values.
    # @param overwrite If True, overwrites last column, doesn't add at end.
    def set_column(self, column_name, column_values, overwrite=False):
        """
        Add a column to the 'end' (with the right most column being the end)
        of the csv file.

        Set overwrite to True if you want to overwrite a column.

        Note, in column_name white space is removed and case is not checked.
        Precondition
        The column_name and column_values cannot have comma's in it.
        """

        # sanity checks
        value_row_count = \
                len(self._attribute_dic[self._title_index_dic.keys()[0]])
        if len(column_values) <> value_row_count:
            msg = 'The number of column values must equal the number of rows.'
            raise DataMissingValuesError, msg

        # check new column name isn't already used, and we aren't overwriting
        if self._attribute_dic.has_key(column_name):
            if not overwrite:
                msg = 'Column name %s already in use!' % column_name
                raise TitleValueError, msg
        else:
            # New title.  Add it to the title index.
            self._title_index_dic[column_name] = len(self._title_index_dic)

        self._attribute_dic[column_name] = column_values

    ##
    # @brief Save the exposure CSV  file.
    # @param file_name If supplied, use this filename, not original.
    def save(self, file_name=None):
        """
        Save the exposure csv file
        """

        if file_name is None:
            file_name = self._file_name

        fd = open(file_name, 'wb')
        writer = csv.writer(fd)

        #Write the title to a cvs file
        line = [None] * len(self._title_index_dic)
        for title in self._title_index_dic.iterkeys():
            line[self._title_index_dic[title]] = title
        writer.writerow(line)

        # Write the values to a cvs file
        value_row_count = \
                len(self._attribute_dic[self._title_index_dic.keys()[0]])
        for row_i in range(value_row_count):
            line = [None] * len(self._title_index_dic)
            for title in self._title_index_dic.iterkeys():
                line[self._title_index_dic[title]] = \
                     self._attribute_dic[title][row_i]
            writer.writerow(line)


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
    
    See underlying function csv2dict for more details.
    """

    X, _ = csv2dict(file_name)

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
# @brief Convert CSV file to a dictionary of arrays.
# @param file_name The path to the file to read.
def csv2array(file_name):
    """
    Convert CSV files of the form:

    time, discharge, velocity
    0.0,  1.2,       0.0
    0.1,  3.2,       1.1
    ...

    to a dictionary of numeric arrays.


    See underlying function csv2dict for more details.
    """

    X, _ = csv2dict(file_name)

    Y = {}
    for key in X.keys():
        Y[key] = num.array([float(x) for x in X[key]])

    return Y


##
# @brief Read a CSV file and convert to a dictionary of {key: column}.
# @param file_name The path to the file to read.
# @param title_check_list List of titles that *must* be columns in the file.
# @return Two dicts: ({key:column}, {title:index}).
# @note WARNING: Values are returned as strings.
def csv2dict(file_name, title_check_list=None):
    """
    Load in the csv as a dictionary, title as key and column info as value.
    Also, create a dictionary, title as key and column index as value,
    to keep track of the column order.

    Two dictionaries are returned.

    WARNING: Values are returned as strings.
             Do this to change a list of strings to a list of floats
                 time = [float(x) for x in time]
    """

    # FIXME(Ole): Consider dealing with files without headers
    # FIXME(Ole): Consider a wrapper automatically converting text fields
    #             to the right type by trying for: int, float, string
    
    attribute_dic = {}
    title_index_dic = {}
    titles_stripped = [] # List of titles

    reader = csv.reader(file(file_name))

    # Read in and manipulate the title info
    titles = reader.next()
    for i, title in enumerate(titles):
        header = title.strip()
        titles_stripped.append(header)
        title_index_dic[header] = i
    title_count = len(titles_stripped)

    # Check required columns
    if title_check_list is not None:
        for title_check in title_check_list:
            if not title_index_dic.has_key(title_check):
                msg = 'Reading error. This row is not present %s' % title_check
                raise IOError, msg

    # Create a dictionary of column values, indexed by column title
    for line in reader:
        n = len(line) # Number of entries
        if n != title_count:
            msg = 'Entry in file %s had %d columns ' % (file_name, n)
            msg += 'although there were %d headers' % title_count
            raise IOError, msg
        for i, value in enumerate(line):
            attribute_dic.setdefault(titles_stripped[i], []).append(value)

    return attribute_dic, title_index_dic


##
# @brief 
# @param filename 
# @param x 
# @param y 
# @param z 
def write_obj(filename, x, y, z):
    """Store x,y,z vectors into filename (obj format).

       Vectors are assumed to have dimension (M,3) where
       M corresponds to the number elements.
       triangles are assumed to be disconnected

       The three numbers in each vector correspond to three vertices,

       e.g. the x coordinate of vertex 1 of element i is in x[i,1]
    """

    import os.path

    root, ext = os.path.splitext(filename)
    if ext == '.obj':
        FN = filename
    else:
        FN = filename + '.obj'

    outfile = open(FN, 'wb')
    outfile.write("# Triangulation as an obj file\n")

    M, N = x.shape
    assert N == 3  #Assuming three vertices per element

    for i in range(M):
        for j in range(N):
            outfile.write("v %f %f %f\n" % (x[i,j], y[i,j], z[i,j]))

    for i in range(M):
        base = i * N
        outfile.write("f %d %d %d\n" % (base+1, base+2, base+3))

    outfile.close()


#########################################################
#Conversion routines
########################################################

##
# @brief Convert SWW data to OBJ data.
# @param basefilename Stem of filename, needs size and extension added.
# @param size The number of lines to write.
def sww2obj(basefilename, size):
    """Convert netcdf based data output to obj
    """

    from Scientific.IO.NetCDF import NetCDFFile

    # Get NetCDF
    FN = create_filename('.', basefilename, 'sww', size)
    log.critical('Reading from %s' % FN)
    fid = NetCDFFile(FN, netcdf_mode_r)  #Open existing file for read

    # Get the variables
    x = fid.variables['x']
    y = fid.variables['y']
    z = fid.variables['elevation']
    time = fid.variables['time']
    stage = fid.variables['stage']

    M = size  #Number of lines
    xx = num.zeros((M,3), num.float)
    yy = num.zeros((M,3), num.float)
    zz = num.zeros((M,3), num.float)

    for i in range(M):
        for j in range(3):
            xx[i,j] = x[i+j*M]
            yy[i,j] = y[i+j*M]
            zz[i,j] = z[i+j*M]

    # Write obj for bathymetry
    FN = create_filename('.', basefilename, 'obj', size)
    write_obj(FN,xx,yy,zz)

    # Now read all the data with variable information, combine with
    # x,y info and store as obj
    for k in range(len(time)):
        t = time[k]
        log.critical('Processing timestep %f' % t)

        for i in range(M):
            for j in range(3):
                zz[i,j] = stage[k,i+j*M]

        #Write obj for variable data
        #FN = create_filename(basefilename, 'obj', size, time=t)
        FN = create_filename('.', basefilename[:5], 'obj', size, time=t)
        write_obj(FN, xx, yy, zz)


##
# @brief 
# @param basefilename Stem of filename, needs size and extension added.
def dat2obj(basefilename):
    """Convert line based data output to obj
    FIXME: Obsolete?
    """

    import glob, os
    from anuga.config import data_dir

    # Get bathymetry and x,y's
    lines = open(data_dir+os.sep+basefilename+'_geometry.dat', 'r').readlines()

    M = len(lines)  #Number of lines
    x = num.zeros((M,3), num.float)
    y = num.zeros((M,3), num.float)
    z = num.zeros((M,3), num.float)

    for i, line in enumerate(lines):
        tokens = line.split()
        values = map(float, tokens)

        for j in range(3):
            x[i,j] = values[j*3]
            y[i,j] = values[j*3+1]
            z[i,j] = values[j*3+2]

    # Write obj for bathymetry
    write_obj(data_dir + os.sep + basefilename + '_geometry', x, y, z)

    # Now read all the data files with variable information, combine with
    # x,y info and store as obj.

    files = glob.glob(data_dir + os.sep + basefilename + '*.dat')
    for filename in files:
        log.critical('Processing %s' % filename)

        lines = open(data_dir + os.sep + filename, 'r').readlines()
        assert len(lines) == M
        root, ext = os.path.splitext(filename)

        # Get time from filename
        i0 = filename.find('_time=')
        if i0 == -1:
            #Skip bathymetry file
            continue

        i0 += 6  #Position where time starts
        i1 = filename.find('.dat')

        if i1 > i0:
            t = float(filename[i0:i1])
        else:
            raise DataTimeError, 'Hmmmm'

        for i, line in enumerate(lines):
            tokens = line.split()
            values = map(float,tokens)

            for j in range(3):
                z[i,j] = values[j]

        # Write obj for variable data
        write_obj(data_dir + os.sep + basefilename + '_time=%.4f' % t, x, y, z)


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
# @brief Convert DEM data  to PTS data.
# @param basename_in Stem of input filename.
# @param basename_out Stem of output filename.
# @param easting_min 
# @param easting_max 
# @param northing_min 
# @param northing_max 
# @param use_cache 
# @param verbose 
# @return 
def dem2pts(basename_in, basename_out=None,
            easting_min=None, easting_max=None,
            northing_min=None, northing_max=None,
            use_cache=False, verbose=False,):
    """Read Digitial Elevation model from the following NetCDF format (.dem)

    Example:

    ncols         3121
    nrows         1800
    xllcorner     722000
    yllcorner     5893000
    cellsize      25
    NODATA_value  -9999
    138.3698 137.4194 136.5062 135.5558 ..........

    Convert to NetCDF pts format which is

    points:  (Nx2) float array
    elevation: N float array
    """

    kwargs = {'basename_out': basename_out,
              'easting_min': easting_min,
              'easting_max': easting_max,
              'northing_min': northing_min,
              'northing_max': northing_max,
              'verbose': verbose}

    if use_cache is True:
        from caching import cache
        result = cache(_dem2pts, basename_in, kwargs,
                       dependencies = [basename_in + '.dem'],
                       verbose = verbose)

    else:
        result = apply(_dem2pts, [basename_in], kwargs)

    return result


##
# @brief 
# @param basename_in 
# @param basename_out 
# @param verbose 
# @param easting_min 
# @param easting_max 
# @param northing_min 
# @param northing_max 
def _dem2pts(basename_in, basename_out=None, verbose=False,
            easting_min=None, easting_max=None,
            northing_min=None, northing_max=None):
    """Read Digitial Elevation model from the following NetCDF format (.dem)

    Internal function. See public function dem2pts for details.
    """

    # FIXME: Can this be written feasibly using write_pts?

    import os
    from Scientific.IO.NetCDF import NetCDFFile

    root = basename_in

    # Get NetCDF
    infile = NetCDFFile(root + '.dem', netcdf_mode_r) 

    if verbose: log.critical('Reading DEM from %s' % (root + '.dem'))

    ncols = infile.ncols[0]
    nrows = infile.nrows[0]
    xllcorner = infile.xllcorner[0]  # Easting of lower left corner
    yllcorner = infile.yllcorner[0]  # Northing of lower left corner
    cellsize = infile.cellsize[0]
    NODATA_value = infile.NODATA_value[0]
    dem_elevation = infile.variables['elevation']

    zone = infile.zone[0]
    false_easting = infile.false_easting[0]
    false_northing = infile.false_northing[0]

    # Text strings
    projection = infile.projection
    datum = infile.datum
    units = infile.units

    # Get output file
    if basename_out == None:
        ptsname = root + '.pts'
    else:
        ptsname = basename_out + '.pts'

    if verbose: log.critical('Store to NetCDF file %s' % ptsname)

    # NetCDF file definition
    outfile = NetCDFFile(ptsname, netcdf_mode_w)

    # Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = 'NetCDF pts format for compact and portable ' \
                          'storage of spatial point data'

    # Assign default values
    if easting_min is None: easting_min = xllcorner
    if easting_max is None: easting_max = xllcorner + ncols*cellsize
    if northing_min is None: northing_min = yllcorner
    if northing_max is None: northing_max = yllcorner + nrows*cellsize

    # Compute offsets to update georeferencing
    easting_offset = xllcorner - easting_min
    northing_offset = yllcorner - northing_min

    # Georeferencing
    outfile.zone = zone
    outfile.xllcorner = easting_min # Easting of lower left corner
    outfile.yllcorner = northing_min # Northing of lower left corner
    outfile.false_easting = false_easting
    outfile.false_northing = false_northing

    outfile.projection = projection
    outfile.datum = datum
    outfile.units = units

    # Grid info (FIXME: probably not going to be used, but heck)
    outfile.ncols = ncols
    outfile.nrows = nrows

    dem_elevation_r = num.reshape(dem_elevation, (nrows, ncols))
    totalnopoints = nrows*ncols

    # Calculating number of NODATA_values for each row in clipped region
    # FIXME: use array operations to do faster
    nn = 0
    k = 0
    i1_0 = 0
    j1_0 = 0
    thisj = 0
    thisi = 0
    for i in range(nrows):
        y = (nrows-i-1)*cellsize + yllcorner
        for j in range(ncols):
            x = j*cellsize + xllcorner
            if easting_min <= x <= easting_max \
               and northing_min <= y <= northing_max:
                thisj = j
                thisi = i
                if dem_elevation_r[i,j] == NODATA_value:
                    nn += 1

                if k == 0:
                    i1_0 = i
                    j1_0 = j

                k += 1

    index1 = j1_0
    index2 = thisj

    # Dimension definitions
    nrows_in_bounding_box = int(round((northing_max-northing_min)/cellsize))
    ncols_in_bounding_box = int(round((easting_max-easting_min)/cellsize))

    clippednopoints = (thisi+1-i1_0)*(thisj+1-j1_0)
    nopoints = clippednopoints-nn

    clipped_dem_elev = dem_elevation_r[i1_0:thisi+1,j1_0:thisj+1]

    if verbose:
        log.critical('There are %d values in the elevation' % totalnopoints)
        log.critical('There are %d values in the clipped elevation'
                     % clippednopoints)
        log.critical('There are %d NODATA_values in the clipped elevation' % nn)

    outfile.createDimension('number_of_points', nopoints)
    outfile.createDimension('number_of_dimensions', 2) #This is 2d data

    # Variable definitions
    outfile.createVariable('points', netcdf_float, ('number_of_points',
                                                    'number_of_dimensions'))
    outfile.createVariable('elevation', netcdf_float, ('number_of_points',))

    # Get handles to the variables
    points = outfile.variables['points']
    elevation = outfile.variables['elevation']

    lenv = index2-index1+1

    # Store data
    global_index = 0
    # for i in range(nrows):
    for i in range(i1_0, thisi+1, 1):
        if verbose and i % ((nrows+10)/10) == 0:
            log.critical('Processing row %d of %d' % (i, nrows))

        lower_index = global_index

        v = dem_elevation_r[i,index1:index2+1]
        no_NODATA = num.sum(v == NODATA_value)
        if no_NODATA > 0:
            newcols = lenv - no_NODATA  # ncols_in_bounding_box - no_NODATA
        else:
            newcols = lenv              # ncols_in_bounding_box

        telev = num.zeros(newcols, num.float)
        tpoints = num.zeros((newcols, 2), num.float)

        local_index = 0

        y = (nrows-i-1)*cellsize + yllcorner
        #for j in range(ncols):
        for j in range(j1_0,index2+1,1):
            x = j*cellsize + xllcorner
            if easting_min <= x <= easting_max \
               and northing_min <= y <= northing_max \
               and dem_elevation_r[i,j] <> NODATA_value:
                tpoints[local_index, :] = [x-easting_min, y-northing_min]
                telev[local_index] = dem_elevation_r[i, j]
                global_index += 1
                local_index += 1

        upper_index = global_index

        if upper_index == lower_index + newcols:
            points[lower_index:upper_index, :] = tpoints
            elevation[lower_index:upper_index] = telev

    assert global_index == nopoints, 'index not equal to number of points'

    infile.close()
    outfile.close()


##
# @brief  Return block of surface lines for each cross section
# @param lines Iterble  of text lines to process.
# @note BROKEN?  UNUSED?
def _read_hecras_cross_sections(lines):
    """Return block of surface lines for each cross section
    Starts with SURFACE LINE,
    Ends with END CROSS-SECTION
    """

    points = []

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


##
# @brief Convert HECRAS (.sdf) file to PTS file.
# @param basename_in Sterm of input filename.
# @param basename_out Sterm of output filename.
# @param verbose True if this function is to be verbose.
def hecras_cross_sections2pts(basename_in,
                              basename_out=None,
                              verbose=False):
    """Read HEC-RAS Elevation datal from the following ASCII format (.sdf)

    Example:

# RAS export file created on Mon 15Aug2005 11:42
# by HEC-RAS Version 3.1.1

BEGIN HEADER:
  UNITS: METRIC
  DTM TYPE: TIN
  DTM: v:\1\cit\perth_topo\river_tin
  STREAM LAYER: c:\local\hecras\21_02_03\up_canning_cent3d.shp
  CROSS-SECTION LAYER: c:\local\hecras\21_02_03\up_can_xs3d.shp
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
    from Scientific.IO.NetCDF import NetCDFFile

    root = basename_in

    # Get ASCII file
    infile = open(root + '.sdf', 'r')

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
    if basename_out == None:
        ptsname = root + '.pts'
    else:
        ptsname = basename_out + '.pts'

    geo_ref = Geo_reference(zone, 0, 0, datum, projection, units)
    geo = Geospatial_data(points, {"elevation":elevation},
                          verbose=verbose, geo_reference=geo_ref)
    geo.export_points_file(ptsname)


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
# @brief 
# @param production_dirs
# @param output_dir
# @param scenario_name
# @param gauges_dir_name
# @param plot_quantity
# @param generate_fig
# @param reportname
# @param surface
# @param time_min
# @param time_max
# @param title_on
# @param verbose
# @param nodes
def get_timeseries(production_dirs, output_dir, scenario_name, gauges_dir_name,
                   plot_quantity, generate_fig=False,
                   reportname=None, surface=False, time_min=None,
                   time_max=None, title_on=False, verbose=True,
                   nodes=None):
    """
    nodes - number of processes used.

    warning - this function has no tests
    """

    if reportname == None:
        report = False
    else:
        report = True

    if nodes is None:
        is_parallel = False
    else:
        is_parallel = True

    # Generate figures
    swwfiles = {}
    if is_parallel is True:
        for i in range(nodes):
            log.critical('Sending node %d of %d' % (i, nodes))
            swwfiles = {}
            if not reportname == None:
                reportname = report_name + '_%s' % i
            for label_id in production_dirs.keys():
                if label_id == 'boundaries':
                    swwfile = best_boundary_sww
                else:
                    file_loc = output_dir + label_id + sep
                    sww_extra = '_P%s_%s' % (i, nodes)
                    swwfile = file_loc + scenario_name + sww_extra + '.sww'
                    log.critical('swwfile %s' % swwfile)
                    swwfiles[swwfile] = label_id

                texname, elev_output = sww2timeseries(swwfiles,
                                              gauges_dir_name,
                                              production_dirs,
                                              report=report,
                                              reportname=reportname,
                                              plot_quantity=plot_quantity,
                                              generate_fig=generate_fig,
                                              surface=surface,
                                              time_min=time_min,
                                              time_max=time_max,
                                              title_on=title_on,
                                              verbose=verbose)
    else:
        for label_id in production_dirs.keys():
            if label_id == 'boundaries':
                log.critical('boundaries')
                file_loc = project.boundaries_in_dir
                swwfile = project.boundaries_dir_name3 + '.sww'
                #  swwfile = boundary_dir_filename
            else:
                file_loc = output_dir + label_id + sep
                swwfile = file_loc + scenario_name + '.sww'
            swwfiles[swwfile] = label_id

        texname, elev_output = sww2timeseries(swwfiles,
                                              gauges_dir_name,
                                              production_dirs,
                                              report=report,
                                              reportname=reportname,
                                              plot_quantity=plot_quantity,
                                              generate_fig=generate_fig,
                                              surface=surface,
                                              time_min=time_min,
                                              time_max=time_max,
                                              title_on=title_on,
                                              verbose=verbose)


##
# @brief Convert SWW file to DEM file (.asc or .ers).
# @param basename_in 
# @param basename_out 
# @param quantity 
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
def sww2dem(basename_in, basename_out=None,
            quantity=None, # defaults to elevation
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
            format='ers',
            block_size=None):
    """Read SWW file and convert to Digitial Elevation model format
    (.asc or .ers)

    Example (ASC):
    ncols         3121
    nrows         1800
    xllcorner     722000
    yllcorner     5893000
    cellsize      25
    NODATA_value  -9999
    138.3698 137.4194 136.5062 135.5558 ..........

    The number of decimal places can be specified by the user to save
    on disk space requirements by specifying in the call to sww2dem.

    Also write accompanying file with same basename_in but extension .prj
    used to fix the UTM zone, datum, false northings and eastings.

    The prj format is assumed to be as

    Projection    UTM
    Zone          56
    Datum         WGS84
    Zunits        NO
    Units         METERS
    Spheroid      WGS84
    Xshift        0.0000000000
    Yshift        10000000.0000000000
    Parameters


    The parameter quantity must be the name of an existing quantity or
    an expression involving existing quantities. The default is
    'elevation'. Quantity is not a list of quantities.

    if timestep (an index) is given, output quantity at that timestep

    if reduction is given and its an index, output quantity at that timestep. If reduction is given
    and is a built in function, use that to reduce quantity over all timesteps.

    datum

    format can be either 'asc' or 'ers'
    block_size - sets the number of slices along the non-time axis to
                 process in one block.
    """

    import sys
    import types

    from anuga.utilities.polygon import inside_polygon, outside_polygon
    from anuga.abstract_2d_finite_volumes.util import \
         apply_expression_to_dictionary

    msg = 'Format must be either asc or ers'
    assert format.lower() in ['asc', 'ers'], msg

    false_easting = 500000
    false_northing = 10000000

    if quantity is None:
        quantity = 'elevation'
    
    if reduction is None:
        reduction = max

    if basename_out is None:
        basename_out = basename_in + '_%s' % quantity

    if quantity_formula.has_key(quantity):
        quantity = quantity_formula[quantity]

    if number_of_decimal_places is None:
        number_of_decimal_places = 3

    if block_size is None:
        block_size = DEFAULT_BLOCK_SIZE

    # Read SWW file
    swwfile = basename_in + '.sww'
    demfile = basename_out + '.' + format

    # Read sww file
    if verbose:
        log.critical('Reading from %s' % swwfile)
        log.critical('Output directory is %s' % basename_out)

    from Scientific.IO.NetCDF import NetCDFFile
    fid = NetCDFFile(swwfile)

    #Get extent and reference
    x = fid.variables['x'][:]
    y = fid.variables['y'][:]
    volumes = fid.variables['volumes'][:]
    if type(reduction) is not types.BuiltinFunctionType:
        times = fid.variables['time'][reduction]
    else:
        times = fid.variables['time'][:]

    number_of_timesteps = fid.dimensions['number_of_timesteps']
    number_of_points = fid.dimensions['number_of_points']

    if origin is None:
        # Get geo_reference
        # sww files don't have to have a geo_ref
        try:
            geo_reference = Geo_reference(NetCDFObject=fid)
        except AttributeError, e:
            geo_reference = Geo_reference() # Default georef object

        xllcorner = geo_reference.get_xllcorner()
        yllcorner = geo_reference.get_yllcorner()
        zone = geo_reference.get_zone()
    else:
        zone = origin[0]
        xllcorner = origin[1]
        yllcorner = origin[2]

    # FIXME: Refactor using code from Interpolation_function.statistics
    # (in interpolate.py)
    # Something like print swwstats(swwname)
    if verbose:
        log.critical('------------------------------------------------')
        log.critical('Statistics of SWW file:')
        log.critical('  Name: %s' % swwfile)
        log.critical('  Reference:')
        log.critical('    Lower left corner: [%f, %f]' % (xllcorner, yllcorner))
        if type(reduction) is not types.BuiltinFunctionType:
            log.critical('    Time: %f' % times)
        else:
            log.critical('    Start time: %f' % fid.starttime[0])
        log.critical('  Extent:')
        log.critical('    x [m] in [%f, %f], len(x) == %d'
                     %(num.min(x), num.max(x), len(x.flat)))
        log.critical('    y [m] in [%f, %f], len(y) == %d'
                     % (num.min(y), num.max(y), len(y.flat)))
        if type(reduction) is not types.BuiltinFunctionType:
            log.critical('    t [s] = %f, len(t) == %d' % (times, 1))
        else:
            log.critical('    t [s] in [%f, %f], len(t) == %d'
                         % (min(times), max(times), len(times)))
        log.critical('  Quantities [SI units]:')
        
        # Comment out for reduced memory consumption
        for name in ['stage', 'xmomentum', 'ymomentum']:
            q = fid.variables[name][:].flatten()
            if type(reduction) is not types.BuiltinFunctionType:
                q = q[reduction*len(x):(reduction+1)*len(x)]
            if verbose: log.critical('    %s in [%f, %f]'
                                     % (name, min(q), max(q)))
        for name in ['elevation']:
            q = fid.variables[name][:].flatten()
            if verbose: log.critical('    %s in [%f, %f]'
                                     % (name, min(q), max(q)))

    # Get the variables in the supplied expression.
    # This may throw a SyntaxError exception.
    var_list = get_vars_in_expression(quantity)

    # Check that we have the required variables in the SWW file.
    missing_vars = []
    for name in var_list:
        try:
            _ = fid.variables[name]
        except:
            missing_vars.append(name)
    if missing_vars:
        msg = ("In expression '%s', variables %s are not in the SWW file '%s'"
               % (quantity, swwfile))
        raise Exception, msg

    # Create result array and start filling, block by block.
    result = num.zeros(number_of_points, num.float)

    for start_slice in xrange(0, number_of_points, block_size):
        # Limit slice size to array end if at last block
        end_slice = min(start_slice + block_size, number_of_points)
        
        # Get slices of all required variables
        q_dict = {}
        for name in var_list:
            # check if variable has time axis
            if len(fid.variables[name].shape) == 2:
                q_dict[name] = fid.variables[name][:,start_slice:end_slice]
            else:       # no time axis
                q_dict[name] = fid.variables[name][start_slice:end_slice]

        # Evaluate expression with quantities found in SWW file
        res = apply_expression_to_dictionary(quantity, q_dict)

        if len(res.shape) == 2:
            new_res = num.zeros(res.shape[1], num.float)
            for k in xrange(res.shape[1]):
                if type(reduction) is not types.BuiltinFunctionType:
                    new_res[k] = res[reduction,k]
                else:
                    new_res[k] = reduction(res[:,k])
            res = new_res

        result[start_slice:end_slice] = res
                                    
    # Post condition: Now q has dimension: number_of_points
    assert len(result.shape) == 1
    assert result.shape[0] == number_of_points

    if verbose:
        log.critical('Processed values for %s are in [%f, %f]'
                     % (quantity, min(result), max(result)))

    # Create grid and update xll/yll corner and x,y
    # Relative extent
    if easting_min is None:
        xmin = min(x)
    else:
        xmin = easting_min - xllcorner

    if easting_max is None:
        xmax = max(x)
    else:
        xmax = easting_max - xllcorner

    if northing_min is None:
        ymin = min(y)
    else:
        ymin = northing_min - yllcorner

    if northing_max is None:
        ymax = max(y)
    else:
        ymax = northing_max - yllcorner

    msg = 'xmax must be greater than or equal to xmin.\n'
    msg += 'I got xmin = %f, xmax = %f' %(xmin, xmax)
    assert xmax >= xmin, msg

    msg = 'ymax must be greater than or equal to xmin.\n'
    msg += 'I got ymin = %f, ymax = %f' %(ymin, ymax)
    assert ymax >= ymin, msg

    if verbose: log.critical('Creating grid')
    ncols = int((xmax-xmin)/cellsize) + 1
    nrows = int((ymax-ymin)/cellsize) + 1

    # New absolute reference and coordinates
    newxllcorner = xmin + xllcorner
    newyllcorner = ymin + yllcorner

    x = x + xllcorner - newxllcorner
    y = y + yllcorner - newyllcorner

    vertex_points = num.concatenate ((x[:,num.newaxis], y[:,num.newaxis]), axis=1)
    assert len(vertex_points.shape) == 2

    grid_points = num.zeros ((ncols*nrows, 2), num.float)

    for i in xrange(nrows):
        if format.lower() == 'asc':
            yg = i * cellsize
        else:
            # this will flip the order of the y values for ers
            yg = (nrows-i) * cellsize

        for j in xrange(ncols):
            xg = j * cellsize
            k = i*ncols + j

            grid_points[k, 0] = xg
            grid_points[k, 1] = yg

    # Interpolate
    from anuga.fit_interpolate.interpolate import Interpolate

    # Remove loners from vertex_points, volumes here
    vertex_points, volumes = remove_lone_verts(vertex_points, volumes)
    # export_mesh_file('monkey.tsh',{'vertices':vertex_points, 'triangles':volumes})
    interp = Interpolate(vertex_points, volumes, verbose = verbose)

    # Interpolate using quantity values
    if verbose: log.critical('Interpolating')
    grid_values = interp.interpolate(result, grid_points).flatten()

    if verbose:
        log.critical('Interpolated values are in [%f, %f]'
                     % (num.min(grid_values), num.max(grid_values)))

    # Assign NODATA_value to all points outside bounding polygon (from interpolation mesh)
    P = interp.mesh.get_boundary_polygon()
    outside_indices = outside_polygon(grid_points, P, closed=True)

    for i in outside_indices:
        grid_values[i] = NODATA_value

    if format.lower() == 'ers':
        # setup ERS header information
        grid_values = num.reshape(grid_values, (nrows, ncols))
        header = {}
        header['datum'] = '"' + datum + '"'
        # FIXME The use of hardwired UTM and zone number needs to be made optional
        # FIXME Also need an automatic test for coordinate type (i.e. EN or LL)
        header['projection'] = '"UTM-' + str(zone) + '"'
        header['coordinatetype'] = 'EN'
        if header['coordinatetype'] == 'LL':
            header['longitude'] = str(newxllcorner)
            header['latitude'] = str(newyllcorner)
        elif header['coordinatetype'] == 'EN':
            header['eastings'] = str(newxllcorner)
            header['northings'] = str(newyllcorner)
        header['nullcellvalue'] = str(NODATA_value)
        header['xdimension'] = str(cellsize)
        header['ydimension'] = str(cellsize)
        header['value'] = '"' + quantity + '"'
        #header['celltype'] = 'IEEE8ByteReal'  #FIXME: Breaks unit test

        #Write
        if verbose: log.critical('Writing %s' % demfile)

        import ermapper_grids

        ermapper_grids.write_ermapper_grid(demfile, grid_values, header)

        fid.close()
    else:
        #Write to Ascii format
        #Write prj file
        prjfile = basename_out + '.prj'

        if verbose: log.critical('Writing %s' % prjfile)
        prjid = open(prjfile, 'w')
        prjid.write('Projection    %s\n' %'UTM')
        prjid.write('Zone          %d\n' %zone)
        prjid.write('Datum         %s\n' %datum)
        prjid.write('Zunits        NO\n')
        prjid.write('Units         METERS\n')
        prjid.write('Spheroid      %s\n' %datum)
        prjid.write('Xshift        %d\n' %false_easting)
        prjid.write('Yshift        %d\n' %false_northing)
        prjid.write('Parameters\n')
        prjid.close()

        if verbose: log.critical('Writing %s' % demfile)

        ascid = open(demfile, 'w')

        ascid.write('ncols         %d\n' %ncols)
        ascid.write('nrows         %d\n' %nrows)
        ascid.write('xllcorner     %d\n' %newxllcorner)
        ascid.write('yllcorner     %d\n' %newyllcorner)
        ascid.write('cellsize      %f\n' %cellsize)
        ascid.write('NODATA_value  %d\n' %NODATA_value)

        #Get bounding polygon from mesh
        #P = interp.mesh.get_boundary_polygon()
        #inside_indices = inside_polygon(grid_points, P)

        # change printoptions so that a long string of zeros in not
        # summarized as [0.0, 0.0, 0.0, ... 0.0, 0.0, 0.0]
        #printoptions = num.get_printoptions()
        #num.set_printoptions(threshold=sys.maxint)

        format = '%.'+'%g' % number_of_decimal_places +'e'
        for i in range(nrows):
            if verbose and i % ((nrows+10)/10) == 0:
                log.critical('Doing row %d of %d' % (i, nrows))

            base_index = (nrows-i-1)*ncols

            slice = grid_values[base_index:base_index+ncols]

            num.savetxt(ascid, slice.reshape(1,ncols), format, ' ' )
            
        
        #Close
        ascid.close()
        fid.close()

        return basename_out

################################################################################
# Obsolete functions - maintained for backwards compatibility
################################################################################

##
# @brief 
# @param basename_in 
# @param basename_out 
# @param quantity 
# @param timestep 
# @param reduction 
# @param cellsize 
# @param verbose 
# @param origin 
# @note OBSOLETE - use sww2dem() instead.
def sww2asc(basename_in, basename_out = None,
            quantity = None,
            reduction = None,
            cellsize = 10,
            verbose = False,
            origin = None):
    
    log.critical('sww2asc will soon be obsoleted - please use sww2dem')
    sww2dem(basename_in,
            basename_out = basename_out,
            quantity = quantity,
            reduction = reduction,
            cellsize = cellsize,
            number_of_decimal_places = number_of_decimal_places,
            verbose = verbose,
            origin = origin,
            datum = 'WGS84',
            format = 'asc')


##
# @brief 
# @param basename_in
# @param basename_out
# @param quantity
# @param timestep
# @param reduction
# @param cellsize
# @param verbose
# @param origin
# @param datum
# @note OBSOLETE - use sww2dem() instead.
def sww2ers(basename_in, basename_out=None,
            quantity=None,
            reduction=None,
            cellsize=10,
            verbose=False,
            origin=None,
            datum='WGS84'):
    log.critical('sww2ers will soon be obsoleted - please use sww2dem')
    sww2dem(basename_in,
            basename_out=basename_out,
            quantity=quantity,
            reduction=reduction,
            cellsize=cellsize,
            number_of_decimal_places=number_of_decimal_places,
            verbose=verbose,
            origin=origin,
            datum=datum,
            format='ers')

################################################################################
# End of obsolete functions
################################################################################


##
# @brief Convert SWW file to PTS file (at selected points).
# @param basename_in Stem name of input SWW file.
# @param basename_out Stem name of output file.
# @param data_points If given, points where quantity is to be computed.
# @param quantity Name (or expression) of existing quantity(s) (def: elevation).
# @param timestep If given, output quantity at that timestep.
# @param reduction If given, reduce quantity by this factor.
# @param NODATA_value The NODATA value (default -9999).
# @param verbose True if this function is to be verbose.
# @param origin ??
def sww2pts(basename_in, basename_out=None,
            data_points=None,
            quantity=None,
            timestep=None,
            reduction=None,
            NODATA_value=-9999,
            verbose=False,
            origin=None):
    """Read SWW file and convert to interpolated values at selected points

    The parameter 'quantity' must be the name of an existing quantity or
    an expression involving existing quantities. The default is 'elevation'.

    if timestep (an index) is given, output quantity at that timestep.

    if reduction is given use that to reduce quantity over all timesteps.

    data_points (Nx2 array) give locations of points where quantity is to 
    be computed.
    """

    import sys
    from anuga.utilities.polygon import inside_polygon, outside_polygon
    from anuga.abstract_2d_finite_volumes.util import \
             apply_expression_to_dictionary
    from anuga.geospatial_data.geospatial_data import Geospatial_data

    if quantity is None:
        quantity = 'elevation'

    if reduction is None:
        reduction = max

    if basename_out is None:
        basename_out = basename_in + '_%s' % quantity

    swwfile = basename_in + '.sww'
    ptsfile = basename_out + '.pts'

    # Read sww file
    if verbose: log.critical('Reading from %s' % swwfile)
    from Scientific.IO.NetCDF import NetCDFFile
    fid = NetCDFFile(swwfile)

    # Get extent and reference
    x = fid.variables['x'][:]
    y = fid.variables['y'][:]
    volumes = fid.variables['volumes'][:]

    number_of_timesteps = fid.dimensions['number_of_timesteps']
    number_of_points = fid.dimensions['number_of_points']
    if origin is None:
        # Get geo_reference
        # sww files don't have to have a geo_ref
        try:
            geo_reference = Geo_reference(NetCDFObject=fid)
        except AttributeError, e:
            geo_reference = Geo_reference() # Default georef object

        xllcorner = geo_reference.get_xllcorner()
        yllcorner = geo_reference.get_yllcorner()
        zone = geo_reference.get_zone()
    else:
        zone = origin[0]
        xllcorner = origin[1]
        yllcorner = origin[2]

    # FIXME: Refactor using code from file_function.statistics
    # Something like print swwstats(swwname)
    if verbose:
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        times = fid.variables['time'][:]
        log.critical('------------------------------------------------')
        log.critical('Statistics of SWW file:')
        log.critical('  Name: %s' % swwfile)
        log.critical('  Reference:')
        log.critical('    Lower left corner: [%f, %f]' % (xllcorner, yllcorner))
        log.critical('    Start time: %f' % fid.starttime[0])
        log.critical('  Extent:')
        log.critical('    x [m] in [%f, %f], len(x) == %d'
                     % (num.min(x), num.max(x), len(x.flat)))
        log.critical('    y [m] in [%f, %f], len(y) == %d'
                     % (num.min(y), num.max(y), len(y.flat)))
        log.critical('    t [s] in [%f, %f], len(t) == %d'
                     % (min(times), max(times), len(times)))
        log.critical('  Quantities [SI units]:')
        for name in ['stage', 'xmomentum', 'ymomentum', 'elevation']:
            q = fid.variables[name][:].flat
            log.critical('    %s in [%f, %f]' % (name, min(q), max(q)))

    # Get quantity and reduce if applicable
    if verbose: log.critical('Processing quantity %s' % quantity)

    # Turn NetCDF objects into numeric arrays
    quantity_dict = {}
    for name in fid.variables.keys():
        quantity_dict[name] = fid.variables[name][:]

    # Convert quantity expression to quantities found in sww file
    q = apply_expression_to_dictionary(quantity, quantity_dict)

    if len(q.shape) == 2:
        # q has a time component and needs to be reduced along
        # the temporal dimension
        if verbose: log.critical('Reducing quantity %s' % quantity)

        q_reduced = num.zeros(number_of_points, num.float)
        for k in range(number_of_points):
            q_reduced[k] = reduction(q[:,k])
        q = q_reduced

    # Post condition: Now q has dimension: number_of_points
    assert len(q.shape) == 1
    assert q.shape[0] == number_of_points

    if verbose:
        log.critical('Processed values for %s are in [%f, %f]'
                     % (quantity, min(q), max(q)))

    # Create grid and update xll/yll corner and x,y
    vertex_points = num.concatenate((x[:, num.newaxis], y[:, num.newaxis]), axis=1)
    assert len(vertex_points.shape) == 2

    # Interpolate
    from anuga.fit_interpolate.interpolate import Interpolate
    interp = Interpolate(vertex_points, volumes, verbose=verbose)

    # Interpolate using quantity values
    if verbose: log.critical('Interpolating')
    interpolated_values = interp.interpolate(q, data_points).flatten()

    if verbose:
        log.critical('Interpolated values are in [%f, %f]'
                     % (num.min(interpolated_values),
                        num.max(interpolated_values)))

    # Assign NODATA_value to all points outside bounding polygon
    # (from interpolation mesh)
    P = interp.mesh.get_boundary_polygon()
    outside_indices = outside_polygon(data_points, P, closed=True)

    for i in outside_indices:
        interpolated_values[i] = NODATA_value

    # Store results
    G = Geospatial_data(data_points=data_points, attributes=interpolated_values)

    G.export_points_file(ptsfile, absolute = True)

    fid.close()


##
# @brief Convert ASC file to DEM file.
# @param basename_in Stem of input filename.
# @param basename_out Stem of output filename.
# @param use_cache ??
# @param verbose True if this function is to be verbose.
# @return 
# @note A PRJ file with same stem basename must exist and is used to fix the
#       UTM zone, datum, false northings and eastings.
def convert_dem_from_ascii2netcdf(basename_in, basename_out=None,
                                  use_cache=False,
                                  verbose=False):
    """Read Digital Elevation model from the following ASCII format (.asc)

    Example:
    ncols         3121
    nrows         1800
    xllcorner     722000
    yllcorner     5893000
    cellsize      25
    NODATA_value  -9999
    138.3698 137.4194 136.5062 135.5558 ..........

    Convert basename_in + '.asc' to NetCDF format (.dem)
    mimicking the ASCII format closely.

    An accompanying file with same basename_in but extension .prj must exist
    and is used to fix the UTM zone, datum, false northings and eastings.

    The prj format is assumed to be as

    Projection    UTM
    Zone          56
    Datum         WGS84
    Zunits        NO
    Units         METERS
    Spheroid      WGS84
    Xshift        0.0000000000
    Yshift        10000000.0000000000
    Parameters
    """

    kwargs = {'basename_out': basename_out, 'verbose': verbose}

    if use_cache is True:
        from caching import cache
        result = cache(_convert_dem_from_ascii2netcdf, basename_in, kwargs,
                       dependencies=[basename_in + '.asc',
                                     basename_in + '.prj'],
                       verbose=verbose)

    else:
        result = apply(_convert_dem_from_ascii2netcdf, [basename_in], kwargs)

    return result


##
# @brief Convert an ASC file to a DEM file.
# @param basename_in Stem of input filename.
# @param basename_out Stem of output filename.
# @param verbose True if this function is to be verbose.
def _convert_dem_from_ascii2netcdf(basename_in, basename_out = None,
                                   verbose = False):
    """Read Digital Elevation model from the following ASCII format (.asc)

    Internal function. See public function convert_dem_from_ascii2netcdf
    for details.
    """

    import os
    from Scientific.IO.NetCDF import NetCDFFile

    root = basename_in

    # Read Meta data
    if verbose: log.critical('Reading METADATA from %s' % (root + '.prj'))

    metadatafile = open(root + '.prj')
    metalines = metadatafile.readlines()
    metadatafile.close()

    L = metalines[0].strip().split()
    assert L[0].strip().lower() == 'projection'
    projection = L[1].strip()                   #TEXT

    L = metalines[1].strip().split()
    assert L[0].strip().lower() == 'zone'
    zone = int(L[1].strip())

    L = metalines[2].strip().split()
    assert L[0].strip().lower() == 'datum'
    datum = L[1].strip()                        #TEXT

    L = metalines[3].strip().split()
    assert L[0].strip().lower() == 'zunits'     #IGNORE
    zunits = L[1].strip()                       #TEXT

    L = metalines[4].strip().split()
    assert L[0].strip().lower() == 'units'
    units = L[1].strip()                        #TEXT

    L = metalines[5].strip().split()
    assert L[0].strip().lower() == 'spheroid'   #IGNORE
    spheroid = L[1].strip()                     #TEXT

    L = metalines[6].strip().split()
    assert L[0].strip().lower() == 'xshift'
    false_easting = float(L[1].strip())

    L = metalines[7].strip().split()
    assert L[0].strip().lower() == 'yshift'
    false_northing = float(L[1].strip())

    #Read DEM data
    datafile = open(basename_in + '.asc')

    if verbose: log.critical('Reading DEM from %s' % (basename_in + '.asc'))

    lines = datafile.readlines()
    datafile.close()

    if verbose: log.critical('Got %d lines' % len(lines))

    ncols = int(lines[0].split()[1].strip())
    nrows = int(lines[1].split()[1].strip())

    # Do cellsize (line 4) before line 2 and 3
    cellsize = float(lines[4].split()[1].strip())

    # Checks suggested by Joaquim Luis
    # Our internal representation of xllcorner
    # and yllcorner is non-standard.
    xref = lines[2].split()
    if xref[0].strip() == 'xllcorner':
        xllcorner = float(xref[1].strip()) # + 0.5*cellsize # Correct offset
    elif xref[0].strip() == 'xllcenter':
        xllcorner = float(xref[1].strip())
    else:
        msg = 'Unknown keyword: %s' % xref[0].strip()
        raise Exception, msg

    yref = lines[3].split()
    if yref[0].strip() == 'yllcorner':
        yllcorner = float(yref[1].strip()) # + 0.5*cellsize # Correct offset
    elif yref[0].strip() == 'yllcenter':
        yllcorner = float(yref[1].strip())
    else:
        msg = 'Unknown keyword: %s' % yref[0].strip()
        raise Exception, msg

    NODATA_value = int(lines[5].split()[1].strip())

    assert len(lines) == nrows + 6

    if basename_out == None:
        netcdfname = root + '.dem'
    else:
        netcdfname = basename_out + '.dem'

    if verbose: log.critical('Store to NetCDF file %s' % netcdfname)

    # NetCDF file definition
    fid = NetCDFFile(netcdfname, netcdf_mode_w)

    #Create new file
    fid.institution = 'Geoscience Australia'
    fid.description = 'NetCDF DEM format for compact and portable storage ' \
                      'of spatial point data'

    fid.ncols = ncols
    fid.nrows = nrows
    fid.xllcorner = xllcorner
    fid.yllcorner = yllcorner
    fid.cellsize = cellsize
    fid.NODATA_value = NODATA_value

    fid.zone = zone
    fid.false_easting = false_easting
    fid.false_northing = false_northing
    fid.projection = projection
    fid.datum = datum
    fid.units = units

    # dimension definitions
    fid.createDimension('number_of_rows', nrows)
    fid.createDimension('number_of_columns', ncols)

    # variable definitions
    fid.createVariable('elevation', netcdf_float, ('number_of_rows',
                                                   'number_of_columns'))

    # Get handles to the variables
    elevation = fid.variables['elevation']

    #Store data
    n = len(lines[6:])
    for i, line in enumerate(lines[6:]):
        fields = line.split()
        if verbose and i % ((n+10)/10) == 0:
            log.critical('Processing row %d of %d' % (i, nrows))
           
        if len(fields) != ncols:
            msg = 'Wrong number of columns in file "%s" line %d\n' % (basename_in + '.asc', i)
            msg += 'I got %d elements, but there should have been %d\n' % (len(fields), ncols)
            raise Exception, msg

        elevation[i, :] = num.array([float(x) for x in fields])

    fid.close()


##
# @brief Convert 'ferret' file to SWW file.
# @param basename_in Stem of input filename.
# @param basename_out Stem of output filename.
# @param verbose True if this function is to be verbose.
# @param minlat 
# @param maxlat 
# @param minlon 
# @param maxlon 
# @param mint 
# @param maxt 
# @param mean_stage 
# @param origin 
# @param zscale 
# @param fail_on_NaN 
# @param NaN_filler 
# @param elevation 
# @param inverted_bathymetry 
def ferret2sww(basename_in, basename_out=None,
               verbose=False,
               minlat=None, maxlat=None,
               minlon=None, maxlon=None,
               mint=None, maxt=None, mean_stage=0,
               origin=None, zscale=1,
               fail_on_NaN=True,
               NaN_filler=0,
               elevation=None,
               inverted_bathymetry=True
               ): #FIXME: Bathymetry should be obtained
                                  #from MOST somehow.
                                  #Alternatively from elsewhere
                                  #or, as a last resort,
                                  #specified here.
                                  #The value of -100 will work
                                  #for the Wollongong tsunami
                                  #scenario but is very hacky
    """Convert MOST and 'Ferret' NetCDF format for wave propagation to
    sww format native to abstract_2d_finite_volumes.

    Specify only basename_in and read files of the form
    basefilename_ha.nc, basefilename_ua.nc, basefilename_va.nc containing
    relative height, x-velocity and y-velocity, respectively.

    Also convert latitude and longitude to UTM. All coordinates are
    assumed to be given in the GDA94 datum.

    min's and max's: If omitted - full extend is used.
    To include a value min may equal it, while max must exceed it.
    Lat and lon are assuemd to be in decimal degrees

    origin is a 3-tuple with geo referenced
    UTM coordinates (zone, easting, northing)

    nc format has values organised as HA[TIME, LATITUDE, LONGITUDE]
    which means that longitude is the fastest
    varying dimension (row major order, so to speak)

    ferret2sww uses grid points as vertices in a triangular grid
    counting vertices from lower left corner upwards, then right
    """

    import os
    from Scientific.IO.NetCDF import NetCDFFile

    precision = num.float

    msg = 'Must use latitudes and longitudes for minlat, maxlon etc'

    if minlat != None:
        assert -90 < minlat < 90 , msg
    if maxlat != None:
        assert -90 < maxlat < 90 , msg
        if minlat != None:
            assert maxlat > minlat
    if minlon != None:
        assert -180 < minlon < 180 , msg
    if maxlon != None:
        assert -180 < maxlon < 180 , msg
        if minlon != None:
            assert maxlon > minlon

    # Get NetCDF data
    if verbose: log.critical('Reading files %s_*.nc' % basename_in)

    file_h = NetCDFFile(basename_in + '_ha.nc', netcdf_mode_r) # Wave amplitude (cm)
    file_u = NetCDFFile(basename_in + '_ua.nc', netcdf_mode_r) # Velocity (x) (cm/s)
    file_v = NetCDFFile(basename_in + '_va.nc', netcdf_mode_r) # Velocity (y) (cm/s)
    file_e = NetCDFFile(basename_in + '_e.nc', netcdf_mode_r)  # Elevation (z) (m)

    if basename_out is None:
        swwname = basename_in + '.sww'
    else:
        swwname = basename_out + '.sww'

    # Get dimensions of file_h
    for dimension in file_h.dimensions.keys():
        if dimension[:3] == 'LON':
            dim_h_longitude = dimension
        if dimension[:3] == 'LAT':
            dim_h_latitude = dimension
        if dimension[:4] == 'TIME':
            dim_h_time = dimension

    times = file_h.variables[dim_h_time]
    latitudes = file_h.variables[dim_h_latitude]
    longitudes = file_h.variables[dim_h_longitude]

    kmin, kmax, lmin, lmax = _get_min_max_indexes(latitudes[:],
                                                  longitudes[:],
                                                  minlat, maxlat,
                                                  minlon, maxlon)
    # get dimensions for file_e
    for dimension in file_e.dimensions.keys():
        if dimension[:3] == 'LON':
            dim_e_longitude = dimension
        if dimension[:3] == 'LAT':
            dim_e_latitude = dimension

    # get dimensions for file_u
    for dimension in file_u.dimensions.keys():
        if dimension[:3] == 'LON':
            dim_u_longitude = dimension
        if dimension[:3] == 'LAT':
            dim_u_latitude = dimension
        if dimension[:4] == 'TIME':
            dim_u_time = dimension

    # get dimensions for file_v
    for dimension in file_v.dimensions.keys():
        if dimension[:3] == 'LON':
            dim_v_longitude = dimension
        if dimension[:3] == 'LAT':
            dim_v_latitude = dimension
        if dimension[:4] == 'TIME':
            dim_v_time = dimension

    # Precision used by most for lat/lon is 4 or 5 decimals
    e_lat = num.around(file_e.variables[dim_e_latitude][:], 5)
    e_lon = num.around(file_e.variables[dim_e_longitude][:], 5)

    # Check that files are compatible
    assert num.allclose(latitudes, file_u.variables[dim_u_latitude])
    assert num.allclose(latitudes, file_v.variables[dim_v_latitude])
    assert num.allclose(latitudes, e_lat)

    assert num.allclose(longitudes, file_u.variables[dim_u_longitude])
    assert num.allclose(longitudes, file_v.variables[dim_v_longitude])
    assert num.allclose(longitudes, e_lon)

    if mint is None:
        jmin = 0
        mint = times[0]
    else:
        jmin = num.searchsorted(times, mint)
        
        # numpy.int32 didn't work in slicing of amplitude below
        jmin = int(jmin)

    if maxt is None:
        jmax = len(times)
        maxt = times[-1]
    else:
        jmax = num.searchsorted(times, maxt)
        
        # numpy.int32 didn't work in slicing of amplitude below
        jmax = int(jmax)        

    kmin, kmax, lmin, lmax = _get_min_max_indexes(latitudes[:],
                                                  longitudes[:],
                                                  minlat, maxlat,
                                                  minlon, maxlon)


    times = times[jmin:jmax]
    latitudes = latitudes[kmin:kmax]
    longitudes = longitudes[lmin:lmax]

    if verbose: log.critical('cropping')

    zname = 'ELEVATION'

    amplitudes = file_h.variables['HA'][jmin:jmax, kmin:kmax, lmin:lmax]
    uspeed = file_u.variables['UA'][jmin:jmax, kmin:kmax, lmin:lmax] #Lon
    vspeed = file_v.variables['VA'][jmin:jmax, kmin:kmax, lmin:lmax] #Lat
    elevations = file_e.variables[zname][kmin:kmax, lmin:lmax]

    #    if latitudes2[0]==latitudes[0] and latitudes2[-1]==latitudes[-1]:
    #        elevations = file_e.variables['ELEVATION'][kmin:kmax, lmin:lmax]
    #    elif latitudes2[0]==latitudes[-1] and latitudes2[-1]==latitudes[0]:
    #        from numpy import asarray
    #        elevations=elevations.tolist()
    #        elevations.reverse()
    #        elevations=asarray(elevations)
    #    else:
    #        from numpy import asarray
    #        elevations=elevations.tolist()
    #        elevations.reverse()
    #        elevations=asarray(elevations)
    #        'print hmmm'

    # Get missing values
    nan_ha = file_h.variables['HA'].missing_value[0]
    nan_ua = file_u.variables['UA'].missing_value[0]
    nan_va = file_v.variables['VA'].missing_value[0]
    if hasattr(file_e.variables[zname],'missing_value'):
        nan_e  = file_e.variables[zname].missing_value[0]
    else:
        nan_e = None

    # Cleanup
    missing = (amplitudes == nan_ha)
    if num.sometrue (missing):
        if fail_on_NaN:
            msg = 'NetCDFFile %s contains missing values' \
                  % basename_in + '_ha.nc'
            raise DataMissingValuesError, msg
        else:
            amplitudes = amplitudes*(missing==0) + missing*NaN_filler

    missing = (uspeed == nan_ua)
    if num.sometrue (missing):
        if fail_on_NaN:
            msg = 'NetCDFFile %s contains missing values' \
                  % basename_in + '_ua.nc'
            raise DataMissingValuesError, msg
        else:
            uspeed = uspeed*(missing==0) + missing*NaN_filler

    missing = (vspeed == nan_va)
    if num.sometrue (missing):
        if fail_on_NaN:
            msg = 'NetCDFFile %s contains missing values' \
                  % basename_in + '_va.nc'
            raise DataMissingValuesError, msg
        else:
            vspeed = vspeed*(missing==0) + missing*NaN_filler

    missing = (elevations == nan_e)
    if num.sometrue (missing):
        if fail_on_NaN:
            msg = 'NetCDFFile %s contains missing values' \
                  % basename_in + '_e.nc'
            raise DataMissingValuesError, msg
        else:
            elevations = elevations*(missing==0) + missing*NaN_filler

    number_of_times = times.shape[0]
    number_of_latitudes = latitudes.shape[0]
    number_of_longitudes = longitudes.shape[0]

    assert amplitudes.shape[0] == number_of_times
    assert amplitudes.shape[1] == number_of_latitudes
    assert amplitudes.shape[2] == number_of_longitudes

    if verbose:
        log.critical('------------------------------------------------')
        log.critical('Statistics:')
        log.critical('  Extent (lat/lon):')
        log.critical('    lat in [%f, %f], len(lat) == %d'
                     % (num.min(latitudes), num.max(latitudes),
                        len(latitudes.flat)))
        log.critical('    lon in [%f, %f], len(lon) == %d'
                     % (num.min(longitudes), num.max(longitudes),
                        len(longitudes.flat)))
        log.critical('    t in [%f, %f], len(t) == %d'
                     % (num.min(times), num.max(times), len(times.flat)))

#        q = amplitudes.flatten()
        name = 'Amplitudes (ha) [cm]'
        log.critical('  %s in [%f, %f]'
                     % (name, num.min(amplitudes), num.max(amplitudes)))

#        q = uspeed.flatten()
        name = 'Speeds (ua) [cm/s]'
        log.critical('  %s in [%f, %f]'
                     % (name, num.min(uspeed), num.max(uspeed)))

#        q = vspeed.flatten()
        name = 'Speeds (va) [cm/s]'
        log.critical('  %s in [%f, %f]'
                     % (name, num.min(vspeed), num.max(vspeed)))

#        q = elevations.flatten()
        name = 'Elevations (e) [m]'
        log.critical('  %s in [%f, %f]'
                     % (name, num.min(elevations), num.max(elevations)))

    # print number_of_latitudes, number_of_longitudes
    number_of_points = number_of_latitudes * number_of_longitudes
    number_of_volumes = (number_of_latitudes-1) * (number_of_longitudes-1) * 2

    file_h.close()
    file_u.close()
    file_v.close()
    file_e.close()

    # NetCDF file definition
    outfile = NetCDFFile(swwname, netcdf_mode_w)

    description = 'Converted from Ferret files: %s, %s, %s, %s' \
                  % (basename_in + '_ha.nc',
                     basename_in + '_ua.nc',
                     basename_in + '_va.nc',
                     basename_in + '_e.nc')

    # Create new file
    starttime = times[0]

    sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])
    sww.store_header(outfile, times, number_of_volumes,
                     number_of_points, description=description,
                     verbose=verbose, sww_precision=netcdf_float)

    # Store
    from anuga.coordinate_transforms.redfearn import redfearn
    x = num.zeros(number_of_points, num.float)  #Easting
    y = num.zeros(number_of_points, num.float)  #Northing

    if verbose: log.critical('Making triangular grid')

    # Check zone boundaries
    refzone, _, _ = redfearn(latitudes[0], longitudes[0])

    vertices = {}
    i = 0
    for k, lat in enumerate(latitudes):       # Y direction
        for l, lon in enumerate(longitudes):  # X direction
            vertices[l,k] = i

            zone, easting, northing = redfearn(lat,lon)

            #msg = 'Zone boundary crossed at longitude =', lon
            #assert zone == refzone, msg
            #print '%7.2f %7.2f %8.2f %8.2f' %(lon, lat, easting, northing)
            x[i] = easting
            y[i] = northing
            i += 1

    #Construct 2 triangles per 'rectangular' element
    volumes = []
    for l in range(number_of_longitudes-1):    # X direction
        for k in range(number_of_latitudes-1): # Y direction
            v1 = vertices[l,k+1]
            v2 = vertices[l,k]
            v3 = vertices[l+1,k+1]
            v4 = vertices[l+1,k]

            volumes.append([v1,v2,v3]) #Upper element
            volumes.append([v4,v3,v2]) #Lower element

    volumes = num.array(volumes, num.int)      #array default#

    if origin is None:
        origin = Geo_reference(refzone, min(x), min(y))
    geo_ref = write_NetCDF_georeference(origin, outfile)

    if elevation is not None:
        z = elevation
    else:
        if inverted_bathymetry:
            z = -1 * elevations
        else:
            z = elevations
    #FIXME: z should be obtained from MOST and passed in here

    #FIXME use the Write_sww instance(sww) to write this info
    z = num.resize(z, outfile.variables['elevation'][:].shape)
    outfile.variables['x'][:] = x - geo_ref.get_xllcorner()
    outfile.variables['y'][:] = y - geo_ref.get_yllcorner()
    #outfile.variables['z'][:] = z             #FIXME HACK for bacwards compat.
    outfile.variables['elevation'][:] = z
    outfile.variables['volumes'][:] = volumes.astype(num.int32) #For Opteron 64

    #Time stepping
    stage = outfile.variables['stage']
    xmomentum = outfile.variables['xmomentum']
    ymomentum = outfile.variables['ymomentum']

    if verbose: log.critical('Converting quantities')

    n = len(times)
    for j in range(n):
        if verbose and j % ((n+10)/10) == 0:
            log.critical('  Doing %d of %d' % (j, n))

        i = 0
        for k in range(number_of_latitudes):      # Y direction
            for l in range(number_of_longitudes): # X direction
                w = zscale * amplitudes[j,k,l] / 100 + mean_stage
                stage[j,i] = w
                h = w - z[i]
                xmomentum[j,i] = uspeed[j,k,l]/100*h
                ymomentum[j,i] = vspeed[j,k,l]/100*h
                i += 1

    #outfile.close()

    #FIXME: Refactor using code from file_function.statistics
    #Something like print swwstats(swwname)
    if verbose:
        x = outfile.variables['x'][:]
        y = outfile.variables['y'][:]
        log.critical('------------------------------------------------')
        log.critical('Statistics of output file:')
        log.critical('  Name: %s' %swwname)
        log.critical('  Reference:')
        log.critical('    Lower left corner: [%f, %f]'
                     % (geo_ref.get_xllcorner(), geo_ref.get_yllcorner()))
        log.critical(' Start time: %f' %starttime)
        log.critical('    Min time: %f' %mint)
        log.critical('    Max time: %f' %maxt)
        log.critical('  Extent:')
        log.critical('    x [m] in [%f, %f], len(x) == %d'
                     % (num.min(x), num.max(x), len(x.flat)))
        log.critical('    y [m] in [%f, %f], len(y) == %d'
                     % (num.min(y), num.max(y), len(y.flat)))
        log.critical('    t [s] in [%f, %f], len(t) == %d'
                     % (min(times), max(times), len(times)))
        log.critical('  Quantities [SI units]:')
        for name in ['stage', 'xmomentum', 'ymomentum', 'elevation']:
            q = outfile.variables[name][:]    # .flatten()
            log.critical('    %s in [%f, %f]' % (name, num.min(q), num.max(q)))

    outfile.close()


##
# @brief Convert time-series text file to TMS file.
# @param filename 
# @param quantity_names 
# @param time_as_seconds 
def timefile2netcdf(filename, quantity_names=None, time_as_seconds=False):
    """Template for converting typical text files with time series to
    NetCDF tms file.

    The file format is assumed to be either two fields separated by a comma:

        time [DD/MM/YY hh:mm:ss], value0 value1 value2 ...

    E.g

      31/08/04 00:00:00, 1.328223 0 0
      31/08/04 00:15:00, 1.292912 0 0

    or time (seconds), value0 value1 value2 ...

      0.0, 1.328223 0 0
      0.1, 1.292912 0 0

    will provide a time dependent function f(t) with three attributes

    filename is assumed to be the rootname with extenisons .txt and .sww
    """

    import time, calendar
    from anuga.config import time_format
    from anuga.utilities.numerical_tools import ensure_numeric

    file_text = filename + '.txt'
    fid = open(file_text)
    line = fid.readline()
    fid.close()

    fields = line.split(',')
    msg = "File %s must have the format 'datetime, value0 value1 value2 ...'" \
          % file_text
    assert len(fields) == 2, msg

    if not time_as_seconds:
        try:
            starttime = calendar.timegm(time.strptime(fields[0], time_format))
        except ValueError:
            msg = 'First field in file %s must be' % file_text
            msg += ' date-time with format %s.\n' % time_format
            msg += 'I got %s instead.' % fields[0]
            raise DataTimeError, msg
    else:
        try:
            starttime = float(fields[0])
        except Error:
            msg = "Bad time format"
            raise DataTimeError, msg

    # Split values
    values = []
    for value in fields[1].split():
        values.append(float(value))

    q = ensure_numeric(values)

    msg = 'ERROR: File must contain at least one independent value'
    assert len(q.shape) == 1, msg

    # Read times proper
    from anuga.config import time_format
    import time, calendar

    fid = open(file_text)
    lines = fid.readlines()
    fid.close()

    N = len(lines)
    d = len(q)

    T = num.zeros(N, num.float)       # Time
    Q = num.zeros((N, d), num.float)  # Values

    for i, line in enumerate(lines):
        fields = line.split(',')
        if not time_as_seconds:
            realtime = calendar.timegm(time.strptime(fields[0], time_format))
        else:
             realtime = float(fields[0])
        T[i] = realtime - starttime

        for j, value in enumerate(fields[1].split()):
            Q[i, j] = float(value)

    msg = 'File %s must list time as a monotonuosly ' % filename
    msg += 'increasing sequence'
    assert num.alltrue(T[1:] - T[:-1] > 0), msg

    #Create NetCDF file
    from Scientific.IO.NetCDF import NetCDFFile

    fid = NetCDFFile(filename + '.tms', netcdf_mode_w)

    fid.institution = 'Geoscience Australia'
    fid.description = 'Time series'

    #Reference point
    #Start time in seconds since the epoch (midnight 1/1/1970)
    #FIXME: Use Georef
    fid.starttime = starttime

    # dimension definitions
    #fid.createDimension('number_of_volumes', self.number_of_volumes)
    #fid.createDimension('number_of_vertices', 3)

    fid.createDimension('number_of_timesteps', len(T))

    fid.createVariable('time', netcdf_float, ('number_of_timesteps',))

    fid.variables['time'][:] = T

    for i in range(Q.shape[1]):
        try:
            name = quantity_names[i]
        except:
            name = 'Attribute%d' % i

        fid.createVariable(name, netcdf_float, ('number_of_timesteps',))
        fid.variables[name][:] = Q[:,i]

    fid.close()


##
# @brief Get the extents of a NetCDF data file.
# @param file_name The path to the SWW file.
# @return A list of x, y, z and stage limits (min, max).
def extent_sww(file_name):
    """Read in an sww file.

    Input:
    file_name - the sww file

    Output:
    A list: [min(x),max(x),min(y),max(y),min(stage.flat),max(stage.flat)]
    """

    from Scientific.IO.NetCDF import NetCDFFile

    #Get NetCDF
    fid = NetCDFFile(file_name, netcdf_mode_r)

    # Get the variables
    x = fid.variables['x'][:]
    y = fid.variables['y'][:]
    stage = fid.variables['stage'][:]

    fid.close()

    return [min(x), max(x), min(y), max(y), num.min(stage), num.max(stage)]


##
# @brief 
# @param filename
# @param boundary
# @param t
# @param fail_if_NaN
# @param NaN_filler
# @param verbose
# @param very_verbose
# @return 
def sww2domain(filename, boundary=None, t=None,
               fail_if_NaN=True, NaN_filler=0,
               verbose=False, very_verbose=False):
    """
    Usage: domain = sww2domain('file.sww',t=time (default = last time in file))

    Boundary is not recommended if domain.smooth is not selected, as it
    uses unique coordinates, but not unique boundaries. This means that
    the boundary file will not be compatable with the coordinates, and will
    give a different final boundary, or crash.
    """
    
    from Scientific.IO.NetCDF import NetCDFFile
    from shallow_water import Domain

    # initialise NaN.
    NaN = 9.969209968386869e+036

    if verbose: log.critical('Reading from %s' % filename)

    fid = NetCDFFile(filename, netcdf_mode_r)    # Open existing file for read
    time = fid.variables['time']       # Timesteps
    if t is None:
        t = time[-1]
    time_interp = get_time_interp(time,t)

    # Get the variables as numeric arrays
    x = fid.variables['x'][:]                   # x-coordinates of vertices
    y = fid.variables['y'][:]                   # y-coordinates of vertices
    elevation = fid.variables['elevation']      # Elevation
    stage = fid.variables['stage']              # Water level
    xmomentum = fid.variables['xmomentum']      # Momentum in the x-direction
    ymomentum = fid.variables['ymomentum']      # Momentum in the y-direction

    starttime = fid.starttime[0]
    volumes = fid.variables['volumes'][:]       # Connectivity
    coordinates = num.transpose(num.asarray([x.tolist(), y.tolist()]))
    # FIXME (Ole): Something like this might be better:
    #                 concatenate((x, y), axis=1)
    # or              concatenate((x[:,num.newaxis], x[:,num.newaxis]), axis=1)

    conserved_quantities = []
    interpolated_quantities = {}
    other_quantities = []

    # get geo_reference
    try:                             # sww files don't have to have a geo_ref
        geo_reference = Geo_reference(NetCDFObject=fid)
    except: # AttributeError, e:
        geo_reference = None

    if verbose: log.critical('    getting quantities')

    for quantity in fid.variables.keys():
        dimensions = fid.variables[quantity].dimensions
        if 'number_of_timesteps' in dimensions:
            conserved_quantities.append(quantity)
            interpolated_quantities[quantity] = \
                  interpolated_quantity(fid.variables[quantity][:], time_interp)
        else:
            other_quantities.append(quantity)

    other_quantities.remove('x')
    other_quantities.remove('y')
    #other_quantities.remove('z')
    other_quantities.remove('volumes')
    try:
        other_quantities.remove('stage_range')
        other_quantities.remove('xmomentum_range')
        other_quantities.remove('ymomentum_range')
        other_quantities.remove('elevation_range')
    except:
        pass

    conserved_quantities.remove('time')

    if verbose: log.critical('    building domain')

    #    From domain.Domain:
    #    domain = Domain(coordinates, volumes,\
    #                    conserved_quantities = conserved_quantities,\
    #                    other_quantities = other_quantities,zone=zone,\
    #                    xllcorner=xllcorner, yllcorner=yllcorner)

    # From shallow_water.Domain:
    coordinates = coordinates.tolist()
    volumes = volumes.tolist()
    # FIXME:should this be in mesh? (peter row)
    if fid.smoothing == 'Yes':
        unique = False
    else:
        unique = True
    if unique:
        coordinates, volumes, boundary = weed(coordinates, volumes,boundary)

      
    
    try:
        domain = Domain(coordinates, volumes, boundary)
    except AssertionError, e:
        fid.close()
        msg = 'Domain could not be created: %s. ' \
              'Perhaps use "fail_if_NaN=False and NaN_filler = ..."' % e
        raise DataDomainError, msg

    if not boundary is None:
        domain.boundary = boundary

    domain.geo_reference = geo_reference

    domain.starttime = float(starttime) + float(t)
    domain.time = 0.0

    for quantity in other_quantities:
        try:
            NaN = fid.variables[quantity].missing_value
        except:
            pass                       # quantity has no missing_value number
        X = fid.variables[quantity][:]
        if very_verbose:
            log.critical('       %s' % str(quantity))
            log.critical('        NaN = %s' % str(NaN))
            log.critical('        max(X)')
            log.critical('       %s' % str(max(X)))
            log.critical('        max(X)==NaN')
            log.critical('       %s' % str(max(X)==NaN))
            log.critical('')
        if max(X) == NaN or min(X) == NaN:
            if fail_if_NaN:
                msg = 'quantity "%s" contains no_data entry' % quantity
                raise DataMissingValuesError, msg
            else:
                data = (X != NaN)
                X = (X*data) + (data==0)*NaN_filler
        if unique:
            X = num.resize(X, (len(X)/3, 3))
        domain.set_quantity(quantity, X)
    #
    for quantity in conserved_quantities:
        try:
            NaN = fid.variables[quantity].missing_value
        except:
            pass                       # quantity has no missing_value number
        X = interpolated_quantities[quantity]
        if very_verbose:
            log.critical('       %s' % str(quantity))
            log.critical('        NaN = %s' % str(NaN))
            log.critical('        max(X)')
            log.critical('       %s' % str(max(X)))
            log.critical('        max(X)==NaN')
            log.critical('       %s' % str(max(X)==NaN))
            log.critical('')
        if max(X) == NaN or min(X) == NaN:
            if fail_if_NaN:
                msg = 'quantity "%s" contains no_data entry' % quantity
                raise DataMissingValuesError, msg
            else:
                data = (X != NaN)
                X = (X*data) + (data==0)*NaN_filler
        if unique:
            X = num.resize(X, (X.shape[0]/3, 3))
        domain.set_quantity(quantity, X)

    fid.close()

    return domain


##
# @brief Interpolate a quantity wrt time.
# @param saved_quantity The quantity to interpolate.
# @param time_interp (index, ratio)
# @return A vector of interpolated values.
def interpolated_quantity(saved_quantity, time_interp):
    '''Given an index and ratio, interpolate quantity with respect to time.'''

    index, ratio = time_interp

    Q = saved_quantity

    if ratio > 0:
        q = (1-ratio)*Q[index] + ratio*Q[index+1]
    else:
        q = Q[index]

    #Return vector of interpolated values
    return q


##
# @brief 
# @parm time 
# @param t 
# @return An (index, ration) tuple.
def get_time_interp(time, t=None):
    #Finds the ratio and index for time interpolation.
    #It is borrowed from previous abstract_2d_finite_volumes code.
    if t is None:
        t=time[-1]
        index = -1
        ratio = 0.
    else:
        T = time
        tau = t
        index=0
        msg = 'Time interval derived from file %s [%s:%s]' \
              % ('FIXMEfilename', T[0], T[-1])
        msg += ' does not match model time: %s' % tau
        if tau < time[0]: raise DataTimeError, msg
        if tau > time[-1]: raise DataTimeError, msg
        while tau > time[index]: index += 1
        while tau < time[index]: index -= 1
        if tau == time[index]:
            #Protect against case where tau == time[-1] (last time)
            # - also works in general when tau == time[i]
            ratio = 0
        else:
            #t is now between index and index+1
            ratio = (tau - time[index])/(time[index+1] - time[index])

    return (index, ratio)


##
# @brief 
# @param coordinates 
# @param volumes 
# @param boundary 
def weed(coordinates, volumes, boundary=None):
    if isinstance(coordinates, num.ndarray):
        coordinates = coordinates.tolist()
    if isinstance(volumes, num.ndarray):
        volumes = volumes.tolist()

    unique = False
    point_dict = {}
    same_point = {}
    for i in range(len(coordinates)):
        point = tuple(coordinates[i])
        if point_dict.has_key(point):
            unique = True
            same_point[i] = point
            #to change all point i references to point j
        else:
            point_dict[point] = i
            same_point[i] = point

    coordinates = []
    i = 0
    for point in point_dict.keys():
        point = tuple(point)
        coordinates.append(list(point))
        point_dict[point] = i
        i += 1

    for volume in volumes:
        for i in range(len(volume)):
            index = volume[i]
            if index > -1:
                volume[i] = point_dict[same_point[index]]

    new_boundary = {}
    if not boundary is None:
        for segment in boundary.keys():
            point0 = point_dict[same_point[segment[0]]]
            point1 = point_dict[same_point[segment[1]]]
            label = boundary[segment]
            #FIXME should the bounday attributes be concaterated
            #('exterior, pond') or replaced ('pond')(peter row)

            if new_boundary.has_key((point0, point1)):
                new_boundary[(point0,point1)] = new_boundary[(point0, point1)]

            elif new_boundary.has_key((point1, point0)):
                new_boundary[(point1,point0)] = new_boundary[(point1, point0)]
            else: new_boundary[(point0, point1)] = label

        boundary = new_boundary

    return coordinates, volumes, boundary


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


##
# @brief 
# @param filename
# @param verbose
def tsh2sww(filename, verbose=False):
    """
    to check if a tsh/msh file 'looks' good.
    """

    if verbose == True: log.critical('Creating domain from %s' % filename)

    domain = pmesh_to_domain_instance(filename, Domain)

    if verbose == True: log.critical("Number of triangles = %s" % len(domain))

    domain.smooth = True
    domain.format = 'sww'   #Native netcdf visualisation format
    file_path, filename = path.split(filename)
    filename, ext = path.splitext(filename)
    domain.set_name(filename)
    domain.reduction = mean

    if verbose == True: log.critical("file_path = %s" % file_path)

    if file_path == "":
        file_path = "."
    domain.set_datadir(file_path)

    if verbose == True:
        log.critical("Output written to %s%s%s.%s"
                     % (domain.get_datadir(), sep, domain.get_name(),
                        domain.format))

    sww = SWW_file(domain)
    sww.store_connectivity()
    sww.store_timestep()


##
# @brief Convert CSIRO ESRI file to an SWW boundary file.
# @param bath_dir 
# @param elevation_dir 
# @param ucur_dir 
# @param vcur_dir 
# @param sww_file 
# @param minlat 
# @param maxlat 
# @param minlon 
# @param maxlon 
# @param zscale 
# @param mean_stage 
# @param fail_on_NaN 
# @param elevation_NaN_filler 
# @param bath_prefix 
# @param elevation_prefix 
# @param verbose 
# @note Also convert latitude and longitude to UTM. All coordinates are
#       assumed to be given in the GDA94 datum.
def asc_csiro2sww(bath_dir,
                  elevation_dir,
                  ucur_dir,
                  vcur_dir,
                  sww_file,
                  minlat=None, maxlat=None,
                  minlon=None, maxlon=None,
                  zscale=1,
                  mean_stage=0,
                  fail_on_NaN=True,
                  elevation_NaN_filler=0,
                  bath_prefix='ba',
                  elevation_prefix='el',
                  verbose=False):
    """
    Produce an sww boundary file, from esri ascii data from CSIRO.

    Also convert latitude and longitude to UTM. All coordinates are
    assumed to be given in the GDA94 datum.

    assume:
    All files are in esri ascii format

    4 types of information
    bathymetry
    elevation
    u velocity
    v velocity

    Assumptions
    The metadata of all the files is the same
    Each type is in a seperate directory
    One bath file with extention .000
    The time period is less than 24hrs and uniform.
    """

    from Scientific.IO.NetCDF import NetCDFFile

    from anuga.coordinate_transforms.redfearn import redfearn

    precision = netcdf_float # So if we want to change the precision its done here

    # go in to the bath dir and load the only file,
    bath_files = os.listdir(bath_dir)
    bath_file = bath_files[0]
    bath_dir_file =  bath_dir + os.sep + bath_file
    bath_metadata, bath_grid =  _read_asc(bath_dir_file)

    #Use the date.time of the bath file as a basis for
    #the start time for other files
    base_start = bath_file[-12:]

    #go into the elevation dir and load the 000 file
    elevation_dir_file = elevation_dir  + os.sep + elevation_prefix \
                         + base_start

    elevation_files = os.listdir(elevation_dir)
    ucur_files = os.listdir(ucur_dir)
    vcur_files = os.listdir(vcur_dir)
    elevation_files.sort()

    # the first elevation file should be the
    # file with the same base name as the bath data
    assert elevation_files[0] == 'el' + base_start

    number_of_latitudes = bath_grid.shape[0]
    number_of_longitudes = bath_grid.shape[1]
    number_of_volumes = (number_of_latitudes-1) * (number_of_longitudes-1) * 2

    longitudes = [bath_metadata['xllcorner'] + x*bath_metadata['cellsize'] \
                  for x in range(number_of_longitudes)]
    latitudes = [bath_metadata['yllcorner'] + y*bath_metadata['cellsize'] \
                 for y in range(number_of_latitudes)]

     # reverse order of lat, so the first lat represents the first grid row
    latitudes.reverse()

    kmin, kmax, lmin, lmax = _get_min_max_indexes(latitudes[:],longitudes[:],
                                                  minlat=minlat, maxlat=maxlat,
                                                  minlon=minlon, maxlon=maxlon)

    bath_grid = bath_grid[kmin:kmax,lmin:lmax]
    latitudes = latitudes[kmin:kmax]
    longitudes = longitudes[lmin:lmax]
    number_of_latitudes = len(latitudes)
    number_of_longitudes = len(longitudes)
    number_of_times = len(os.listdir(elevation_dir))
    number_of_points = number_of_latitudes * number_of_longitudes
    number_of_volumes = (number_of_latitudes-1) * (number_of_longitudes-1) * 2

    #Work out the times
    if len(elevation_files) > 1:
        # Assume: The time period is less than 24hrs.
        time_period = (int(elevation_files[1][-3:]) \
                       - int(elevation_files[0][-3:])) * 60*60
        times = [x*time_period for x in range(len(elevation_files))]
    else:
        times = [0.0]

    if verbose:
        log.critical('------------------------------------------------')
        log.critical('Statistics:')
        log.critical('  Extent (lat/lon):')
        log.critical('    lat in [%f, %f], len(lat) == %d'
                     % (min(latitudes), max(latitudes), len(latitudes)))
        log.critical('    lon in [%f, %f], len(lon) == %d'
                     % (min(longitudes), max(longitudes), len(longitudes)))
        log.critical('    t in [%f, %f], len(t) == %d'
                     % (min(times), max(times), len(times)))

    ######### WRITE THE SWW FILE #############

    # NetCDF file definition
    outfile = NetCDFFile(sww_file, netcdf_mode_w)

    #Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = 'Converted from XXX'

    #For sww compatibility
    outfile.smoothing = 'Yes'
    outfile.order = 1

    #Start time in seconds since the epoch (midnight 1/1/1970)
    outfile.starttime = starttime = times[0]

    # dimension definitions
    outfile.createDimension('number_of_volumes', number_of_volumes)
    outfile.createDimension('number_of_vertices', 3)
    outfile.createDimension('number_of_points', number_of_points)
    outfile.createDimension('number_of_timesteps', number_of_times)

    # variable definitions
    outfile.createVariable('x', precision, ('number_of_points',))
    outfile.createVariable('y', precision, ('number_of_points',))
    outfile.createVariable('elevation', precision, ('number_of_points',))

    #FIXME: Backwards compatibility
    #outfile.createVariable('z', precision, ('number_of_points',))
    #################################

    outfile.createVariable('volumes', netcdf_int, ('number_of_volumes',
                                                   'number_of_vertices'))

    outfile.createVariable('time', precision, ('number_of_timesteps',))

    outfile.createVariable('stage', precision, ('number_of_timesteps',
                                                'number_of_points'))

    outfile.createVariable('xmomentum', precision, ('number_of_timesteps',
                                                    'number_of_points'))

    outfile.createVariable('ymomentum', precision, ('number_of_timesteps',
                                                    'number_of_points'))

    #Store
    from anuga.coordinate_transforms.redfearn import redfearn

    x = num.zeros(number_of_points, num.float)  #Easting
    y = num.zeros(number_of_points, num.float)  #Northing

    if verbose: log.critical('Making triangular grid')

    #Get zone of 1st point.
    refzone, _, _ = redfearn(latitudes[0], longitudes[0])

    vertices = {}
    i = 0
    for k, lat in enumerate(latitudes):
        for l, lon in enumerate(longitudes):
            vertices[l,k] = i

            zone, easting, northing = redfearn(lat, lon)

            #msg = 'Zone boundary crossed at longitude =', lon
            #assert zone == refzone, msg
            #print '%7.2f %7.2f %8.2f %8.2f' %(lon, lat, easting, northing)
            x[i] = easting
            y[i] = northing
            i += 1

    #Construct 2 triangles per 'rectangular' element
    volumes = []
    for l in range(number_of_longitudes-1):    #X direction
        for k in range(number_of_latitudes-1): #Y direction
            v1 = vertices[l,k+1]
            v2 = vertices[l,k]
            v3 = vertices[l+1,k+1]
            v4 = vertices[l+1,k]

            #Note, this is different to the ferrit2sww code
            #since the order of the lats is reversed.
            volumes.append([v1,v3,v2]) #Upper element
            volumes.append([v4,v2,v3]) #Lower element

    volumes = num.array(volumes, num.int)      #array default#

    geo_ref = Geo_reference(refzone, min(x), min(y))
    geo_ref.write_NetCDF(outfile)

    # This will put the geo ref in the middle
    #geo_ref = Geo_reference(refzone, (max(x)+min(x))/2., (max(x)+min(y))/2.)

    if verbose:
        log.critical('------------------------------------------------')
        log.critical('More Statistics:')
        log.critical('  Extent (/lon):')
        log.critical('    x in [%f, %f], len(lat) == %d'
                     % (min(x), max(x), len(x)))
        log.critical('    y in [%f, %f], len(lon) == %d'
                     % (min(y), max(y), len(y)))
        log.critical('geo_ref: ', geo_ref)

    z = num.resize(bath_grid,outfile.variables['elevation'][:].shape)
    outfile.variables['x'][:] = x - geo_ref.get_xllcorner()
    outfile.variables['y'][:] = y - geo_ref.get_yllcorner()
    # FIXME (Ole): Remove once viewer has been recompiled and changed
    #              to use elevation instead of z
    #outfile.variables['z'][:] = z
    outfile.variables['elevation'][:] = z
    outfile.variables['volumes'][:] = volumes.astype(num.int32) # On Opteron 64

    stage = outfile.variables['stage']
    xmomentum = outfile.variables['xmomentum']
    ymomentum = outfile.variables['ymomentum']

    outfile.variables['time'][:] = times   #Store time relative

    if verbose: log.critical('Converting quantities')

    n = number_of_times
    for j in range(number_of_times):
        # load in files
        elevation_meta, elevation_grid = \
            _read_asc(elevation_dir + os.sep + elevation_files[j])

        _, u_momentum_grid = _read_asc(ucur_dir + os.sep + ucur_files[j])
        _, v_momentum_grid = _read_asc(vcur_dir + os.sep + vcur_files[j])

        #cut matrix to desired size
        elevation_grid = elevation_grid[kmin:kmax,lmin:lmax]
        u_momentum_grid = u_momentum_grid[kmin:kmax,lmin:lmax]
        v_momentum_grid = v_momentum_grid[kmin:kmax,lmin:lmax]

        # handle missing values
        missing = (elevation_grid == elevation_meta['NODATA_value'])
        if num.sometrue (missing):
            if fail_on_NaN:
                msg = 'File %s contains missing values' \
                      % (elevation_files[j])
                raise DataMissingValuesError, msg
            else:
                elevation_grid = elevation_grid*(missing==0) \
                                 + missing*elevation_NaN_filler

        if verbose and j % ((n+10)/10) == 0: log.critical('  Doing %d of %d'
                                                          % (j, n))

        i = 0
        for k in range(number_of_latitudes):      #Y direction
            for l in range(number_of_longitudes): #X direction
                w = zscale*elevation_grid[k,l] + mean_stage
                stage[j,i] = w
                h = w - z[i]
                xmomentum[j,i] = u_momentum_grid[k,l]*h
                ymomentum[j,i] = v_momentum_grid[k,l]*h
                i += 1

    outfile.close()


##
# @brief Return max&min indexes (for slicing) of area specified.
# @param latitudes_ref ??
# @param longitudes_ref ??
# @param minlat Minimum latitude of specified area.
# @param maxlat Maximum latitude of specified area.
# @param minlon Minimum longitude of specified area.
# @param maxlon Maximum longitude of specified area.
# @return Tuple (lat_min_index, lat_max_index, lon_min_index, lon_max_index)
def _get_min_max_indexes(latitudes_ref,longitudes_ref,
                         minlat=None, maxlat=None,
                         minlon=None, maxlon=None):
    """
    Return max, min indexes (for slicing) of the lat and long arrays to cover
    the area specified with min/max lat/long.

    Think of the latitudes and longitudes describing a 2d surface.
    The area returned is, if possible, just big enough to cover the
    inputed max/min area. (This will not be possible if the max/min area
    has a section outside of the latitudes/longitudes area.)

    asset  longitudes are sorted,
    long - from low to high (west to east, eg 148 - 151)
    assert latitudes are sorted, ascending or decending
    """

    latitudes = latitudes_ref[:]
    longitudes = longitudes_ref[:]

    latitudes = ensure_numeric(latitudes)
    longitudes = ensure_numeric(longitudes)

    assert num.allclose(num.sort(longitudes), longitudes)

    #print latitudes[0],longitudes[0]
    #print len(latitudes),len(longitudes)
    #print latitudes[len(latitudes)-1],longitudes[len(longitudes)-1]

    lat_ascending = True
    if not num.allclose(num.sort(latitudes), latitudes):
        lat_ascending = False
        # reverse order of lat, so it's in ascending order
        latitudes = latitudes[::-1]
        assert num.allclose(num.sort(latitudes), latitudes)

    largest_lat_index = len(latitudes)-1

    #Cut out a smaller extent.
    if minlat == None:
        lat_min_index = 0
    else:
        lat_min_index = num.searchsorted(latitudes, minlat)-1
        if lat_min_index <0:
            lat_min_index = 0

    if maxlat == None:
        lat_max_index = largest_lat_index #len(latitudes)
    else:
        lat_max_index = num.searchsorted(latitudes, maxlat)
        if lat_max_index > largest_lat_index:
            lat_max_index = largest_lat_index

    if minlon == None:
        lon_min_index = 0
    else:
        lon_min_index = num.searchsorted(longitudes, minlon)-1
        if lon_min_index <0:
            lon_min_index = 0

    if maxlon == None:
        lon_max_index = len(longitudes)
    else:
        lon_max_index = num.searchsorted(longitudes, maxlon)

    # Reversing the indexes, if the lat array is decending
    if lat_ascending is False:
        lat_min_index, lat_max_index = largest_lat_index - lat_max_index, \
                                       largest_lat_index - lat_min_index
    lat_max_index = lat_max_index + 1 # taking into account how slicing works
    lon_max_index = lon_max_index + 1 # taking into account how slicing works

    return lat_min_index, lat_max_index, lon_min_index, lon_max_index


##
# @brief Read an ASC file.
# @parem filename Path to the file to read.
# @param verbose True if this function is to be verbose.
def _read_asc(filename, verbose=False):
    """Read esri file from the following ASCII format (.asc)

    Example:

    ncols         3121
    nrows         1800
    xllcorner     722000
    yllcorner     5893000
    cellsize      25
    NODATA_value  -9999
    138.3698 137.4194 136.5062 135.5558 ..........
    """

    datafile = open(filename)

    if verbose: log.critical('Reading DEM from %s' % filename)

    lines = datafile.readlines()
    datafile.close()

    if verbose: log.critical('Got %d lines' % len(lines))

    ncols = int(lines.pop(0).split()[1].strip())
    nrows = int(lines.pop(0).split()[1].strip())
    xllcorner = float(lines.pop(0).split()[1].strip())
    yllcorner = float(lines.pop(0).split()[1].strip())
    cellsize = float(lines.pop(0).split()[1].strip())
    NODATA_value = float(lines.pop(0).split()[1].strip())

    assert len(lines) == nrows

    #Store data
    grid = []

    n = len(lines)
    for i, line in enumerate(lines):
        cells = line.split()
        assert len(cells) == ncols
        grid.append(num.array([float(x) for x in cells]))
    grid = num.array(grid)

    return {'xllcorner':xllcorner,
            'yllcorner':yllcorner,
            'cellsize':cellsize,
            'NODATA_value':NODATA_value}, grid


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
# @brief Convert URS file to SWW file.
# @param basename_in Stem of the input filename.
# @param basename_out Stem of the output filename.
# @param verbose True if this function is to be verbose.
# @param remove_nc_files 
# @param minlat Sets extent of area to be used.  If not supplied, full extent.
# @param maxlat Sets extent of area to be used.  If not supplied, full extent.
# @param minlon Sets extent of area to be used.  If not supplied, full extent.
# @param maxlon Sets extent of area to be used.  If not supplied, full extent.
# @param mint 
# @param maxt 
# @param mean_stage 
# @param origin A 3-tuple with geo referenced UTM coordinates
# @param zscale 
# @param fail_on_NaN 
# @param NaN_filler 
# @param elevation 
# @note Also convert latitude and longitude to UTM. All coordinates are
#       assumed to be given in the GDA94 datum.
def urs2sww(basename_in='o', basename_out=None, verbose=False,
            remove_nc_files=True,
            minlat=None, maxlat=None,
            minlon=None, maxlon=None,
            mint=None, maxt=None,
            mean_stage=0,
            origin=None,
            zscale=1,
            fail_on_NaN=True,
            NaN_filler=0,
            elevation=None):
    """Convert a URS file to an SWW file.
    Convert URS C binary format for wave propagation to
    sww format native to abstract_2d_finite_volumes.

    Specify only basename_in and read files of the form
    basefilename-z-mux2, basefilename-e-mux2 and
    basefilename-n-mux2 containing relative height,
    x-velocity and y-velocity, respectively.

    Also convert latitude and longitude to UTM. All coordinates are
    assumed to be given in the GDA94 datum. The latitude and longitude
    information is for  a grid.

    min's and max's: If omitted - full extend is used.
    To include a value min may equal it, while max must exceed it.
    Lat and lon are assumed to be in decimal degrees.
    NOTE: minlon is the most east boundary.

    origin is a 3-tuple with geo referenced
    UTM coordinates (zone, easting, northing)
    It will be the origin of the sww file. This shouldn't be used,
    since all of anuga should be able to handle an arbitary origin.

    URS C binary format has data orgainised as TIME, LONGITUDE, LATITUDE
    which means that latitude is the fastest
    varying dimension (row major order, so to speak)

    In URS C binary the latitudes and longitudes are in assending order.
    """

    if basename_out == None:
        basename_out = basename_in

    files_out = urs2nc(basename_in, basename_out)

    ferret2sww(basename_out,
               minlat=minlat,
               maxlat=maxlat,
               minlon=minlon,
               maxlon=maxlon,
               mint=mint,
               maxt=maxt,
               mean_stage=mean_stage,
               origin=origin,
               zscale=zscale,
               fail_on_NaN=fail_on_NaN,
               NaN_filler=NaN_filler,
               inverted_bathymetry=True,
               verbose=verbose)
    
    if remove_nc_files:
        for file_out in files_out:
            os.remove(file_out)


##
# @brief Convert 3 URS files back to 4 NC files.
# @param basename_in Stem of the input filenames.
# @param basename_outStem of the output filenames.
# @note The name of the urs file names must be:
#          [basename_in]-z-mux
#          [basename_in]-e-mux
#          [basename_in]-n-mux
def urs2nc(basename_in='o', basename_out='urs'):
    """Convert the 3 urs files to 4 nc files.

    The name of the urs file names must be;
    [basename_in]-z-mux
    [basename_in]-e-mux
    [basename_in]-n-mux
    """

    files_in = [basename_in + WAVEHEIGHT_MUX_LABEL,
                basename_in + EAST_VELOCITY_LABEL,
                basename_in + NORTH_VELOCITY_LABEL]
    files_out = [basename_out + '_ha.nc',
                 basename_out + '_ua.nc',
                 basename_out + '_va.nc']
    quantities = ['HA', 'UA', 'VA']

    #if os.access(files_in[0]+'.mux', os.F_OK) == 0 :
    for i, file_name in enumerate(files_in):
        if os.access(file_name, os.F_OK) == 0:
            if os.access(file_name + '.mux', os.F_OK) == 0 :
                msg = 'File %s does not exist or is not accessible' % file_name
                raise IOError, msg
            else:
               files_in[i] += '.mux'
               log.critical("file_name %s" % file_name)

    hashed_elevation = None
    for file_in, file_out, quantity in map(None, files_in,
                                           files_out,
                                           quantities):
        lonlatdep, lon, lat, depth = _binary_c2nc(file_in,
                                                  file_out,
                                                  quantity)
        if hashed_elevation == None:
            elevation_file = basename_out + '_e.nc'
            write_elevation_nc(elevation_file,
                               lon,
                               lat,
                               depth)
            hashed_elevation = myhash(lonlatdep)
        else:
            msg = "The elevation information in the mux files is inconsistent"
            assert hashed_elevation == myhash(lonlatdep), msg

    files_out.append(elevation_file)

    return files_out


##
# @brief Convert a quantity URS file to a NetCDF file.
# @param file_in Path to input URS file.
# @param file_out Path to the output file.
# @param quantity Name of the quantity to be written to the output file.
# @return A tuple (lonlatdep, lon, lat, depth).
def _binary_c2nc(file_in, file_out, quantity):
    """Reads in a quantity urs file and writes a quantity nc file.
    Additionally, returns the depth and lat, long info,
    so it can be written to a file.
    """

    columns = 3                           # long, lat , depth
    mux_file = open(file_in, 'rb')

    # Number of points/stations
    (points_num,) = unpack('i', mux_file.read(4))

    # nt, int - Number of time steps
    (time_step_count,) = unpack('i', mux_file.read(4))

    #dt, float - time step, seconds
    (time_step,) = unpack('f', mux_file.read(4))

    msg = "Bad data in the mux file."
    if points_num < 0:
        mux_file.close()
        raise ANUGAError, msg
    if time_step_count < 0:
        mux_file.close()
        raise ANUGAError, msg
    if time_step < 0:
        mux_file.close()
        raise ANUGAError, msg

    lonlatdep = p_array.array('f')
    lonlatdep.read(mux_file, columns * points_num)
    lonlatdep = num.array(lonlatdep, dtype=num.float)
    lonlatdep = num.reshape(lonlatdep, (points_num, columns))

    lon, lat, depth = lon_lat2grid(lonlatdep)
    lon_sorted = list(lon)
    lon_sorted.sort()

    if not num.alltrue(lon == lon_sorted):
        msg = "Longitudes in mux file are not in ascending order"
        raise IOError, msg

    lat_sorted = list(lat)
    lat_sorted.sort()

    nc_file = Write_nc(quantity,
                       file_out,
                       time_step_count,
                       time_step,
                       lon,
                       lat)

    for i in range(time_step_count):
        #Read in a time slice from mux file
        hz_p_array = p_array.array('f')
        hz_p_array.read(mux_file, points_num)
        hz_p = num.array(hz_p_array, dtype=num.float)
        hz_p = num.reshape(hz_p, (len(lon), len(lat)))
        hz_p = num.transpose(hz_p)  # mux has lat varying fastest, nc has long v.f.

        #write time slice to nc file
        nc_file.store_timestep(hz_p)

    mux_file.close()
    nc_file.close()

    return lonlatdep, lon, lat, depth


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

################################################################################
# URS UNGRIDDED 2 SWW
################################################################################

### PRODUCING THE POINTS NEEDED FILE ###

# Ones used for FESA 2007 results
#LL_LAT = -50.0
#LL_LONG = 80.0
#GRID_SPACING = 1.0/60.0
#LAT_AMOUNT = 4800
#LONG_AMOUNT = 3600


##
# @brief 
# @param file_name 
# @param boundary_polygon 
# @param zone 
# @param ll_lat 
# @param ll_long 
# @param grid_spacing 
# @param lat_amount 
# @param long_amount 
# @param isSouthernHemisphere 
# @param export_csv 
# @param use_cache 
# @param verbose True if this function is to be verbose.
# @return 
def URS_points_needed_to_file(file_name, boundary_polygon, zone,
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

    geo = URS_points_needed(boundary_polygon, zone, ll_lat, ll_long,
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


##
# @brief 
# @param boundary_polygon
# @param zone
# @param ll_lat
# @param ll_long
# @param grid_spacing
# @param lat_amount
# @param long_amount
# @param isSouthHemisphere
# @param use_cache
# @param verbose
def URS_points_needed(boundary_polygon, zone, ll_lat,
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
            raise msg

        geo = cache(_URS_points_needed,
                    args, kwargs,
                    verbose=verbose,
                    compression=False)
    else:
        geo = apply(_URS_points_needed, args, kwargs)

    return geo


##
# @brief 
# @param boundary_polygon An iterable of points that describe a polygon.
# @param zone
# @param ll_lat Lower left latitude, in decimal degrees
# @param ll_long Lower left longitude, in decimal degrees
# @param grid_spacing Grid spacing in decimal degrees.
# @param lat_amount
# @param long_amount
# @param isSouthHemisphere
def _URS_points_needed(boundary_polygon,
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
        raise ValueError, msg

    # Warning there is no info in geospatial saying the hemisphere of
    # these points.  There should be.
    geo = Geospatial_data(data_points=list(lat_long_set),
                          points_are_lats_longs=True)

    return geo


##
# @brief Get the points that are needed to interpolate any point a a segment.
# @param seg Two points in the UTM.
# @param ll_lat Lower left latitude, in decimal degrees
# @param ll_long Lower left longitude, in decimal degrees
# @param grid_spacing 
# @param lat_amount 
# @param long_amount 
# @param zone 
# @param isSouthHemisphere 
# @return A list of points.
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

    first_row = (min_long - ll_long) / grid_spacing

    # To round up
    first_row_long = int(round(first_row + 0.5))

    last_row = (max_long - ll_long) / grid_spacing # round down
    last_row_long = int(round(last_row))

    first_row = (min_lat - ll_lat) / grid_spacing
    # To round up
    first_row_lat = int(round(first_row + 0.5))

    last_row = (max_lat - ll_lat) / grid_spacing # round down
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


##
# @brief 
# @param lat
# @param long
# @param seg Two points in UTM.
# @param max_distance
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
    try:
        d = abs((x2_1)*(y1-y0)-(x1-x0)*(y2_1))/sqrt( \
            (x2_1)*(x2_1)+(y2_1)*(y2_1))
    except ZeroDivisionError:
        if sqrt((x2_1)*(x2_1)+(y2_1)*(y2_1)) == 0 \
           and abs((x2_1)*(y1-y0)-(x1-x0)*(y2_1)) == 0:
            return True
        else:
            return False

    return d <= max_distance

################################################################################
# CONVERTING UNGRIDDED URS DATA TO AN SWW FILE
################################################################################

WAVEHEIGHT_MUX_LABEL = '-z-mux'
EAST_VELOCITY_LABEL =  '-e-mux'
NORTH_VELOCITY_LABEL =  '-n-mux'

##
# @brief Convert URS file(s) (wave prop) to an SWW file.
# @param basename_in Stem of the input filenames.
# @param basename_out Path to the output SWW file.
# @param verbose True if this function is to be verbose.
# @param mint
# @param maxt
# @param mean_stage
# @param origin Tuple with geo-ref UTM coordinates (zone, easting, northing).
# @param hole_points_UTM
# @param zscale
# @note Also convert latitude and longitude to UTM. All coordinates are
#       assumed to be given in the GDA94 datum.
# @note Input filename stem has suffixes '-z-mux', '-e-mux' and '-n-mux'
#       added for relative height, x-velocity and y-velocity respectively.
def urs_ungridded2sww(basename_in='o', basename_out=None, verbose=False,
                      mint=None, maxt=None,
                      mean_stage=0,
                      origin=None,
                      hole_points_UTM=None,
                      zscale=1):
    """
    Convert URS C binary format for wave propagation to
    sww format native to abstract_2d_finite_volumes.

    Specify only basename_in and read files of the form
    basefilename-z-mux, basefilename-e-mux and
    basefilename-n-mux containing relative height,
    x-velocity and y-velocity, respectively.

    Also convert latitude and longitude to UTM. All coordinates are
    assumed to be given in the GDA94 datum. The latitude and longitude
    information is assumed ungridded grid.

    min's and max's: If omitted - full extend is used.
    To include a value min ans max may equal it.
    Lat and lon are assumed to be in decimal degrees.

    origin is a 3-tuple with geo referenced
    UTM coordinates (zone, easting, northing)
    It will be the origin of the sww file. This shouldn't be used,
    since all of anuga should be able to handle an arbitary origin.
    The mux point info is NOT relative to this origin.

    URS C binary format has data organised as TIME, LONGITUDE, LATITUDE
    which means that latitude is the fastest
    varying dimension (row major order, so to speak)

    In URS C binary the latitudes and longitudes are in assending order.

    Note, interpolations of the resulting sww file will be different
    from results of urs2sww.  This is due to the interpolation
    function used, and the different grid structure between urs2sww
    and this function.

    Interpolating data that has an underlying gridded source can
    easily end up with different values, depending on the underlying
    mesh.

    consider these 4 points
    50  -50

    0     0

    The grid can be
     -
    |\|   A
     -
     or;
      -
     |/|  B
      -

    If a point is just below the center of the midpoint, it will have a
    +ve value in grid A and a -ve value in grid B.
    """

    from anuga.mesh_engine.mesh_engine import NoTrianglesError
    from anuga.pmesh.mesh import Mesh

    files_in = [basename_in + WAVEHEIGHT_MUX_LABEL,
                basename_in + EAST_VELOCITY_LABEL,
                basename_in + NORTH_VELOCITY_LABEL]
    quantities = ['HA','UA','VA']

    # instantiate urs_points of the three mux files.
    mux = {}
    for quantity, file in map(None, quantities, files_in):
        mux[quantity] = Urs_points(file)

    # Could check that the depth is the same. (hashing)

    # handle to a mux file to do depth stuff
    a_mux = mux[quantities[0]]

    # Convert to utm
    lat = a_mux.lonlatdep[:,1]
    long = a_mux.lonlatdep[:,0]
    points_utm, zone = convert_from_latlon_to_utm(latitudes=lat,
                                                  longitudes=long)

    elevation = a_mux.lonlatdep[:,2] * -1

    # grid (create a mesh from the selected points)
    # This mesh has a problem.  Triangles are streched over ungridded areas.
    # If these areas could be described as holes in pmesh, that would be great.

    # I can't just get the user to selection a point in the middle.
    # A boundary is needed around these points.
    # But if the zone of points is obvious enough auto-segment should do
    # a good boundary.
    mesh = Mesh()
    mesh.add_vertices(points_utm)
    mesh.auto_segment(smooth_indents=True, expand_pinch=True)

    # To try and avoid alpha shape 'hugging' too much
    mesh.auto_segment(mesh.shape.get_alpha() * 1.1)
    if hole_points_UTM is not None:
        point = ensure_absolute(hole_points_UTM)
        mesh.add_hole(point[0], point[1])

    try:
        mesh.generate_mesh(minimum_triangle_angle=0.0, verbose=False)
    except NoTrianglesError:
        # This is a bit of a hack, going in and changing the data structure.
        mesh.holes = []
        mesh.generate_mesh(minimum_triangle_angle=0.0, verbose=False)

    mesh_dic = mesh.Mesh2MeshList()

    #mesh.export_mesh_file(basename_in + '_168.tsh')
    #import sys; sys.exit()
    # These are the times of the mux file
    mux_times = []
    for i in range(a_mux.time_step_count):
        mux_times.append(a_mux.time_step * i)
    (mux_times_start_i, mux_times_fin_i) = mux2sww_time(mux_times, mint, maxt)
    times = mux_times[mux_times_start_i:mux_times_fin_i]

    if mux_times_start_i == mux_times_fin_i:
        # Close the mux files
        for quantity, file in map(None, quantities, files_in):
            mux[quantity].close()
        msg = "Due to mint and maxt there's no time info in the boundary SWW."
        raise Exception, msg

    # If this raise is removed there is currently no downstream errors

    points_utm=ensure_numeric(points_utm)
    assert num.alltrue(ensure_numeric(mesh_dic['generatedpointlist'])
                       == ensure_numeric(points_utm))

    volumes = mesh_dic['generatedtrianglelist']

    # Write sww intro and grid stuff.
    if basename_out is None:
        swwname = basename_in + '.sww'
    else:
        swwname = basename_out + '.sww'

    if verbose: log.critical('Output to %s' % swwname)

    outfile = NetCDFFile(swwname, netcdf_mode_w)

    # For a different way of doing this, check out tsh2sww
    # work out sww_times and the index range this covers
    sww = Write_sww(['elevation'], ['stage', 'xmomentum', 'ymomentum'])
    sww.store_header(outfile, times, len(volumes), len(points_utm),
                     verbose=verbose, sww_precision=netcdf_float)
    outfile.mean_stage = mean_stage
    outfile.zscale = zscale

    sww.store_triangulation(outfile, points_utm, volumes,
                            zone,  
                            new_origin=origin,
                            verbose=verbose)
    sww.store_static_quantities(outfile, elevation=elevation)

    if verbose: log.critical('Converting quantities')

    # Read in a time slice from each mux file and write it to the SWW file
    j = 0
    for ha, ua, va in map(None, mux['HA'], mux['UA'], mux['VA']):
        if j >= mux_times_start_i and j < mux_times_fin_i:
            stage = zscale*ha + mean_stage
            h = stage - elevation
            xmomentum = ua*h
            ymomentum = -1 * va * h # -1 since in mux files south is positive.
            sww.store_quantities(outfile,
                                 slice_index=j-mux_times_start_i,
                                 verbose=verbose,
                                 stage=stage,
                                 xmomentum=xmomentum,
                                 ymomentum=ymomentum,
                                 sww_precision=num.float)
        j += 1

    if verbose: sww.verbose_quantities(outfile)

    outfile.close()

    # Do some conversions while writing the sww file


################################################################################
# READ MUX2 FILES line of points
################################################################################

WAVEHEIGHT_MUX2_LABEL = '-z-mux2'
EAST_VELOCITY_MUX2_LABEL = '-e-mux2'
NORTH_VELOCITY_MUX2_LABEL = '-n-mux2'

##
# @brief 
# @param filenames List of mux2 format input filenames.
# @param weights Weights associated with each source file.
# @param permutation The gauge numbers for which data is extracted.
# @param verbose True if this function is to be verbose.
# @return (times, latitudes, longitudes, elevation, quantity, starttime)
def read_mux2_py(filenames,
                 weights=None,
                 permutation=None,
                 verbose=False):
    """Access the mux files specified in the filenames list. Combine the
       data found therin as a weighted linear sum as specifed by the weights.
       If permutation is None or empty extract timeseries data for all gauges
       within the files.

       Input:
           filenames:   List of filenames specifiying the file containing the
                        timeseries data (mux2 format) for each source
           weights:     Weighs associated with each source
                        (defaults to 1 for each source)
           permutation: Specifies the gauge numbers that for which data is to be
                        extracted
    """

    from urs_ext import read_mux2

    numSrc = len(filenames)

    file_params = -1 * num.ones(3, num.float)                    # [nsta,dt,nt]

    # Convert verbose to int C flag
    if verbose:
        verbose=1
    else:
        verbose=0

    if weights is None:
        weights = num.ones(numSrc)

    if permutation is None:
        permutation = ensure_numeric([], num.float)

    # Call underlying C implementation urs2sts_ext.c
    data = read_mux2(numSrc, filenames, weights, file_params,
                     permutation, verbose)

    msg = 'File parameter values were not read in correctly from c file'
    assert len(num.compress(file_params > 0, file_params)) != 0, msg

    msg = 'The number of stations specifed in the c array and in the file ' \
          'are inconsistent'
    assert file_params[0] >= len(permutation), msg

    msg = 'The number of stations returned is inconsistent with ' \
          'the requested number'
    assert len(permutation) == 0 or len(permutation) == data.shape[0], msg

    nsta = int(file_params[0])
    msg = 'Must have at least one station'
    assert nsta > 0, msg

    dt = file_params[1]
    msg = 'Must have a postive timestep'
    assert dt > 0, msg

    nt = int(file_params[2])
    msg = 'Must have at least one gauge value'
    assert nt > 0, msg

    OFFSET = 5 # Number of site parameters p passed back with data
               # p = [geolat,geolon,depth,start_tstep,finish_tstep]

    # FIXME (Ole): What is the relationship with params and data.shape ?
    # It looks as if the following asserts should pass but they don't always
    #
    #msg = 'nt = %d, data.shape[1] == %d' %(nt, data.shape[1])
    #assert nt == data.shape[1] - OFFSET, msg
    #
    #msg = 'nsta = %d, data.shape[0] == %d' %(nsta, data.shape[0])
    #assert nsta == data.shape[0], msg

    # Number of stations in ordering file
    number_of_selected_stations = data.shape[0]

    # Index where data ends and parameters begin
    parameters_index = data.shape[1] - OFFSET

    times = dt * num.arange(parameters_index)
    latitudes = num.zeros(number_of_selected_stations, num.float)
    longitudes = num.zeros(number_of_selected_stations, num.float)
    elevation = num.zeros(number_of_selected_stations, num.float)
    quantity = num.zeros((number_of_selected_stations, parameters_index), num.float)

    starttime = 1e16
    for i in range(number_of_selected_stations):
        quantity[i][:] = data[i][:parameters_index]
        latitudes[i] = data[i][parameters_index]
        longitudes[i] = data[i][parameters_index+1]
        elevation[i] = -data[i][parameters_index+2]
        first_time_step = data[i][parameters_index+3]
        starttime = min(dt*first_time_step, starttime)

    return times, latitudes, longitudes, elevation, quantity, starttime


##
# @brief ??
# @param mux_times ??
# @param mint ??
# @param maxt ??
# @return ??
def mux2sww_time(mux_times, mint, maxt):
    """
    """

    if mint == None:
        mux_times_start_i = 0
    else:
        mux_times_start_i = num.searchsorted(mux_times, mint)

    if maxt == None:
        mux_times_fin_i = len(mux_times)
    else:
        maxt += 0.5 # so if you specify a time where there is
                    # data that time will be included
        mux_times_fin_i = num.searchsorted(mux_times, maxt)

    return mux_times_start_i, mux_times_fin_i


##
# @brief Convert a URS (mux2, wave propagation) file to an STS file.
# @param basename_in String (or list) of source file stems.
# @param basename_out Stem of output STS file path.
# @param weights
# @param verbose True if this function is to be verbose.
# @param origin Tuple with with geo-ref UTM coords (zone, easting, northing).
# @param zone
# @param mean_stage
# @param zscale
# @param ordering_filename Path of a file specifying which mux2 gauge points are
#                          to be stored.
# @note Also convert latitude and longitude to UTM. All coordinates are
#       assumed to be given in the GDA94 datum.
def urs2sts(basename_in, basename_out=None,
            weights=None,
            verbose=False,
            origin=None,
            zone=None,
            central_meridian=None,            
            mean_stage=0.0,
            zscale=1.0,
            ordering_filename=None):
    """Convert URS mux2 format for wave propagation to sts format

    Also convert latitude and longitude to UTM. All coordinates are
    assumed to be given in the GDA94 datum

    origin is a 3-tuple with geo referenced
    UTM coordinates (zone, easting, northing)

    inputs:

    basename_in: list of source file prefixes

        These are combined with the extensions:
        WAVEHEIGHT_MUX2_LABEL = '-z-mux2' for stage
        EAST_VELOCITY_MUX2_LABEL = '-e-mux2' xmomentum
        NORTH_VELOCITY_MUX2_LABEL = '-n-mux2' and ymomentum

        to create a 2D list of mux2 file. The rows are associated with each
        quantity and must have the above extensions
        the columns are the list of file prefixes.

    ordering: a .txt file name specifying which mux2 gauge points are
              to be stored. This is indicated by the index of the gauge
              in the ordering file.

              ordering file format:
              1st line:    'index,longitude,latitude\n'
              other lines: index,longitude,latitude

              If ordering is None or ordering file is empty then
               all points are taken in the order they
              appear in the mux2 file.


    output:
      basename_out: name of sts file in which mux2 data is stored.
      
      
      
    NOTE: South is positive in mux files so sign of y-component of velocity is reverted
    """

    import os
    from Scientific.IO.NetCDF import NetCDFFile
    from types import ListType,StringType
    from operator import __and__

    if not isinstance(basename_in, ListType):
        if verbose: log.critical('Reading single source')
        basename_in = [basename_in]

    # This is the value used in the mux file format to indicate NAN data
    # FIXME (Ole): This should be changed everywhere to IEEE NAN when
    #              we upgrade to Numpy
    NODATA = 99

    # Check that basename is a list of strings
    if not reduce(__and__, map(lambda z:isinstance(z,StringType), basename_in)):
        msg= 'basename_in must be a string or list of strings'
        raise Exception, msg

    # Find the number of sources to be used
    numSrc = len(basename_in)

    # A weight must be specified for each source
    if weights is None:
        # Default is equal weighting
        weights = num.ones(numSrc, num.float) / numSrc
    else:
        weights = ensure_numeric(weights)
        msg = 'When combining multiple sources must specify a weight for ' \
              'mux2 source file'
        assert len(weights) == numSrc, msg

    if verbose: log.critical('Weights used in urs2sts: %s' % str(weights))

    # Check output filename
    if basename_out is None:
        msg = 'STS filename must be specified as basename_out ' \
              'in function urs2sts'
        raise Exception, msg

    if basename_out.endswith('.sts'):
        stsname = basename_out
    else:
        stsname = basename_out + '.sts'

    # Create input filenames from basenames and check their existence
    files_in = [[], [], []]
    for files in basename_in:
        files_in[0].append(files + WAVEHEIGHT_MUX2_LABEL),
        files_in[1].append(files + EAST_VELOCITY_MUX2_LABEL)
        files_in[2].append(files + NORTH_VELOCITY_MUX2_LABEL)

    quantities = ['HA','UA','VA'] # Quantity names used in the MUX2 format
    for i in range(len(quantities)):
        for file_in in files_in[i]:
            if (os.access(file_in, os.R_OK) == 0):
                msg = 'File %s does not exist or is not accessible' % file_in
                raise IOError, msg

    # Establish permutation array
    if ordering_filename is not None:
        if verbose is True: log.critical('Reading ordering file %s'
                                         % ordering_filename)

        # Read ordering file
        try:
            fid = open(ordering_filename, 'r')
            file_header = fid.readline().split(',')
            ordering_lines = fid.readlines()
            fid.close()
        except:
            msg = 'Cannot open %s' % ordering_filename
            raise Exception, msg

        reference_header = 'index, longitude, latitude\n'
        reference_header_split = reference_header.split(',')
        for i in range(3):
            if not file_header[i].strip() == reference_header_split[i].strip():
                msg = 'File must contain header: ' + reference_header
                raise Exception, msg

        if len(ordering_lines) < 2:
            msg = 'File must contain at least two points'
            raise Exception, msg

        permutation = [int(line.split(',')[0]) for line in ordering_lines]
        permutation = ensure_numeric(permutation)
    else:
        permutation = None

    # Read MUX2 files
    if (verbose): log.critical('reading mux2 file')

    mux={}
    for i, quantity in enumerate(quantities):
        # For each quantity read the associated list of source mux2 file with
        # extention associated with that quantity

        times, latitudes, longitudes, elevation, mux[quantity], starttime \
            = read_mux2_py(files_in[i], weights, permutation, verbose=verbose)

        # Check that all quantities have consistent time and space information
        if quantity != quantities[0]:
            msg = '%s, %s and %s have inconsistent gauge data' \
                  % (files_in[0], files_in[1], files_in[2])
            assert num.allclose(times, times_old), msg
            assert num.allclose(latitudes, latitudes_old), msg
            assert num.allclose(longitudes, longitudes_old), msg
            assert num.allclose(elevation, elevation_old), msg
            assert num.allclose(starttime, starttime_old), msg
        times_old = times
        latitudes_old = latitudes
        longitudes_old = longitudes
        elevation_old = elevation
        starttime_old = starttime

        # Self check - can be removed to improve speed
        #ref_longitudes = [float(line.split(',')[1]) for line in ordering_lines]
        #ref_latitudes = [float(line.split(',')[2]) for line in ordering_lines]
        #
        #msg = 'Longitudes specified in ordering file do not match those ' \
        #      'found in mux files. ' \
        #      'I got %s instead of %s (only beginning shown)' \
        #      % (str(longitudes[:10]) + '...',
        #         str(ref_longitudes[:10]) + '...')
        #assert allclose(longitudes, ref_longitudes), msg
        #
        #msg = 'Latitudes specified in ordering file do not match those ' \
        #      'found in mux files. '
        #      'I got %s instead of %s (only beginning shown)' \
        #      % (str(latitudes[:10]) + '...',
        #         str(ref_latitudes[:10]) + '...')
        #assert allclose(latitudes, ref_latitudes), msg

    # Store timeseries in STS file
    msg = 'File is empty and or clipped region not in file region'
    assert len(latitudes > 0), msg

    number_of_points = latitudes.shape[0]      # Number of stations retrieved
    number_of_times = times.shape[0]           # Number of timesteps
    number_of_latitudes = latitudes.shape[0]   # Number latitudes
    number_of_longitudes = longitudes.shape[0] # Number longitudes

    # The permutation vector of contains original indices
    # as given in ordering file or None in which case points
    # are assigned the trivial indices enumerating them from
    # 0 to number_of_points-1
    if permutation is None:
        permutation = num.arange(number_of_points, dtype=num.int)

    # NetCDF file definition
    outfile = NetCDFFile(stsname, netcdf_mode_w)

    description = 'Converted from URS mux2 files: %s' % basename_in

    # Create new file
    sts = Write_sts()
    sts.store_header(outfile,
                     times+starttime,
                     number_of_points,
                     description=description,
                     verbose=verbose,
                     sts_precision=netcdf_float)

    # Store
    from anuga.coordinate_transforms.redfearn import redfearn

    x = num.zeros(number_of_points, num.float)  # Easting
    y = num.zeros(number_of_points, num.float)  # Northing

    # Check zone boundaries
    if zone is None:
        refzone, _, _ = redfearn(latitudes[0], longitudes[0],
                                 central_meridian=central_meridian)
    else:
        refzone = zone

    old_zone = refzone

    for i in range(number_of_points):
        computed_zone, easting, northing = redfearn(latitudes[i], longitudes[i],
                                                    zone=zone,
                                                    central_meridian=central_meridian)
        x[i] = easting
        y[i] = northing
        if computed_zone != refzone:
            msg = 'All sts gauges need to be in the same zone. \n'
            msg += 'offending gauge:Zone %d,%.4f, %4f\n' \
                   % (computed_zone, easting, northing)
            msg += 'previous gauge:Zone %d,%.4f, %4f' \
                   % (old_zone, old_easting, old_northing)
            raise Exception, msg
        old_zone = computed_zone
        old_easting = easting
        old_northing = northing

    if origin is None:
        origin = Geo_reference(refzone, min(x), min(y))
    geo_ref = write_NetCDF_georeference(origin, outfile)

    elevation = num.resize(elevation, outfile.variables['elevation'][:].shape)
    outfile.variables['permutation'][:] = permutation.astype(num.int32) # Opteron 64
    outfile.variables['x'][:] = x - geo_ref.get_xllcorner()
    outfile.variables['y'][:] = y - geo_ref.get_yllcorner()
    outfile.variables['elevation'][:] = elevation

    stage = outfile.variables['stage']
    xmomentum = outfile.variables['xmomentum']
    ymomentum = outfile.variables['ymomentum']

    if verbose: log.critical('Converting quantities')

    for j in range(len(times)):
        for i in range(number_of_points):
            ha = mux['HA'][i,j]
            ua = mux['UA'][i,j]
            va = mux['VA'][i,j]
            if ha == NODATA:
                if verbose:
                    msg = 'Setting nodata value %d to 0 at time = %f, ' \
                          'point = %d' % (ha, times[j], i)
                    log.critical(msg)
                ha = 0.0
                ua = 0.0
                va = 0.0

            w = zscale*ha + mean_stage
            h = w - elevation[i]
            stage[j,i] = w

            xmomentum[j,i] = ua * h
            ymomentum[j,i] = -va * h # South is positive in mux files


    outfile.close()


##
# @brief Create a list of points defining a boundary from an STS file.
# @param stsname Stem of path to the STS file to read.
# @return A list of boundary points.
def create_sts_boundary(stsname):
    """Create a list of points defining a boundary from an STS file.

    Create boundary segments from .sts file. Points can be stored in
    arbitrary order within the .sts file. The order in which the .sts points
    make up the boundary are given in order.txt file

    FIXME:
    Point coordinates are stored in relative eastings and northings.
    But boundary is produced in absolute coordinates
    """

    try:
        fid = NetCDFFile(stsname + '.sts', netcdf_mode_r)
    except:
        msg = 'Cannot open %s' % stsname + '.sts'
        raise msg

    xllcorner = fid.xllcorner[0]
    yllcorner = fid.yllcorner[0]

    #Points stored in sts file are normalised to [xllcorner,yllcorner] but
    #we cannot assume that boundary polygon will be. At least the
    #additional points specified by the user after this function is called
    x = fid.variables['x'][:] + xllcorner
    y = fid.variables['y'][:] + yllcorner

    x = num.reshape(x, (len(x),1))
    y = num.reshape(y, (len(y),1))
    sts_points = num.concatenate((x,y), axis=1)

    return sts_points.tolist()


# @brief A class to write an SWW file.
class Write_sww:
    from anuga.shallow_water.shallow_water_domain import Domain

    RANGE = '_range'
    EXTREMA = ':extrema'

    ##
    # brief Instantiate the SWW writer class.
    def __init__(self, static_quantities, dynamic_quantities):
        """Initialise Write_sww with two list af quantity names: 
        
        static_quantities (e.g. elevation or friction): 
            Stored once at the beginning of the simulation in a 1D array
            of length number_of_points   
        dynamic_quantities (e.g stage):
            Stored every timestep in a 2D array with 
            dimensions number_of_points X number_of_timesteps        
        
        """
        self.static_quantities = static_quantities   
        self.dynamic_quantities = dynamic_quantities


    ##
    # @brief Store a header in the SWW file.
    # @param outfile Open handle to the file that will be written.
    # @param times A list of time slices *or* a start time.
    # @param number_of_volumes The number of triangles.
    # @param number_of_points The number of points.
    # @param description The internal file description string.
    # @param smoothing True if smoothing is to be used.
    # @param order 
    # @param sww_precision Data type of the quantity written (netcdf constant)
    # @param verbose True if this function is to be verbose.
    # @note If 'times' is a list, the info will be made relative.
    def store_header(self,
                     outfile,
                     times,
                     number_of_volumes,
                     number_of_points,
                     description='Generated by ANUGA',
                     smoothing=True,
                     order=1,
                     sww_precision=netcdf_float32,
                     verbose=False):
        """Write an SWW file header.

        outfile - the open file that will be written
        times - A list of the time slice times OR a start time
        Note, if a list is given the info will be made relative.
        number_of_volumes - the number of triangles
        """

        outfile.institution = 'Geoscience Australia'
        outfile.description = description

        # For sww compatibility
        if smoothing is True:
            # Smoothing to be depreciated
            outfile.smoothing = 'Yes'
            outfile.vertices_are_stored_uniquely = 'False'
        else:
            # Smoothing to be depreciated
            outfile.smoothing = 'No'
            outfile.vertices_are_stored_uniquely = 'True'
        outfile.order = order

        try:
            revision_number = get_revision_number()
        except:
            revision_number = None
        # Allow None to be stored as a string
        outfile.revision_number = str(revision_number)

        # This is being used to seperate one number from a list.
        # what it is actually doing is sorting lists from numeric arrays.
        if isinstance(times, (list, num.ndarray)):
            number_of_times = len(times)
            times = ensure_numeric(times)
            if number_of_times == 0:
                starttime = 0
            else:
                starttime = times[0]
                times = times - starttime  #Store relative times
        else:
            number_of_times = 0
            starttime = times


        outfile.starttime = starttime

        # dimension definitions
        outfile.createDimension('number_of_volumes', number_of_volumes)
        outfile.createDimension('number_of_vertices', 3)
        outfile.createDimension('numbers_in_range', 2)

        if smoothing is True:
            outfile.createDimension('number_of_points', number_of_points)
            # FIXME(Ole): This will cause sww files for parallel domains to
            # have ghost nodes stored (but not used by triangles).
            # To clean this up, we have to change get_vertex_values and
            # friends in quantity.py (but I can't be bothered right now)
        else:
            outfile.createDimension('number_of_points', 3*number_of_volumes)

        outfile.createDimension('number_of_timesteps', number_of_times)

        # variable definitions
        outfile.createVariable('x', sww_precision, ('number_of_points',))
        outfile.createVariable('y', sww_precision, ('number_of_points',))

        outfile.createVariable('volumes', netcdf_int, ('number_of_volumes',
                                                       'number_of_vertices'))

        # Doing sww_precision instead of Float gives cast errors.
        outfile.createVariable('time', netcdf_float,
                               ('number_of_timesteps',))

                               
        for q in self.static_quantities:
            
            outfile.createVariable(q, sww_precision,
                                   ('number_of_points',))
            
            outfile.createVariable(q + Write_sww.RANGE, sww_precision,
                                   ('numbers_in_range',))
                                   
            # Initialise ranges with small and large sentinels.
            # If this was in pure Python we could have used None sensibly
            outfile.variables[q+Write_sww.RANGE][0] = max_float  # Min
            outfile.variables[q+Write_sww.RANGE][1] = -max_float # Max

        #if 'elevation' in self.static_quantities:    
        #    # FIXME: Backwards compat - get rid of z once old view has retired
        #    outfile.createVariable('z', sww_precision,
        #                           ('number_of_points',))
                               
        for q in self.dynamic_quantities:
            outfile.createVariable(q, sww_precision, ('number_of_timesteps',
                                                      'number_of_points'))
            outfile.createVariable(q + Write_sww.RANGE, sww_precision,
                                   ('numbers_in_range',))

            # Initialise ranges with small and large sentinels.
            # If this was in pure Python we could have used None sensibly
            outfile.variables[q+Write_sww.RANGE][0] = max_float  # Min
            outfile.variables[q+Write_sww.RANGE][1] = -max_float # Max

        if isinstance(times, (list, num.ndarray)):
            outfile.variables['time'][:] = times    # Store time relative

        if verbose:
            log.critical('------------------------------------------------')
            log.critical('Statistics:')
            log.critical('    t in [%f, %f], len(t) == %d'
                         % (num.min(times), num.max(times), len(times.flat)))

    ##
    # @brief Store triangulation data in the underlying file.
    # @param outfile Open handle to underlying file.
    # @param points_utm List or array of points in UTM.
    # @param volumes 
    # @param zone 
    # @param new_origin georeference that the points can be set to.
    # @param points_georeference The georeference of the points_utm.
    # @param verbose True if this function is to be verbose.
    def store_triangulation(self,
                            outfile,
                            points_utm,
                            volumes,
                            zone=None, 
                            new_origin=None,
                            points_georeference=None, 
                            verbose=False):
        """
        new_origin - qa georeference that the points can be set to. (Maybe
        do this before calling this function.)

        points_utm - currently a list or array of the points in UTM.
        points_georeference - the georeference of the points_utm

        How about passing new_origin and current_origin.
        If you get both, do a convertion from the old to the new.

        If you only get new_origin, the points are absolute,
        convert to relative

        if you only get the current_origin the points are relative, store
        as relative.

        if you get no georefs create a new georef based on the minimums of
        points_utm.  (Another option would be to default to absolute)

        Yes, and this is done in another part of the code.
        Probably geospatial.

        If you don't supply either geo_refs, then supply a zone. If not
        the default zone will be used.

        precon:
            header has been called.
        """

        number_of_points = len(points_utm)
        volumes = num.array(volumes)
        points_utm = num.array(points_utm)

        # Given the two geo_refs and the points, do the stuff
        # described in the method header
        # if this is needed else where, pull out as a function
        points_georeference = ensure_geo_reference(points_georeference)
        new_origin = ensure_geo_reference(new_origin)
        if new_origin is None and points_georeference is not None:
            points = points_utm
            geo_ref = points_georeference
        else:
            if new_origin is None:
                new_origin = Geo_reference(zone, min(points_utm[:,0]),
                                                 min(points_utm[:,1]))
            points = new_origin.change_points_geo_ref(points_utm,
                                                      points_georeference)
            geo_ref = new_origin

        # At this stage I need a georef and points
        # the points are relative to the georef
        geo_ref.write_NetCDF(outfile)

        # This will put the geo ref in the middle
        #geo_ref = Geo_reference(refzone,(max(x)+min(x))/2.0,(max(x)+min(y))/2.)

        x =  points[:,0]
        y =  points[:,1]

        if verbose:
            log.critical('------------------------------------------------')
            log.critical('More Statistics:')
            log.critical('  Extent (/lon):')
            log.critical('    x in [%f, %f], len(lat) == %d'
                         % (min(x), max(x), len(x)))
            log.critical('    y in [%f, %f], len(lon) == %d'
                         % (min(y), max(y), len(y)))
            #log.critical('    z in [%f, %f], len(z) == %d'
            #             % (min(elevation), max(elevation), len(elevation)))
            log.critical('geo_ref: %s' % str(geo_ref))
            log.critical('------------------------------------------------')

        outfile.variables['x'][:] = points[:,0] #- geo_ref.get_xllcorner()
        outfile.variables['y'][:] = points[:,1] #- geo_ref.get_yllcorner()
        outfile.variables['volumes'][:] = volumes.astype(num.int32) #On Opteron 64



    # @brief Write the static quantity data to the underlying file.
    # @param outfile Handle to open underlying file.
    # @param sww_precision Format of quantity data to write (default Float32).
    # @param verbose True if this function is to be verbose.
    # @param **quant
    def store_static_quantities(self, 
                                outfile, 
                                sww_precision=num.float32,
                                verbose=False, 
                                **quant):
        """
        Write the static quantity info.

        **quant is extra keyword arguments passed in. These must be
          the numpy arrays to be stored in the sww file at each timestep.

        The argument sww_precision allows for storing as either 
        * single precision (default): num.float32
        * double precision: num.float64 or num.float 

        Precondition:
            store_triangulation and
            store_header have been called.
        """

        # The dictionary quant must contain numpy arrays for each name.
        # These will typically be quantities from Domain such as friction 
        #
        # Arrays not listed in static_quantitiues will be ignored, silently.
        #
        # This method will also write the ranges for each quantity, 
        # e.g. stage_range, xmomentum_range and ymomentum_range
        for q in self.static_quantities:
            if not quant.has_key(q):
                msg = 'Values for quantity %s was not specified in ' % q
                msg += 'store_quantities so they cannot be stored.'
                raise NewQuantity, msg
            else:
                q_values = ensure_numeric(quant[q])
                
                x = q_values.astype(sww_precision)
                outfile.variables[q][:] = x
        
                # This populates the _range values
                outfile.variables[q + Write_sww.RANGE][0] = num.min(x)
                outfile.variables[q + Write_sww.RANGE][1] = num.max(x)
                    
        # FIXME: Hack for backwards compatibility with old viewer
        #if 'elevation' in self.static_quantities:
        #    outfile.variables['z'][:] = outfile.variables['elevation'][:]

                    
                    
        
        
    ##
    # @brief Write the quantity data to the underlying file.
    # @param outfile Handle to open underlying file.
    # @param sww_precision Format of quantity data to write (default Float32).
    # @param slice_index
    # @param time
    # @param verbose True if this function is to be verbose.
    # @param **quant
    def store_quantities(self, 
                         outfile, 
                         sww_precision=num.float32,
                         slice_index=None,
                         time=None,
                         verbose=False, 
                         **quant):
        """
        Write the quantity info at each timestep.

        **quant is extra keyword arguments passed in. These must be
          the numpy arrays to be stored in the sww file at each timestep.

        if the time array is already been built, use the slice_index
        to specify the index.

        Otherwise, use time to increase the time dimension

        Maybe make this general, but the viewer assumes these quantities,
        so maybe we don't want it general - unless the viewer is general
        
        The argument sww_precision allows for storing as either 
        * single precision (default): num.float32
        * double precision: num.float64 or num.float 

        Precondition:
            store_triangulation and
            store_header have been called.
        """

        if time is not None:
            file_time = outfile.variables['time']
            slice_index = len(file_time)
            file_time[slice_index] = time
        else:
            slice_index = int(slice_index) # Has to be cast in case it was numpy.int    

        # Write the named dynamic quantities
        # The dictionary quant must contain numpy arrays for each name.
        # These will typically be the conserved quantities from Domain 
        # (Typically stage,  xmomentum, ymomentum).
        #
        # Arrays not listed in dynamic_quantitiues will be ignored, silently.
        #
        # This method will also write the ranges for each quantity, 
        # e.g. stage_range, xmomentum_range and ymomentum_range
        for q in self.dynamic_quantities:
            if not quant.has_key(q):
                msg = 'Values for quantity %s was not specified in ' % q
                msg += 'store_quantities so they cannot be stored.'
                raise NewQuantity, msg
            else:
                q_values = ensure_numeric(quant[q])
                
                x = q_values.astype(sww_precision)
                outfile.variables[q][slice_index] = x
                    
        
                # This updates the _range values
                q_range = outfile.variables[q + Write_sww.RANGE][:]
                q_values_min = num.min(q_values)
                if q_values_min < q_range[0]:
                    outfile.variables[q + Write_sww.RANGE][0] = q_values_min
                q_values_max = num.max(q_values)
                if q_values_max > q_range[1]:
                    outfile.variables[q + Write_sww.RANGE][1] = q_values_max

    ##
    # @brief Print the quantities in the underlying file.
    # @param outfile UNUSED.
    def verbose_quantities(self, outfile):
        log.critical('------------------------------------------------')
        log.critical('More Statistics:')
        for q in self.dynamic_quantities:
            log.critical('  %s in [%f, %f]'
                         % (q, outfile.variables[q+Write_sww.RANGE][0],
                            outfile.variables[q+Write_sww.RANGE][1]))
        log.critical('------------------------------------------------')


##
# @brief Obsolete?
# @param outfile
# @param has
# @param uas
# @param vas
# @param elevation
# @param mean_stage
# @param zscale
# @param verbose
def obsolete_write_sww_time_slices(outfile, has, uas, vas, elevation,
                                   mean_stage=0, zscale=1,
                                   verbose=False):
    #Time stepping
    stage = outfile.variables['stage']
    xmomentum = outfile.variables['xmomentum']
    ymomentum = outfile.variables['ymomentum']

    n = len(has)
    j = 0
    for ha, ua, va in map(None, has, uas, vas):
        if verbose and j % ((n+10)/10) == 0: log.critical('  Doing %d of %d'
                                                          % (j, n))
        w = zscale*ha + mean_stage
        stage[j] = w
        h = w - elevation
        xmomentum[j] = ua * h
        ymomentum[j] = -1 * va * h  # -1 since in mux files south is positive.
        j += 1


##
# @brief Convert a set of URS files to a text file.
# @param basename_in Stem path to the 3 URS files.
# @param location_index ??
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
    for quantity, file in map(None, quantities, files_in):
        mux[quantity] = Urs_points(file)

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
                  str(points_utm[li][01]) + '\n')

    # the non-time dependent stuff
    #Title
    fid.write('location_index' + d + 'lat' + d + 'long' + d +
              'Easting' + d + 'Northing' + d + 'depth m' + '\n')
    i = 0
    for depth, point_utm, lat, long in map(None, depths, points_utm,
                                           latitudes, longitudes):

        fid.write(str(i) + d + str(lat) + d + str(long) + d +
                  str(point_utm[0]) + d + str(point_utm[01]) + d +
                  str(depth) + '\n')
        i += 1

    #Time dependent
    if location_index is not None:
        time_step = a_mux.time_step
        i = 0
        #Title
        fid.write('time' + d + 'HA depth m' + d + 'UA momentum East x m/sec' +
                  d + 'VA momentum North y m/sec' + '\n')
        for HA, UA, VA in map(None, mux['HA'], mux['UA'], mux['VA']):
            fid.write(str(i*time_step) + d + str(HA[location_index]) + d +
                      str(UA[location_index]) + d +
                      str(VA[location_index]) + '\n')
            i += 1




##
# @brief A class to write STS files.
class Write_sts:
    sts_quantities = ['stage','xmomentum','ymomentum']
    RANGE = '_range'
    EXTREMA = ':extrema'

    ##
    # @brief Instantiate this STS writer class.
    def __init__(self):
        pass

    ##
    # @brief Write a header to the underlying data file.
    # @param outfile Handle to open file to write.
    # @param times A list of the time slice times *or* a start time.
    # @param number_of_points The number of URS gauge sites.
    # @param description Description string to write into the STS file.
    # @param sts_precision Format of data to write (netcdf constant ONLY).
    # @param verbose True if this function is to be verbose.
    # @note If 'times' is a list, the info will be made relative.
    def store_header(self,
                     outfile,
                     times,
                     number_of_points,
                     description='Converted from URS mux2 format',
                     sts_precision=netcdf_float32,
                     verbose=False):
        """
        outfile - the name of the file that will be written
        times - A list of the time slice times OR a start time
        Note, if a list is given the info will be made relative.
        number_of_points - the number of urs gauges sites
        """

        outfile.institution = 'Geoscience Australia'
        outfile.description = description

        try:
            revision_number = get_revision_number()
        except:
            revision_number = None

        # Allow None to be stored as a string
        outfile.revision_number = str(revision_number)

        # Start time in seconds since the epoch (midnight 1/1/1970)
        # This is being used to seperate one number from a list.
        # what it is actually doing is sorting lists from numeric arrays.
        if isinstance(times, (list, num.ndarray)):
            number_of_times = len(times)
            times = ensure_numeric(times)
            if number_of_times == 0:
                starttime = 0
            else:
                starttime = times[0]
                times = times - starttime  #Store relative times
        else:
            number_of_times = 0
            starttime = times

        outfile.starttime = starttime

        # Dimension definitions
        outfile.createDimension('number_of_points', number_of_points)
        outfile.createDimension('number_of_timesteps', number_of_times)
        outfile.createDimension('numbers_in_range', 2)

        # Variable definitions
        outfile.createVariable('permutation', netcdf_int, ('number_of_points',))
        outfile.createVariable('x', sts_precision, ('number_of_points',))
        outfile.createVariable('y', sts_precision, ('number_of_points',))
        outfile.createVariable('elevation',sts_precision, ('number_of_points',))

        q = 'elevation'
        outfile.createVariable(q + Write_sts.RANGE, sts_precision,
                               ('numbers_in_range',))

        # Initialise ranges with small and large sentinels.
        # If this was in pure Python we could have used None sensibly
        outfile.variables[q + Write_sts.RANGE][0] = max_float  # Min
        outfile.variables[q + Write_sts.RANGE][1] = -max_float # Max

        # Doing sts_precision instead of Float gives cast errors.
        outfile.createVariable('time', netcdf_float, ('number_of_timesteps',))

        for q in Write_sts.sts_quantities:
            outfile.createVariable(q, sts_precision, ('number_of_timesteps',
                                                      'number_of_points'))
            outfile.createVariable(q + Write_sts.RANGE, sts_precision,
                                   ('numbers_in_range',))
            # Initialise ranges with small and large sentinels.
            # If this was in pure Python we could have used None sensibly
            outfile.variables[q + Write_sts.RANGE][0] = max_float  # Min
            outfile.variables[q + Write_sts.RANGE][1] = -max_float # Max

        if isinstance(times, (list, num.ndarray)):
            outfile.variables['time'][:] = times    #Store time relative

        if verbose:
            log.critical('------------------------------------------------')
            log.critical('Statistics:')
            log.critical('    t in [%f, %f], len(t) == %d'
                         % (num.min(times), num.max(times), len(times.flat)))

    ##
    # @brief 
    # @param outfile 
    # @param points_utm 
    # @param elevation 
    # @param zone 
    # @param new_origin 
    # @param points_georeference 
    # @param verbose True if this function is to be verbose.
    def store_points(self,
                     outfile,
                     points_utm,
                     elevation, zone=None, new_origin=None,
                     points_georeference=None, verbose=False):

        """
        points_utm - currently a list or array of the points in UTM.
        points_georeference - the georeference of the points_utm

        How about passing new_origin and current_origin.
        If you get both, do a convertion from the old to the new.

        If you only get new_origin, the points are absolute,
        convert to relative

        if you only get the current_origin the points are relative, store
        as relative.

        if you get no georefs create a new georef based on the minimums of
        points_utm.  (Another option would be to default to absolute)

        Yes, and this is done in another part of the code.
        Probably geospatial.

        If you don't supply either geo_refs, then supply a zone. If not
        the default zone will be used.

        precondition:
             header has been called.
        """

        number_of_points = len(points_utm)
        points_utm = num.array(points_utm)

        # given the two geo_refs and the points, do the stuff
        # described in the method header
        points_georeference = ensure_geo_reference(points_georeference)
        new_origin = ensure_geo_reference(new_origin)

        if new_origin is None and points_georeference is not None:
            points = points_utm
            geo_ref = points_georeference
        else:
            if new_origin is None:
                new_origin = Geo_reference(zone, min(points_utm[:,0]),
                                                 min(points_utm[:,1]))
            points = new_origin.change_points_geo_ref(points_utm,
                                                      points_georeference)
            geo_ref = new_origin

        # At this stage I need a georef and points
        # the points are relative to the georef
        geo_ref.write_NetCDF(outfile)

        x = points[:,0]
        y = points[:,1]
        z = outfile.variables['elevation'][:]

        if verbose:
            log.critical('------------------------------------------------')
            log.critical('More Statistics:')
            log.critical('  Extent (/lon):')
            log.critical('    x in [%f, %f], len(lat) == %d'
                         % (min(x), max(x), len(x)))
            log.critical('    y in [%f, %f], len(lon) == %d'
                         % (min(y), max(y), len(y)))
            log.critical('    z in [%f, %f], len(z) == %d'
                         % (min(elevation), max(elevation), len(elevation)))
            log.critical('geo_ref: %s' % str(geo_ref))
            log.critical('------------------------------------------------')

        z = resize(bath_grid,outfile.variables['elevation'][:].shape)
        outfile.variables['x'][:] = points[:,0] #- geo_ref.get_xllcorner()
        outfile.variables['y'][:] = points[:,1] #- geo_ref.get_yllcorner()
        #outfile.variables['z'][:] = elevation
        outfile.variables['elevation'][:] = elevation  #FIXME HACK4

        # This updates the _range values
        q = 'elevation'
        outfile.variables[q + Write_sts.RANGE][0] = min(elevation)
        outfile.variables[q + Write_sts.RANGE][1] = max(elevation)

    ##
    # @brief Store quantity data in the underlying file.
    # @param outfile
    # @param sts_precision
    # @param slice_index
    # @param time
    # @param verboseTrue if this function is to be verbose.
    # @param **quant Extra keyword args.
    def store_quantities(self, outfile, sts_precision=num.float32,
                         slice_index=None, time=None,
                         verbose=False, **quant):
        """Write the quantity info.

        **quant is extra keyword arguments passed in. These must be
          the sts quantities, currently; stage.

        if the time array is already been built, use the slice_index
        to specify the index.

        Otherwise, use time to increase the time dimension

        Maybe make this general, but the viewer assumes these quantities,
        so maybe we don't want it general - unless the viewer is general

        precondition:
            triangulation and header have been called.
        """

        if time is not None:
            file_time = outfile.variables['time']
            slice_index = len(file_time)
            file_time[slice_index] = time

        # Write the conserved quantities from Domain.
        # Typically stage, xmomentum, ymomentum
        # other quantities will be ignored, silently.
        # Also write the ranges: stage_range
        for q in Write_sts.sts_quantities:
            if not quant.has_key(q):
                msg = 'STS file can not write quantity %s' % q
                raise NewQuantity, msg
            else:
                q_values = quant[q]
                outfile.variables[q][slice_index] = \
                                q_values.astype(sts_precision)

                # This updates the _range values
                q_range = outfile.variables[q + Write_sts.RANGE][:]
                q_values_min = num.min(q_values)
                if q_values_min < q_range[0]:
                    outfile.variables[q + Write_sts.RANGE][0] = q_values_min
                q_values_max = num.max(q_values)
                if q_values_max > q_range[1]:
                    outfile.variables[q + Write_sts.RANGE][1] = q_values_max


##
# @brief 
class Urs_points:
    """
    Read the info in URS mux files.

    for the quantities here's a correlation between the file names and
    what they mean;
    z-mux is height above sea level, m
    e-mux is velocity is Eastern direction, m/s
    n-mux is velocity is Northern direction, m/s
    """

    ##
    # @brief Initialize this instance of Urs_points.
    # @param urs_file Path to the underlying data file.
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
            raise ANUGAError, msg
        if self.time_step_count < 0:
            mux_file.close()
            raise ANUGAError, msg
        if self.time_step < 0:
            mux_file.close()
            raise ANUGAError, msg

        # The depth is in meters, and it is the distance from the ocean
        # to the sea bottom.
        lonlatdep = p_array.array('f')
        lonlatdep.read(mux_file, columns * self.points_num)
        lonlatdep = num.array(lonlatdep, dtype=num.float)
        lonlatdep = num.reshape(lonlatdep, (self.points_num, columns))
        self.lonlatdep = lonlatdep

        self.mux_file = mux_file
        # check this array

    ##
    # @brief Allow iteration over quantity data wrt time.
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

    ##
    # @brief 
    def next(self):
        if self.time_step_count == self.iter_time_step:
            self.close()
            raise StopIteration

        #Read in a time slice from mux file
        hz_p_array = p_array.array('f')
        hz_p_array.read(self.mux_file, self.points_num)
        hz_p = num.array(hz_p_array, dtype=num.float)
        self.iter_time_step += 1

        return hz_p

    ##
    # @brief Close the mux file.
    def close(self):
        self.mux_file.close()

################################################################################
# END URS UNGRIDDED 2 SWW
################################################################################

##
# @brief Store screen output and errors to a file.
# @param dir_name Path to directory for output files (default '.').
# @param myid 
# @param numprocs 
# @param extra_info 
# @param verbose True if this function is to be verbose.
def start_screen_catcher(dir_name=None, myid='', numprocs='', extra_info='',
                         verbose=True):
    """
    Used to store screen output and errors to file, if run on multiple
    processes each processor will have its own output and error file.

    extra_info - is used as a string that can identify outputs with another
    string eg. '_other'

    FIXME: Would be good if you could suppress all the screen output and
    only save it to file... however it seems a bit tricky as this capture
    techique response to sys.stdout and by this time it is already printed out.
    """

    import sys

    if dir_name == None:
        dir_name = getcwd()

    if access(dir_name, W_OK) == 0:
        if verbose: log.critical('Making directory %s' % dir_name)
        mkdir (dir_name,0777)

    if myid != '':
        myid = '_' + str(myid)
    if numprocs != '':
        numprocs = '_' + str(numprocs)
    if extra_info != '':
        extra_info = '_' + str(extra_info)

    screen_output_name = join(dir_name, "screen_output%s%s%s.txt" %
                                        (myid, numprocs, extra_info))
    screen_error_name = join(dir_name, "screen_error%s%s%s.txt" %
                                       (myid, numprocs, extra_info))

    if verbose: log.critical('Starting ScreenCatcher, all output will be '
                             'stored in %s' % screen_output_name)

    # used to catch screen output to file
    sys.stdout = Screen_Catcher(screen_output_name)
    sys.stderr = Screen_Catcher(screen_error_name)


##
# @brief A class to catch stdout and stderr and write to files.
class Screen_Catcher:
    """this simply catches the screen output and stores it to file defined by
    start_screen_catcher (above)
    """

    ##
    # @brief Initialize this instance of Screen_Catcher.
    # @param filename The path to the file to write to.
    def __init__(self, filename):
        self.filename = filename
        if exists(self.filename)is True:
            log.critical('Old existing file "%s" has been deleted'
                         % self.filename)
            remove(self.filename)

    ##
    # @brief Write output to the file.
    # @param stuff The string to write.
    def write(self, stuff):
        fid = open(self.filename, 'a')
        fid.write(stuff)
        fid.close()


##
# @brief Copy a file to a directory, and optionally append another file to it.
# @param dir_name Target directory.
# @param filename Path to file to copy to directory 'dir_name'.
# @param filename2 Optional path to file to append to copied file.
# @param verbose True if this function is to be verbose.
# @note Allow filenames to be either a string or sequence of strings.
def copy_code_files(dir_name, filename1, filename2=None, verbose=False):
    """Copies "filename1" and "filename2" to "dir_name".

    Each 'filename' may be a string or list of filename strings.

    Filenames must be absolute pathnames
    """

    ##
    # @brief copies a file or sequence to destination directory.
    # @param dest The destination directory to copy to.
    # @param file A filename string or sequence of filename strings.
    def copy_file_or_sequence(dest, file):
        if hasattr(file, '__iter__'):
            for f in file:
                shutil.copy(f, dir_name)
                if verbose:
                    log.critical('File %s copied' % f)
        else:
            shutil.copy(file, dir_name)
            if verbose:
                log.critical('File %s copied' % file)

    # check we have a destination directory, create if necessary
    if not os.path.isdir(dir_name):
        if verbose:
            log.critical('Make directory %s' % dir_name)
        mkdir(dir_name, 0777)

    if verbose:
        log.critical('Output directory: %s' % dir_name)        

    copy_file_or_sequence(dir_name, filename1)

    if not filename2 is None:
        copy_file_or_sequence(dir_name, filename2)


##
# @brief Get data from a text file.
# @param filename Path to file to read.
# @param separator_value String to split header line with
# @return (header_fields, data), header_fields is a list of fields from header,
#         data is an array (N columns x M lines) of data from the file.
def get_data_from_file(filename, separator_value=','):
    """
    Read in data information from file and

    Returns:
        header_fields, a string? of the first line separated
        by the 'separator_value'

        data, an array (N data columns X M lines) in the file
        excluding the header

    NOTE: won't deal with columns with different lengths and there must be
          no blank lines at the end.
    """

    fid = open(filename)
    lines = fid.readlines()
    fid.close()

    header_line = lines[0]
    header_fields = header_line.split(separator_value)

    # array to store data, number in there is to allow float...
    # i'm sure there is a better way!
    data = num.array([], dtype=num.float)
    data = num.resize(data, ((len(lines)-1), len(header_fields)))

    array_number = 0
    line_number = 1
    while line_number < len(lines):
        for i in range(len(header_fields)):
            #this get line below the header, explaining the +1
            #and also the line_number can be used as the array index
            fields = lines[line_number].split(separator_value)
            #assign to array
            data[array_number,i] = float(fields[i])

        line_number = line_number + 1
        array_number = array_number + 1

    return header_fields, data


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

################################################################################
# Functions to obtain diagnostics from sww files
################################################################################

##
# @brief Get mesh and quantity data from an SWW file.
# @param filename Path to data file to read.
# @param quantities UNUSED!
# @param verbose True if this function is to be verbose.
# @return (mesh, quantities, time) where mesh is the mesh data, quantities is
#         a dictionary of {name: value}, and time is the time vector.
# @note Quantities extracted: 'elevation', 'stage', 'xmomentum' and 'ymomentum'
def get_mesh_and_quantities_from_file(filename,
                                      quantities=None,
                                      verbose=False):
    """Get and rebuild mesh structure and associated quantities from sww file

    Input:
        filename - Name os sww file
        quantities - Names of quantities to load

    Output:
        mesh - instance of class Interpolate
               (including mesh and interpolation functionality)
        quantities - arrays with quantity values at each mesh node
        time - vector of stored timesteps

    This function is used by e.g.:
        get_interpolated_quantities_at_polyline_midpoints
    """

    # FIXME (Ole): Maybe refactor filefunction using this more fundamental code.

    import types
    from Scientific.IO.NetCDF import NetCDFFile
    from shallow_water import Domain
    from anuga.abstract_2d_finite_volumes.neighbour_mesh import Mesh

    if verbose: log.critical('Reading from %s' % filename)

    fid = NetCDFFile(filename, netcdf_mode_r)    # Open existing file for read
    time = fid.variables['time'][:]    # Time vector
    time += fid.starttime[0]

    # Get the variables as numeric arrays
    x = fid.variables['x'][:]                   # x-coordinates of nodes
    y = fid.variables['y'][:]                   # y-coordinates of nodes
    elevation = fid.variables['elevation'][:]   # Elevation
    stage = fid.variables['stage'][:]           # Water level
    xmomentum = fid.variables['xmomentum'][:]   # Momentum in the x-direction
    ymomentum = fid.variables['ymomentum'][:]   # Momentum in the y-direction

    # Mesh (nodes (Mx2), triangles (Nx3))
    nodes = num.concatenate((x[:,num.newaxis], y[:,num.newaxis]), axis=1)
    triangles = fid.variables['volumes'][:]

    # Get geo_reference
    try:
        geo_reference = Geo_reference(NetCDFObject=fid)
    except: #AttributeError, e:
        # Sww files don't have to have a geo_ref
        geo_reference = None

    if verbose: log.critical('    building mesh from sww file %s' % filename)

    boundary = None

    #FIXME (Peter Row): Should this be in mesh?
    if fid.smoothing != 'Yes':
        nodes = nodes.tolist()
        triangles = triangles.tolist()
        nodes, triangles, boundary = weed(nodes, triangles, boundary)

    try:
        mesh = Mesh(nodes, triangles, boundary, geo_reference=geo_reference)
    except AssertionError, e:
        fid.close()
        msg = 'Domain could not be created: %s. "' % e
        raise DataDomainError, msg

    quantities = {}
    quantities['elevation'] = elevation
    quantities['stage'] = stage
    quantities['xmomentum'] = xmomentum
    quantities['ymomentum'] = ymomentum

    fid.close()

    return mesh, quantities, time


##
# @brief Get values for quantities interpolated to polyline midpoints from SWW.
# @param filename Path to file to read.
# @param quantity_names Quantity names to get.
# @param polyline Representation of desired cross-section.
# @param verbose True if this function is to be verbose.
# @return (segments, i_func) where segments is a list of Triangle_intersection
#         instances and i_func is an instance of Interpolation_function.
# @note For 'polyline' assume absolute UTM coordinates.
def get_interpolated_quantities_at_polyline_midpoints(filename,
                                                      quantity_names=None,
                                                      polyline=None,
                                                      verbose=False):
    """Get values for quantities interpolated to polyline midpoints from SWW

    Input:
        filename - Name of sww file
        quantity_names - Names of quantities to load
        polyline: Representation of desired cross section - it may contain
                  multiple sections allowing for complex shapes. Assume
                  absolute UTM coordinates.
                  Format [[x0, y0], [x1, y1], ...]

    Output:
        segments: list of instances of class Triangle_intersection
        interpolation_function: Instance of class Interpolation_function


      This function is used by get_flow_through_cross_section and
      get_energy_through_cross_section
    """

    from anuga.fit_interpolate.interpolate import Interpolation_function

    # Get mesh and quantities from sww file
    X = get_mesh_and_quantities_from_file(filename,
                                          quantities=quantity_names,
                                          verbose=verbose)
    mesh, quantities, time = X

    # Find all intersections and associated triangles.
    segments = mesh.get_intersecting_segments(polyline, verbose=verbose)

    # Get midpoints
    interpolation_points = segment_midpoints(segments)

    # Interpolate
    if verbose:
        log.critical('Interpolating - total number of interpolation points = %d'
                     % len(interpolation_points))

    I = Interpolation_function(time,
                               quantities,
                               quantity_names=quantity_names,
                               vertex_coordinates=mesh.nodes,
                               triangles=mesh.triangles,
                               interpolation_points=interpolation_points,
                               verbose=verbose)

    return segments, I


##
# @brief Obtain flow (m^3/s) perpendicular to specified cross section.
# @param filename Path to file to read.
# @param polyline Representation of desired cross-section.
# @param verbose Trie if this function is to be verbose.
# @return (time, Q) where time and Q are lists of time and flow respectively.
def get_flow_through_cross_section(filename, polyline, verbose=False):
    """Obtain flow (m^3/s) perpendicular to specified cross section.

    Inputs:
        filename: Name of sww file
        polyline: Representation of desired cross section - it may contain
                  multiple sections allowing for complex shapes. Assume
                  absolute UTM coordinates.
                  Format [[x0, y0], [x1, y1], ...]

    Output:
        time: All stored times in sww file
        Q: Hydrograph of total flow across given segments for all stored times.

    The normal flow is computed for each triangle intersected by the polyline
    and added up.  Multiple segments at different angles are specified the
    normal flows may partially cancel each other.

    The typical usage of this function would be to get flow through a channel,
    and the polyline would then be a cross section perpendicular to the flow.
    """

    quantity_names =['elevation',
                     'stage',
                     'xmomentum',
                     'ymomentum']

    # Get values for quantities at each midpoint of poly line from sww file
    X = get_interpolated_quantities_at_polyline_midpoints(filename,
                                                          quantity_names=\
                                                              quantity_names,
                                                          polyline=polyline,
                                                          verbose=verbose)
    segments, interpolation_function = X

    # Get vectors for time and interpolation_points
    time = interpolation_function.time
    interpolation_points = interpolation_function.interpolation_points

    if verbose: log.critical('Computing hydrograph')

    # Compute hydrograph
    Q = []
    for t in time:
        total_flow = 0
        for i in range(len(interpolation_points)):
            elevation, stage, uh, vh = interpolation_function(t, point_id=i)
            normal = segments[i].normal

            # Inner product of momentum vector with segment normal [m^2/s]
            normal_momentum = uh*normal[0] + vh*normal[1]

            # Flow across this segment [m^3/s]
            segment_flow = normal_momentum * segments[i].length

            # Accumulate
            total_flow += segment_flow

        # Store flow at this timestep
        Q.append(total_flow)


    return time, Q


##
# @brief Get average energy across a cross-section.
# @param filename Path to file of interest.
# @param polyline Representation of desired cross-section.
# @param kind Select energy to compute: 'specific' or 'total'.
# @param verbose True if this function is to be verbose.
# @return (time, E) where time and E are lists of timestep and energy.
def get_energy_through_cross_section(filename,
                                     polyline,
                                     kind='total',
                                     verbose=False):
    """Obtain average energy head [m] across specified cross section.

    Inputs:
        polyline: Representation of desired cross section - it may contain
                  multiple sections allowing for complex shapes. Assume
                  absolute UTM coordinates.
                  Format [[x0, y0], [x1, y1], ...]
        kind:     Select which energy to compute.
                  Options are 'specific' and 'total' (default)

    Output:
        E: Average energy [m] across given segments for all stored times.

    The average velocity is computed for each triangle intersected by
    the polyline and averaged weighted by segment lengths.

    The typical usage of this function would be to get average energy of
    flow in a channel, and the polyline would then be a cross section
    perpendicular to the flow.

    #FIXME (Ole) - need name for this energy reflecting that its dimension
    is [m].
    """

    from anuga.config import g, epsilon, velocity_protection as h0

    quantity_names =['elevation',
                     'stage',
                     'xmomentum',
                     'ymomentum']

    # Get values for quantities at each midpoint of poly line from sww file
    X = get_interpolated_quantities_at_polyline_midpoints(filename,
                                                          quantity_names=\
                                                              quantity_names,
                                                          polyline=polyline,
                                                          verbose=verbose)
    segments, interpolation_function = X

    # Get vectors for time and interpolation_points
    time = interpolation_function.time
    interpolation_points = interpolation_function.interpolation_points

    if verbose: log.critical('Computing %s energy' % kind)

    # Compute total length of polyline for use with weighted averages
    total_line_length = 0.0
    for segment in segments:
        total_line_length += segment.length

    # Compute energy
    E = []
    for t in time:
        average_energy = 0.0
        for i, p in enumerate(interpolation_points):
            elevation, stage, uh, vh = interpolation_function(t, point_id=i)

            # Depth
            h = depth = stage-elevation

            # Average velocity across this segment
            if h > epsilon:
                # Use protection against degenerate velocities
                u = uh / (h + h0/h)
                v = vh / (h + h0/h)
            else:
                u = v = 0.0

            speed_squared = u*u + v*v
            kinetic_energy = 0.5 * speed_squared / g

            if kind == 'specific':
                segment_energy = depth + kinetic_energy
            elif kind == 'total':
                segment_energy = stage + kinetic_energy
            else:
                msg = 'Energy kind must be either "specific" or "total". '
                msg += 'I got %s' % kind

            # Add to weighted average
            weigth = segments[i].length / total_line_length
            average_energy += segment_energy * weigth

        # Store energy at this timestep
        E.append(average_energy)

    return time, E


##
# @brief Return highest elevation where depth > 0.
# @param filename Path to SWW file of interest.
# @param polygon If specified resrict to points inside this polygon.
# @param time_interval If specified resrict to within the time specified.
# @param verbose True if this function is  to be verbose.
def get_maximum_inundation_elevation(filename,
                                     polygon=None,
                                     time_interval=None,
                                     verbose=False):
    """Return highest elevation where depth > 0

    Usage:
    max_runup = get_maximum_inundation_elevation(filename,
                                                 polygon=None,
                                                 time_interval=None,
                                                 verbose=False)

    filename is a NetCDF sww file containing ANUGA model output.
    Optional arguments polygon and time_interval restricts the maximum
    runup calculation
    to a points that lie within the specified polygon and time interval.

    If no inundation is found within polygon and time_interval the return value
    is None signifying "No Runup" or "Everything is dry".

    See general function get_maximum_inundation_data for details.
    """

    runup, _ = get_maximum_inundation_data(filename,
                                           polygon=polygon,
                                           time_interval=time_interval,
                                           verbose=verbose)
    return runup


##
# @brief Return location of highest elevation where h > 0
# @param filename Path to SWW file to read.
# @param polygon If specified resrict to points inside this polygon.
# @param time_interval If specified resrict to within the time specified.
# @param verbose True if this function is  to be verbose.
def get_maximum_inundation_location(filename,
                                    polygon=None,
                                    time_interval=None,
                                    verbose=False):
    """Return location of highest elevation where h > 0

    Usage:
    max_runup_location = get_maximum_inundation_location(filename,
                                                         polygon=None,
                                                         time_interval=None,
                                                         verbose=False)

    filename is a NetCDF sww file containing ANUGA model output.
    Optional arguments polygon and time_interval restricts the maximum
    runup calculation
    to a points that lie within the specified polygon and time interval.

    If no inundation is found within polygon and time_interval the return value
    is None signifying "No Runup" or "Everything is dry".

    See general function get_maximum_inundation_data for details.
    """

    _, max_loc = get_maximum_inundation_data(filename,
                                             polygon=polygon,
                                             time_interval=time_interval,
                                             verbose=verbose)
    return max_loc


##
# @brief Compute maximum run up height from SWW file.
# @param filename Path to SWW file to read.
# @param polygon If specified resrict to points inside this polygon.
# @param time_interval If specified resrict to within the time specified.
# @param use_centroid_values 
# @param verbose True if this function is to be verbose.
# @return (maximal_runup, maximal_runup_location)
def get_maximum_inundation_data(filename, polygon=None, time_interval=None,
                                use_centroid_values=False,
                                verbose=False):
    """Compute maximum run up height from sww file.

    Usage:
    runup, location = get_maximum_inundation_data(filename,
                                                  polygon=None,
                                                  time_interval=None,
                                                  verbose=False)

    Algorithm is as in get_maximum_inundation_elevation from
    shallow_water_domain except that this function works with the sww file and
    computes the maximal runup height over multiple timesteps.

    Optional arguments polygon and time_interval restricts the maximum runup
    calculation to a points that lie within the specified polygon and time
    interval.

    Polygon is assumed to be in (absolute) UTM coordinates in the same zone
    as domain.

    If no inundation is found within polygon and time_interval the return value
    is None signifying "No Runup" or "Everything is dry".
    """

    # We are using nodal values here as that is what is stored in sww files.

    # Water depth below which it is considered to be 0 in the model
    # FIXME (Ole): Allow this to be specified as a keyword argument as well

    from anuga.utilities.polygon import inside_polygon
    from anuga.config import minimum_allowed_height
    from Scientific.IO.NetCDF import NetCDFFile

    dir, base = os.path.split(filename)

    iterate_over = get_all_swwfiles(dir, base)

    # Read sww file
    if verbose: log.critical('Reading from %s' % filename)
    # FIXME: Use general swwstats (when done)

    maximal_runup = None
    maximal_runup_location = None

    for file, swwfile in enumerate (iterate_over):
        # Read sww file
        filename = join(dir, swwfile+'.sww')

        if verbose: log.critical('Reading from %s' % filename)
        # FIXME: Use general swwstats (when done)

        fid = NetCDFFile(filename)

        # Get geo_reference
        # sww files don't have to have a geo_ref
        try:
            geo_reference = Geo_reference(NetCDFObject=fid)
        except AttributeError, e:
            geo_reference = Geo_reference() # Default georef object

        xllcorner = geo_reference.get_xllcorner()
        yllcorner = geo_reference.get_yllcorner()
        zone = geo_reference.get_zone()

        # Get extent
        volumes = fid.variables['volumes'][:]
        x = fid.variables['x'][:] + xllcorner
        y = fid.variables['y'][:] + yllcorner

        # Get the relevant quantities (Convert from single precison)
        elevation = num.array(fid.variables['elevation'][:], num.float)
        stage = num.array(fid.variables['stage'][:], num.float)

        # Here's where one could convert nodal information to centroid
        # information but is probably something we need to write in C.
        # Here's a Python thought which is NOT finished!!!
        if use_centroid_values is True:
            x = get_centroid_values(x, volumes)
            y = get_centroid_values(y, volumes)
            elevation = get_centroid_values(elevation, volumes)

        # Spatial restriction
        if polygon is not None:
            msg = 'polygon must be a sequence of points.'
            assert len(polygon[0]) == 2, msg
            # FIXME (Ole): Make a generic polygon input check in polygon.py
            # and call it here
            points = num.ascontiguousarray(num.concatenate((x[:,num.newaxis],
                                                            y[:,num.newaxis]),
                                                            axis=1))
            point_indices = inside_polygon(points, polygon)

            # Restrict quantities to polygon
            elevation = num.take(elevation, point_indices, axis=0)
            stage = num.take(stage, point_indices, axis=1)

            # Get info for location of maximal runup
            points_in_polygon = num.take(points, point_indices, axis=0)

            x = points_in_polygon[:,0]
            y = points_in_polygon[:,1]
        else:
            # Take all points
            point_indices = num.arange(len(x))

        # Temporal restriction
        time = fid.variables['time'][:]
        all_timeindices = num.arange(len(time))
        if time_interval is not None:
            msg = 'time_interval must be a sequence of length 2.'
            assert len(time_interval) == 2, msg
            msg = 'time_interval %s must not be decreasing.' % time_interval
            assert time_interval[1] >= time_interval[0], msg
            msg = 'Specified time interval [%.8f:%.8f] ' % tuple(time_interval)
            msg += 'must does not match model time interval: [%.8f, %.8f]\n' \
                   % (time[0], time[-1])
            if time_interval[1] < time[0]: raise ValueError(msg)
            if time_interval[0] > time[-1]: raise ValueError(msg)

            # Take time indices corresponding to interval (& is bitwise AND)
            timesteps = num.compress((time_interval[0] <= time) \
                                     & (time <= time_interval[1]),
                                     all_timeindices)

            msg = 'time_interval %s did not include any model timesteps.' \
                  % time_interval
            assert not num.alltrue(timesteps == 0), msg
        else:
            # Take them all
            timesteps = all_timeindices

        fid.close()

        # Compute maximal runup for each timestep
        #maximal_runup = None
        #maximal_runup_location = None
        #maximal_runups = [None]
        #maximal_runup_locations = [None]

        for i in timesteps:
            if use_centroid_values is True:
                stage_i = get_centroid_values(stage[i,:], volumes)
            else:
                stage_i = stage[i,:]

            depth = stage_i - elevation

            # Get wet nodes i.e. nodes with depth>0 within given region
            # and timesteps
            wet_nodes = num.compress(depth > minimum_allowed_height,
                                     num.arange(len(depth)))

            if num.alltrue(wet_nodes == 0):
                runup = None
            else:
                # Find maximum elevation among wet nodes
                wet_elevation = num.take(elevation, wet_nodes, axis=0)
                runup_index = num.argmax(wet_elevation)
                runup = max(wet_elevation)
                assert wet_elevation[runup_index] == runup       # Must be True

            if runup > maximal_runup:
                maximal_runup = runup      # works even if maximal_runup is None

                # Record location
                wet_x = num.take(x, wet_nodes, axis=0)
                wet_y = num.take(y, wet_nodes, axis=0)
                maximal_runup_location = [wet_x[runup_index],wet_y[runup_index]]

    return maximal_runup, maximal_runup_location


##
# @brief Find all SWW files in a directory with given stem name.
# @param look_in_dir The directory to look in.
# @param base_name The file stem name.
# @param verbose True if this function is to be verbose.
# @return A list of found filename strings.
# @note Will accept 'base_name' with or without '.sww' extension.
# @note If no files found, raises IOError exception.
def get_all_swwfiles(look_in_dir='', base_name='', verbose=False):
    '''
    Finds all the sww files in a "look_in_dir" which contains a "base_name".
    will accept base_name with or without the extension ".sww"

    Returns: a list of strings

    Usage:     iterate_over = get_all_swwfiles(dir, name)
    then
               for swwfile in iterate_over:
                   do stuff

    Check "export_grids" and "get_maximum_inundation_data" for examples
    '''

    # plus tests the extension
    name, extension = os.path.splitext(base_name)

    if extension != '' and extension != '.sww':
        msg = 'file %s%s must be a NetCDF sww file!' % (base_name, extension)
        raise IOError, msg

    if look_in_dir == "":
        look_in_dir = "."                                   # Unix compatibility

    dir_ls = os.listdir(look_in_dir)
    iterate_over = [x[:-4] for x in dir_ls if name in x and x[-4:] == '.sww']
    if len(iterate_over) == 0:
        msg = 'No files of the base name %s' % name
        raise IOError, msg

    if verbose: log.critical('iterate over %s' % iterate_over)

    return iterate_over


##
# @brief Find all files in a directory that contain a string and have extension.
# @param look_in_dir Path to the directory to look in.
# @param base_name Stem filename of the file(s) of interest.
# @param extension Extension of the files to look for.
# @param verbose True if this function is to be verbose.
# @return A list of found filename strings.
# @note If no files found, raises IOError exception.
def get_all_files_with_extension(look_in_dir='',
                                 base_name='',
                                 extension='.sww',
                                 verbose=False):
    '''Find all files in a directory with given stem name.
    Finds all the sww files in a "look_in_dir" which contains a "base_name".

    Returns: a list of strings

    Usage:     iterate_over = get_all_swwfiles(dir, name)
    then
               for swwfile in iterate_over:
                   do stuff

    Check "export_grids" and "get_maximum_inundation_data" for examples
    '''

    # plus tests the extension
    name, ext = os.path.splitext(base_name)

    if ext != '' and ext != extension:
        msg = 'base_name %s must be a file with %s extension!' \
              % (base_name, extension)
        raise IOError, msg

    if look_in_dir == "":
        look_in_dir = "."                               # Unix compatibility

    dir_ls = os.listdir(look_in_dir)
    iterate_over = [x[:-4] for x in dir_ls if name in x and x[-4:] == extension]

    if len(iterate_over) == 0:
        msg = 'No files of the base name %s in %s' % (name, look_in_dir)
        raise IOError, msg

    if verbose: log.critical('iterate over %s' % iterate_over)

    return iterate_over


##
# @brief Find all files in a directory that contain a given string.
# @param look_in_dir Path to the directory to look in.
# @param base_name String that files must contain.
# @param verbose True if this function is to be verbose.
def get_all_directories_with_name(look_in_dir='', base_name='', verbose=False):
    '''
    Finds all the directories in a "look_in_dir" which contains a "base_name".

    Returns: a list of strings

    Usage:     iterate_over = get_all_directories_with_name(dir, name)
    then:      for swwfile in iterate_over:
                   do stuff

    Check "export_grids" and "get_maximum_inundation_data" for examples
    '''

    if look_in_dir == "":
        look_in_dir = "."                                  # Unix compatibility

    dir_ls = os.listdir(look_in_dir)
    iterate_over = [x for x in dir_ls if base_name in x]

    if len(iterate_over) == 0:
        msg = 'No files of the base name %s' % base_name
        raise IOError, msg

    if verbose: log.critical('iterate over %s' % iterate_over)

    return iterate_over


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
