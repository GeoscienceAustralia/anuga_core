"""Class Geospatial_data

Manipulation of locations on the planet and associated attributes.
"""



from sys import maxsize
from os import access, F_OK, R_OK, remove
#from types import DictType
from warnings import warn
#from string import lower
from copy import deepcopy
import copy

try:
    from exceptions import Exception
except:
    pass

from anuga.file.netcdf import NetCDFFile
import numpy as num
from numpy.random import randint, seed

from anuga.coordinate_transforms.lat_long_UTM_conversion import UTMtoLL
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.coordinate_transforms.geo_reference import Geo_reference, \
    TitleError, DEFAULT_ZONE, ensure_geo_reference, write_NetCDF_georeference
from anuga.coordinate_transforms.redfearn import convert_from_latlon_to_utm
from anuga.utilities.system_tools import clean_line
from anuga.anuga_exceptions import ANUGAError
from anuga.config import points_file_block_line_size as MAX_READ_LINES
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
from anuga.config import netcdf_float
import anuga.utilities.log as log


DEFAULT_ATTRIBUTE = 'elevation'


class Geospatial_data(object):

    def __init__(self,
                 data_points=None,  # this can also be a points file name
                 attributes=None,
                 geo_reference=None,
                 default_attribute_name=None,
                 file_name=None,
                 latitudes=None,
                 longitudes=None,
                 points_are_lats_longs=False,
                 max_read_lines=None,
                 load_file_now=True,
                 verbose=False):
        """Create instance from data points and associated attributes

        data_points: x,y coordinates in meters. Type must be either a
        sequence of 2-tuples or an Mx2 numeric array of floats.  A file name
        with extension .txt, .cvs, .xya or .pts can also be passed in here.

        attributes: Associated values for each data point. The type
        must be either a list or an array of length M or a dictionary
        of lists (or arrays) of length M. In the latter case the keys
        in the dictionary represent the attribute names, in the former
        the attribute will get the default name "elevation".

        geo_reference: Object representing the origin of the data
        points. It contains UTM zone, easting and northing and data
        points are assumed to be relative to this origin.
        If geo_reference is None, the default geo ref object is used.

        default_attribute_name: Name of default attribute to be used with
        get_attribute_values. The idea is that the dataset can be
        equipped with information about which attribute to return.
        If None, the default is the "first"

        latitudes, longitudes: Vectors of latitudes and longitudes,
        used to specify location instead of points.

        points_are_lats_longs: Set this as true if the points are actually
        lats and longs, not UTM

        max_read_lines: The number of rows read into memory when using
        blocking to read a file.

        load_file_now:  If true the file is automatically loaded
        into the geospatial instance. Used when blocking.

        file_name: Name of input netCDF file or  text file (.xya, .txt). An netCDF file must
        have dimensions "points" etc.
        
        Text file is a comma seperated file with x, y and attribute
        data.

        The first line has the titles of the columns.  The first two
        column titles are checked to see if they start with lat or
        long (not case sensitive).  If so the data is assumed to be
        latitude and longitude, in decimal format and converted to
        UTM.  Otherwise the first two columns are assumed to be the x
        and y, and the title names acually used are ignored.


        The format for a text file is:
            1st line:     [column names]
            other lines:  x y [attributes]

            for example:
            x, y, elevation, friction
            0.6, 0.7, 4.9, 0.3
            1.9, 2.8, 5, 0.3
            2.7, 2.4, 5.2, 0.3

        The first two columns have to be x, y or lat, long
        coordinates.


        The format for a Points dictionary is:
          ['pointlist'] a 2 column array describing points. 1st column x,
          2nd column y.
          ['attributelist'], a dictionary of 1D arrays, representing
          attribute values at the point.  The dictionary key is the attribute
          header.
          ['geo_reference'] a Geo_refernece object. Use if the point
          information is relative. This is optional.
            eg
            dic['pointlist'] = [[1.0,2.0],[3.0,5.0]]
            dic['attributelist']['elevation'] = [[7.0,5.0]

        verbose:
        """

        if isinstance(data_points, str):
            # assume data_points is really a file name
            file_name = data_points

        self.set_verbose(verbose)
        self.file_name = file_name

        if geo_reference is None:
            self.geo_reference = Geo_reference()  # Use default geo_reference
        else:
            self.geo_reference = geo_reference

        if max_read_lines is None:
            self.max_read_lines = int(MAX_READ_LINES)
        else:
            self.max_read_lines = int(max_read_lines)

        if file_name is None:
            if (latitudes is not None or longitudes is not None
                    or points_are_lats_longs):
                data_points, geo_reference = \
                    _set_using_lat_long(latitudes=latitudes,
                                        longitudes=longitudes,
                                        geo_reference=geo_reference,
                                        data_points=data_points,
                                        points_are_lats_longs=points_are_lats_longs)
            self.check_data_points(data_points)
            self.set_attributes(attributes)
            self.set_geo_reference(geo_reference)
            self.set_default_attribute_name(default_attribute_name)
        elif load_file_now is True:
            # watch for case where file name and points,
            # attributes etc are provided!!
            # if file name then all provided info will be removed!

            if verbose is True:
                if file_name is not None:
                    log.critical('Geospatial_data: Loading Geospatial data from file: %s'
                                 % file_name)

            self.import_points_file(file_name, verbose=verbose)

            self.check_data_points(self.data_points)
            self.set_attributes(self.attributes)
            self.set_geo_reference(self.geo_reference)
            self.set_default_attribute_name(default_attribute_name)

        if verbose is True:
            if file_name is not None:
                log.critical('Geospatial_data: Created from file: %s'
                             % file_name)
                if load_file_now is False:
                    log.critical('Data will be loaded blockwise on demand')

                    if file_name.endswith('csv') or file_name.endswith('txt'):
                        pass
                        # This message was misleading.
                        # FIXME (Ole): Are we blocking here or not?
                        # print 'ASCII formats are not that great for '
                        # print 'blockwise reading. Consider storing this'
                        # print 'data as a pts NetCDF format'

    def __len__(self):
        return len(self.data_points)

    def __repr__(self):
        return str(self.get_data_points(absolute=True))

    def check_data_points(self, data_points):
        """Checks data points"""

        if data_points is None:
            self.data_points = None
            msg = 'There is no data or file provided!'
            raise ValueError(msg)
        else:
            self.data_points = ensure_numeric(data_points)
            if not (0,) == self.data_points.shape:
                assert len(self.data_points.shape) == 2
                assert self.data_points.shape[1] == 2

    def set_attributes(self, attributes):
        """Check and assign attributes dictionary"""

        if attributes is None:
            self.attributes = None
            return

        if not isinstance(attributes, dict):
            # Convert single attribute into dictionary
            attributes = {DEFAULT_ATTRIBUTE: attributes}

        # Check input attributes
        for key in list(attributes.keys()):
            try:
                attributes[key] = ensure_numeric(attributes[key])
            except:
                msg = ("Attribute '%s' (%s) could not be converted to a"
                       "numeric vector" % (str(key), str(attributes[key])))
                raise Exception(msg)

        self.attributes = attributes

    def set_geo_reference(self, geo_reference):
        """Set the georeference of geospatial data.

        It can also be used to change the georeference and will ensure that
        the absolute coordinate values are unchanged.
        """

        if geo_reference is None:
            # Use default - points are in absolute coordinates
            geo_reference = Geo_reference()

        # Allow for tuple (zone, xllcorner, yllcorner)
        geo_reference = ensure_geo_reference(geo_reference)

        if not isinstance(geo_reference, Geo_reference):
            # FIXME (Ole): This exception will be raised even
            # if geo_reference is None. Is that the intent Duncan?
            msg = ('Argument geo_reference must be a valid Geo_reference '
                   'object or None.')
            raise (msg)

        # If a geo_reference already exists, change the point data according to
        # the new geo reference
        if self.geo_reference is not None:
            self.data_points = self.get_data_points(
                geo_reference=geo_reference)

        self.geo_reference = geo_reference

    def set_default_attribute_name(self, default_attribute_name):
        self.default_attribute_name = default_attribute_name

    def set_verbose(self, verbose=False):
        if verbose in [False, True]:
            self.verbose = verbose
        else:
            msg = 'Illegal value: %s' % str(verbose)
            raise Exception(msg)

    def clip(self, polygon, closed=True, verbose=False):
        """Clip geospatial data by a polygon

        Input
          polygon - Either a list of points, an Nx2 array or
                    a Geospatial data object.
          closed - (optional) determine whether points on boundary should be
          regarded as belonging to the polygon (closed = True)
          or not (closed = False). Default is True.

        Output
          New geospatial data object representing points inside
          specified polygon.


        Note - this method is non-destructive and leaves the data in 'self'
        unchanged
        """

        from anuga.geometry.polygon import inside_polygon

        if isinstance(polygon, Geospatial_data):
            # Polygon is an object - extract points
            polygon = polygon.get_data_points()

        points = self.get_data_points()
        inside_indices = inside_polygon(points, polygon, closed, verbose)

        clipped_G = self.get_sample(inside_indices)

        return clipped_G

    def clip_outside(self, polygon, closed=True, verbose=False):
        """Clip geospatial date by a polygon, keeping data OUTSIDE of polygon

        Input
          polygon - Either a list of points, an Nx2 array or
                    a Geospatial data object.
          closed - (optional) determine whether points on boundary should be
          regarded as belonging to the polygon (closed = True)
          or not (closed = False). Default is True.

        Output
          Geospatial data object representing point OUTSIDE specified polygon
        """

        from anuga.geometry.polygon import outside_polygon

        if isinstance(polygon, Geospatial_data):
            # Polygon is an object - extract points
            polygon = polygon.get_data_points()

        points = self.get_data_points()
        outside_indices = outside_polygon(points, polygon, closed, verbose)

        clipped_G = self.get_sample(outside_indices)

        return clipped_G

    def get_geo_reference(self):
        return self.geo_reference

    def get_data_points(self,
                        absolute=True,
                        geo_reference=None,
                        as_lat_long=False,
                        isSouthHemisphere=True):
        """Get coordinates for all data points as an Nx2 array

        If absolute is False returned coordinates are relative to the
        internal georeference's xll and yll corners, otherwise
        absolute UTM coordinates are returned.

        If a geo_reference is passed the points are returned relative
        to that geo_reference.

        isSH (isSouthHemisphere) is only used when getting data
        points "as_lat_long" is True and if FALSE will return lats and
        longs valid for the Northern Hemisphere.

        Default: absolute is True.
        """

        if as_lat_long is True:
            msg = "Points need a zone to be converted into lats and longs"
            assert self.geo_reference is not None, msg
            zone = self.geo_reference.get_zone()
            assert self.geo_reference.get_zone() is not DEFAULT_ZONE, msg
            lats_longs = []
            for point in self.get_data_points(True):
                # UTMtoLL(northing, easting, zone,
                lat_calced, long_calced = UTMtoLL(point[1], point[0],
                                                  zone, isSouthHemisphere)
                lats_longs.append((lat_calced, long_calced))  # to hash
            return lats_longs

        if absolute is True and geo_reference is None:
            return self.geo_reference.get_absolute(self.data_points)
        elif geo_reference is not None:
            return geo_reference.change_points_geo_ref(self.data_points,
                                                       self.geo_reference)
        else:
            # If absolute is False
            return self.data_points

    def get_attributes(self, attribute_name=None):
        """Return values for one named attribute.

        If attribute_name is None, default_attribute_name is used
        """

        if attribute_name is None:
            if self.default_attribute_name is not None:
                attribute_name = self.default_attribute_name
            else:
                attribute_name = list(self.attributes.keys())[0]
                # above line takes the first one from keys

        if self.verbose is True:
            log.critical('Geospatial_data: Using attribute %s' %
                         attribute_name)
            log.critical('Geospatial_data: Available attributes: %s' %
                         (list(self.attributes.keys())))

        msg = 'Attribute name %s does not exist in data set' % attribute_name
        assert attribute_name in self.attributes, msg

        return self.attributes[attribute_name]

    def get_all_attributes(self):
        """Return values for all attributes.
        The return value is either None or a dictionary (possibly empty).
        """

        return self.attributes

    def __add__(self, other):
        """Returns the addition of 2 geospatial objects,
        objects are concatencated to the end of each other

        NOTE: doesn't add if objects contain different
        attributes

        FIXME(Ole): I am not sure this is a good idea, but leave it for now until all tests pass
        Always return absolute points (i.e. xllcorner = yllcorner = 0)
        This also means, that if you add None to the object,
        it will be turned into absolute coordinates

        other can be None in which case nothing is added to self.
        """

        # find objects zone and checks if the same
        geo_ref1 = self.get_geo_reference()
        zone1 = geo_ref1.get_zone()

        if other is not None:
            geo_ref2 = other.get_geo_reference()
            zone2 = geo_ref2.get_zone()
            geo_ref1.reconcile_zones(geo_ref2)
            new_points = num.concatenate((self.get_data_points(absolute=True),
                                          other.get_data_points(absolute=True)),
                                         axis=0)

            # Concatenate attributes if any
            if self.attributes is None:
                if other.attributes is not None:
                    msg = ('Geospatial data must have the same '
                           'attributes to allow addition.')
                    raise Exception(msg)

                new_attributes = None
            else:
                new_attributes = {}
                for x in list(self.attributes.keys()):
                    if x in other.attributes:
                        attrib1 = self.attributes[x]
                        attrib2 = other.attributes[x]
                        new_attributes[x] = num.concatenate((attrib1, attrib2),
                                                            axis=0)  # ??default#
                    else:
                        msg = ('Geospatial data must have the same '
                               'attributes to allow addition.')
                        raise Exception(msg)
        else:
            new_points = self.get_data_points(absolute=True)
            new_attributes = self.attributes

        # Instantiate new data object and return absolute coordinates
        # FIXME (Ole): This does not make sense - need to revisit.
        new_geo_ref = Geo_reference(geo_ref1.zone, 0, 0)

        return Geospatial_data(data_points=new_points,
                               attributes=new_attributes,
                               geo_reference=new_geo_ref)


    def __radd__(self, other):
        """Handle cases like None + Geospatial_data(...)"""

        return self + other

################################################################################
#  IMPORT/EXPORT POINTS FILES
################################################################################

    def import_points_file(self, file_name, delimiter=None, verbose=False):
        """ load an .txt, .xya, .csv or .pts file

        Note: will throw an IOError/SyntaxError if it can't load the file.
        Catch these!

        Post condition: self.attributes dictionary has been set
        """

        if access(file_name, F_OK) == 0:
            msg = 'File %s does not exist or is not accessible' % file_name
            raise IOError(msg)

        attributes = {}
        if file_name[-4:] == ".pts":
            try:
                data_points, attributes, geo_reference = \
                    _read_pts_file(file_name, verbose)
            except IOError as e:
                msg = 'Could not open file %s ' % file_name
                raise IOError(msg)
        elif file_name[-4:] == ".txt" or file_name[-4:] == ".csv" or file_name[-4:] == ".xya":
            try:
                data_points, attributes, geo_reference = \
                    _read_csv_file(file_name, verbose)
            except IOError as e:
                # This should only be if a file is not found
                msg = ('Could not open file %s. Check the file location.'
                       % file_name)
                raise IOError(msg)
            except SyntaxError as e:
                # This should only be if there is a format error
                msg = ('Problem with format of file %s.\n%s'
                       % (file_name, Error_message['IOError']))
                raise SyntaxError(msg)
        else:
            msg = 'Extension %s is unknown' % file_name[-4:]
            raise IOError(msg)

        self.data_points = data_points
        self.attributes = attributes
        self.geo_reference = geo_reference

    def export_points_file(self, file_name, absolute=True,
                           as_lat_long=False, isSouthHemisphere=True):
        """write a points file as a text (.csv) or binary (.pts) file

        file_name is the file name, including the extension
        The point_dict is defined at the top of this file.

        If absolute is True data the xll and yll are added to the points value
        and the xll and yll of the geo_reference are set to 0.

        If absolute is False data points at returned as relative to the xll
        and yll and geo_reference remains uneffected

        isSouthHemisphere: is only used when getting data
        points "as_lat_long" is True and if FALSE will return lats and
        longs valid for the Northern Hemisphere.
        """

        if (file_name[-4:] == ".pts"):
            if absolute is True:
                geo_ref = deepcopy(self.geo_reference)
                geo_ref.xllcorner = 0
                geo_ref.yllcorner = 0
                _write_pts_file(file_name,
                                self.get_data_points(absolute),
                                self.get_all_attributes(),
                                geo_ref)
            else:
                _write_pts_file(file_name,
                                self.get_data_points(absolute),
                                self.get_all_attributes(),
                                self.get_geo_reference())
        elif file_name[-4:] == ".txt" or file_name[-4:] == ".csv" or file_name[-4:] == ".xya":
            msg = "ERROR: trying to write a .txt file with relative data."
            assert absolute, msg
            _write_csv_file(file_name,
                            self.get_data_points(absolute=True,
                                                 as_lat_long=as_lat_long,
                                                 isSouthHemisphere=isSouthHemisphere),
                            self.get_all_attributes(),
                            as_lat_long=as_lat_long)
        elif file_name[-4:] == ".urs":
            msg = "ERROR: Can not write a .urs file as a relative file."
            assert absolute, msg
            _write_urs_file(file_name,
                            self.get_data_points(as_lat_long=True,
                                                 isSouthHemisphere=isSouthHemisphere))
        else:
            msg = 'Unknown file type %s ' % file_name
            raise IOError(msg)

    def get_sample(self, indices):
        """ Returns a object which is a subset of the original
        and the data points and attributes in this new object refer to
        the indices provided

        Input
            indices- a list of integers that represent the new object
        Output
            New geospatial data object representing points specified by
            the indices
        """

        # FIXME: add the geo_reference to this
        points = self.get_data_points()
        sampled_points = num.take(points, indices, axis=0)

        attributes = self.get_all_attributes()

        sampled_attributes = {}
        if attributes is not None:
            for key, att in list(attributes.items()):
                sampled_attributes[key] = num.take(att, indices, axis=0)

        return Geospatial_data(sampled_points, sampled_attributes)

    def split(self, factor=0.5, seed_num=None, verbose=False):
        """Returns two geospatial_data object, first is the size of the 'factor'
        smaller the original and the second is the remainder. The two
        new objects are disjoint sets of each other.

        Points of the two new object have selected RANDOMLY.

        This method create two lists of indices which are passed into
        get_sample.  The lists are created using random numbers, and
        they are unique sets eg.  total_list(1,2,3,4,5,6,7,8,9)
        random_list(1,3,6,7,9) and remainder_list(0,2,4,5,8)

        Input -  the factor which to split the object, if 0.1 then 10% of the
                 together object will be returned

        Output - two geospatial_data objects that are disjoint sets of the
                 original
        """

        i = 0
        self_size = len(self)
        random_list = []
        remainder_list = []
        new_size = round(factor * self_size)

        # Find unique random numbers
        if verbose:
            log.critical("make unique random number list "
                         "and get indices")

        total = num.array(list(range(self_size)), int)  # array default#
        total_list = total.tolist()

        if verbose:
            log.critical("total list len=%d" % len(total_list))

        # There will be repeated random numbers however will not be a
        # problem as they are being 'pop'ed out of array so if there
        # are two numbers the same they will pop different indicies,
        # still basically random
        # create list of non-unquie random numbers
        if verbose:
            log.critical("create random numbers list %s long"
                         % str(new_size))

        # Set seed if provided, mainly important for unit test!
        # plus recalcule seed when no seed provided.
        if seed_num is not None:
            seed(seed_num)
        else:
            seed()

        random_num = randint(0, self_size-1, (int(new_size),))
        random_num = random_num.tolist()

        # need to sort and reverse so the pop() works correctly
        random_num.sort()
        random_num.reverse()

        if verbose:
            log.critical("make random number list and get indices")

        j = 0
        k = 1
        remainder_list = total_list[:]

        # pops array index (random_num) from remainder_list
        # (which starts as the total_list and appends to random_list)
        random_num_len = len(random_num)
        for i in random_num:
            random_list.append(remainder_list.pop(i))
            j += 1
            # prints progress
            if verbose and round(random_num_len/10*k) == j:
                log.critical('(%s/%s)' % (j, random_num_len))
                k += 1

        # FIXME: move to tests, it might take a long time
        # then create an array of random length between 500 and 1000,
        # and use a random factor between 0 and 1
        # setup for assertion
        test_total = random_list[:]
        test_total.extend(remainder_list)
        test_total.sort()
        msg = ('The two random lists made from the original list when added '
               'together DO NOT equal the original list')
        assert total_list == test_total, msg

        # Get new samples
        if verbose:
            log.critical("get values of indices for random list")
        G1 = self.get_sample(random_list)
        if verbose:
            log.critical("get values of indices for "
                         "opposite of random list")
        G2 = self.get_sample(remainder_list)

        return G1, G2

    def __iter__(self):
        """Read in the header, number_of_points and save the
        file pointer position
        """

        # FIXME - what to do if the file isn't there

        # FIXME (Ole): Shouldn't this go into the constructor?
        # This method acts like the constructor when blocking.
        # ... and shouldn't it be called block_size?
        #
        if self.max_read_lines is None:
            self.max_read_lines = int(MAX_READ_LINES)

        if self.file_name[-4:] == ".pts":
            # See if the file is there.  Throw a QUIET IO error if it isn't
            fd = open(self.file_name, 'r')
            fd.close()

            # Throws prints to screen if file not present
            self.fid = NetCDFFile(self.file_name, netcdf_mode_r)

            (self.blocking_georef,
             self.blocking_keys,
             self.number_of_points) = _read_pts_file_header(self.fid,
                                                            self.verbose)
            self.start_row = 0
            self.last_row = self.number_of_points
            self.show_verbose = 0
            self.verbose_block_size = (self.last_row + 10)//10
            self.block_number = 0
            self.number_of_blocks = int(
                self.number_of_points/self.max_read_lines)
            # This computes the number of full blocks. The last block may be
            # smaller and won't be included in this estimate.

            if self.verbose is True:
                log.critical('Geospatial_data: Reading %d points (in %d block(s)) from file %s. '
                             % (self.number_of_points, self.number_of_blocks+1,
                                self.file_name))
                log.critical('Geospatial_data: Each block consists of %d data points'
                             % self.max_read_lines)
        else:
            # Assume the file is a csv file
            file_pointer = open(self.file_name)
            self.header, self.file_pointer = _read_csv_file_header(
                file_pointer)
            self.blocking_georef = None  # Used for reconciling zones

        return self

    def __next__(self):
        """read a block, instanciate a new geospatial and return it"""

        if self.file_name[-4:] == ".pts":
            if self.start_row == self.last_row:
                # Read the end of the file last iteration
                # Remove blocking attributes
                self.fid.close()
                del self.max_read_lines
                del self.blocking_georef
                del self.last_row
                del self.start_row
                del self.blocking_keys
                del self.fid
                raise StopIteration
            # Make sure we have an integer index
            fin_row = self.start_row + self.max_read_lines
            if fin_row > self.last_row:
                fin_row = self.last_row

            if self.verbose is True:
                if (self.show_verbose >= self.start_row
                        and self.show_verbose < fin_row):

                    log.critical('\nGeospatial_data: Reading block %d (points %d to %d) out of %d'
                                 % (self.block_number, self.start_row,
                                    fin_row, self.number_of_blocks))

                    self.show_verbose += max(self.max_read_lines,
                                             self.verbose_block_size)

            # Read next block
            pointlist, att_dict, = _read_pts_file_blocking(self.fid,
                                                           self.start_row,
                                                           fin_row,
                                                           self.blocking_keys)

            geo = Geospatial_data(pointlist, att_dict, self.blocking_georef)
            self.start_row = fin_row

            self.block_number += 1
        else:
            # Assume the file is a csv file
            try:
                (pointlist,
                 att_dict,
                 geo_ref,
                 self.file_pointer) = _read_csv_file_blocking(self.file_pointer,
                                                              self.header[:],
                                                              max_read_lines=self.max_read_lines,
                                                              verbose=self.verbose)

                # Check that the zones haven't changed.
                if geo_ref is not None:
                    geo_ref.reconcile_zones(self.blocking_georef)
                    self.blocking_georef = geo_ref
                elif self.blocking_georef is not None:
                    msg = ('Geo reference given, then not given.'
                           ' This should not happen.')
                    raise ValueError(msg)
                geo = Geospatial_data(pointlist, att_dict, geo_ref)
            except StopIteration:
                self.file_pointer.close()
                del self.header
                del self.file_pointer
                raise StopIteration
            except ANUGAError:
                self.file_pointer.close()
                del self.header
                del self.file_pointer
                raise
            except SyntaxError:
                self.file_pointer.close()
                del self.header
                del self.file_pointer
                # This should only be if there is a format error
                msg = ('Could not open file %s.\n%s'
                       % (self.file_name, Error_message['IOError']))
                raise SyntaxError(msg)
        return geo


##################### Error messages ###########
Error_message = {}
Em = Error_message
Em['IOError'] = ('NOTE: The format for a comma separated .txt/.csv file is:\n'
                 '        1st line:     [column names]\n'
                 '        other lines:  [x value], [y value], [attributes]\n'
                 '\n'
                 '           for example:\n'
                 '           x, y, elevation, friction\n'
                 '           0.6, 0.7, 4.9, 0.3\n'
                 '           1.9, 2.8, 5, 0.3\n'
                 '           2.7, 2.4, 5.2, 0.3\n'
                 '\n'
                 'The first two columns are assumed to be x, y coordinates.\n'
                 'The attribute values must be numeric.\n')


def _set_using_lat_long(latitudes,
                        longitudes,
                        geo_reference,
                        data_points,
                        points_are_lats_longs):
    """If the points has lat long info, assume it is in (lat, long) order."""

    if geo_reference is not None:
        msg = ('A georeference is specified yet latitude and longitude '
               'are also specified!')
        raise ValueError(msg)

    if data_points is not None and not points_are_lats_longs:
        msg = ('Data points are specified yet latitude and longitude are '
               'also specified.')
        raise ValueError(msg)

    if points_are_lats_longs:
        if data_points is None:
            msg = "Data points are not specified."
            raise ValueError(msg)
        lats_longs = ensure_numeric(data_points)
        latitudes = num.ravel(lats_longs[:, 0:1])
        longitudes = num.ravel(lats_longs[:, 1:])

    if latitudes is None and longitudes is None:
        msg = "Latitudes and Longitudes are not specified."
        raise ValueError(msg)

    if latitudes is None:
        msg = "Longitudes are specified yet latitudes aren't."
        raise ValueError(msg)

    if longitudes is None:
        msg = "Latitudes are specified yet longitudes aren't."
        raise ValueError(msg)

    data_points, zone = convert_from_latlon_to_utm(latitudes=latitudes,
                                                   longitudes=longitudes)
    return data_points, Geo_reference(zone=zone)


def _read_pts_file(file_name, verbose=False):
    """Read .pts NetCDF file

    Return a (dict_points, dict_attribute, geo_ref)
    eg
    dict['points'] = [[1.0,2.0],[3.0,5.0]]
    dict['attributelist']['elevation'] = [[7.0,5.0]]
    """

    if verbose:
        log.critical('Geospatial_data: Reading %s' % file_name)

    # See if the file is there.  Throw a QUIET IO error if it isn't
    fd = open(file_name, 'r')
    fd.close()

    # Throws prints to screen if file not present
    fid = NetCDFFile(file_name, netcdf_mode_r)

    pointlist = fid.variables['points'][:]
    keys = list(fid.variables.keys())

    if verbose:
        log.critical('Geospatial_data: Got %d variables: %s' %
                     (len(keys), keys))

    try:
        keys.remove('points')
    except IOError as e:
        fid.close()
        msg = "Expected keyword 'points' but could not find it"
        raise IOError(msg)

    attributes = {}
    for key in keys:
        if verbose:
            log.critical("Geospatial_data: Reading attribute '%s'" % key)

        if not (key == 'points'):
            attributes[key] = fid.variables[key][:]

    try:
        geo_reference = Geo_reference(NetCDFObject=fid)
    except AttributeError as e:
        geo_reference = None

    fid.close()

    if verbose:
        log.critical("Geospatial_data: %g data points" % len(pointlist))

    return pointlist, attributes, geo_reference


def _read_csv_file(file_name, verbose=False):
    """Read .csv file

    Return a dic of array of points, and dic of array of attribute
    eg
    dic['points'] = [[1.0,2.0],[3.0,5.0]]
    dic['attributelist']['elevation'] = [[7.0,5.0]]
    """

    file_pointer = open(file_name)
    header, file_pointer = _read_csv_file_header(file_pointer)
    try:
        (pointlist,
         att_dict,
         geo_ref,
         file_pointer) = _read_csv_file_blocking(file_pointer,
                                                 header,
                                                 max_read_lines=1e30)
        # If the file is bigger that this, block..
        # FIXME (Ole) What's up here?
    except ANUGAError:
        file_pointer.close()
        raise

    file_pointer.close()

    return pointlist, att_dict, geo_ref


CSV_DELIMITER = ','


def _read_csv_file_header(file_pointer,
                          delimiter=CSV_DELIMITER,
                          verbose=False):
    """Read the header of a .csv file
    Return a list of the header names
    """

    line = file_pointer.readline()
    header = clean_line(line, delimiter)

    return header, file_pointer


def _read_csv_file_blocking(file_pointer,
                            header,
                            delimiter=CSV_DELIMITER,
                            max_read_lines=MAX_READ_LINES,
                            verbose=False):
    """Read the body of a .csv file.
    header: The list header of the csv file, with the x and y labels.
    """

    points = []
    pointattributes = []
    att_dict = {}

    # This is to remove the x and y headers.
    header = header[:]
    try:
        x_header = header.pop(0)
        y_header = header.pop(0)
    except IndexError:
        # if there are not two columns this will occur.
        # eg if it is a space seperated file
        raise SyntaxError

    read_lines = 0
    while read_lines < max_read_lines:
        line = file_pointer.readline()
        numbers = clean_line(line, delimiter)
        if line == '\n':
            continue
        if len(numbers) <= 1:
            break
        if line[0] == '#':
            continue

        read_lines += 1

        try:
            x = float(numbers[0])
            y = float(numbers[1])
            points.append([x, y])
            numbers.pop(0)
            numbers.pop(0)
            if len(header) != len(numbers):
                file_pointer.close()
                msg = ('File load error. '
                       'There might be a problem with the file header.')
                raise SyntaxError(msg)
            for i, n in enumerate(numbers):
                n.strip()
                if n != '\n' and n != '':
                    att_dict.setdefault(header[i], []).append(float(n))
        except ValueError:
            raise SyntaxError

    if points == []:
        raise StopIteration

    pointlist = num.array(points, float)
    for key in list(att_dict.keys()):
        att_dict[key] = num.array(att_dict[key], float)

    # Do stuff here so the info is in lat's and longs
    geo_ref = None
    x_header = x_header[:3].lower()
    y_header = y_header[:3].lower()
    if (x_header == 'lon' or x_header == 'lat') \
       and (y_header == 'lon' or y_header == 'lat'):
        if x_header == 'lon':
            longitudes = num.ravel(pointlist[:, 0:1])
            latitudes = num.ravel(pointlist[:, 1:])
        else:
            latitudes = num.ravel(pointlist[:, 0:1])
            longitudes = num.ravel(pointlist[:, 1:])

        pointlist, geo_ref = _set_using_lat_long(latitudes,
                                                 longitudes,
                                                 geo_reference=None,
                                                 data_points=None,
                                                 points_are_lats_longs=False)

    return pointlist, att_dict, geo_ref, file_pointer


def _read_pts_file_header(fid, verbose=False):
    '''Read the geo_reference and number_of_points from a .pts file'''

    keys = list(fid.variables.keys())
    try:
        keys.remove('points')
    except IOError as e:
        fid.close()
        msg = "Expected keyword 'points' but could not find it."
        raise IOError(msg)

    if verbose:
        log.critical('Got %d variables: %s' % (len(keys), keys))

    try:
        geo_reference = Geo_reference(NetCDFObject=fid)
    except AttributeError as e:
        geo_reference = None

    try:  # netcdf4
        number_of_points = len(fid.dimensions['number_of_points'])
    except:  # scientific python
        number_of_points = fid.dimensions['number_of_points']

    return geo_reference, keys, number_of_points


def _read_pts_file_blocking(fid, start_row, fin_row, keys):
    '''Read the body of a .pts file.'''

    pointlist = num.array(fid.variables['points'][start_row:fin_row])

    attributes = {}
    for key in keys:
        attributes[key] = num.array(fid.variables[key][start_row:fin_row])

    return pointlist, attributes


def _write_pts_file(file_name,
                    write_data_points,
                    write_attributes=None,
                    write_geo_reference=None):
    """Write .pts NetCDF file

    NOTE: Below might not be valid ask Duncan : NB 5/2006

    WARNING: This function mangles the point_atts data structure
    # F??ME: (DSG)This format has issues.
    # There can't be an attribute called points
    # consider format change
    # method changed by NB not sure if above statement is correct

    should create new test for this
    legal_keys = ['pointlist', 'attributelist', 'geo_reference']
    for key in point_atts.keys():
        msg = 'Key %s is illegal. Valid keys are %s' %(key, legal_keys)
        assert key in legal_keys, msg
    """

    # NetCDF file definition
    outfile = NetCDFFile(file_name, netcdf_mode_w)

    # Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = ('NetCDF format for compact and portable storage '
                           'of spatial point data')

    # Dimension definitions
    shape = write_data_points.shape[0]
    outfile.createDimension('number_of_points', shape)
    outfile.createDimension('number_of_dimensions', 2)  # This is 2d data

    # Variable definition
    outfile.createVariable('points', netcdf_float,
                           ('number_of_points', 'number_of_dimensions'))

    # create variables
    outfile.variables['points'][:] = write_data_points

    if write_attributes is not None:
        for key in list(write_attributes.keys()):
            outfile.createVariable(key, netcdf_float, ('number_of_points',))
            outfile.variables[key][:] = write_attributes[key]

    if write_geo_reference is not None:
        write_NetCDF_georeference(write_geo_reference, outfile)

    outfile.close()


def _write_csv_file(file_name,
                    write_data_points,
                    write_attributes=None,
                    as_lat_long=False,
                    delimiter=','):
    """Write a .csv file."""

    points = write_data_points
    pointattributes = write_attributes

    fd = open(file_name, 'w')

    if as_lat_long:
        titlelist = "latitude" + delimiter + "longitude" + delimiter
    else:
        titlelist = "x" + delimiter + "y" + delimiter

    if pointattributes is not None:
        for title in list(pointattributes.keys()):
            titlelist = titlelist + title + delimiter
        titlelist = titlelist[0:-len(delimiter)]  # remove the last delimiter

    fd.write(titlelist + "\n")

    # <x/lat> <y/long> [attributes]
    for i, vert in enumerate(points):
        if pointattributes is not None:
            attlist = ","
            for att in list(pointattributes.keys()):
                attlist = attlist + str(pointattributes[att][i]) + delimiter
            attlist = attlist[0:-len(delimiter)]  # remove the last delimiter
            attlist.strip()
        else:
            attlist = ''

        fd.write(str(vert[0]) + delimiter + str(vert[1]) + attlist + "\n")

    fd.close()


def _write_urs_file(file_name, points, delimiter=' '):
    """Write a URS format file.
    export a file, file_name, with the urs format
    the data points are in lats and longs
    """

    fd = open(file_name, 'w')

    # first line is # points
    fd.write(str(len(points)) + "\n")

    # <lat> <long> <id#>
    for i, vert in enumerate(points):
        fd.write(str(round(vert[0], 7)) + delimiter +
                 str(round(vert[1], 7)) + delimiter + str(i) + "\n")

    fd.close()


def _point_atts2array(point_atts):
    point_atts['pointlist'] = num.array(point_atts['pointlist'], float)

    for key in list(point_atts['attributelist'].keys()):
        point_atts['attributelist'][key] = \
            num.array(point_atts['attributelist'][key], float)

    return point_atts


def geospatial_data2points_dictionary(geospatial_data):
    """Convert geospatial data to points_dictionary"""

    points_dictionary = {}
    points_dictionary['pointlist'] = geospatial_data.data_points

    points_dictionary['attributelist'] = {}

    for attribute_name in list(geospatial_data.attributes.keys()):
        val = geospatial_data.attributes[attribute_name]
        points_dictionary['attributelist'][attribute_name] = val

    points_dictionary['geo_reference'] = geospatial_data.geo_reference

    return points_dictionary


def points_dictionary2geospatial_data(points_dictionary):
    """Convert points_dictionary to geospatial data object"""

    msg = "Points dictionary must have key 'pointlist'"
    assert 'pointlist' in points_dictionary, msg

    msg = "Points dictionary must have key 'attributelist'"
    assert 'attributelist' in points_dictionary, msg

    if 'geo_reference' in points_dictionary:
        geo = points_dictionary['geo_reference']
    else:
        geo = None

    return Geospatial_data(points_dictionary['pointlist'],
                           points_dictionary['attributelist'],
                           geo_reference=geo)


def ensure_absolute(points, geo_reference=None):
    """Ensure that points are in absolute coordinates.

    This function inputs several formats and
    outputs one format. - a numeric array of absolute points.

    Input formats are;
      points: List or numeric array of coordinate pairs [xi, eta] of
              points or geospatial object or points file name

    mesh_origin: A geo_reference object or 3-tuples consisting of
                 UTM zone, easting and northing.
                 If specified vertex coordinates are assumed to be
                 relative to their respective origins.
    """

    # Input check
    if isinstance(points, str):
        # It's a string - assume it is a point file
        points = Geospatial_data(file_name=points)

    if isinstance(points, Geospatial_data):
        points = points.get_data_points(absolute=True)
        msg = 'Use a Geospatial_data object or a mesh origin, not both.'
        assert geo_reference is None, msg
    else:
        points = ensure_numeric(copy.copy(points), float)

    # Sort of geo_reference and convert points
    if geo_reference is None:
        geo = None    # Geo_reference()
    else:
        if isinstance(geo_reference, Geo_reference):
            geo = geo_reference
        else:
            geo = Geo_reference(geo_reference[0],
                                geo_reference[1],
                                geo_reference[2])
        points = geo.get_absolute(points)

    return points


def ensure_geospatial(points, geo_reference=None):
    """Convert various data formats to a geospatial_data instance.

    Inputed formats are;
    points:      List or numeric array of coordinate pairs [xi, eta] of
                 points or geospatial object

    mesh_origin: A geo_reference object or 3-tuples consisting of
                 UTM zone, easting and northing.
                 If specified vertex coordinates are assumed to be
                 relative to their respective origins.
    """

    # Input check
    if isinstance(points, Geospatial_data):
        msg = "Use a Geospatial_data object or a mesh origin, not both."
        assert geo_reference is None, msg
        return points
    else:
        # List or numeric array of absolute points
        points = ensure_numeric(points, float)

    # Sort out geo reference
    if geo_reference is None:
        geo = None
    else:
        if isinstance(geo_reference, Geo_reference):
            geo = geo_reference
        else:
            geo = Geo_reference(geo_reference[0],
                                geo_reference[1],
                                geo_reference[2])

    # Create Geospatial_data object with appropriate geo reference and return
    points = Geospatial_data(data_points=points, geo_reference=geo)

    return points


def find_optimal_smoothing_parameter(data_file,
                                     alpha_list=None,
                                     mesh_file=None,
                                     boundary_poly=None,
                                     mesh_resolution=100000,
                                     north_boundary=None,
                                     south_boundary=None,
                                     east_boundary=None,
                                     west_boundary=None,
                                     plot_name='all_alphas',
                                     split_factor=0.1,
                                     seed_num=None,
                                     cache=False,
                                     verbose=False):
    """Removes a small random sample of points from 'data_file'.
    Then creates models with different alpha values from 'alpha_list' and
    cross validates the predicted value to the previously removed point data.
    Returns the alpha value which has the smallest covariance.

    data_file: must not contain points outside the boundaries defined
               and it must be either a pts, txt or csv file.

    alpha_list: the alpha values to test in a single list

    mesh_file: name of the created mesh file or if passed in will read it.
               NOTE, if there is a mesh file mesh_resolution,
               north_boundary, south... etc will be ignored.

    mesh_resolution: the maximum area size for a triangle

    north_boundary... west_boundary: the value of the boundary

    plot_name: the name for the plot contain the results

    seed_num: the seed to the random number generator

    USAGE:
        value, alpha = find_optimal_smoothing_parameter(data_file=fileName,
                                             alpha_list=[0.0001, 0.01, 1],
                                             mesh_file=None,
                                             mesh_resolution=3,
                                             north_boundary=5,
                                             south_boundary=-5,
                                             east_boundary=5,
                                             west_boundary=-5,
                                             plot_name='all_alphas',
                                             seed_num=100000,
                                             verbose=False)

    OUTPUT: returns the minumum normalised covalance calculate AND the
           alpha that created it. PLUS writes a plot of the results

    NOTE: code will not work if the data_file extent is greater than the
    boundary_polygon or any of the boundaries, eg north_boundary...west_boundary
    """

    from anuga.shallow_water.shallow_water_domain import Domain
    from anuga.geospatial_data.geospatial_data import Geospatial_data
    from anuga.pmesh.mesh_interface import create_mesh_from_regions
    from anuga.utilities.numerical_tools import cov
    from anuga.geometry.polygon import is_inside_polygon
    from anuga.fit_interpolate.benchmark_least_squares import mem_usage

    attribute_smoothed = 'elevation'

    if mesh_file is None:
        if verbose:
            log.critical("building mesh")
        mesh_file = 'temp.msh'

        if (north_boundary is None or south_boundary is None
                or east_boundary is None or west_boundary is None):
            no_boundary = True
        else:
            no_boundary = False

        if no_boundary is True:
            msg = 'All boundaries must be defined'
            raise Exception(msg)

        poly_topo = [[east_boundary, south_boundary],
                     [east_boundary, north_boundary],
                     [west_boundary, north_boundary],
                     [west_boundary, south_boundary]]

        create_mesh_from_regions(poly_topo,
                                 boundary_tags={'back': [2],
                                                'side': [1, 3],
                                                'ocean': [0]},
                                 maximum_triangle_area=mesh_resolution,
                                 filename=mesh_file,
                                 use_cache=cache,
                                 verbose=verbose)

    else:  # if mesh file provided
        # test mesh file exists?
        if verbose:
            "reading from file: %s" % mesh_file
        if access(mesh_file, F_OK) == 0:
            msg = "file %s doesn't exist!" % mesh_file
            raise IOError(msg)

    # split topo data
    if verbose:
        log.critical('Reading elevation file: %s' % data_file)
    G = Geospatial_data(file_name=data_file)
    if verbose:
        log.critical('Start split')
    G_small, G_other = G.split(split_factor, seed_num, verbose=verbose)
    if verbose:
        log.critical('Finish split')
    points = G_small.get_data_points()

    if verbose:
        log.critical("Number of points in sample to compare: %d"
                     % len(points))

    if alpha_list is None:
        alphas = [0.001, 0.01, 100]
        # alphas = [0.000001, 0.00001, 0.0001, 0.001, 0.01,
        #          0.1, 1.0, 10.0, 100.0,1000.0,10000.0]
    else:
        alphas = alpha_list

    # creates array with columns 1 and 2 are x, y. column 3 is elevation
    # 4 onwards is the elevation_predicted using the alpha, which will
    # be compared later against the real removed data
    data = num.array([], dtype=float)

    data = num.resize(data, (len(points), 3+len(alphas)))

    # gets relative point from sample
    data[:, 0] = points[:, 0]
    data[:, 1] = points[:, 1]
    elevation_sample = G_small.get_attributes(
        attribute_name=attribute_smoothed)
    data[:, 2] = elevation_sample

    normal_cov = num.array(num.zeros([len(alphas), 2]), dtype=float)

    if verbose:
        log.critical('Setup computational domains with '
                     'different alphas')

    for i, alpha in enumerate(alphas):
        # add G_other data to domains with different alphas
        if verbose:
            log.critical('Calculating domain and mesh for Alpha=%s'
                         % str(alpha))
        domain = Domain(mesh_file, use_cache=cache, verbose=verbose)
        if verbose:
            log.critical(domain.statistics())
        domain.set_quantity(attribute_smoothed,
                            geospatial_data=G_other,
                            use_cache=cache,
                            verbose=verbose,
                            alpha=alpha)

        # Convert points to geospatial data for use with get_values below
        points_geo = Geospatial_data(points, domain.geo_reference)

        # returns the predicted elevation of the points that were "split" out
        # of the original data set for one particular alpha
        if verbose:
            log.critical('Get predicted elevation for location '
                         'to be compared')
        elevation_predicted = \
            domain.quantities[attribute_smoothed].\
            get_values(interpolation_points=points_geo)

        # add predicted elevation to array that starts with x, y, z...
        data[:, i+3] = elevation_predicted

        sample_cov = cov(elevation_sample)
        ele_cov = cov(elevation_sample - elevation_predicted)
        normal_cov[i, :] = [alpha, ele_cov/sample_cov]

        if verbose:
            log.critical('Covariance for alpha %s=%s'
                         % (normal_cov[i][0], normal_cov[i][1]))
            log.critical('--------------------------------------------')

    normal_cov0 = normal_cov[:, 0]
    normal_cov_new = num.take(normal_cov, num.argsort(normal_cov0), axis=0)

    if plot_name is not None:
        from pylab import savefig, semilogx, loglog

        semilogx(normal_cov_new[:, 0], normal_cov_new[:, 1])
        loglog(normal_cov_new[:, 0], normal_cov_new[:, 1])
        savefig(plot_name, dpi=300)

    if mesh_file == 'temp.msh':
        remove(mesh_file)

    if verbose:
        log.critical('Final results:')
        for i, alpha in enumerate(alphas):
            log.critical('covariance for alpha %s = %s '
                         % (normal_cov[i][0], normal_cov[i][1]))
        log.critical('Optimal alpha is: %s '
                     % normal_cov_new[(num.argmin(normal_cov_new, axis=0))[1], 0])

    # covariance and optimal alpha
    return (min(normal_cov_new[:, 1]),
            normal_cov_new[(num.argmin(normal_cov_new, axis=0))[1], 0])


def old_find_optimal_smoothing_parameter(data_file,
                                         alpha_list=None,
                                         mesh_file=None,
                                         boundary_poly=None,
                                         mesh_resolution=100000,
                                         north_boundary=None,
                                         south_boundary=None,
                                         east_boundary=None,
                                         west_boundary=None,
                                         plot_name='all_alphas',
                                         split_factor=0.1,
                                         seed_num=None,
                                         cache=False,
                                         verbose=False):
    """
    data_file: must not contain points outside the boundaries defined
               and it either a pts, txt or csv file.

    alpha_list: the alpha values to test in a single list

    mesh_file: name of the created mesh file or if passed in will read it.
               NOTE, if there is a mesh file mesh_resolution,
               north_boundary, south... etc will be ignored.

    mesh_resolution: the maximum area size for a triangle

    north_boundary... west_boundary: the value of the boundary

    plot_name: the name for the plot contain the results

    seed_num: the seed to the random number generator

    USAGE:
        value, alpha = find_optimal_smoothing_parameter(data_file=fileName,
                                             alpha_list=[0.0001, 0.01, 1],
                                             mesh_file=None,
                                             mesh_resolution=3,
                                             north_boundary=5,
                                             south_boundary=-5,
                                             east_boundary=5,
                                             west_boundary=-5,
                                             plot_name='all_alphas',
                                             seed_num=100000,
                                             verbose=False)

    OUTPUT: returns the minumum normalised covalance calculate AND the
            alpha that created it. PLUS writes a plot of the results

    NOTE: code will not work if the data_file extend is greater than the
          boundary_polygon or the north_boundary...west_boundary
    """

    from anuga.shallow_water.shallow_water_domain import Domain
    from anuga.geospatial_data.geospatial_data import Geospatial_data
    from anuga.pmesh.mesh_interface import create_mesh_from_regions
    from anuga.utilities.numerical_tools import cov
    from anuga.geometry.polygon import is_inside_polygon
    from anuga.fit_interpolate.benchmark_least_squares import mem_usage

    attribute_smoothed = 'elevation'

    if mesh_file is None:
        mesh_file = 'temp.msh'

        if (north_boundary is None or south_boundary is None
                or east_boundary is None or west_boundary is None):
            no_boundary = True
        else:
            no_boundary = False

        if no_boundary is True:
            msg = 'All boundaries must be defined'
            raise Exception(msg)

        poly_topo = [[east_boundary, south_boundary],
                     [east_boundary, north_boundary],
                     [west_boundary, north_boundary],
                     [west_boundary, south_boundary]]

        create_mesh_from_regions(poly_topo,
                                 boundary_tags={'back': [2],
                                                'side': [1, 3],
                                                'ocean': [0]},
                                 maximum_triangle_area=mesh_resolution,
                                 filename=mesh_file,
                                 use_cache=cache,
                                 verbose=verbose)

    else:  # if mesh file provided
        # test mesh file exists?
        if access(mesh_file, F_OK) == 0:
            msg = "file %s doesn't exist!" % mesh_file
            raise IOError(msg)

    # split topo data
    G = Geospatial_data(file_name=data_file)
    if verbose:
        log.critical('start split')
    G_small, G_other = G.split(split_factor, seed_num, verbose=verbose)
    if verbose:
        log.critical('finish split')
    points = G_small.get_data_points()

    if verbose:
        log.critical("Number of points in sample to compare: %d"
                     % len(points))

    if alpha_list is None:
        alphas = [0.001, 0.01, 100]
        # alphas = [0.000001, 0.00001, 0.0001, 0.001, 0.01,
        #          0.1, 1.0, 10.0, 100.0,1000.0,10000.0]
    else:
        alphas = alpha_list

    domains = {}

    if verbose:
        log.critical('Setup computational domains with '
                     'different alphas')

    for alpha in alphas:
        # add G_other data to domains with different alphas
        if verbose:
            log.critical('Calculating domain and mesh for Alpha = %s' %
                         str(alpha))
        domain = Domain(mesh_file, use_cache=cache, verbose=verbose)
        if verbose:
            log.critical(domain.statistics())
        domain.set_quantity(attribute_smoothed,
                            geospatial_data=G_other,
                            use_cache=cache,
                            verbose=verbose,
                            alpha=alpha)
        domains[alpha] = domain

    # creates array with columns 1 and 2 are x, y. column 3 is elevation
    # 4 onwards is the elevation_predicted using the alpha, which will
    # be compared later against the real removed data
    data = num.array([], dtype=float)

    data = num.resize(data, (len(points), 3+len(alphas)))

    # gets relative point from sample
    data[:, 0] = points[:, 0]
    data[:, 1] = points[:, 1]
    elevation_sample = G_small.get_attributes(
        attribute_name=attribute_smoothed)
    data[:, 2] = elevation_sample

    normal_cov = num.array(num.zeros([len(alphas), 2]), dtype=float)

    if verbose:
        log.critical('Determine difference between predicted results '
                     'and actual data')

    for i, alpha in enumerate(domains):
        if verbose:
            print('Alpha =', alpha)

        points_geo = domains[alpha].geo_reference.change_points_geo_ref(points)
        # returns the predicted elevation of the points that were "split" out
        # of the original data set for one particular alpha
        elevation_predicted = \
            domains[alpha].quantities[attribute_smoothed].\
            get_values(interpolation_points=points_geo)

        # add predicted elevation to array that starts with x, y, z...
        data[:, i+3] = elevation_predicted

        sample_cov = cov(elevation_sample)
        ele_cov = cov(elevation_sample - elevation_predicted)
        normal_cov[i, :] = [alpha, ele_cov/sample_cov]
        log.critical('memory usage during compare: %s' % str(mem_usage()))
        if verbose:
            log.critical('cov %s = %s'
                         % (normal_cov[i][0], normal_cov[i][1]))

    normal_cov0 = normal_cov[:, 0]
    normal_cov_new = num.take(normal_cov, num.argsort(normal_cov0), axis=0)

    if plot_name is not None:
        from pylab import savefig, semilogx, loglog

        semilogx(normal_cov_new[:, 0], normal_cov_new[:, 1])
        loglog(normal_cov_new[:, 0], normal_cov_new[:, 1])
        savefig(plot_name, dpi=300)
    if mesh_file == 'temp.msh':
        remove(mesh_file)

    return (min(normal_cov_new[:, 1]),
            normal_cov_new[(num.argmin(normal_cov_new, axis=0))[1], 0])


if __name__ == "__main__":
    pass
