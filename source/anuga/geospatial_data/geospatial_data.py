"""Class Geospatial_data - Manipulation of locations on the planet and 
associated attributes.

"""
from sys import maxint
from os import access, F_OK, R_OK
from types import DictType
from warnings import warn
from string import lower
from Numeric import concatenate, array, Float, shape, reshape, ravel, take, \
                        size, shape
from random import randint
from copy import deepcopy

#from MA import tolist

from Scientific.IO.NetCDF import NetCDFFile    
from anuga.coordinate_transforms.lat_long_UTM_conversion import UTMtoLL    
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.coordinate_transforms.geo_reference import Geo_reference, \
     TitleError, DEFAULT_ZONE, ensure_geo_reference, write_NetCDF_georeference
from anuga.coordinate_transforms.redfearn import convert_from_latlon_to_utm
from anuga.utilities.anuga_exceptions import ANUGAError
from anuga.config import points_file_block_line_size as MAX_READ_LINES    

class Geospatial_data:

    def __init__(self,
                 data_points=None, # this can also be a points file name
                 attributes=None,
                 geo_reference=None,
                 default_attribute_name=None,
                 file_name=None,
                 delimiter=None,
                 latitudes=None,
                 longitudes=None,
                 points_are_lats_longs=False,
                 max_read_lines=None,
                 load_file_now=True,
                 verbose=False):

        
        """
        Create instance from data points and associated attributes

        data_points: x,y coordinates in meters. Type must be either a
        sequence of 2-tuples or an Mx2 Numeric array of floats.  A file name
        can also be passed in here.

        attributes: Associated values for each data point. The type
        must be either a list or an array of length M or a dictionary
        of lists (or arrays) of length M. In the latter case the keys
        in the dictionary represent the attribute names, in the former
        the attribute will get the default name "attribute".
        
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
        
        file_name: Name of input netCDF file or .txt file. netCDF file must 
        have dimensions "points" etc.
        .txt file is a comma seperated file with x, y and attribute
        data.
        
        The first line has the titles of the columns.  The first two
        column titles are checked to see if they start with lat or
        long (not case sensitive).  If so the data is assumed to be
        latitude and longitude, in decimal format and converted to
        UTM.  Otherwise the first two columns are assumed to be the x
        and y, and the title names acually used are ignored.
 
        
        The format for a .txt file is:
            1st line:     [column names]
            other lines:  x y [attributes]

            for example:
            x, y, elevation, friction
            0.6, 0.7, 4.9, 0.3
            1.9, 2.8, 5, 0.3
            2.7, 2.4, 5.2, 0.3

        The first two columns are always  assumed to be x, y
        coordinates.
     
        An issue with the xya format is that the attribute column order
        is not be controlled.  The info is stored in a dictionary and it's
        written in an order dependent on the hash order
        
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
                
        delimiter: is the file delimiter that will be used when 
            importing a .xya file, which is being phased out.
            
        verbose:

          
        """

        if isinstance(data_points, basestring):
            # assume data point is really a file name
            file_name = data_points

        self.set_verbose(verbose)
        self.geo_reference=None #create the attribute 
        self.file_name = file_name
        self.max_read_lines = max_read_lines

        if delimiter is not None:
            msg = 'Specifying delimiters will be removed.'
            msg = 'Text file format is moving to comma seperated .txt files.'
            warn(msg, DeprecationWarning) 
        if file_name is None:
            if delimiter is not None:
                msg = 'No file specified yet a delimiter is provided!'
                raise ValueError, msg
            file_name = None
            if latitudes is not None or longitudes is not None or \
                   points_are_lats_longs:
                data_points, geo_reference =  \
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
            self.import_points_file(file_name, delimiter, verbose)
                
            self.check_data_points(self.data_points)
            self.set_attributes(self.attributes) 
            self.set_geo_reference(self.geo_reference)
            self.set_default_attribute_name(default_attribute_name)

        #Why?    
        #assert self.attributes is None or isinstance(self.attributes, DictType)
        #This is a hassle when blocking, so I've removed it.
        

    def __len__(self):
        return len(self.data_points)

    def __repr__(self):
        return str(self.get_data_points(absolute=True))
    
    
    def check_data_points(self, data_points):
        """Checks data points
        """
        
        if data_points is None:
            self.data_points = None
            msg = 'There is no data or file provided!'
            raise ValueError, msg
            
        else:
            self.data_points = ensure_numeric(data_points)
            #print "self.data_points.shape",self.data_points.shape
            if not (0,) == self.data_points.shape:
                assert len(self.data_points.shape) == 2
                assert self.data_points.shape[1] == 2

    def set_attributes(self, attributes):
        """Check and assign attributes dictionary
        """
        
        if attributes is None:
            self.attributes = None
            return
        
        if not isinstance(attributes, DictType):
            #Convert single attribute into dictionary
            attributes = {'attribute': attributes}

        #Check input attributes    
        for key in attributes.keys():
            try:
                attributes[key] = ensure_numeric(attributes[key])
            except:
                msg = 'Attribute %s could not be converted' %key
                msg += 'to a numeric vector'
                raise msg

        self.attributes = attributes    


    def set_geo_reference(self, geo_reference):
        """
        Set's the georeference of geospatial.
        It can also be used to change the georeference
        """
        from anuga.coordinate_transforms.geo_reference import Geo_reference

        if geo_reference is None:
            geo_reference = Geo_reference() # Use default
        geo_reference = ensure_geo_reference(geo_reference)
        if not isinstance(geo_reference, Geo_reference):
            msg = 'Argument geo_reference must be a valid Geo_reference \n'
            msg += 'object or None.'
            raise msg

        # if a geo ref already exists, change the point data to
        # represent the new geo-ref
        if  self.geo_reference is not None:
            #FIXME: Maybe put out a warning here...
            self.data_points = self.get_data_points \
                               (geo_reference=geo_reference)
            
        self.geo_reference = geo_reference


    def set_default_attribute_name(self, default_attribute_name):
        self.default_attribute_name = default_attribute_name

    def set_verbose(self, verbose=False):
        if verbose in [False, True]:
            self.verbose = verbose
        else:
            msg = 'Illegal value: %s' %str(verbose)
            raise Exception, msg

    def clip(self, polygon, closed=True):
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
        
        """

        from anuga.utilities.polygon import inside_polygon

        if isinstance(polygon, Geospatial_data):
            # Polygon is an object - extract points
            polygon = polygon.get_data_points()

        points = self.get_data_points()    
        inside_indices = inside_polygon(points, polygon, closed)

        clipped_G = self.get_sample(inside_indices)
#        clipped_points = take(points, inside_indices)

        # Clip all attributes
#        attributes = self.get_all_attributes()

#        clipped_attributes = {}
#        if attributes is not None:
#            for key, att in attributes.items():
#                clipped_attributes[key] = take(att, inside_indices)

#        return Geospatial_data(clipped_points, clipped_attributes)
        return clipped_G
        
        
    def clip_outside(self, polygon, closed=True):
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

        from anuga.utilities.polygon import outside_polygon

        if isinstance(polygon, Geospatial_data):
            # Polygon is an object - extract points
            polygon = polygon.get_data_points()

        points = self.get_data_points()    
        outside_indices = outside_polygon(points, polygon, closed)

        clipped_G = self.get_sample(outside_indices)

#        clipped_points = take(points, outside_indices)

        # Clip all attributes
#        attributes = self.get_all_attributes()

#        clipped_attributes = {}
#        if attributes is not None:
#            for key, att in attributes.items():
#                clipped_attributes[key] = take(att, outside_indices)

#        return Geospatial_data(clipped_points, clipped_attributes)
        return clipped_G

    
    def get_geo_reference(self):
        return self.geo_reference
       
    def get_data_points(self, absolute=True, geo_reference=None,
                        as_lat_long=False):
        """Get coordinates for all data points as an Nx2 array

        If absolute is False returned coordinates are relative to the
        internal georeference's xll and yll corners, otherwise
        absolute UTM coordinates are returned.

        If a geo_reference is passed the points are returned relative
        to that geo_reference.

        Default: absolute is True.
        """
        if as_lat_long is True:
            msg = "Points need a zone to be converted into lats and longs"
            assert self.geo_reference is not None, msg
            zone = self.geo_reference.get_zone()
            assert self.geo_reference.get_zone() is not DEFAULT_ZONE, msg
            lats_longs = []
            for point in self.get_data_points(True):
                ### UTMtoLL(northing, easting, zone,
                lat_calced, long_calced = UTMtoLL(point[1],point[0], zone)
                lats_longs.append((lat_calced, long_calced)) # to hash
            return lats_longs
            
        if absolute is True and geo_reference is None:
            return self.geo_reference.get_absolute(self.data_points)
        elif geo_reference is not None:
            return geo_reference.change_points_geo_ref \
                               (self.data_points, 
                                self.geo_reference)
        else:
            return self.data_points
        
    
    def get_attributes(self, attribute_name=None):
        """Return values for one named attribute.

        If attribute_name is None, default_attribute_name is used
        """

        if attribute_name is None:
            if self.default_attribute_name is not None:
                attribute_name = self.default_attribute_name
            else:
                attribute_name = self.attributes.keys()[0] 
                # above line takes the first one from keys
        
        if self.verbose is True:
            print 'Using attribute %s' %attribute_name
            print 'Available attributes: %s' %(self.attributes.keys())        

        msg = 'Attribute name %s does not exist in data set' %attribute_name
        assert self.attributes.has_key(attribute_name), msg

        return self.attributes[attribute_name]

    def get_all_attributes(self):
        """
        Return values for all attributes.
        The return value is either None or a dictionary (possibly empty).
        """

        return self.attributes

    def __add__(self, other):
        """
        Returns the addition of 2 geospatical objects,
        objects are concatencated to the end of each other
            
        NOTE: doesn't add if objects contain different 
        attributes  
        
        Always return relative points!
        """

        # find objects zone and checks if the same
        geo_ref1 = self.get_geo_reference()
        zone1 = geo_ref1.get_zone()
        
        geo_ref2 = other.get_geo_reference()
        zone2 = geo_ref2.get_zone()

        geo_ref1.reconcile_zones(geo_ref2)


        # sets xll and yll as the smallest from self and other
        # FIXME (Duncan and Ole): use lower left corner derived from
        # absolute coordinates
        if self.geo_reference.xllcorner <= other.geo_reference.xllcorner:
            xll = self.geo_reference.xllcorner
        else:
            xll = other.geo_reference.xllcorner

        if self.geo_reference.yllcorner <= other.geo_reference.yllcorner:
            yll = self.geo_reference.yllcorner
        else:
            yll = other.geo_reference.yllcorner
        new_geo_ref = Geo_reference(geo_ref1.get_zone(), xll, yll)

        xll = yll = 0. 
        
        relative_points1 = self.get_data_points(absolute = False)
        relative_points2 = other.get_data_points(absolute = False)

        
        new_relative_points1 = new_geo_ref.\
                               change_points_geo_ref(relative_points1,
                                                     geo_ref1)
        new_relative_points2 = new_geo_ref.\
                               change_points_geo_ref(relative_points2,
                                                     geo_ref2)
        
        # Now both point sets are relative to new_geo_ref and
        # zones have been reconciled

        # Concatenate points
        new_points = concatenate((new_relative_points1,
                                  new_relative_points2),
                                  axis = 0)
      
        # Concatenate attributes if any
        if self.attributes is None:
            if other.attributes is not None:
                msg = 'Both geospatial_data objects must have the same \n'
                msg += 'attributes to allow addition.'
                raise Exception, msg
            
            new_attributes = None
        else:    
            new_attributes = {}
            for x in self.attributes.keys():
                if other.attributes.has_key(x):

                    attrib1 = self.attributes[x]
                    attrib2 = other.attributes[x]
                    new_attributes[x] = concatenate((attrib1, attrib2))

                else:
                    msg = 'Both geospatial_data objects must have the same \n'
                    msg += 'attributes to allow addition.'
                    raise Exception, msg

        # Instantiate new data object and return    
        return Geospatial_data(new_points,
                               new_attributes,
                               new_geo_ref)
    
    ###
    #  IMPORT/EXPORT POINTS FILES
    ###

    def import_points_file(self, file_name, delimiter=None, verbose=False):
        """ load an .xya or .pts file
        Note: will throw an IOError if it can't load the file.
        Catch these!

        Post condition: self.attributes dictionary has been set
        """
        
        if access(file_name, F_OK) == 0 :
            msg = 'File %s does not exist or is not accessible' %file_name
            raise IOError, msg
        
        attributes = {}
        if file_name[-4:]== ".xya":
            msg = 'Text file format is moving to comma seperated .txt files.'
            warn(msg, DeprecationWarning) 
            try:
                if delimiter == None:
                    try:
                        fd = open(file_name)
                        data_points, attributes, geo_reference =\
                                     _read_xya_file(fd, ',')
                    except TitleError:
                        fd.close()
                        fd = open(file_name)
                        data_points, attributes, geo_reference =\
                                     _read_xya_file(fd, ' ')
                else:
                    fd = open(file_name)
                    data_points, attributes, geo_reference =\
                                 _read_xya_file(fd, delimiter)
                fd.close()
            except (IndexError,ValueError,SyntaxError):
                fd.close()    
                msg = 'Could not open file %s ' %file_name
                raise IOError, msg
            except IOError, e:
                fd.close()  
                # Catch this to add an error message
                msg = 'Could not open file or incorrect file format %s:%s'\
                      %(file_name, e)
                raise IOError, msg
                
        elif file_name[-4:]== ".pts":
            try:
                data_points, attributes, geo_reference =\
                             _read_pts_file(file_name, verbose)
            except IOError, e:    
                msg = 'Could not open file %s ' %file_name
                raise IOError, msg  
        
        elif file_name[-4:]== ".txt" or file_name[-4:]== ".csv":
            #let's do ticket#116 stuff
            #
            try:
                data_points, attributes, geo_reference =\
                             _read_csv_file(file_name, verbose)
            except IOError, e:    
                msg = 'Could not open file %s ' %file_name
                raise IOError, msg        
        else:      
            msg = 'Extension %s is unknown' %file_name[-4:]
            raise IOError, msg
#        print'in import data_points', data_points
#        print'in import attributes', attributes
#        print'in import data_points', geo_reference
        self.data_points = data_points
        self.attributes = attributes
        self.geo_reference = geo_reference
    
#        return all_data
    
    def export_points_file(self, file_name, absolute=True, as_lat_long=False):
        
        """
        write a points file, file_name, as a text (.xya) or binary (.pts) file
        file_name is the file name, including the extension
        The point_dict is defined at the top of this file.
        
        If absolute is True data the xll and yll are added to the points value 
        and the xll and yll of the geo_reference are set to 0.
        
        If absolute is False data points at returned as relative to the xll 
        and yll and geo_reference remains uneffected
        """

        if absolute is False and file_name[-4:] == ".xya":
            msg = 'The text file values must be absolute.   '
            msg += 'Text file format is moving to comma seperated .txt files.'
            warn(msg, DeprecationWarning) 

        if (file_name[-4:] == ".xya"):
            msg = '.xya format is deprecated.  Please use .txt.'
            warn(msg, DeprecationWarning)
            if absolute is True:     
                geo_ref = deepcopy(self.geo_reference)
                geo_ref.xllcorner = 0
                geo_ref.yllcorner = 0    
                _write_xya_file(file_name,
                                self.get_data_points(absolute=True), 
                                self.get_all_attributes(),
                                geo_ref)
            else:
                _write_xya_file(file_name,
                                self.get_data_points(absolute=False), 
                                self.get_all_attributes(),
                                self.get_geo_reference())
                                    
        elif (file_name[-4:] == ".pts"):
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
                
        elif file_name[-4:] == ".txt" or file_name[-4:] == ".csv":
            msg = "ERROR: trying to write a .txt file with relative data."
            assert absolute, msg
            _write_csv_file(file_name,
                            self.get_data_points(absolute=True,
                                                 as_lat_long=as_lat_long), 
                            self.get_all_attributes(),
                            as_lat_long=as_lat_long)
              
        elif file_name[-4:] == ".urs" :
            msg = "ERROR: Can not write a .urs file as a relative file."
            assert absolute, msg
            _write_urs_file(file_name,
                            self.get_data_points(as_lat_long=True))
                                                          
        else:
            msg = 'Unknown file type %s ' %file_name
            raise IOError, msg 
        
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
        #FIXME: add the geo_reference to this
        
        points = self.get_data_points()
        sampled_points = take(points, indices)

        attributes = self.get_all_attributes()

        sampled_attributes = {}
        if attributes is not None:
            for key, att in attributes.items():
                sampled_attributes[key] = take(att, indices)

        return Geospatial_data(sampled_points, sampled_attributes)
    
    
    def split(self, factor=0.5, verbose=False):
        """Returns two geospatial_data object, first is size of the 'factor'
        smaller the original and the second is the remainer. The two new 
        object are disjoin set of each other. 
        
        Points of the two new object have selected RANDOMLY. 
        AND if factor is a decimal it will round (2.25 to 2 and 2.5 to 3)
        
        
        Input - the factor which to split the object, if 0.1 then 10% of the
            object will be returned
        
        Output - two geospatial_data objects that are disjoint sets of the 
            original
            """
        
        i=0
        self_size = len(self)
        random_list = []
        remainder_list = []
        new_size = round(factor*self_size)
   #     print'Split original %s by %s' %(self_size, factor)
   #     print'New samples are %s and %s in size' %(int(round(factor*self_size)),int(self_size-new_size))
        
        #find unique random numbers
        if verbose: print "make unique random number list and get indices"
        while i < new_size:
            random_num = randint(0,self_size-1)
            if random_num not in random_list:
                random_list.append(random_num)
                i=i+1

        #Make list of opposite to random_list
        if verbose: print "make list of opposite to random list"
        for i in range(0,self_size,1):
            remainder_list.append(i)

        #remove random list from remainder_list to get correct remainder_list
        #need to sort and reverse so the pop() works correctly
        random_list.sort()
        random_list.reverse()
        if verbose: print "get indices of opposite to random list"
        for i in random_list:
            remainder_list.pop(i)
            if verbose:
                if ((i/100)==(float(i)/100)): print "reached: ",i
            
        #get new samples
        if verbose: print "get values of indices for random list"
        G1 = self.get_sample(random_list)
        if verbose: print "get values of indices for opposite of random list"
        G2 = self.get_sample(remainder_list)

        return G1, G2

    def __iter__(self):
        """
        read in the header and save the file pointer position
        """

        from Scientific.IO.NetCDF import NetCDFFile
        
        #FIXME - what to do if the file isn't there

        if self.file_name[-4:]== ".xya":
            #let's just read it all
            pass
        elif self.file_name[-4:]== ".pts":
            
            # see if the file is there.  Throw a QUIET IO error if it isn't
            fd = open(self.file_name,'r')
            fd.close()
    
            #throws prints to screen if file not present
            self.fid = NetCDFFile(self.file_name, 'r')
            
            self.blocking_georef, self.blocking_keys, self.last_row = \
                     _read_pts_file_header(self.fid, self.verbose)
            self.start_row=0
            self.show_verbose = 0
            self.verbose_block_size =  (self.last_row +10) /10
        else:
            file_pointer = open(self.file_name)
            self.header, self.file_pointer = \
                         _read_csv_file_header(file_pointer)
            self.blocking_georef = None # Used for reconciling zones
        if self.max_read_lines is None:
            self.max_read_lines = MAX_READ_LINES
        return self
    
    def next(self):
        # read a block, instanciate a new geospatial and return it
        if self.file_name[-4:]== ".xya" :
            if not hasattr(self,'finished_reading') or \
                   self.finished_reading is False:
                #let's just read it all
                geo = Geospatial_data(self.file_name)
                self.finished_reading = True
            else:
                raise StopIteration
                self.finished_reading = False
                
        elif self.file_name[-4:]== ".pts":
            if self.start_row == self.last_row:
                # read the end of the file last iteration
                # remove blocking attributes
                self.fid.close()
                del self.max_read_lines
                del self.blocking_georef
                del self.last_row
                del self.start_row
                del self.blocking_keys
                del self.fid
                raise StopIteration
            fin_row = self.start_row + self.max_read_lines
            if fin_row > self.last_row:
                fin_row = self.last_row

                
            
            if self.verbose is True:
                if self.show_verbose >= self.start_row and \
                                       self.show_verbose < fin_row:
                 print 'Doing %d of %d' %(self.start_row, self.last_row+10)
                 self.show_verbose += max(self.max_read_lines,
                                          self.verbose_block_size)
            #call stuff
            pointlist, att_dict, = \
                   _read_pts_file_blocking( self.fid,
                                            self.start_row,
                                            fin_row,
                                            self.blocking_keys
                                            ) 
            geo = Geospatial_data(pointlist, att_dict, self.blocking_georef)
            self.start_row = fin_row
            
        else:
            try:
                pointlist, att_dict, geo_ref, self.file_pointer = \
                   _read_csv_file_blocking( self.file_pointer,
                                            self.header[:],
                                            max_read_lines=self.max_read_lines,
                                            verbose=self.verbose)

                # Check that the zones haven't changed.
                if geo_ref is not None:
                    geo_ref.reconcile_zones(self.blocking_georef)
                    self.blocking_georef = geo_ref
                elif self.blocking_georef is not None:
                    
                    msg = 'Geo reference given, then not given.'
                    msg += ' This should not happen.' 
                    raise ValueError, msg
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
        return geo

def _set_using_lat_long(latitudes,
                        longitudes,
                        geo_reference,
                        data_points,
                        points_are_lats_longs):
    """
    if the points has lat long info, assume it is in (lat, long) order.
    """
    
    if geo_reference is not None:
        msg = """A georeference is specified yet latitude and longitude
        are also specified!"""
        raise ValueError, msg
    
    if data_points is not None and not points_are_lats_longs:
        msg = """Data points are specified yet latitude and
        longitude are also specified!"""
        raise ValueError, msg
    
    if points_are_lats_longs:
        if data_points is None:
            msg = """Data points are not specified !"""
            raise ValueError, msg
        lats_longs = ensure_numeric(data_points)
        latitudes = ravel(lats_longs[:,0:1])
        longitudes = ravel(lats_longs[:,1:])
        
    if latitudes is None and longitudes is None:
        msg = """Latitudes and Longitudes are not."""
        raise ValueError, msg
    
    if latitudes is None:
        msg = """Longitudes are specified yet latitudes aren't!"""
        raise ValueError, msg
    
    if longitudes is None:
        msg = """Latitudes are specified yet longitudes aren't!"""
        raise ValueError, msg
    
    data_points, zone  = convert_from_latlon_to_utm(latitudes=latitudes,
                                                    longitudes=longitudes)
    return data_points, Geo_reference(zone=zone)
    
def _read_pts_file(file_name, verbose=False):
    """Read .pts NetCDF file
    
    Return a dic of array of points, and dic of array of attribute
    eg
    dic['points'] = [[1.0,2.0],[3.0,5.0]]
    dic['attributelist']['elevation'] = [[7.0,5.0]
    """    

    from Scientific.IO.NetCDF import NetCDFFile
    
    if verbose: print 'Reading ', file_name
    
        
    # see if the file is there.  Throw a QUIET IO error if it isn't
    fd = open(file_name,'r')
    fd.close()
    
    #throws prints to screen if file not present
    fid = NetCDFFile(file_name, 'r') 
    
    pointlist = array(fid.variables['points'])
    keys = fid.variables.keys()
    if verbose: print 'Got %d variables: %s' %(len(keys), keys)
    try:
        keys.remove('points')
    except IOError, e:       
        fid.close()    
        msg = 'Expected keyword "points" but could not find it'
        raise IOError, msg
    
    attributes = {}
    for key in keys:
        if verbose: print "reading attribute '%s'" %key
            
        attributes[key] = array(fid.variables[key])
    
    
    try:
        geo_reference = Geo_reference(NetCDFObject=fid)
    except AttributeError, e:
        geo_reference = None
    
    fid.close()
    
    return pointlist, attributes, geo_reference


def _read_csv_file(file_name, verbose=False):
    """Read .csv file
    
    Return a dic of array of points, and dic of array of attribute
    eg
    dic['points'] = [[1.0,2.0],[3.0,5.0]]
    dic['attributelist']['elevation'] = [[7.0,5.0]
    """
    
    #from anuga.shallow_water.data_manager import Exposure_csv
    #csv =Exposure_csv(file_name)
    
    file_pointer = open(file_name)
    header, file_pointer = _read_csv_file_header(file_pointer)
    try:
        pointlist, att_dict, geo_ref, file_pointer = \
                   _read_csv_file_blocking( \
                file_pointer,
                header,
                max_read_lines=1e30) #If the file is bigger that this, block..
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

def _read_csv_file_blocking(file_pointer, header,
                            delimiter=CSV_DELIMITER,
                            max_read_lines=MAX_READ_LINES,
                            verbose=False):
    

    """
    Read the body of a .csv file.
    header: The list header of the csv file, with the x and y labels.
    """
    points = []
    pointattributes = []
    att_dict = {}

    #This is to remove the x and y headers.
    header = header[:]
    x_header = header.pop(0)
    y_header = header.pop(0)
    
    read_lines = 0
    while read_lines<max_read_lines:
        line = file_pointer.readline()
        #print "line",line
        numbers = clean_line(line,delimiter)
        if len(numbers) <= 1:
            break
        if line[0] == '#':
            continue
        read_lines += 1
        if True: # remove.. #if numbers != []:
            try:
                x = float(numbers[0])
                y = float(numbers[1])
                points.append([x,y])
                numbers.pop(0)
                numbers.pop(0)
                if len(header) != len(numbers):
                    
                    file_pointer.close() 
                    # It might not be a problem with the header
                    #raise TitleAmountError
                    msg = "File load error.  There might be a problem with the file header"
                    raise IOError, msg
                for i,num in enumerate(numbers):
                    num.strip()
                    if num != '\n' and num != '':
                        #attributes.append(float(num))
                        att_dict.setdefault(header[i],[]).append(float(num))
            #except IOError:           
            except ValueError:
                raise SyntaxError
    if points == []:
        raise StopIteration
        
    
    pointlist = array(points).astype(Float)
    for key in att_dict.keys():
        att_dict[key] = array(att_dict[key]).astype(Float)

    #Do stuff here so the info is in lat's and longs
    geo_ref = None
    x_header = lower(x_header[:3])
    y_header = lower(y_header[:3])
    if (x_header == 'lon' or  x_header == 'lat') and \
       (y_header == 'lon' or  y_header == 'lat'):
        if x_header == 'lon':
            longitudes = ravel(pointlist[:,0:1])
            latitudes = ravel(pointlist[:,1:])
        else:
            latitudes = ravel(pointlist[:,0:1])
            longitudes = ravel(pointlist[:,1:])
        
        pointlist, geo_ref = _set_using_lat_long(latitudes,
                                                 longitudes,
                                                 geo_reference=None,
                                                 data_points=None,
                                                 points_are_lats_longs=False)
    return pointlist, att_dict, geo_ref, file_pointer 

def _read_pts_file_header(fid, verbose=False):

    """
    Read the geo_reference of a .pts file
    """
    
    keys = fid.variables.keys()
    try:
        keys.remove('points')
    except IOError, e:       
        fid.close()    
        msg = 'Expected keyword "points" but could not find it'
        raise IOError, msg
    if verbose: print 'Got %d variables: %s' %(len(keys), keys)
    
    try:
        geo_reference = Geo_reference(NetCDFObject=fid)
    except AttributeError, e:
        geo_reference = None

    return geo_reference, keys, fid.dimensions['number_of_points']

def _read_pts_file_blocking(fid, start_row, fin_row, keys):
                            #verbose=False):
    

    """
    Read the body of a .csv file.
    header: The list header of the csv file, with the x and y labels.
    """
    
    pointlist = array(fid.variables['points'][start_row:fin_row])
    
    attributes = {}
    for key in keys:
        attributes[key] = array(fid.variables[key][start_row:fin_row])

    return pointlist, attributes
    
    
def _read_xya_file(fd, delimiter):
    points = []
    pointattributes = []
    title = fd.readline()
    att_names = clean_line(title,delimiter)
    att_dict = {}
    line = fd.readline()
    numbers = clean_line(line,delimiter)
    
    while len(numbers) > 1 and line[0] <> '#':
        if numbers != []:
            try:
                x = float(numbers[0])
                y = float(numbers[1])
                points.append([x,y])
                numbers.pop(0)
                numbers.pop(0)
                if len(att_names) != len(numbers):
                    fd.close()
                    # It might not be a problem with the title
                    #raise TitleAmountError
                    raise IOError
                for i,num in enumerate(numbers):
                    num.strip()
                    if num != '\n' and num != '':
                        #attributes.append(float(num))
                        att_dict.setdefault(att_names[i],[]).append(float(num))
            except ValueError:
                raise SyntaxError
        line = fd.readline()
        numbers = clean_line(line,delimiter)
    
    if line == '':
        geo_reference = None
    else:
        geo_reference = Geo_reference(ASCIIFile=fd,read_title=line)
        
    
    pointlist = array(points).astype(Float)
    for key in att_dict.keys():
        att_dict[key] = array(att_dict[key]).astype(Float)
    
    return pointlist, att_dict, geo_reference

def _write_pts_file(file_name,
                    write_data_points,
                    write_attributes=None, 
                    write_geo_reference=None):
    """
    Write .pts NetCDF file   

    NOTE: Below might not be valid ask Duncan : NB 5/2006
    
    WARNING: This function mangles the point_atts data structure
    #F??ME: (DSG)This format has issues.
    # There can't be an attribute called points 
    # consider format change
    # method changed by NB not sure if above statement is correct
    
    should create new test for this
    legal_keys = ['pointlist', 'attributelist', 'geo_reference']
    for key in point_atts.keys():
        msg = 'Key %s is illegal. Valid keys are %s' %(key, legal_keys) 
        assert key in legal_keys, msg
    """    
    from Scientific.IO.NetCDF import NetCDFFile
    # NetCDF file definition
    outfile = NetCDFFile(file_name, 'w')
    
    #Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = 'NetCDF format for compact and portable storage ' +\
                          'of spatial point data'
    
    # dimension definitions
    shape = write_data_points.shape[0]
    outfile.createDimension('number_of_points', shape)  
    outfile.createDimension('number_of_dimensions', 2) #This is 2d data
    
    # variable definition
    outfile.createVariable('points', Float, ('number_of_points',
                                             'number_of_dimensions'))

    #create variables  
    outfile.variables['points'][:] = write_data_points #.astype(Float32)

    if write_attributes is not None:
        for key in write_attributes.keys():
            outfile.createVariable(key, Float, ('number_of_points',))
            outfile.variables[key][:] = write_attributes[key] #.astype(Float32)
        
    if write_geo_reference is not None:
        write_NetCDF_georeference(write_geo_reference, outfile)
        
    outfile.close() 
  


def _write_xya_file(file_name,
                    write_data_points,
                    write_attributes=None, 
                    write_geo_reference=None, 
                    delimiter=','):
    """
    export a file, file_name, with the xya format
    
    """
    points = write_data_points 
    pointattributes = write_attributes
    
    fd = open(file_name,'w')
    titlelist = ""
    if pointattributes is not None:    
        for title in pointattributes.keys():
            titlelist = titlelist + title + delimiter
        titlelist = titlelist[0:-len(delimiter)] # remove the last delimiter
    fd.write(titlelist+"\n")
    
    #<vertex #> <x> <y> [attributes]
    for i, vert in enumerate( points):


        if pointattributes is not None:            
            attlist = ","
            for att in pointattributes.keys():
                attlist = attlist + str(pointattributes[att][i])+ delimiter
            attlist = attlist[0:-len(delimiter)] # remove the last delimiter
            attlist.strip()
        else:
            attlist = ''

        fd.write(str(vert[0]) + delimiter +
                 str(vert[1]) + attlist + "\n")

    if  write_geo_reference is not None:
        write_geo_reference = ensure_geo_reference(write_geo_reference)
        write_geo_reference.write_ASCII(fd)
    fd.close()


def _write_csv_file(file_name,
                    write_data_points,
                    write_attributes=None,
                    as_lat_long=False,
                    delimiter=','):
    """
    export a file, file_name, with the xya format
    
    """
    points = write_data_points 
    pointattributes = write_attributes
    
    fd = open(file_name,'w')
    if as_lat_long:
        titlelist = "latitude" + delimiter + "longitude"  + delimiter
    else:
        titlelist = "x" + delimiter + "y"  + delimiter
    if pointattributes is not None:    
        for title in pointattributes.keys():
            titlelist = titlelist + title + delimiter
        titlelist = titlelist[0:-len(delimiter)] # remove the last delimiter
    fd.write(titlelist+"\n")
    
    # <x/lat> <y/long> [attributes]
    for i, vert in enumerate( points):


        if pointattributes is not None:            
            attlist = ","
            for att in pointattributes.keys():
                attlist = attlist + str(pointattributes[att][i])+ delimiter
            attlist = attlist[0:-len(delimiter)] # remove the last delimiter
            attlist.strip()
        else:
            attlist = ''

        fd.write(str(vert[0]) + delimiter +
                 str(vert[1]) + attlist + "\n")

    fd.close()
  
def _write_urs_file(file_name,
                    points,
                    delimiter=' '):
    """
    export a file, file_name, with the urs format
    the data points are in lats and longs
    
    """    
    fd = open(file_name,'w')   
    fd.write(str(len(points))+"\n")   
    # <lat> <long> <id#>
    for i, vert in enumerate( points):
        fd.write(str(round(vert[0],7)) + delimiter + \
                 str(round(vert[1],7)) + delimiter +str(i)+ "\n")
    fd.close()
    
def _point_atts2array(point_atts):
    point_atts['pointlist'] = array(point_atts['pointlist']).astype(Float)
    
    for key in point_atts['attributelist'].keys():
        point_atts['attributelist'][key]=\
                array(point_atts['attributelist'][key]).astype(Float)
    return point_atts

   


def geospatial_data2points_dictionary(geospatial_data):
    """Convert geospatial data to points_dictionary
    """

    points_dictionary = {}
    points_dictionary['pointlist'] = geospatial_data.data_points

    points_dictionary['attributelist'] = {}

    for attribute_name in geospatial_data.attributes.keys():
        val = geospatial_data.attributes[attribute_name]
        points_dictionary['attributelist'][attribute_name] = val

    points_dictionary['geo_reference'] = geospatial_data.geo_reference

    return points_dictionary

    
def points_dictionary2geospatial_data(points_dictionary):
    """Convert points_dictionary to geospatial data object
    """

    msg = 'Points dictionary must have key pointlist' 
    assert points_dictionary.has_key('pointlist'), msg

    msg = 'Points dictionary must have key attributelist'     
    assert points_dictionary.has_key('attributelist'), msg        

    if points_dictionary.has_key('geo_reference'):
        geo = points_dictionary['geo_reference']
    else:
        geo = None
    
    return Geospatial_data(points_dictionary['pointlist'],
                           points_dictionary['attributelist'],
                           geo_reference = geo)

def clean_line(line,delimiter):      
    """Remove whitespace
    """
    #print ">%s" %line
    line = line.strip()
    #print "stripped>%s" %line
    numbers = line.split(delimiter)
    i = len(numbers) - 1
    while i >= 0:
        if numbers[i] == '':
            numbers.pop(i)
        else:
            numbers[i] = numbers[i].strip()
        
        i += -1
    #for num in numbers:
    #    print "num>%s<" %num
    return numbers
            
def ensure_absolute(points, geo_reference=None):
    """
    This function inputs several formats and
    outputs one format. - a numeric array of absolute points.

    Inputed formats are;
    points: List or numeric array of coordinate pairs [xi, eta] of
              points or geospatial object or points file name

    mesh_origin: A geo_reference object or 3-tuples consisting of
                 UTM zone, easting and northing.
                 If specified vertex coordinates are assumed to be
                 relative to their respective origins.
    """
    if isinstance(points,type('')):
        #It's a string
        #assume it is a point file
        points = Geospatial_data(file_name = points)
        
    if isinstance(points,Geospatial_data):
        points = points.get_data_points( \
                absolute = True)
        msg = "Use a Geospatial_data object or a mesh origin. Not both."
        assert geo_reference == None, msg
            
    else:
        points = ensure_numeric(points, Float)
    if geo_reference is None:
        geo = None #Geo_reference()
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
    """
    This function inputs several formats and
    outputs one format. - a geospatial_data instance.

    Inputed formats are;
    points: List or numeric array of coordinate pairs [xi, eta] of
              points or geospatial object

    mesh_origin: A geo_reference object or 3-tuples consisting of
                 UTM zone, easting and northing.
                 If specified vertex coordinates are assumed to be
                 relative to their respective origins.
    """
    if isinstance(points,Geospatial_data):
        msg = "Use a Geospatial_data object or a mesh origin. Not both."
        assert geo_reference == None, msg
        return points    
    else:
        points = ensure_numeric(points, Float)
    if geo_reference is None:
        geo = None
    else:
        if isinstance(geo_reference, Geo_reference):
            geo = geo_reference
        else:
            geo = Geo_reference(geo_reference[0],
                                geo_reference[1],
                                geo_reference[2])
    points = Geospatial_data(data_points=points, geo_reference=geo)        
    return points

#def file2xya(filename):
    
#    G = Geospatial_data(filename)
#    G.export_points_file()

             

         
if __name__ == "__main__":
    pass
    
