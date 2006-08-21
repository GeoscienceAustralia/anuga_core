"""datamanager.py - input output for AnuGA


This module takes care of reading and writing datafiles such as topograhies,
model output, etc


Formats used within AnuGA:

.sww: Netcdf format for storing model output f(t,x,y)
.tms: Netcdf format for storing time series f(t)

.xya: ASCII format for storing arbitrary points and associated attributes
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
TSH:          Triangular meshes (e.g. created from pmesh)
NC            Model outputs for use as boundary conditions (e.g from MOST)


AUTOMATICALLY CREATED FILES:

ASC, PRJ  ->  DEM  ->  PTS: Conversion of DEM's to native pts file

NC -> SWW: Conversion of MOST bundary files to boundary sww

PTS + TSH -> TSH with elevation: Least squares fit

TSH -> SWW:  Conversion of TSH to sww viewable using Swollen

TSH + Boundary SWW -> SWW: Simluation using pyvolution

"""

import exceptions
class TitleValueError(exceptions.Exception): pass
class DataMissingValuesError(exceptions.Exception): pass
class DataFileNotOpenError(exceptions.Exception): pass
class DataTimeError(exceptions.Exception): pass
class DataDomainError(exceptions.Exception): pass



import csv
import os

from Numeric import concatenate, array, Float, Int, Int32, resize, sometrue, \
     searchsorted, zeros, allclose, around


from coordinate_transforms.geo_reference import Geo_reference
from geospatial_data.geospatial_data import Geospatial_data


def make_filename(s):
    """Transform argument string into a suitable filename
    """

    s = s.strip()
    s = s.replace(' ', '_')
    s = s.replace('(', '')
    s = s.replace(')', '')
    s = s.replace('__', '_')

    return s


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

    import os, sys
    import os.path

    if sys.platform in ['nt', 'dos', 'win32', 'what else?']:
        unix = 0
    else:
        unix = 1


    if path[-1] != os.sep:
        path = path + os.sep  # Add separator for directories

    path = os.path.expanduser(path) # Expand ~ or ~user in pathname
    if not (os.access(path,os.R_OK and os.W_OK) or path == ''):
        try:
            exitcode=os.mkdir(path)

            # Change access rights if possible
            #
            if unix:
                exitcode=os.system('chmod 775 '+path)
            else:
                pass  # FIXME: What about acces rights under Windows?

            if verbose: print 'MESSAGE: Directory', path, 'created.'

        except:
            print 'WARNING: Directory', path, 'could not be created.'
            if unix:
                path = '/tmp/'
            else:
                path = 'C:'

            print 'Using directory %s instead' %path

    return(path)



def del_dir(path):
    """Recursively delete directory path and all its contents
    """

    import os

    if os.path.isdir(path):
        for file in os.listdir(path):
            X = os.path.join(path, file)


            if os.path.isdir(X) and not os.path.islink(X):
                del_dir(X)
            else:
                try:
                    os.remove(X)
                except:
                    print "Could not remove file %s" %X

        os.rmdir(path)



def create_filename(datadir, filename, format, size=None, time=None):

    import os
    #from config import data_dir

    FN = check_dir(datadir) + filename

    if size is not None:
        FN += '_size%d' %size

    if time is not None:
        FN += '_time%.2f' %time

    FN += '.' + format
    return FN


def get_files(datadir, filename, format, size):
    """Get all file (names) with gven name, size and format
    """

    import glob

    import os
    #from config import data_dir

    dir = check_dir(datadir)

    pattern = dir + os.sep + filename + '_size=%d*.%s' %(size, format)
    return glob.glob(pattern)



#Generic class for storing output to e.g. visualisation or checkpointing
class Data_format:
    """Generic interface to data formats
    """


    def __init__(self, domain, extension, mode = 'w'):
        assert mode in ['r', 'w', 'a'], '''Mode %s must be either:''' %mode +\
                                        '''   'w' (write)'''+\
                                        '''   'r' (read)''' +\
                                        '''   'a' (append)'''

        #Create filename
        #self.filename = create_filename(domain.get_datadir(),
        #domain.get_name(), extension, len(domain))


        self.filename = create_filename(domain.get_datadir(),
            domain.get_name(), extension)

        #print 'F', self.filename
        self.timestep = 0
        self.number_of_volumes = len(domain)
        self.domain = domain


        #FIXME: Should we have a general set_precision function?



#Class for storing output to e.g. visualisation
class Data_format_sww(Data_format):
    """Interface to native NetCDF format (.sww) for storing model output

    There are two kinds of data

    1: Constant data: Vertex coordinates and field values. Stored once
    2: Variable data: Conserved quantities. Stored once per timestep.

    All data is assumed to reside at vertex locations.
    """


    def __init__(self, domain, mode = 'w',\
                 max_size = 2000000000,
                 recursion = False):
        from Scientific.IO.NetCDF import NetCDFFile
        from Numeric import Int, Float, Float32

        self.precision = Float32 #Use single precision for quantities
        if hasattr(domain, 'max_size'):
            self.max_size = domain.max_size #file size max is 2Gig
        else:
            self.max_size = max_size
        self.recursion = recursion
        self.mode = mode

        Data_format.__init__(self, domain, 'sww', mode)


        # NetCDF file definition
        fid = NetCDFFile(self.filename, mode)

        if mode == 'w':

            #Create new file
            fid.institution = 'Geoscience Australia'
            fid.description = 'Output from pyvolution suitable for plotting'

            if domain.smooth:
                fid.smoothing = 'Yes'
            else:
                fid.smoothing = 'No'

            fid.order = domain.default_order

            if hasattr(domain, 'texture'):
                fid.texture = domain.texture
        #else:
            #    fid.texture = 'None'

            #Reference point
            #Start time in seconds since the epoch (midnight 1/1/1970)
            #FIXME: Use Georef
            fid.starttime = domain.starttime

            # dimension definitions
            fid.createDimension('number_of_volumes', self.number_of_volumes)
            fid.createDimension('number_of_vertices', 3)

            if domain.smooth is True:
                fid.createDimension('number_of_points', len(domain.vertexlist))
            else:
                fid.createDimension('number_of_points', 3*self.number_of_volumes)

            fid.createDimension('number_of_timesteps', None) #extensible

            # variable definitions
            fid.createVariable('x', self.precision, ('number_of_points',))
            fid.createVariable('y', self.precision, ('number_of_points',))
            fid.createVariable('elevation', self.precision, ('number_of_points',))
            if domain.geo_reference is not None:
                domain.geo_reference.write_NetCDF(fid)

            #FIXME: Backwards compatibility
            fid.createVariable('z', self.precision, ('number_of_points',))
            #################################

            fid.createVariable('volumes', Int, ('number_of_volumes',
                                                'number_of_vertices'))

            fid.createVariable('time', Float,  # Always use full precision lest two timesteps 
	                                       # close to each other may appear as the same step 
                               ('number_of_timesteps',))

            fid.createVariable('stage', self.precision,
                               ('number_of_timesteps',
                                'number_of_points'))

            fid.createVariable('xmomentum', self.precision,
                               ('number_of_timesteps',
                                'number_of_points'))

            fid.createVariable('ymomentum', self.precision,
                               ('number_of_timesteps',
                                'number_of_points'))

        #Close
        fid.close()


    def store_connectivity(self):
        """Specialisation of store_connectivity for net CDF format

        Writes x,y,z coordinates of triangles constituting
        the bed elevation.
        """

        from Scientific.IO.NetCDF import NetCDFFile

        from Numeric import concatenate, Int

        domain = self.domain

        #Get NetCDF
        fid = NetCDFFile(self.filename, 'a')  #Open existing file for append

        # Get the variables
        x = fid.variables['x']
        y = fid.variables['y']
        z = fid.variables['elevation']

        volumes = fid.variables['volumes']

        # Get X, Y and bed elevation Z
        Q = domain.quantities['elevation']
        X,Y,Z,V = Q.get_vertex_values(xy=True,
                      precision = self.precision)



        x[:] = X.astype(self.precision)
        y[:] = Y.astype(self.precision)
        z[:] = Z.astype(self.precision)

        #FIXME: Backwards compatibility
        z = fid.variables['z']
        z[:] = Z.astype(self.precision)
        ################################

        volumes[:] = V.astype(volumes.typecode())

        #Close
        fid.close()

    def store_timestep(self, names):
        """Store time and named quantities to file
        """
        from Scientific.IO.NetCDF import NetCDFFile
        import types
        from time import sleep
        from os import stat

        minimum_allowed_depth = 0.001
        #minimum_allowed_depth = 0.0  #FIXME pass in or read from domain
        from Numeric import choose

        #Get NetCDF
        retries = 0
        file_open = False
        while not file_open and retries < 10:
            try:
                fid = NetCDFFile(self.filename, 'a') #Open existing file
            except IOError:
                #This could happen if someone was reading the file.
                #In that case, wait a while and try again
                msg = 'Warning (store_timestep): File %s could not be opened'\
                      %self.filename
                msg += ' - trying step %s again' %self.domain.time
                print msg
                retries += 1
                sleep(1)
            else:
                file_open = True

        if not file_open:
            msg = 'File %s could not be opened for append' %self.filename
            raise DataFileNotOpenError, msg



        #Check to see if the file is already too big:
        time = fid.variables['time']
        i = len(time)+1
        file_size = stat(self.filename)[6]
        file_size_increase =  file_size/i
        if file_size + file_size_increase > self.max_size*(2**self.recursion):
            #in order to get the file name and start time correct,
            #I change the domian.filename and domain.starttime.
            #This is the only way to do this without changing
            #other modules (I think).

            #write a filename addon that won't break swollens reader
            #(10.sww is bad)
            filename_ext = '_time_%s'%self.domain.time
            filename_ext = filename_ext.replace('.', '_')
            #remember the old filename, then give domain a
            #name with the extension
            old_domain_filename = self.domain.filename
            if not self.recursion:
                self.domain.filename = self.domain.filename+filename_ext

            #change the domain starttime to the current time
            old_domain_starttime = self.domain.starttime
            self.domain.starttime = self.domain.time

            #build a new data_structure.
            next_data_structure=\
                Data_format_sww(self.domain, mode=self.mode,\
                                max_size = self.max_size,\
                                recursion = self.recursion+1)
            if not self.recursion:
                print '    file_size = %s'%file_size
                print '    saving file to %s'%next_data_structure.filename
            #set up the new data_structure
            self.domain.writer = next_data_structure

            #FIXME - could be cleaner to use domain.store_timestep etc.
            next_data_structure.store_connectivity()
            next_data_structure.store_timestep(names)
            fid.sync()
            fid.close()

            #restore the old starttime and filename
            self.domain.starttime = old_domain_starttime
            self.domain.filename = old_domain_filename
        else:
            self.recursion = False
            domain = self.domain

            # Get the variables
            time = fid.variables['time']
            stage = fid.variables['stage']
            xmomentum = fid.variables['xmomentum']
            ymomentum = fid.variables['ymomentum']
            i = len(time)

            #Store time
            time[i] = self.domain.time


            if type(names) not in [types.ListType, types.TupleType]:
                names = [names]

            for name in names:
                # Get quantity
                Q = domain.quantities[name]
                A,V = Q.get_vertex_values(xy = False,
                                          precision = self.precision)

                #FIXME: Make this general (see below)
                if name == 'stage':
                    z = fid.variables['elevation']
                    #print z[:]
                    #print A-z[:]
                    A = choose( A-z[:] >= minimum_allowed_depth, (z[:], A))
                    stage[i,:] = A.astype(self.precision)
                elif name == 'xmomentum':
                    xmomentum[i,:] = A.astype(self.precision)
                elif name == 'ymomentum':
                    ymomentum[i,:] = A.astype(self.precision)

                #As in....
            #eval( name + '[i,:] = A.astype(self.precision)' )
            #FIXME: But we need a UNIT test for that before refactoring



            #Flush and close
            fid.sync()
            fid.close()



#Class for handling checkpoints data
class Data_format_cpt(Data_format):
    """Interface to native NetCDF format (.cpt)
    """


    def __init__(self, domain, mode = 'w'):
        from Scientific.IO.NetCDF import NetCDFFile
        from Numeric import Int, Float, Float

        self.precision = Float #Use full precision

        Data_format.__init__(self, domain, 'sww', mode)


        # NetCDF file definition
        fid = NetCDFFile(self.filename, mode)

        if mode == 'w':
            #Create new file
            fid.institution = 'Geoscience Australia'
            fid.description = 'Checkpoint data'
            #fid.smooth = domain.smooth
            fid.order = domain.default_order

            # dimension definitions
            fid.createDimension('number_of_volumes', self.number_of_volumes)
            fid.createDimension('number_of_vertices', 3)

            #Store info at all vertices (no smoothing)
            fid.createDimension('number_of_points', 3*self.number_of_volumes)
            fid.createDimension('number_of_timesteps', None) #extensible

            # variable definitions

            #Mesh
            fid.createVariable('x', self.precision, ('number_of_points',))
            fid.createVariable('y', self.precision, ('number_of_points',))


            fid.createVariable('volumes', Int, ('number_of_volumes',
                                                'number_of_vertices'))

            fid.createVariable('time', self.precision,
                               ('number_of_timesteps',))

            #Allocate space for all quantities
            for name in domain.quantities.keys():
                fid.createVariable(name, self.precision,
                                   ('number_of_timesteps',
                                    'number_of_points'))

        #Close
        fid.close()


    def store_checkpoint(self):
        """
        Write x,y coordinates of triangles.
        Write connectivity (
        constituting
        the bed elevation.
        """

        from Scientific.IO.NetCDF import NetCDFFile

        from Numeric import concatenate

        domain = self.domain

        #Get NetCDF
        fid = NetCDFFile(self.filename, 'a')  #Open existing file for append

        # Get the variables
        x = fid.variables['x']
        y = fid.variables['y']

        volumes = fid.variables['volumes']

        # Get X, Y and bed elevation Z
        Q = domain.quantities['elevation']
        X,Y,Z,V = Q.get_vertex_values(xy=True,
                      precision = self.precision)



        x[:] = X.astype(self.precision)
        y[:] = Y.astype(self.precision)
        z[:] = Z.astype(self.precision)

        volumes[:] = V

        #Close
        fid.close()


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
                fid = NetCDFFile(self.filename, 'a') #Open existing file
            except IOError:
                #This could happen if someone was reading the file.
                #In that case, wait a while and try again
                msg = 'Warning (store_timestep): File %s could not be opened'\
                  %self.filename
                msg += ' - trying again'
                print msg
                retries += 1
                sleep(1)
            else:
                file_open = True

        if not file_open:
            msg = 'File %s could not be opened for append' %self.filename
            raise DataFileNotOPenError, msg


        domain = self.domain

        # Get the variables
        time = fid.variables['time']
        stage = fid.variables['stage']
        i = len(time)

        #Store stage
        time[i] = self.domain.time

        # Get quantity
        Q = domain.quantities[name]
        A,V = Q.get_vertex_values(xy=False,
                  precision = self.precision)

        stage[i,:] = A.astype(self.precision)

        #Flush and close
        fid.sync()
        fid.close()





#Function for storing xya output
#FIXME Not done yet for this version
#This is obsolete.  Use geo_spatial_data instead
class Data_format_xya(Data_format):
    """Generic interface to data formats
    """


    def __init__(self, domain, mode = 'w'):
        from Scientific.IO.NetCDF import NetCDFFile
        from Numeric import Int, Float, Float32

        self.precision = Float32 #Use single precision

        Data_format.__init__(self, domain, 'xya', mode)



    #FIXME -This is the old xya format
    def store_all(self):
        """Specialisation of store all for xya format

        Writes x,y,z coordinates of triangles constituting
        the bed elevation.
        """

        from Numeric import concatenate

        domain = self.domain

        fd = open(self.filename, 'w')


        if domain.smooth is True:
            number_of_points =  len(domain.vertexlist)
        else:
            number_of_points = 3*self.number_of_volumes

        numVertAttrib = 3 #Three attributes is what is assumed by the xya format

        fd.write(str(number_of_points) + " " + str(numVertAttrib) +\
                 " # <vertex #> <x> <y> [attributes]" + "\n")


        # Get X, Y, bed elevation and friction (index=0,1)
        X,Y,A,V = domain.get_vertex_values(xy=True, value_array='field_values',
                                           indices = (0,1), precision = self.precision)

        bed_eles = A[:,0]
        fricts = A[:,1]

        # Get stage (index=0)
        B,V = domain.get_vertex_values(xy=False, value_array='conserved_quantities',
                                       indices = (0,), precision = self.precision)

        stages = B[:,0]

        #<vertex #> <x> <y> [attributes]
        for x, y, bed_ele, stage, frict in map(None, X, Y, bed_eles,
                                               stages, fricts):

            s = '%.6f %.6f %.6f %.6f %.6f\n' %(x, y, bed_ele, stage, frict)
            fd.write(s)

        #close
        fd.close()


    def store_timestep(self, t, V0, V1, V2):
        """Store time, water heights (and momentums) to file
        """
        pass


#### NBED national exposure database
    
LAT_TITLE = 'LATITUDE'
LONG_TITLE = 'LONGITUDE'
X_TITLE = 'x'
Y_TITLE = 'y'
class Exposure_csv:
    def __init__(self,file_name, latitude_title=LAT_TITLE,
                 longitude_title=LONG_TITLE, is_x_y_locations=None,
                 x_title=X_TITLE, y_title=Y_TITLE,
                 refine_polygon=None):
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
        self._load_exposure_csv(self._file_name)
        try:
            lats = self._attribute_dic[latitude_title]
            longs = self._attribute_dic[longitude_title]
            
        except KeyError:
            # maybe a warning..
            #Let's see if this works..
            if False != is_x_y_locations:
                is_x_y_locations = True
            pass
        else:
            self._geospatial = Geospatial_data(latitudes = lats,
                 longitudes = longs)

        if is_x_y_locations is True:
            if self._geospatial is not None:
                pass #fixme throw an error
            try:
                xs = self._attribute_dic[x_title]
                ys = self._attribute_dic[y_title]
                points = [[float(i),float(j)] for i,j in map(None,xs,ys)]
            except KeyError:
                # maybe a warning..
                raise TitleValueError, msg
            else:
                self._geospatial = Geospatial_data(data_points=points)
            
        # create a list of points that are in the refining_polygon
        # described by a list of indexes representing the points

    def __cmp__(self, other):
        #print "self._attribute_dic",self._attribute_dic
        #print "other._attribute_dic",other._attribute_dic
        #print "self._title_index_dic", self._title_index_dic
        #print "other._title_index_dic", other._title_index_dic
        
        #check that a is an instance of this class
        if isinstance(self, type(other)):
            result = cmp(self._attribute_dic, other._attribute_dic)
            if result <>0:
                return result
            # The order of the columns is important. Therefore.. 
            result = cmp(self._title_index_dic, other._title_index_dic)
            if result <>0:
                return result
            for self_ls, other_ls in map(None,self._attribute_dic, \
                    other._attribute_dic):
                result = cmp(self._attribute_dic[self_ls],
                             other._attribute_dic[other_ls])
                if result <>0:
                    return result
            return 0
        else:
            return 1
    
    def _load_exposure_csv(self, file_name):
        """
        Load in the csv as a dic, title as key and column info as value, .
        Also, create a dic, title as key and column index as value. 
        """
        #
        attribute_dic = {}
        title_index_dic = {}
        titles_stripped = [] # list of titles
        
        reader = csv.reader(file(file_name))

        # Read in and manipulate the title info
        titles = reader.next()
        for i,title in enumerate(titles):
            titles_stripped.append(title.strip())
            title_index_dic[title.strip()] = i
        title_count = len(titles_stripped)       
        #print "title_index_dic",title_index_dic

        
        #create a dic of colum values, indexed by column title
        for line in reader:
            if len(line) <> title_count:
                raise IOError #FIXME make this nicer
            for i, value in enumerate(line):
                attribute_dic.setdefault(titles_stripped[i],[]).append(value)
            
        return attribute_dic, title_index_dic

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
            msg = 'Therer is  no column called %s!' %column_name
            raise TitleValueError, msg
        return self._attribute_dic[column_name]


    def get_value(self, value_column_name,
                  known_column_name,
                  known_values,
                  use_refind_polygon=False):
        """
        Do linear interpolation on the known_colum, using the known_value,
        to return a value of the column_value_name.
        """
        pass


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

    def set_column(self, column_name, column_values, overwrite=False):
        """
        Add a column to the 'end' (with the right most column being the end)
        of the csv file.

        Set overwrite to True if you want to overwrite a column.
        
        Note, in column_name white space is removed and case is not checked.
        Precondition
        The column_name and column_values cannot have comma's in it.
        """
        
        value_row_count = \
            len(self._attribute_dic[self._title_index_dic.keys()[0]])
        if len(column_values) <> value_row_count: 
            msg = 'The number of column values must equal the number of rows.'
            raise DataMissingValuesError, msg
        
        if self._attribute_dic.has_key(column_name):
            if not overwrite:
                msg = 'Column name %s already in use!' %column_name
                raise TitleValueError, msg
        else:
            # New title.  Add it to the title index.
            self._title_index_dic[column_name] = len(self._title_index_dic)
        self._attribute_dic[column_name] = column_values
        #print "self._title_index_dic[column_name]",self._title_index_dic 

    def save(self, file_name=None):
        """
        Save the exposure csv file
        """
        if file_name is None:
            file_name = self._file_name
        
        fd = open(file_name,'wb')
        writer = csv.writer(fd)
        
        #Write the title to a cvs file
        line = [None]* len(self._title_index_dic)
        for title in self._title_index_dic.iterkeys():
            line[self._title_index_dic[title]]= title
        writer.writerow(line)
        
        # Write the values to a cvs file
        value_row_count = \
            len(self._attribute_dic[self._title_index_dic.keys()[0]])
        for row_i in range(value_row_count):
            line = [None]* len(self._title_index_dic)
            for title in self._title_index_dic.iterkeys():
                line[self._title_index_dic[title]]= \
                     self._attribute_dic[title][row_i]
            writer.writerow(line)


#Auxiliary
def write_obj(filename,x,y,z):
    """Store x,y,z vectors into filename (obj format)
       Vectors are assumed to have dimension (M,3) where
       M corresponds to the number elements.
       triangles are assumed to be disconnected

       The three numbers in each vector correspond to three vertices,

       e.g. the x coordinate of vertex 1 of element i is in x[i,1]

    """
    #print 'Writing obj to %s' % filename

    import os.path

    root, ext = os.path.splitext(filename)
    if ext == '.obj':
        FN = filename
    else:
        FN = filename + '.obj'


    outfile = open(FN, 'wb')
    outfile.write("# Triangulation as an obj file\n")

    M, N = x.shape
    assert N==3  #Assuming three vertices per element

    for i in range(M):
        for j in range(N):
            outfile.write("v %f %f %f\n" % (x[i,j],y[i,j],z[i,j]))

    for i in range(M):
        base = i*N
        outfile.write("f %d %d %d\n" % (base+1,base+2,base+3))

    outfile.close()


#########################################################
#Conversion routines
########################################################

def sww2obj(basefilename, size):
    """Convert netcdf based data output to obj
    """
    from Scientific.IO.NetCDF import NetCDFFile

    from Numeric import Float, zeros

    #Get NetCDF
    FN = create_filename('.', basefilename, 'sww', size)
    print 'Reading from ', FN
    fid = NetCDFFile(FN, 'r')  #Open existing file for read


    # Get the variables
    x = fid.variables['x']
    y = fid.variables['y']
    z = fid.variables['elevation']
    time = fid.variables['time']
    stage = fid.variables['stage']

    M = size  #Number of lines
    xx = zeros((M,3), Float)
    yy = zeros((M,3), Float)
    zz = zeros((M,3), Float)

    for i in range(M):
        for j in range(3):
            xx[i,j] = x[i+j*M]
            yy[i,j] = y[i+j*M]
            zz[i,j] = z[i+j*M]

    #Write obj for bathymetry
    FN = create_filename('.', basefilename, 'obj', size)
    write_obj(FN,xx,yy,zz)


    #Now read all the data with variable information, combine with
    #x,y info and store as obj

    for k in range(len(time)):
        t = time[k]
        print 'Processing timestep %f' %t

        for i in range(M):
            for j in range(3):
                zz[i,j] = stage[k,i+j*M]


        #Write obj for variable data
        #FN = create_filename(basefilename, 'obj', size, time=t)
        FN = create_filename('.', basefilename[:5], 'obj', size, time=t)
        write_obj(FN,xx,yy,zz)


def dat2obj(basefilename):
    """Convert line based data output to obj
    FIXME: Obsolete?
    """

    import glob, os
    from config import data_dir


    #Get bathymetry and x,y's
    lines = open(data_dir+os.sep+basefilename+'_geometry.dat', 'r').readlines()

    from Numeric import zeros, Float

    M = len(lines)  #Number of lines
    x = zeros((M,3), Float)
    y = zeros((M,3), Float)
    z = zeros((M,3), Float)

    ##i = 0
    for i, line in enumerate(lines):
        tokens = line.split()
        values = map(float,tokens)

        for j in range(3):
            x[i,j] = values[j*3]
            y[i,j] = values[j*3+1]
            z[i,j] = values[j*3+2]

        ##i += 1


    #Write obj for bathymetry
    write_obj(data_dir+os.sep+basefilename+'_geometry',x,y,z)


    #Now read all the data files with variable information, combine with
    #x,y info
    #and store as obj

    files = glob.glob(data_dir+os.sep+basefilename+'*.dat')

    for filename in files:
        print 'Processing %s' % filename

        lines = open(data_dir+os.sep+filename,'r').readlines()
        assert len(lines) == M
        root, ext = os.path.splitext(filename)

        #Get time from filename
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



        ##i = 0
        for i, line in enumerate(lines):
            tokens = line.split()
            values = map(float,tokens)

            for j in range(3):
                z[i,j] = values[j]

            ##i += 1

        #Write obj for variable data
        write_obj(data_dir+os.sep+basefilename+'_time=%.4f' %t,x,y,z)


def filter_netcdf(filename1, filename2, first=0, last=None, step = 1):
    """Read netcdf filename1, pick timesteps first:step:last and save to
    nettcdf file filename2
    """
    from Scientific.IO.NetCDF import NetCDFFile

    #Get NetCDF
    infile = NetCDFFile(filename1, 'r')  #Open existing file for read
    outfile = NetCDFFile(filename2, 'w')  #Open new file


    #Copy dimensions
    for d in infile.dimensions:
        outfile.createDimension(d, infile.dimensions[d])

    for name in infile.variables:
        var = infile.variables[name]
        outfile.createVariable(name, var.typecode(), var.dimensions)


    #Copy the static variables
    for name in infile.variables:
        if name == 'time' or name == 'stage':
            pass
        else:
            #Copy
            outfile.variables[name][:] = infile.variables[name][:]

    #Copy selected timesteps
    time = infile.variables['time']
    stage = infile.variables['stage']

    newtime = outfile.variables['time']
    newstage = outfile.variables['stage']

    if last is None:
        last = len(time)

    selection = range(first, last, step)
    for i, j in enumerate(selection):
        print 'Copying timestep %d of %d (%f)' %(j, last-first, time[j])
        newtime[i] = time[j]
        newstage[i,:] = stage[j,:]

    #Close
    infile.close()
    outfile.close()


#Get data objects
def get_dataobject(domain, mode='w'):
    """Return instance of class of given format using filename
    """

    cls = eval('Data_format_%s' %domain.format)
    return cls(domain, mode)

#FIXME move into geospatial.  There should only be one method that
# reads xya, and writes pts.
def xya2pts(basename_in, basename_out=None, verbose=False,
            #easting_min=None, easting_max=None,
            #northing_min=None, northing_max=None,
            stride = 1,
            attribute_name = 'elevation',
            z_func = None):
    """Read points data from ascii (.xya)

    Example:

              x(m)        y(m)        z(m)
         0.00000e+00  0.00000e+00  1.3535000e-01
         0.00000e+00  1.40000e-02  1.3535000e-01



    Convert to NetCDF pts format which is

    points:  (Nx2) Float array
    elevation: N Float array

    Only lines that contain three numeric values are processed

    If z_func is specified, it will be applied to the third field
    """

    import os
    #from Scientific.IO.NetCDF import NetCDFFile
    from Numeric import Float, arrayrange, concatenate

    root, ext = os.path.splitext(basename_in)

    if ext == '': ext = '.xya'

    #Get NetCDF
    infile = open(root + ext, 'r')  #Open existing xya file for read

    if verbose: print 'Reading xya points from %s' %(root + ext)

    points = []
    attribute = []
    for i, line in enumerate(infile.readlines()):

        if i % stride != 0: continue

        fields = line.split()

        try:
            assert len(fields) == 3
        except:
            print 'WARNING: Line %d doesn\'t have 3 elements: %s' %(i, line)

        try:
            x = float( fields[0] )
            y = float( fields[1] )
            z = float( fields[2] )
        except:
            continue

        points.append( [x, y] )

        if callable(z_func):
            attribute.append(z_func(z))
        else:
            attribute.append(z)


    #Get output file
    if basename_out == None:
        ptsname = root + '.pts'
    else:
        ptsname = basename_out + '.pts'

    if verbose: print 'Store to NetCDF file %s' %ptsname
    write_ptsfile(ptsname, points, attribute, attribute_name)



######Obsoleted by export_points in load_mesh
def write_ptsfile(ptsname, points, attribute, attribute_name = None,
                  zone=None, xllcorner=None, yllcorner=None):
    """Write points and associated attribute to pts (NetCDF) format
    """

    print 'WARNING: write_ptsfile is obsolete. Use export_points from load_mesh.loadASCII instead.'

    from Numeric import Float

    if attribute_name is None:
        attribute_name = 'attribute'


    from Scientific.IO.NetCDF import NetCDFFile

    # NetCDF file definition
    outfile = NetCDFFile(ptsname, 'w')


    #Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = 'NetCDF pts format for compact and '\
                          'portable storage of spatial point data'


    #Georeferencing
    from coordinate_transforms.geo_reference import Geo_reference
    if zone is None:
        assert xllcorner is None, 'xllcorner must be None'
        assert yllcorner is None, 'yllcorner must be None'
        Geo_reference().write_NetCDF(outfile)
    else:
        Geo_reference(zone, xllcorner, yllcorner).write_NetCDF(outfile)



    outfile.createDimension('number_of_points', len(points))
    outfile.createDimension('number_of_dimensions', 2) #This is 2d data

    # variable definitions
    outfile.createVariable('points', Float, ('number_of_points',
                                             'number_of_dimensions'))
    outfile.createVariable(attribute_name, Float, ('number_of_points',))

    # Get handles to the variables
    nc_points = outfile.variables['points']
    nc_attribute = outfile.variables[attribute_name]

    #Store data
    nc_points[:, :] = points
    nc_attribute[:] = attribute

    outfile.close()


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

    points:  (Nx2) Float array
    elevation: N Float array
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


def _dem2pts(basename_in, basename_out=None, verbose=False,
            easting_min=None, easting_max=None,
            northing_min=None, northing_max=None):
    """Read Digitial Elevation model from the following NetCDF format (.dem)

    Internal function. See public function dem2pts for details.
    """

    #FIXME: Can this be written feasibly using write_pts?

    import os
    from Scientific.IO.NetCDF import NetCDFFile
    from Numeric import Float, zeros, reshape, sum

    root = basename_in

    #Get NetCDF
    infile = NetCDFFile(root + '.dem', 'r')  #Open existing netcdf file for read

    if verbose: print 'Reading DEM from %s' %(root + '.dem')

    ncols = infile.ncols[0]
    nrows = infile.nrows[0]
    xllcorner = infile.xllcorner[0]  #Easting of lower left corner
    yllcorner = infile.yllcorner[0]  #Northing of lower left corner
    cellsize = infile.cellsize[0]
    NODATA_value = infile.NODATA_value[0]
    dem_elevation = infile.variables['elevation']

    zone = infile.zone[0]
    false_easting = infile.false_easting[0]
    false_northing = infile.false_northing[0]

    #Text strings
    projection = infile.projection
    datum = infile.datum
    units = infile.units


    #Get output file
    if basename_out == None:
        ptsname = root + '.pts'
    else:
        ptsname = basename_out + '.pts'

    if verbose: print 'Store to NetCDF file %s' %ptsname
    # NetCDF file definition
    outfile = NetCDFFile(ptsname, 'w')

    #Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = 'NetCDF pts format for compact and portable storage ' +\
                      'of spatial point data'
    #assign default values
    if easting_min is None: easting_min = xllcorner
    if easting_max is None: easting_max = xllcorner + ncols*cellsize
    if northing_min is None: northing_min = yllcorner
    if northing_max is None: northing_max = yllcorner + nrows*cellsize

    #compute offsets to update georeferencing
    easting_offset = xllcorner - easting_min
    northing_offset = yllcorner - northing_min

    #Georeferencing
    outfile.zone = zone
    outfile.xllcorner = easting_min #Easting of lower left corner
    outfile.yllcorner = northing_min #Northing of lower left corner
    outfile.false_easting = false_easting
    outfile.false_northing = false_northing

    outfile.projection = projection
    outfile.datum = datum
    outfile.units = units


    #Grid info (FIXME: probably not going to be used, but heck)
    outfile.ncols = ncols
    outfile.nrows = nrows

    dem_elevation_r = reshape(dem_elevation, (nrows, ncols))
    totalnopoints = nrows*ncols

    # calculating number of NODATA_values for each row in clipped region
    #FIXME: use array operations to do faster
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
            if easting_min <= x <= easting_max and \
               northing_min <= y <= northing_max:
                thisj = j
                thisi = i
                if dem_elevation_r[i,j] == NODATA_value: nn += 1

                if k == 0:
                    i1_0 = i
                    j1_0 = j
                k += 1

    index1 = j1_0
    index2 = thisj

    # dimension definitions
    nrows_in_bounding_box = int(round((northing_max-northing_min)/cellsize))
    ncols_in_bounding_box = int(round((easting_max-easting_min)/cellsize))

    clippednopoints = (thisi+1-i1_0)*(thisj+1-j1_0)
    nopoints = clippednopoints-nn

    clipped_dem_elev = dem_elevation_r[i1_0:thisi+1,j1_0:thisj+1]

    if verbose and nn > 0:
        print 'There are %d values in the elevation' %totalnopoints
        print 'There are %d values in the clipped elevation' %clippednopoints
        print 'There are %d NODATA_values in the clipped elevation' %nn

    outfile.createDimension('number_of_points', nopoints)
    outfile.createDimension('number_of_dimensions', 2) #This is 2d data

    # variable definitions
    outfile.createVariable('points', Float, ('number_of_points',
                                             'number_of_dimensions'))
    outfile.createVariable('elevation', Float, ('number_of_points',))

    # Get handles to the variables
    points = outfile.variables['points']
    elevation = outfile.variables['elevation']

    lenv = index2-index1+1
    #Store data
    global_index = 0
    #for i in range(nrows):
    for i in range(i1_0,thisi+1,1):
        if verbose and i%((nrows+10)/10)==0:
            print 'Processing row %d of %d' %(i, nrows)

        lower_index = global_index

        v = dem_elevation_r[i,index1:index2+1]
        no_NODATA = sum(v == NODATA_value)
        if no_NODATA > 0:
            newcols = lenv - no_NODATA #ncols_in_bounding_box - no_NODATA
        else:
            newcols = lenv #ncols_in_bounding_box

        telev = zeros(newcols, Float)
        tpoints = zeros((newcols, 2), Float)

        local_index = 0

        y = (nrows-i-1)*cellsize + yllcorner
        #for j in range(ncols):
        for j in range(j1_0,index2+1,1):

            x = j*cellsize + xllcorner
            if easting_min <= x <= easting_max and \
               northing_min <= y <= northing_max and \
               dem_elevation_r[i,j] <> NODATA_value:
                tpoints[local_index, :] = [x-easting_min,y-northing_min]
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

    points:  (Nx2) Float array
    elevation: N Float array
    """

    #FIXME: Can this be written feasibly using write_pts?

    import os
    from Scientific.IO.NetCDF import NetCDFFile
    from Numeric import Float, zeros, reshape

    root = basename_in

    #Get ASCII file
    infile = open(root + '.sdf', 'r')  #Open SDF file for read

    if verbose: print 'Reading DEM from %s' %(root + '.sdf')

    lines = infile.readlines()
    infile.close()

    if verbose: print 'Converting to pts format'

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


    #Now read all points
    points = []
    elevation = []
    for j, entries in enumerate(_read_hecras_cross_sections(lines[i:])):
        for k, entry in enumerate(entries):
            points.append(entry[:2])
            elevation.append(entry[2])


    msg = 'Actual #number_of_cross_sections == %d, Reported as %d'\
          %(j+1, number_of_cross_sections)
    assert j+1 == number_of_cross_sections, msg

    #Get output file
    if basename_out == None:
        ptsname = root + '.pts'
    else:
        ptsname = basename_out + '.pts'

    #FIXME (DSG-ON): use loadASCII export_points_file
    if verbose: print 'Store to NetCDF file %s' %ptsname
    # NetCDF file definition
    outfile = NetCDFFile(ptsname, 'w')

    #Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = 'NetCDF pts format for compact and portable ' +\
                          'storage of spatial point data derived from HEC-RAS'

    #Georeferencing
    outfile.zone = zone
    outfile.xllcorner = 0.0
    outfile.yllcorner = 0.0
    outfile.false_easting = 500000    #FIXME: Use defaults from e.g. config.py
    outfile.false_northing = 1000000

    outfile.projection = projection
    outfile.datum = datum
    outfile.units = units


    outfile.createDimension('number_of_points', len(points))
    outfile.createDimension('number_of_dimensions', 2) #This is 2d data

    # variable definitions
    outfile.createVariable('points', Float, ('number_of_points',
                                             'number_of_dimensions'))
    outfile.createVariable('elevation', Float, ('number_of_points',))

    # Get handles to the variables
    outfile.variables['points'][:] = points
    outfile.variables['elevation'][:] = elevation

    outfile.close()



def sww2dem(basename_in, basename_out = None,
            quantity = None,
            timestep = None,
            reduction = None,
            cellsize = 10,
            NODATA_value = -9999,
            easting_min = None,
            easting_max = None,
            northing_min = None,
            northing_max = None,
            verbose = False,
            origin = None,
            datum = 'WGS84',
        format = 'ers'):

    """Read SWW file and convert to Digitial Elevation model format (.asc or .ers)

    Example (ASC):

    ncols         3121
    nrows         1800
    xllcorner     722000
    yllcorner     5893000
    cellsize      25
    NODATA_value  -9999
    138.3698 137.4194 136.5062 135.5558 ..........

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
    'elevation'.

    if timestep (an index) is given, output quantity at that timestep

    if reduction is given use that to reduce quantity over all timesteps.

    datum

    format can be either 'asc' or 'ers'
    """

    import sys
    from Numeric import array, Float, concatenate, NewAxis, zeros, reshape, sometrue
    from Numeric import array2string

    from utilities.polygon import inside_polygon, outside_polygon, separate_points_by_polygon
    from util import apply_expression_to_dictionary

    msg = 'Format must be either asc or ers'
    assert format.lower() in ['asc', 'ers'], msg


    false_easting = 500000
    false_northing = 10000000

    if quantity is None:
        quantity = 'elevation'

    if reduction is None:
        reduction = max

    if basename_out is None:
        basename_out = basename_in + '_%s' %quantity

    swwfile = basename_in + '.sww'
    demfile = basename_out + '.' + format
    # Note the use of a .ers extension is optional (write_ermapper_grid will
    # deal with either option

    #Read sww file
    if verbose: print 'Reading from %s' %swwfile
    from Scientific.IO.NetCDF import NetCDFFile
    fid = NetCDFFile(swwfile)

    #Get extent and reference
    x = fid.variables['x'][:]
    y = fid.variables['y'][:]
    volumes = fid.variables['volumes'][:]

    number_of_timesteps = fid.dimensions['number_of_timesteps']
    number_of_points = fid.dimensions['number_of_points']
    if origin is None:

        #Get geo_reference
        #sww files don't have to have a geo_ref
        try:
            geo_reference = Geo_reference(NetCDFObject=fid)
        except AttributeError, e:
            geo_reference = Geo_reference() #Default georef object

        xllcorner = geo_reference.get_xllcorner()
        yllcorner = geo_reference.get_yllcorner()
        zone = geo_reference.get_zone()
    else:
        zone = origin[0]
        xllcorner = origin[1]
        yllcorner = origin[2]



    #FIXME: Refactor using code from file_function.statistics
    #Something like print swwstats(swwname)
    if verbose:
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        times = fid.variables['time'][:]
        print '------------------------------------------------'
        print 'Statistics of SWW file:'
        print '  Name: %s' %swwfile
        print '  Reference:'
        print '    Lower left corner: [%f, %f]'\
              %(xllcorner, yllcorner)
        print '    Start time: %f' %fid.starttime[0]
        print '  Extent:'
        print '    x [m] in [%f, %f], len(x) == %d'\
              %(min(x.flat), max(x.flat), len(x.flat))
        print '    y [m] in [%f, %f], len(y) == %d'\
              %(min(y.flat), max(y.flat), len(y.flat))
        print '    t [s] in [%f, %f], len(t) == %d'\
              %(min(times), max(times), len(times))
        print '  Quantities [SI units]:'
        for name in ['stage', 'xmomentum', 'ymomentum', 'elevation']:
            q = fid.variables[name][:].flat
            print '    %s in [%f, %f]' %(name, min(q), max(q))





    # Get quantity and reduce if applicable
    if verbose: print 'Processing quantity %s' %quantity

    # Turn NetCDF objects into Numeric arrays
    quantity_dict = {}
    for name in fid.variables.keys():
        quantity_dict[name] = fid.variables[name][:]



    # Convert quantity expression to quantities found in sww file    
    q = apply_expression_to_dictionary(quantity, quantity_dict)



    if len(q.shape) == 2:
        #q has a time component and needs to be reduced along
        #the temporal dimension
        if verbose: print 'Reducing quantity %s' %quantity
        q_reduced = zeros( number_of_points, Float )

        for k in range(number_of_points):
            q_reduced[k] = reduction( q[:,k] )

        q = q_reduced

    #Post condition: Now q has dimension: number_of_points
    assert len(q.shape) == 1
    assert q.shape[0] == number_of_points


    if verbose:
        print 'Processed values for %s are in [%f, %f]' %(quantity, min(q), max(q))


    #Create grid and update xll/yll corner and x,y

    #Relative extent
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



    if verbose: print 'Creating grid'
    ncols = int((xmax-xmin)/cellsize)+1
    nrows = int((ymax-ymin)/cellsize)+1


    #New absolute reference and coordinates
    newxllcorner = xmin+xllcorner
    newyllcorner = ymin+yllcorner

    x = x+xllcorner-newxllcorner
    y = y+yllcorner-newyllcorner

    vertex_points = concatenate ((x[:, NewAxis] ,y[:, NewAxis]), axis = 1)
    assert len(vertex_points.shape) == 2



    grid_points = zeros ( (ncols*nrows, 2), Float )


    for i in xrange(nrows):
        if format.lower() == 'asc':
            yg = i*cellsize
        else:
        #this will flip the order of the y values for ers
            yg = (nrows-i)*cellsize

        for j in xrange(ncols):
            xg = j*cellsize
            k = i*ncols + j

            grid_points[k,0] = xg
            grid_points[k,1] = yg

    #Interpolate
    #from least_squares import Interpolation
    from fit_interpolate.interpolate import Interpolate


    interp = Interpolate(vertex_points, volumes, verbose = verbose)



    #Interpolate using quantity values
    if verbose: print 'Interpolating'
    grid_values = interp.interpolate(q, grid_points).flat


    if verbose:
        print 'Interpolated values are in [%f, %f]' %(min(grid_values),
                                                      max(grid_values))

    #Assign NODATA_value to all points outside bounding polygon (from interpolation mesh)
    P = interp.mesh.get_boundary_polygon()
    outside_indices = outside_polygon(grid_points, P, closed=True)

    for i in outside_indices:
        grid_values[i] = NODATA_value




    if format.lower() == 'ers':
        # setup ERS header information
        grid_values = reshape(grid_values,(nrows, ncols))
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
        if verbose: print 'Writing %s' %demfile
        import ermapper_grids
        ermapper_grids.write_ermapper_grid(demfile, grid_values, header)

        fid.close()
    else:
        #Write to Ascii format

        #Write prj file
        prjfile = basename_out + '.prj'

        if verbose: print 'Writing %s' %prjfile
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



        if verbose: print 'Writing %s' %demfile

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

        for i in range(nrows):
            if verbose and i%((nrows+10)/10)==0:
                print 'Doing row %d of %d' %(i, nrows)

            base_index = (nrows-i-1)*ncols

            slice = grid_values[base_index:base_index+ncols]
            s = array2string(slice, max_line_width=sys.maxint)
            ascid.write(s[1:-1] + '\n')


            #print
            #for j in range(ncols):
            #    index = base_index+j#
            #    print grid_values[index],
            #    ascid.write('%f ' %grid_values[index])
            #ascid.write('\n')

        #Close
        ascid.close()
        fid.close()

#Backwards compatibility
def sww2asc(basename_in, basename_out = None,
            quantity = None,
            timestep = None,
            reduction = None,
            cellsize = 10,
            verbose = False,
            origin = None):
    print 'sww2asc will soon be obsoleted - please use sww2dem'
    sww2dem(basename_in,
            basename_out = basename_out,
            quantity = quantity,
            timestep = timestep,
            reduction = reduction,
            cellsize = cellsize,
            verbose = verbose,
            origin = origin,
        datum = 'WGS84',
        format = 'asc')

def sww2ers(basename_in, basename_out = None,
            quantity = None,
            timestep = None,
            reduction = None,
            cellsize = 10,
            verbose = False,
            origin = None,
            datum = 'WGS84'):
    print 'sww2ers will soon be obsoleted - please use sww2dem'
    sww2dem(basename_in,
            basename_out = basename_out,
            quantity = quantity,
            timestep = timestep,
            reduction = reduction,
            cellsize = cellsize,
            verbose = verbose,
            origin = origin,
            datum = datum,
            format = 'ers')
################################# END COMPATIBILITY ##############



def sww2pts(basename_in, basename_out=None,
            data_points=None,
            quantity=None,
            timestep=None,
            reduction=None,
            NODATA_value=-9999,
            verbose=False,
            origin=None):
            #datum = 'WGS84')


    """Read SWW file and convert to interpolated values at selected points

    The parameter quantity' must be the name of an existing quantity or
    an expression involving existing quantities. The default is
    'elevation'.

    if timestep (an index) is given, output quantity at that timestep

    if reduction is given use that to reduce quantity over all timesteps.

    data_points (Nx2 array) give locations of points where quantity is to be computed.
    
    """

    import sys
    from Numeric import array, Float, concatenate, NewAxis, zeros, reshape, sometrue
    from Numeric import array2string

    from utilities.polygon import inside_polygon, outside_polygon, separate_points_by_polygon
    from util import apply_expression_to_dictionary

    from geospatial_data import Geospatial_data

    if quantity is None:
        quantity = 'elevation'

    if reduction is None:
        reduction = max

    if basename_out is None:
        basename_out = basename_in + '_%s' %quantity

    swwfile = basename_in + '.sww'
    ptsfile = basename_out + '.pts'

    # Read sww file
    if verbose: print 'Reading from %s' %swwfile
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
            geo_reference = Geo_reference() #Default georef object

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
        print '------------------------------------------------'
        print 'Statistics of SWW file:'
        print '  Name: %s' %swwfile
        print '  Reference:'
        print '    Lower left corner: [%f, %f]'\
              %(xllcorner, yllcorner)
        print '    Start time: %f' %fid.starttime[0]
        print '  Extent:'
        print '    x [m] in [%f, %f], len(x) == %d'\
              %(min(x.flat), max(x.flat), len(x.flat))
        print '    y [m] in [%f, %f], len(y) == %d'\
              %(min(y.flat), max(y.flat), len(y.flat))
        print '    t [s] in [%f, %f], len(t) == %d'\
              %(min(times), max(times), len(times))
        print '  Quantities [SI units]:'
        for name in ['stage', 'xmomentum', 'ymomentum', 'elevation']:
            q = fid.variables[name][:].flat
            print '    %s in [%f, %f]' %(name, min(q), max(q))



    # Get quantity and reduce if applicable
    if verbose: print 'Processing quantity %s' %quantity

    # Turn NetCDF objects into Numeric arrays
    quantity_dict = {}
    for name in fid.variables.keys():
        quantity_dict[name] = fid.variables[name][:]



    # Convert quantity expression to quantities found in sww file    
    q = apply_expression_to_dictionary(quantity, quantity_dict)



    if len(q.shape) == 2:
        # q has a time component and needs to be reduced along
        # the temporal dimension
        if verbose: print 'Reducing quantity %s' %quantity
        q_reduced = zeros( number_of_points, Float )

        for k in range(number_of_points):
            q_reduced[k] = reduction( q[:,k] )

        q = q_reduced

    # Post condition: Now q has dimension: number_of_points
    assert len(q.shape) == 1
    assert q.shape[0] == number_of_points


    if verbose:
        print 'Processed values for %s are in [%f, %f]' %(quantity, min(q), max(q))


    # Create grid and update xll/yll corner and x,y
    vertex_points = concatenate ((x[:, NewAxis] ,y[:, NewAxis]), axis = 1)
    assert len(vertex_points.shape) == 2

    # Interpolate
    from fit_interpolate.interpolate import Interpolate
    interp = Interpolate(vertex_points, volumes, verbose = verbose)

    # Interpolate using quantity values
    if verbose: print 'Interpolating'
    interpolated_values = interp.interpolate(q, data_points).flat


    if verbose:
        print 'Interpolated values are in [%f, %f]' %(min(interpolated_values),
                                                      max(interpolated_values))

    # Assign NODATA_value to all points outside bounding polygon (from interpolation mesh)
    P = interp.mesh.get_boundary_polygon()
    outside_indices = outside_polygon(data_points, P, closed=True)

    for i in outside_indices:
        interpolated_values[i] = NODATA_value


    # Store results    
    G = Geospatial_data(data_points=data_points,
                        attributes=interpolated_values)

    G.export_points_file(ptsfile, absolute = True)

    fid.close()


def convert_dem_from_ascii2netcdf(basename_in, basename_out = None,
                                  use_cache = False,
                                  verbose = False):
    """Read Digitial Elevation model from the following ASCII format (.asc)

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
                       dependencies = [basename_in + '.asc'],
                       verbose = verbose)

    else:
        result = apply(_convert_dem_from_ascii2netcdf, [basename_in], kwargs)

    return result






def _convert_dem_from_ascii2netcdf(basename_in, basename_out = None,
                                  verbose = False):
    """Read Digitial Elevation model from the following ASCII format (.asc)

    Internal function. See public function convert_dem_from_ascii2netcdf for details.
    """

    import os
    from Scientific.IO.NetCDF import NetCDFFile
    from Numeric import Float, array

    #root, ext = os.path.splitext(basename_in)
    root = basename_in

    ###########################################
    # Read Meta data
    if verbose: print 'Reading METADATA from %s' %root + '.prj'
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

    #print false_easting, false_northing, zone, datum


    ###########################################
    #Read DEM data

    datafile = open(basename_in + '.asc')

    if verbose: print 'Reading DEM from %s' %(basename_in + '.asc')
    lines = datafile.readlines()
    datafile.close()

    if verbose: print 'Got', len(lines), ' lines'

    ncols = int(lines[0].split()[1].strip())
    nrows = int(lines[1].split()[1].strip())
    xllcorner = float(lines[2].split()[1].strip())
    yllcorner = float(lines[3].split()[1].strip())
    cellsize = float(lines[4].split()[1].strip())
    NODATA_value = int(lines[5].split()[1].strip())

    assert len(lines) == nrows + 6


    ##########################################


    if basename_out == None:
        netcdfname = root + '.dem'
    else:
        netcdfname = basename_out + '.dem'

    if verbose: print 'Store to NetCDF file %s' %netcdfname
    # NetCDF file definition
    fid = NetCDFFile(netcdfname, 'w')

    #Create new file
    fid.institution = 'Geoscience Australia'
    fid.description = 'NetCDF DEM format for compact and portable storage ' +\
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
    fid.createVariable('elevation', Float, ('number_of_rows',
                                            'number_of_columns'))

    # Get handles to the variables
    elevation = fid.variables['elevation']

    #Store data
    n = len(lines[6:])
    for i, line in enumerate(lines[6:]):
        fields = line.split()
        if verbose and i%((n+10)/10)==0:
            print 'Processing row %d of %d' %(i, nrows)

        elevation[i, :] = array([float(x) for x in fields])

    fid.close()





def ferret2sww(basename_in, basename_out = None,
               verbose = False,
               minlat = None, maxlat = None,
               minlon = None, maxlon = None,
               mint = None, maxt = None, mean_stage = 0,
               origin = None, zscale = 1,
               fail_on_NaN = True,
               NaN_filler = 0,
               elevation = None,
               inverted_bathymetry = False
               ): #FIXME: Bathymetry should be obtained
                                  #from MOST somehow.
                                  #Alternatively from elsewhere
                                  #or, as a last resort,
                                  #specified here.
                                  #The value of -100 will work
                                  #for the Wollongong tsunami
                                  #scenario but is very hacky
    """Convert MOST and 'Ferret' NetCDF format for wave propagation to
    sww format native to pyvolution.

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
    from Numeric import Float, Int, Int32, searchsorted, zeros, array
    from Numeric import allclose, around

    precision = Float

    msg = 'Must use latitudes and longitudes for minlat, maxlon etc'

    if minlat != None:
        assert -90 < minlat < 90 , msg
    if maxlat != None:
        assert -90 < maxlat < 90 , msg
    if minlon != None:
        assert -180 < minlon < 180 , msg
    if maxlon != None:
        assert -180 < maxlon < 180 , msg


    #Get NetCDF data
    if verbose: print 'Reading files %s_*.nc' %basename_in
    file_h = NetCDFFile(basename_in + '_ha.nc', 'r') #Wave amplitude (cm)
    file_u = NetCDFFile(basename_in + '_ua.nc', 'r') #Velocity (x) (cm/s)
    file_v = NetCDFFile(basename_in + '_va.nc', 'r') #Velocity (y) (cm/s)
    file_e = NetCDFFile(basename_in + '_e.nc', 'r')  #Elevation (z) (m)

    if basename_out is None:
        swwname = basename_in + '.sww'
    else:
        swwname = basename_out + '.sww'

#    Get dimensions of file_h
    for dimension in file_h.dimensions.keys():
        if dimension[:3] == 'LON':
            dim_h_longitude = dimension
        if dimension[:3] == 'LAT':
            dim_h_latitude = dimension
        if dimension[:4] == 'TIME':
            dim_h_time = dimension

#    print 'long:', dim_h_longitude
#    print 'lats:', dim_h_latitude
#    print 'times:', dim_h_time

    times = file_h.variables[dim_h_time]
    latitudes = file_h.variables[dim_h_latitude]
    longitudes = file_h.variables[dim_h_longitude]
    
# get dimensions for file_e
    for dimension in file_e.dimensions.keys():
        if dimension[:3] == 'LON':
            dim_e_longitude = dimension
        if dimension[:3] == 'LAT':
            dim_e_latitude = dimension
        if dimension[:4] == 'TIME':
            dim_e_time = dimension

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


    #Precision used by most for lat/lon is 4 or 5 decimals
    e_lat = around(file_e.variables[dim_e_latitude][:], 5)
    e_lon = around(file_e.variables[dim_e_longitude][:], 5)

    #Check that files are compatible
    assert allclose(latitudes, file_u.variables[dim_u_latitude])
    assert allclose(latitudes, file_v.variables[dim_v_latitude])
    assert allclose(latitudes, e_lat)

    assert allclose(longitudes, file_u.variables[dim_u_longitude])
    assert allclose(longitudes, file_v.variables[dim_v_longitude])
    assert allclose(longitudes, e_lon)

    if mint == None:
        jmin = 0
    else:
        jmin = searchsorted(times, mint)

    if maxt == None:
        jmax=len(times)
    else:
        jmax = searchsorted(times, maxt)

    if minlat == None:
        kmin=0
    else:
        kmin = searchsorted(latitudes, minlat)

    if maxlat == None:
        kmax = len(latitudes)
    else:
        kmax = searchsorted(latitudes, maxlat)

    if minlon == None:
        lmin=0
    else:
        lmin = searchsorted(longitudes, minlon)

    if maxlon == None:
        lmax = len(longitudes)
    else:
        lmax = searchsorted(longitudes, maxlon)

#    print' j', jmin, jmax
    times = times[jmin:jmax]
    latitudes = latitudes[kmin:kmax]
    longitudes = longitudes[lmin:lmax]


    if verbose: print 'cropping'
    zname = 'ELEVATION'

    amplitudes = file_h.variables['HA'][jmin:jmax, kmin:kmax, lmin:lmax]
    uspeed = file_u.variables['UA'][jmin:jmax, kmin:kmax, lmin:lmax] #Lon
    vspeed = file_v.variables['VA'][jmin:jmax, kmin:kmax, lmin:lmax] #Lat
    elevations = file_e.variables[zname][kmin:kmax, lmin:lmax]

    #    if latitudes2[0]==latitudes[0] and latitudes2[-1]==latitudes[-1]:
    #        elevations = file_e.variables['ELEVATION'][kmin:kmax, lmin:lmax]
    #    elif latitudes2[0]==latitudes[-1] and latitudes2[-1]==latitudes[0]:
    #        from Numeric import asarray
    #        elevations=elevations.tolist()
    #        elevations.reverse()
    #        elevations=asarray(elevations)
    #    else:
    #        from Numeric import asarray
    #        elevations=elevations.tolist()
    #        elevations.reverse()
    #        elevations=asarray(elevations)
    #        'print hmmm'



    #Get missing values
    nan_ha = file_h.variables['HA'].missing_value[0]
    nan_ua = file_u.variables['UA'].missing_value[0]
    nan_va = file_v.variables['VA'].missing_value[0]
    if hasattr(file_e.variables[zname],'missing_value'):
        nan_e  = file_e.variables[zname].missing_value[0]
    else:
        nan_e = None

    #Cleanup
    from Numeric import sometrue

    missing = (amplitudes == nan_ha)
    if sometrue (missing):
        if fail_on_NaN:
            msg = 'NetCDFFile %s contains missing values'\
                  %(basename_in+'_ha.nc')
            raise DataMissingValuesError, msg
        else:
            amplitudes = amplitudes*(missing==0) + missing*NaN_filler

    missing = (uspeed == nan_ua)
    if sometrue (missing):
        if fail_on_NaN:
            msg = 'NetCDFFile %s contains missing values'\
                  %(basename_in+'_ua.nc')
            raise DataMissingValuesError, msg
        else:
            uspeed = uspeed*(missing==0) + missing*NaN_filler

    missing = (vspeed == nan_va)
    if sometrue (missing):
        if fail_on_NaN:
            msg = 'NetCDFFile %s contains missing values'\
                  %(basename_in+'_va.nc')
            raise DataMissingValuesError, msg
        else:
            vspeed = vspeed*(missing==0) + missing*NaN_filler


    missing = (elevations == nan_e)
    if sometrue (missing):
        if fail_on_NaN:
            msg = 'NetCDFFile %s contains missing values'\
                  %(basename_in+'_e.nc')
            raise DataMissingValuesError, msg
        else:
            elevations = elevations*(missing==0) + missing*NaN_filler

    #######



    number_of_times = times.shape[0]
    number_of_latitudes = latitudes.shape[0]
    number_of_longitudes = longitudes.shape[0]

    assert amplitudes.shape[0] == number_of_times
    assert amplitudes.shape[1] == number_of_latitudes
    assert amplitudes.shape[2] == number_of_longitudes

    if verbose:
        print '------------------------------------------------'
        print 'Statistics:'
        print '  Extent (lat/lon):'
        print '    lat in [%f, %f], len(lat) == %d'\
              %(min(latitudes.flat), max(latitudes.flat),
                len(latitudes.flat))
        print '    lon in [%f, %f], len(lon) == %d'\
              %(min(longitudes.flat), max(longitudes.flat),
                len(longitudes.flat))
        print '    t in [%f, %f], len(t) == %d'\
              %(min(times.flat), max(times.flat), len(times.flat))

        q = amplitudes.flat
        name = 'Amplitudes (ha) [cm]'
        print '  %s in [%f, %f]' %(name, min(q), max(q))

        q = uspeed.flat
        name = 'Speeds (ua) [cm/s]'
        print '  %s in [%f, %f]' %(name, min(q), max(q))

        q = vspeed.flat
        name = 'Speeds (va) [cm/s]'
        print '  %s in [%f, %f]' %(name, min(q), max(q))

        q = elevations.flat
        name = 'Elevations (e) [m]'
        print '  %s in [%f, %f]' %(name, min(q), max(q))


    #print number_of_latitudes, number_of_longitudes
    number_of_points = number_of_latitudes*number_of_longitudes
    number_of_volumes = (number_of_latitudes-1)*(number_of_longitudes-1)*2


    file_h.close()
    file_u.close()
    file_v.close()
    file_e.close()


    # NetCDF file definition
    outfile = NetCDFFile(swwname, 'w')

    #Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = 'Converted from Ferret files: %s, %s, %s, %s'\
                          %(basename_in + '_ha.nc',
                            basename_in + '_ua.nc',
                            basename_in + '_va.nc',
                            basename_in + '_e.nc')


    #For sww compatibility
    outfile.smoothing = 'Yes'
    outfile.order = 1

    #Start time in seconds since the epoch (midnight 1/1/1970)
    outfile.starttime = starttime = times[0]
    times = times - starttime  #Store relative times

    # dimension definitions
    outfile.createDimension('number_of_volumes', number_of_volumes)

    outfile.createDimension('number_of_vertices', 3)
    outfile.createDimension('number_of_points', number_of_points)


    #outfile.createDimension('number_of_timesteps', len(times))
    outfile.createDimension('number_of_timesteps', len(times))

    # variable definitions
    outfile.createVariable('x', precision, ('number_of_points',))
    outfile.createVariable('y', precision, ('number_of_points',))
    outfile.createVariable('elevation', precision, ('number_of_points',))

    #FIXME: Backwards compatibility
    outfile.createVariable('z', precision, ('number_of_points',))
    #################################

    outfile.createVariable('volumes', Int, ('number_of_volumes',
                                            'number_of_vertices'))

    outfile.createVariable('time', precision,
                           ('number_of_timesteps',))

    outfile.createVariable('stage', precision,
                           ('number_of_timesteps',
                            'number_of_points'))

    outfile.createVariable('xmomentum', precision,
                           ('number_of_timesteps',
                            'number_of_points'))

    outfile.createVariable('ymomentum', precision,
                           ('number_of_timesteps',
                            'number_of_points'))


    #Store
    from coordinate_transforms.redfearn import redfearn
    x = zeros(number_of_points, Float)  #Easting
    y = zeros(number_of_points, Float)  #Northing


    if verbose: print 'Making triangular grid'
    #Check zone boundaries
    refzone, _, _ = redfearn(latitudes[0],longitudes[0])

    vertices = {}
    i = 0
    for k, lat in enumerate(latitudes):       #Y direction
        for l, lon in enumerate(longitudes):  #X direction

            vertices[l,k] = i

            zone, easting, northing = redfearn(lat,lon)

            msg = 'Zone boundary crossed at longitude =', lon
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

            volumes.append([v1,v2,v3]) #Upper element
            volumes.append([v4,v3,v2]) #Lower element

    volumes = array(volumes)

    if origin == None:
        zone = refzone
        xllcorner = min(x)
        yllcorner = min(y)
    else:
        zone = origin[0]
        xllcorner = origin[1]
        yllcorner = origin[2]


    outfile.xllcorner = xllcorner
    outfile.yllcorner = yllcorner
    outfile.zone = zone


    if elevation is not None:
        z = elevation
    else:
        if inverted_bathymetry:
            z = -1*elevations
        else:
            z = elevations
    #FIXME: z should be obtained from MOST and passed in here

    from Numeric import resize
    z = resize(z,outfile.variables['z'][:].shape)
    outfile.variables['x'][:] = x - xllcorner
    outfile.variables['y'][:] = y - yllcorner
    outfile.variables['z'][:] = z             #FIXME HACK
    outfile.variables['elevation'][:] = z
    outfile.variables['time'][:] = times   #Store time relative
    outfile.variables['volumes'][:] = volumes.astype(Int32) #On Opteron 64



    #Time stepping
    stage = outfile.variables['stage']
    xmomentum = outfile.variables['xmomentum']
    ymomentum = outfile.variables['ymomentum']

    if verbose: print 'Converting quantities'
    n = len(times)
    for j in range(n):
        if verbose and j%((n+10)/10)==0: print '  Doing %d of %d' %(j, n)
        i = 0
        for k in range(number_of_latitudes):      #Y direction
            for l in range(number_of_longitudes): #X direction
                w = zscale*amplitudes[j,k,l]/100 + mean_stage
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
        times = outfile.variables['time'][:]
        print '------------------------------------------------'
        print 'Statistics of output file:'
        print '  Name: %s' %swwname
        print '  Reference:'
        print '    Lower left corner: [%f, %f]'\
              %(xllcorner, yllcorner)
        print '    Start time: %f' %starttime
        print '  Extent:'
        print '    x [m] in [%f, %f], len(x) == %d'\
              %(min(x.flat), max(x.flat), len(x.flat))
        print '    y [m] in [%f, %f], len(y) == %d'\
              %(min(y.flat), max(y.flat), len(y.flat))
        print '    t [s] in [%f, %f], len(t) == %d'\
              %(min(times), max(times), len(times))
        print '  Quantities [SI units]:'
        for name in ['stage', 'xmomentum', 'ymomentum', 'elevation']:
            q = outfile.variables[name][:].flat
            print '    %s in [%f, %f]' %(name, min(q), max(q))



    outfile.close()





def timefile2netcdf(filename, quantity_names = None):
    """Template for converting typical text files with time series to
    NetCDF tms file.


    The file format is assumed to be either two fields separated by a comma:

        time [DD/MM/YY hh:mm:ss], value0 value1 value2 ...

    E.g

      31/08/04 00:00:00, 1.328223 0 0
      31/08/04 00:15:00, 1.292912 0 0

    will provide a time dependent function f(t) with three attributes

    filename is assumed to be the rootname with extenisons .txt and .sww
    """

    import time, calendar
    from Numeric import array
    from config import time_format
    from utilities.numerical_tools import ensure_numeric



    fid = open(filename + '.txt')
    line = fid.readline()
    fid.close()

    fields = line.split(',')
    msg = 'File %s must have the format date, value0 value1 value2 ...'
    assert len(fields) == 2, msg

    try:
        starttime = calendar.timegm(time.strptime(fields[0], time_format))
    except ValueError:
        msg = 'First field in file %s must be' %filename
        msg += ' date-time with format %s.\n' %time_format
        msg += 'I got %s instead.' %fields[0]
        raise DataTimeError, msg


    #Split values
    values = []
    for value in fields[1].split():
        values.append(float(value))

    q = ensure_numeric(values)

    msg = 'ERROR: File must contain at least one independent value'
    assert len(q.shape) == 1, msg



    #Read times proper
    from Numeric import zeros, Float, alltrue
    from config import time_format
    import time, calendar

    fid = open(filename + '.txt')
    lines = fid.readlines()
    fid.close()

    N = len(lines)
    d = len(q)

    T = zeros(N, Float)       #Time
    Q = zeros((N, d), Float)  #Values

    for i, line in enumerate(lines):
        fields = line.split(',')
        realtime = calendar.timegm(time.strptime(fields[0], time_format))

        T[i] = realtime - starttime

        for j, value in enumerate(fields[1].split()):
            Q[i, j] = float(value)

    msg = 'File %s must list time as a monotonuosly ' %filename
    msg += 'increasing sequence'
    assert alltrue( T[1:] - T[:-1] > 0 ), msg


    #Create NetCDF file
    from Scientific.IO.NetCDF import NetCDFFile

    fid = NetCDFFile(filename + '.tms', 'w')


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

    fid.createVariable('time', Float, ('number_of_timesteps',))

    fid.variables['time'][:] = T

    for i in range(Q.shape[1]):
        try:
            name = quantity_names[i]
        except:
            name = 'Attribute%d'%i

        fid.createVariable(name, Float, ('number_of_timesteps',))
        fid.variables[name][:] = Q[:,i]

    fid.close()


def extent_sww(file_name):
    """
    Read in an sww file.

    Input;
    file_name - the sww file

    Output;
    z - Vector of bed elevation
    volumes - Array.  Each row has 3 values, representing
    the vertices that define the volume
    time - Vector of the times where there is stage information
    stage - array with respect to time and vertices (x,y)
    """


    from Scientific.IO.NetCDF import NetCDFFile

    #Check contents
    #Get NetCDF
    fid = NetCDFFile(file_name, 'r')

    # Get the variables
    x = fid.variables['x'][:]
    y = fid.variables['y'][:]
    stage = fid.variables['stage'][:]
    #print "stage",stage
    #print "stage.shap",stage.shape
    #print "min(stage.flat), mpythonax(stage.flat)",min(stage.flat), max(stage.flat)
    #print "min(stage)",min(stage)

    fid.close()

    return [min(x),max(x),min(y),max(y),min(stage.flat),max(stage.flat)]


def sww2domain(filename,boundary=None,t=None,\
               fail_if_NaN=True,NaN_filler=0\
               ,verbose = False,very_verbose = False):
    """
    Usage: domain = sww2domain('file.sww',t=time (default = last time in file))

    Boundary is not recommended if domain.smooth is not selected, as it
    uses unique coordinates, but not unique boundaries. This means that
    the boundary file will not be compatable with the coordinates, and will
    give a different final boundary, or crash.
    """
    NaN=9.969209968386869e+036
    #initialise NaN.

    from Scientific.IO.NetCDF import NetCDFFile
    from shallow_water import Domain
    from Numeric import asarray, transpose, resize

    if verbose: print 'Reading from ', filename
    fid = NetCDFFile(filename, 'r')    #Open existing file for read
    time = fid.variables['time']       #Timesteps
    if t is None:
        t = time[-1]
    time_interp = get_time_interp(time,t)

    # Get the variables as Numeric arrays
    x = fid.variables['x'][:]             #x-coordinates of vertices
    y = fid.variables['y'][:]             #y-coordinates of vertices
    elevation = fid.variables['elevation']     #Elevation
    stage = fid.variables['stage']     #Water level
    xmomentum = fid.variables['xmomentum']   #Momentum in the x-direction
    ymomentum = fid.variables['ymomentum']   #Momentum in the y-direction

    starttime = fid.starttime[0]
    volumes = fid.variables['volumes'][:] #Connectivity
    coordinates=transpose(asarray([x.tolist(),y.tolist()]))

    conserved_quantities = []
    interpolated_quantities = {}
    other_quantities = []

    # get geo_reference
    #sww files don't have to have a geo_ref
    try:
        geo_reference = Geo_reference(NetCDFObject=fid)
    except: #AttributeError, e:
        geo_reference = None

    if verbose: print '    getting quantities'
    for quantity in fid.variables.keys():
        dimensions = fid.variables[quantity].dimensions
        if 'number_of_timesteps' in dimensions:
            conserved_quantities.append(quantity)
            interpolated_quantities[quantity]=\
                  interpolated_quantity(fid.variables[quantity][:],time_interp)
        else: other_quantities.append(quantity)

    other_quantities.remove('x')
    other_quantities.remove('y')
    other_quantities.remove('z')
    other_quantities.remove('volumes')

    conserved_quantities.remove('time')

    if verbose: print '    building domain'
    #    From domain.Domain:
    #    domain = Domain(coordinates, volumes,\
    #                    conserved_quantities = conserved_quantities,\
    #                    other_quantities = other_quantities,zone=zone,\
    #                    xllcorner=xllcorner, yllcorner=yllcorner)

    #   From shallow_water.Domain:
    coordinates=coordinates.tolist()
    volumes=volumes.tolist()
    #FIXME:should this be in mesh?(peter row)
    if fid.smoothing == 'Yes': unique = False
    else: unique = True
    if unique:
        coordinates,volumes,boundary=weed(coordinates,volumes,boundary)


    try:
        domain = Domain(coordinates, volumes, boundary)
    except AssertionError, e:
        fid.close()
        msg = 'Domain could not be created: %s. Perhaps use "fail_if_NaN=False and NaN_filler = ..."' %e
        raise DataDomainError, msg

    if not boundary is None:
        domain.boundary = boundary

    domain.geo_reference = geo_reference

    domain.starttime=float(starttime)+float(t)
    domain.time=0.0

    for quantity in other_quantities:
        try:
            NaN = fid.variables[quantity].missing_value
        except:
            pass #quantity has no missing_value number
        X = fid.variables[quantity][:]
        if very_verbose:
            print '       ',quantity
            print '        NaN =',NaN
            print '        max(X)'
            print '       ',max(X)
            print '        max(X)==NaN'
            print '       ',max(X)==NaN
            print ''
        if (max(X)==NaN) or (min(X)==NaN):
            if fail_if_NaN:
                msg = 'quantity "%s" contains no_data entry'%quantity
                raise DataMissingValuesError, msg
            else:
                data = (X<>NaN)
                X = (X*data)+(data==0)*NaN_filler
        if unique:
            X = resize(X,(len(X)/3,3))
        domain.set_quantity(quantity,X)
    #
    for quantity in conserved_quantities:
        try:
            NaN = fid.variables[quantity].missing_value
        except:
            pass #quantity has no missing_value number
        X = interpolated_quantities[quantity]
        if very_verbose:
            print '       ',quantity
            print '        NaN =',NaN
            print '        max(X)'
            print '       ',max(X)
            print '        max(X)==NaN'
            print '       ',max(X)==NaN
            print ''
        if (max(X)==NaN) or (min(X)==NaN):
            if fail_if_NaN:
                msg = 'quantity "%s" contains no_data entry'%quantity
                raise DataMissingValuesError, msg
            else:
                data = (X<>NaN)
                X = (X*data)+(data==0)*NaN_filler
        if unique:
            X = resize(X,(X.shape[0]/3,3))
        domain.set_quantity(quantity,X)

    fid.close()
    return domain

def interpolated_quantity(saved_quantity,time_interp):

    #given an index and ratio, interpolate quantity with respect to time.
    index,ratio = time_interp
    Q = saved_quantity
    if ratio > 0:
        q = (1-ratio)*Q[index]+ ratio*Q[index+1]
    else:
        q = Q[index]
    #Return vector of interpolated values
    return q

def get_time_interp(time,t=None):
    #Finds the ratio and index for time interpolation.
    #It is borrowed from previous pyvolution code.
    if t is None:
        t=time[-1]
        index = -1
        ratio = 0.
    else:
        T = time
        tau = t
        index=0
        msg = 'Time interval derived from file %s [%s:%s]'\
            %('FIXMEfilename', T[0], T[-1])
        msg += ' does not match model time: %s' %tau
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
    return (index,ratio)


def weed(coordinates,volumes,boundary = None):
    if type(coordinates)=='array':
        coordinates = coordinates.tolist()
    if type(volumes)=='array':
        volumes = volumes.tolist()

    unique = False
    point_dict = {}
    same_point = {}
    for i in range(len(coordinates)):
        point = tuple(coordinates[i])
        if point_dict.has_key(point):
            unique = True
            same_point[i]=point
            #to change all point i references to point j
        else:
            point_dict[point]=i
            same_point[i]=point

    coordinates = []
    i = 0
    for point in point_dict.keys():
        point = tuple(point)
        coordinates.append(list(point))
        point_dict[point]=i
        i+=1


    for volume in volumes:
        for i in range(len(volume)):
            index = volume[i]
            if index>-1:
                volume[i]=point_dict[same_point[index]]

    new_boundary = {}
    if not boundary is None:
        for segment in boundary.keys():
            point0 = point_dict[same_point[segment[0]]]
            point1 = point_dict[same_point[segment[1]]]
            label = boundary[segment]
            #FIXME should the bounday attributes be concaterated
            #('exterior, pond') or replaced ('pond')(peter row)

            if new_boundary.has_key((point0,point1)):
                new_boundary[(point0,point1)]=new_boundary[(point0,point1)]#\
                                              #+','+label

            elif new_boundary.has_key((point1,point0)):
                new_boundary[(point1,point0)]=new_boundary[(point1,point0)]#\
                                              #+','+label
            else: new_boundary[(point0,point1)]=label

        boundary = new_boundary

    return coordinates,volumes,boundary


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
    from Numeric import Float, zeros, sum, reshape, equal

    root = basename_in
    inname = root + '.dem'

    #Open existing netcdf file to read
    infile = NetCDFFile(inname, 'r')
    if verbose: print 'Reading DEM from %s' %inname

    #Read metadata
    ncols = infile.ncols[0]
    nrows = infile.nrows[0]
    xllcorner = infile.xllcorner[0]
    yllcorner = infile.yllcorner[0]
    cellsize = infile.cellsize[0]
    NODATA_value = infile.NODATA_value[0]
    zone = infile.zone[0]
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

    if verbose: print 'Write decimated NetCDF file to %s' %outname

    #Determine some dimensions for decimated grid
    (nrows_stencil, ncols_stencil) = stencil.shape
    x_offset = ncols_stencil / 2
    y_offset = nrows_stencil / 2
    cellsize_ratio = int(cellsize_new / cellsize)
    ncols_new = 1 + (ncols - ncols_stencil) / cellsize_ratio
    nrows_new = 1 + (nrows - nrows_stencil) / cellsize_ratio

    #Open netcdf file for output
    outfile = NetCDFFile(outname, 'w')

    #Create new file
    outfile.institution = 'Geoscience Australia'
    outfile.description = 'NetCDF DEM format for compact and portable storage ' +\
                           'of spatial point data'
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
    outfile.createDimension('number_of_points', nrows_new*ncols_new)

    # variable definition
    outfile.createVariable('elevation', Float, ('number_of_points',))

    # Get handle to the variable
    elevation = outfile.variables['elevation']

    dem_elevation_r = reshape(dem_elevation, (nrows, ncols))

    #Store data
    global_index = 0
    for i in range(nrows_new):
        if verbose: print 'Processing row %d of %d' %(i, nrows_new)
        lower_index = global_index
        telev =  zeros(ncols_new, Float)
        local_index = 0
        trow = i * cellsize_ratio

        for j in range(ncols_new):
            tcol = j * cellsize_ratio
            tmp = dem_elevation_r[trow:trow+nrows_stencil, tcol:tcol+ncols_stencil]

            #if dem contains 1 or more NODATA_values set value in
            #decimated dem to NODATA_value, else compute decimated
            #value using stencil
            if sum(sum(equal(tmp, NODATA_value))) > 0:
                telev[local_index] = NODATA_value
            else:
                telev[local_index] = sum(sum(tmp * stencil))

            global_index += 1
            local_index += 1

        upper_index = global_index

        elevation[lower_index:upper_index] = telev

    assert global_index == nrows_new*ncols_new, 'index not equal to number of points'

    infile.close()
    outfile.close()




def tsh2sww(filename, verbose=False): #test_tsh2sww
    """
    to check if a tsh/msh file 'looks' good.
    """

    from shallow_water import Domain
    from pmesh2domain import pmesh_to_domain_instance
    import time, os
    from data_manager import get_dataobject
    from os import sep, path
    from utilities.numerical_tools import mean

    if verbose == True:print 'Creating domain from', filename
    domain = pmesh_to_domain_instance(filename, Domain)
    if verbose == True:print "Number of triangles = ", len(domain)

    domain.smooth = True
    domain.format = 'sww'   #Native netcdf visualisation format
    file_path, filename = path.split(filename)
    filename, ext = path.splitext(filename)
    domain.filename = filename
    domain.reduction = mean
    if verbose == True:print "file_path",file_path
    if file_path == "":file_path = "."
    domain.set_datadir(file_path)

    if verbose == True:
        print "Output written to " + domain.get_datadir() + sep + \
              domain.filename + "." + domain.format
    sww = get_dataobject(domain)
    sww.store_connectivity()
    sww.store_timestep('stage')


def asc_csiro2sww(bath_dir,
                  elevation_dir,
                  ucur_dir,
                  vcur_dir,
                  sww_file,
                  minlat = None, maxlat = None,
                  minlon = None, maxlon = None,
                  zscale=1,
                  mean_stage = 0,
                  fail_on_NaN = True,
                  elevation_NaN_filler = 0,
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

    from coordinate_transforms.redfearn import redfearn

    precision = Float # So if we want to change the precision its done here

    # go in to the bath dir and load the only file,
    bath_files = os.listdir(bath_dir)
    #print "bath_files",bath_files

    #fixme: make more general?
    bath_file = bath_files[0]
    bath_dir_file =  bath_dir + os.sep + bath_file
    bath_metadata,bath_grid =  _read_asc(bath_dir_file)
    #print "bath_metadata",bath_metadata
    #print "bath_grid",bath_grid

    #Use the date.time of the bath file as a basis for
    #the start time for other files
    base_start = bath_file[-12:]

    #go into the elevation dir and load the 000 file
    elevation_dir_file = elevation_dir  + os.sep + elevation_prefix \
                         + base_start
    #elevation_metadata,elevation_grid =  _read_asc(elevation_dir_file)
    #print "elevation_dir_file",elevation_dir_file
    #print "elevation_grid", elevation_grid

    elevation_files = os.listdir(elevation_dir)
    ucur_files = os.listdir(ucur_dir)
    vcur_files = os.listdir(vcur_dir)

    # the first elevation file should be the
    # file with the same base name as the bath data
    #print "elevation_files[0]",elevation_files[0]
    #print "'el' + base_start",'el' + base_start
    assert elevation_files[0] == 'el' + base_start

    #assert bath_metadata == elevation_metadata



    number_of_latitudes = bath_grid.shape[0]
    number_of_longitudes = bath_grid.shape[1]
    #number_of_times = len(os.listdir(elevation_dir))
    #number_of_points = number_of_latitudes*number_of_longitudes
    number_of_volumes = (number_of_latitudes-1)*(number_of_longitudes-1)*2

    longitudes = [bath_metadata['xllcorner']+x*bath_metadata['cellsize'] \
                  for x in range(number_of_longitudes)]
    latitudes = [bath_metadata['yllcorner']+y*bath_metadata['cellsize'] \
                 for y in range(number_of_latitudes)]

     # reverse order of lat, so the fist lat represents the first grid row
    latitudes.reverse()

    #print "latitudes - before _get_min_max_indexes",latitudes
    kmin, kmax, lmin, lmax = _get_min_max_indexes(latitudes,longitudes,
                                                 minlat=minlat, maxlat=maxlat,
                                                 minlon=minlon, maxlon=maxlon)


    bath_grid = bath_grid[kmin:kmax,lmin:lmax]
    #print "bath_grid",bath_grid
    latitudes = latitudes[kmin:kmax]
    longitudes = longitudes[lmin:lmax]
    number_of_latitudes = len(latitudes)
    number_of_longitudes = len(longitudes)
    number_of_times = len(os.listdir(elevation_dir))
    number_of_points = number_of_latitudes*number_of_longitudes
    number_of_volumes = (number_of_latitudes-1)*(number_of_longitudes-1)*2
    #print "number_of_points",number_of_points

    #Work out the times
    if len(elevation_files) > 1:
        # Assume: The time period is less than 24hrs.
        time_period = (int(elevation_files[1][-3:]) - \
                      int(elevation_files[0][-3:]))*60*60
        times = [x*time_period for x in range(len(elevation_files))]
    else:
        times = [0.0]
    #print "times", times
    #print "number_of_latitudes", number_of_latitudes
    #print "number_of_longitudes", number_of_longitudes
    #print "number_of_times", number_of_times
    #print "latitudes", latitudes
    #print "longitudes", longitudes


    if verbose:
        print '------------------------------------------------'
        print 'Statistics:'
        print '  Extent (lat/lon):'
        print '    lat in [%f, %f], len(lat) == %d'\
              %(min(latitudes), max(latitudes),
                len(latitudes))
        print '    lon in [%f, %f], len(lon) == %d'\
              %(min(longitudes), max(longitudes),
                len(longitudes))
        print '    t in [%f, %f], len(t) == %d'\
              %(min(times), max(times), len(times))

    ######### WRITE THE SWW FILE #############
    # NetCDF file definition
    outfile = NetCDFFile(sww_file, 'w')

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
    outfile.createVariable('z', precision, ('number_of_points',))
    #################################

    outfile.createVariable('volumes', Int, ('number_of_volumes',
                                            'number_of_vertices'))

    outfile.createVariable('time', precision,
                           ('number_of_timesteps',))

    outfile.createVariable('stage', precision,
                           ('number_of_timesteps',
                            'number_of_points'))

    outfile.createVariable('xmomentum', precision,
                           ('number_of_timesteps',
                            'number_of_points'))

    outfile.createVariable('ymomentum', precision,
                           ('number_of_timesteps',
                            'number_of_points'))

    #Store
    from coordinate_transforms.redfearn import redfearn
    x = zeros(number_of_points, Float)  #Easting
    y = zeros(number_of_points, Float)  #Northing

    if verbose: print 'Making triangular grid'
    #Get zone of 1st point.
    refzone, _, _ = redfearn(latitudes[0],longitudes[0])

    vertices = {}
    i = 0
    for k, lat in enumerate(latitudes):
        for l, lon in enumerate(longitudes):

            vertices[l,k] = i

            zone, easting, northing = redfearn(lat,lon)

            msg = 'Zone boundary crossed at longitude =', lon
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

    volumes = array(volumes)

    geo_ref = Geo_reference(refzone,min(x),min(y))
    geo_ref.write_NetCDF(outfile)

    # This will put the geo ref in the middle
    #geo_ref = Geo_reference(refzone,(max(x)+min(x))/2.0,(max(x)+min(y))/2.)


    if verbose:
        print '------------------------------------------------'
        print 'More Statistics:'
        print '  Extent (/lon):'
        print '    x in [%f, %f], len(lat) == %d'\
              %(min(x), max(x),
                len(x))
        print '    y in [%f, %f], len(lon) == %d'\
              %(min(y), max(y),
                len(y))
        print 'geo_ref: ',geo_ref

    z = resize(bath_grid,outfile.variables['z'][:].shape)
    outfile.variables['x'][:] = x - geo_ref.get_xllcorner()
    outfile.variables['y'][:] = y - geo_ref.get_yllcorner()
    outfile.variables['z'][:] = z
    outfile.variables['elevation'][:] = z  #FIXME HACK
    outfile.variables['volumes'][:] = volumes.astype(Int32) #On Opteron 64

    # do this to create an ok sww file.
    #outfile.variables['time'][:] = [0]   #Store time relative
    #outfile.variables['stage'] = z
    # put the next line up in the code after outfile.order = 1
    #number_of_times = 1

    stage = outfile.variables['stage']
    xmomentum = outfile.variables['xmomentum']
    ymomentum = outfile.variables['ymomentum']

    outfile.variables['time'][:] = times   #Store time relative

    if verbose: print 'Converting quantities'
    n = number_of_times
    for j in range(number_of_times):
        # load in files
        elevation_meta, elevation_grid = \
           _read_asc(elevation_dir + os.sep + elevation_files[j])

        _, u_momentum_grid =  _read_asc(ucur_dir + os.sep + ucur_files[j])
        _, v_momentum_grid =  _read_asc(vcur_dir + os.sep + vcur_files[j])

        #print "elevation_grid",elevation_grid
        #cut matrix to desired size
        elevation_grid = elevation_grid[kmin:kmax,lmin:lmax]
        u_momentum_grid = u_momentum_grid[kmin:kmax,lmin:lmax]
        v_momentum_grid = v_momentum_grid[kmin:kmax,lmin:lmax]
        #print "elevation_grid",elevation_grid
        # handle missing values
        missing = (elevation_grid == elevation_meta['NODATA_value'])
        if sometrue (missing):
            if fail_on_NaN:
                msg = 'File %s contains missing values'\
                      %(elevation_files[j])
                raise DataMissingValuesError, msg
            else:
                elevation_grid = elevation_grid*(missing==0) + \
                                 missing*elevation_NaN_filler


        if verbose and j%((n+10)/10)==0: print '  Doing %d of %d' %(j, n)
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

def _get_min_max_indexes(latitudes,longitudes,
                        minlat=None, maxlat=None,
                        minlon=None, maxlon=None):
    """
    return max, min indexes (for slicing) of the lat and long arrays to cover the area
    specified with min/max lat/long

    Think of the latitudes and longitudes describing a 2d surface.
    The area returned is, if possible, just big enough to cover the
    inputed max/min area. (This will not be possible if the max/min area
    has a section outside of the latitudes/longitudes area.)

    assume latitudess & longitudes are sorted,
    long - from low to high (west to east, eg 148 - 151)
    lat - from high to low (north to south, eg -35 - -38)
    """

    # reverse order of lat, so it's in ascending order
    latitudes.reverse()
    largest_lat_index = len(latitudes)-1
    #Cut out a smaller extent.
    if minlat == None:
        lat_min_index = 0
    else:
        lat_min_index = searchsorted(latitudes, minlat)-1
        if lat_min_index <0:
            lat_min_index = 0


    if maxlat == None:
        lat_max_index = largest_lat_index #len(latitudes)
    else:
        lat_max_index = searchsorted(latitudes, maxlat)
        if lat_max_index > largest_lat_index:
            lat_max_index = largest_lat_index

    if minlon == None:
        lon_min_index = 0
    else:
        lon_min_index = searchsorted(longitudes, minlon)-1
        if lon_min_index <0:
            lon_min_index = 0

    if maxlon == None:
        lon_max_index = len(longitudes)
    else:
        lon_max_index = searchsorted(longitudes, maxlon)

    #Take into account that the latitude list was reversed
    latitudes.reverse() # Python passes by reference, need to swap back
    lat_min_index, lat_max_index = largest_lat_index - lat_max_index , \
                                   largest_lat_index - lat_min_index
    lat_max_index = lat_max_index + 1 # taking into account how slicing works
    lon_max_index = lon_max_index + 1 # taking into account how slicing works

    return lat_min_index, lat_max_index, lon_min_index, lon_max_index


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

    if verbose: print 'Reading DEM from %s' %(filename)
    lines = datafile.readlines()
    datafile.close()

    if verbose: print 'Got', len(lines), ' lines'

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
        grid.append(array([float(x) for x in cells]))
    grid = array(grid)

    return {'xllcorner':xllcorner,
            'yllcorner':yllcorner,
            'cellsize':cellsize,
            'NODATA_value':NODATA_value}, grid


# FIXME (Ole): Is this doing anything at all? 
def sww2timeseries(swwfile,
                   gauge_filename,
                   gauge_data_outname,
                   quantity = None,
                   time_min = None,
                   time_max = None,
                   verbose = False):

    """Read SWW file and extract time series for prescribed quantities at
    gauge locations.

    The gauge locations are defined in gauge_filename. This file should be
    in the form: gaugename, easting, northing, and should be stored in a
    .csv or .xya file.

    Time series data at the gauges can be written to file (for example,
    Benfield requested raw data) with the default being no data written.

    The parameter quantity must be the name of an existing quantity or
    an expression involving existing quantities. The default is
    'depth'.

    The user can define a list of quantities. The possibilities are
    the conserved quantitues of the shallow water wave equation and other
    quantities which can be derived from those, i.e.
    ['depth', 'xmomentum', 'ymomentum', 'momentum', 'velocity', 'bearing'].

    Momentum is the absolute momentum, sqrt(xmomentum^2 + ymomentum^2).
    Note, units of momentum are m^2/s and depth is m.

    Velocity is absolute momentum divided by depth. (Confirming correct units:
    vel = abs mom / depth = (m^2/s)/m = m/s.)

    Bearing returns the angle of the velocity vector from North.

    If time_min and time_max is given, output plotted in that time range.
    The default is to plot the entire range of the time evident in sww file.

    The export graphic format in 'png' and will be stored in the same
    directory as the input file.
    """

    if quantity is None: quantity = 'depth'

    # extract gauge locations from gauge file

    # extract all quantities from swwfile (f = file_function)
    if time_min is None: time_min = min(f.get_time())
    if time_max is None: time_max = max(f.get_time())

    # loop through appropriate range of time
    # plot prescribed quantities and export data if requested

    #if gauge_data_outname is None:

    return

#-------------------------------------------------------------
if __name__ == "__main__":
    pass

