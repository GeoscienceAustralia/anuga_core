"""

"""


#FIXME: Ensure that all attributes of a georef are treated everywhere
#and unit test

from builtins import str
from builtins import object
import sys
import copy

from anuga.utilities.numerical_tools import ensure_numeric
from anuga.anuga_exceptions import ANUGAError, TitleError, \
                                             ParsingError, ShapeError
from anuga.config import netcdf_float, netcdf_int, netcdf_float32
import anuga.utilities.log as log

import numpy as num


DEFAULT_ZONE = -1  # This signifies that simulation isn't located within a UTM framework.
                   # This is the case for hypothetical simulations or something relative
                   # to an arbitrary origin (e.g. the corner of a wavetank).

DEFAULT_PROJECTION = 'UTM'
DEFAULT_DATUM = 'wgs84'
DEFAULT_UNITS = 'm'
DEFAULT_FALSE_EASTING = 500000
DEFAULT_FALSE_NORTHING = 10000000    
DEFAULT_HEMISPHERE = 'undefined'

TITLE = '#geo reference' + "\n" # this title is referred to in the test format

class Geo_reference(object):
    """
    Attributes of the Geo_reference class:
        .zone           The UTM zone (default is -1)
        .hemisphere     southern, northern or undefined hemisphere
        .false_easting  ??
        .false_northing ??
        .datum          The Datum used (default is wgs84)
        .projection     The projection used (default is 'UTM')
        .units          The units of measure used (default metres)
        .xllcorner      The X coord of origin (default is 0.0 wrt UTM grid)
        .yllcorner      The y coord of origin (default is 0.0 wrt UTM grid)
        .is_absolute    ??
    

    """

    def __init__(self,
                 zone=None,
                 xllcorner=0.0,
                 yllcorner=0.0,
                 datum=DEFAULT_DATUM,
                 projection=DEFAULT_PROJECTION,
                 units=DEFAULT_UNITS,
                 false_easting=DEFAULT_FALSE_EASTING,
                 false_northing=DEFAULT_FALSE_NORTHING,
                 hemisphere=DEFAULT_HEMISPHERE,
                 NetCDFObject=None,
                 ASCIIFile=None,
                 read_title=None):
        """
        zone            the UTM zone.
        hemisphere      southern or northern hemisphere
        xllcorner       X coord of origin of georef.
        yllcorner       Y coord of origin of georef.
        datum           ??
        projection      the projection used (default UTM).
        units           units used in measuring distance (default m).
        false_easting   ??
        false_northing  ??
        NetCDFObject    NetCDF file *handle* to write to.
        ASCIIFile       ASCII text file *handle* to write to.
        read_title      title of the georeference text.

        If the function that calls this has already read the title line,
        it can't unread it, so this info has to be passed.
        If you know of a way to unread this info, then tell us.

        Note, the text file only saves a sub set of the info the
        points file does.  Currently the info not written in text
        must be the default info, since ANUGA assumes it isn't
        changing.
        """

        if zone is None:
            zone = DEFAULT_ZONE

        self.set_zone(zone)
        self.set_hemisphere(hemisphere)
        self.false_easting = int(false_easting)
        self.false_northing = int(false_northing)
        self.datum = datum
        self.projection = projection
    
        self.units = units
        self.xllcorner = float(xllcorner)
        self.yllcorner = float(yllcorner)
            
        if NetCDFObject is not None:
            self.read_NetCDF(NetCDFObject)

        if ASCIIFile is not None:
            self.read_ASCII(ASCIIFile, read_title=read_title)
            
        # Set flag for absolute points (used by get_absolute)
        # FIXME (Ole): It would be more robust to always use get_absolute()
        self.absolute = num.allclose([self.xllcorner, self.yllcorner], 0)
        
    def __eq__(self, other):

        # FIXME (Ole): Can this be automatically done for all attributes?
        # Anyway, it is arranged like this so one can step through and find out
        # why two objects might not be equal
        
        equal = True
        if self.false_easting != other.false_easting: equal = False
        if self.false_northing != other.false_northing: equal = False
        if self.datum != other.datum: equal = False
        if self.projection != other.projection: equal = False
        if self.zone != other.zone: equal = False
        if self.hemisphere != other.hemisphere: equal = False
        if self.units != other.units: equal = False
        if self.xllcorner != other.xllcorner: equal = False
        if self.yllcorner != other.yllcorner: equal = False
        if self.absolute != other.absolute: equal = False

        return(equal)

    def get_xllcorner(self):
        """Get the X coordinate of the origin of this georef."""
        return self.xllcorner

    def get_yllcorner(self):
        """Get the Y coordinate of the origin of this georef."""

        return self.yllcorner

    def set_zone(self, zone):
        """set zone as an integer in [1,60] or -1."""

        zone = int(zone)

        assert (zone == -1 or (zone >= 1 and zone <= 60)), f'zone {zone} not valid.'

        self.zone = zone

    def get_zone(self):
        """Get the zone of this georef."""

        return self.zone

    def get_hemisphere(self):
        """Check if this georef has a defined hemisphere."""

        return self.hemisphere

    def set_hemisphere(self, hemisphere):

        msg = f"'{hemisphere}' not corresponding to allowed hemisphere values 'southern', 'northern' or 'undefined'" 
        assert hemisphere in ['southern', 'northern', 'undefined'], msg

        self.hemisphere=str(hemisphere)

    def write_NetCDF(self, outfile):
        """Write georef attributes to an open NetCDF file.

        outfile  handle to open NetCDF file
        """

        outfile.xllcorner = self.xllcorner
        outfile.yllcorner = self.yllcorner
        outfile.zone = self.zone
        outfile.hemisphere = self.hemisphere

        outfile.false_easting = self.false_easting
        outfile.false_northing = self.false_northing

        outfile.datum = self.datum
        outfile.projection = self.projection
        outfile.units = self.units

    def read_NetCDF(self, infile):
        """Set georef attributes from open NetCDF file.

        infile Handle to open NetCDF file
        """

        self.xllcorner = float(infile.xllcorner)
        self.yllcorner = float(infile.yllcorner)
        self.zone = int(infile.zone)
        try:
            self.hemisphere = str(infile.hemisphere)
        except:
            self.hemisphere = DEFAULT_HEMISPHERE

        self.false_easting = int(infile.false_easting)
        self.false_northing = int(infile.false_northing)

        self.datum = infile.datum
        self.projection = infile.projection
        self.units = infile.units

        # Set flag for absolute points (used by get_absolute)
        # FIXME (Ole): It would be more robust to always use get_absolute()        
        self.absolute = num.allclose([self.xllcorner, self.yllcorner], 0)
        
        if self.false_easting != DEFAULT_FALSE_EASTING:
            log.critical("WARNING: False easting of %f specified."
                         % self.false_easting)
            log.critical("Default false easting is %f." % DEFAULT_FALSE_EASTING)
            log.critical("ANUGA does not correct for differences in "
                         "False Eastings.")

        if self.false_northing != DEFAULT_FALSE_NORTHING:
            log.critical("WARNING: False northing of %f specified."
                         % self.false_northing)
            log.critical("Default false northing is %f."
                         % DEFAULT_FALSE_NORTHING)
            log.critical("ANUGA does not correct for differences in "
                         "False Northings.")

        if self.datum.upper() != DEFAULT_DATUM.upper():
            log.critical("WARNING: Datum of %s specified." % self.datum)
            log.critical("Default Datum is %s." % DEFAULT_DATUM)
            log.critical("ANUGA does not correct for differences in datums.")

        if self.projection.upper() != DEFAULT_PROJECTION.upper():
            log.critical("WARNING: Projection of %s specified."
                         % self.projection)
            log.critical("Default Projection is %s." % DEFAULT_PROJECTION)
            log.critical("ANUGA does not correct for differences in "
                         "Projection.")

        if self.units.upper() != DEFAULT_UNITS.upper():
            log.critical("WARNING: Units of %s specified." % self.units)
            log.critical("Default units is %s." % DEFAULT_UNITS)
            log.critical("ANUGA does not correct for differences in units.")



################################################################################
# ASCII files with geo-refs are currently not used
################################################################################

    def write_ASCII(self, fd):
        """Write georef attriutes to an open text file.

        fd  handle to open text file
        """

        fd.write(TITLE)
        fd.write(str(self.zone) + "\n")
        fd.write(str(self.xllcorner) + "\n")
        fd.write(str(self.yllcorner) + "\n")

    def read_ASCII(self, fd, read_title=None):
        """Set georef attribtes from open text file.

        fd  handle to open text file
        """

        try:
            if read_title is None:
                read_title = fd.readline()     # remove the title line
            if read_title[0:2].upper() != TITLE[0:2].upper():
                msg = ('File error.  Expecting line: %s.  Got this line: %s'
                       % (TITLE, read_title))
                raise TitleError(msg)
            self.zone = int(fd.readline())
            self.xllcorner = float(fd.readline())
            self.yllcorner = float(fd.readline())
        except SyntaxError:
            msg = 'File error.  Got syntax error while parsing geo reference'
            raise ParsingError(msg)

        # Fix some assertion failures
        if isinstance(self.zone, num.ndarray) and self.zone.shape == ():
            self.zone = self.zone[0]
        if (isinstance(self.xllcorner, num.ndarray) and
                self.xllcorner.shape == ()):
            self.xllcorner = self.xllcorner[0]
        if (isinstance(self.yllcorner, num.ndarray) and
                self.yllcorner.shape == ()):
            self.yllcorner = self.yllcorner[0]

        assert isinstance(self.xllcorner, float)
        assert isinstance(self.yllcorner, float)
        assert isinstance(self.zone, int)

################################################################################

    def change_points_geo_ref(self, points, points_geo_ref=None):
        """Change points to be absolute wrt new georef 'points_geo_ref'.

        points          the points to change
        points_geo_ref  the new georef to make points absolute wrt

        Returns the changed points data.
        If the points do not have a georef, assume 'absolute' values.
        """

        import copy
       
        # remember if we got a list
        is_list = isinstance(points, list)

        points = ensure_numeric(points, float)

        # sanity checks	
        if len(points.shape) == 1:
            #One point has been passed
            msg = 'Single point must have two elements'
            assert len(points) == 2, msg
            points = num.reshape(points, (1,2))

        msg = 'Points array must be two dimensional.\n'
        msg += 'I got %d dimensions' %len(points.shape)
        assert len(points.shape) == 2, msg

        msg = 'Input must be an N x 2 array or list of (x,y) values. '
        msg += 'I got an %d x %d array' %points.shape    
        assert points.shape[1] == 2, msg                

        # FIXME (Ole): Could also check if zone, xllcorner, yllcorner 
        # are identical in the two geo refs.    
        if points_geo_ref is not self:
            # If georeferences are different
            points = copy.copy(points) # Don't destroy input                    
            if not points_geo_ref is None:
                # Convert points to absolute coordinates
                points[:,0] += points_geo_ref.xllcorner 
                points[:,1] += points_geo_ref.yllcorner 
        
            # Make points relative to primary geo reference
            points[:,0] -= self.xllcorner 
            points[:,1] -= self.yllcorner

        if is_list:
            points = points.tolist()
            
        return points

    def is_absolute(self):
        """Test if points in georef are absolute.

        Return True if xllcorner==yllcorner==0 indicating that points in
        question are absolute.
        """
        
        # FIXME(Ole): It is unfortunate that decision about whether points
        # are absolute or not lies with the georeference object. Ross pointed this out.
        # Moreover, this little function is responsible for a large fraction of the time
        # using in data fitting (something in like 40 - 50%.
        # This was due to the repeated calls to allclose.
        # With the flag method fitting is much faster (18 Mar 2009).

        # FIXME(Ole): HACK to be able to reuse data already cached (18 Mar 2009). 
        # Remove at some point
        if not hasattr(self, 'absolute'):
            self.absolute = num.allclose([self.xllcorner, self.yllcorner], 0)
            
        # Return absolute flag    
        return self.absolute

    def get_absolute(self, points):
        """Given a set of points geo referenced to this instance, return the
        points as absolute values.
        """

        # remember if we got a list
        is_list = isinstance(points, list)

        points = ensure_numeric(points, float)
        if len(points.shape) == 1:
            # One point has been passed
            msg = 'Single point must have two elements'
            if not len(points) == 2:
                raise ShapeError(msg)    


        msg = 'Input must be an N x 2 array or list of (x,y) values. '
        msg += 'I got an %d x %d array' %points.shape    
        if not points.shape[1] == 2:
            raise ShapeError(msg)    
            
        
        # Add geo ref to points
        if not self.is_absolute():
            points = copy.copy(points) # Don't destroy input                    
            points[:,0] += self.xllcorner 
            points[:,1] += self.yllcorner

        
        if is_list:
            points = points.tolist()
             
        return points

    def get_relative(self, points):
        """Convert points to relative measurement.

        points Points to convert to relative measurements

        Returns a set of points relative to the geo_reference instance.

        This is the inverse of get_absolute().
        """

        # remember if we got a list
        is_list = isinstance(points, list)

        points = ensure_numeric(points, float)
        if len(points.shape) == 1:
            #One point has been passed
            msg = 'Single point must have two elements'
            if not len(points) == 2:
                raise ShapeError(msg)    

        if not points.shape[1] == 2:
            msg = ('Input must be an N x 2 array or list of (x,y) values. '
                   'I got an %d x %d array' % points.shape)
            raise ShapeError(msg)    

        # Subtract geo ref from points
        if not self.is_absolute():
            points = copy.copy(points) # Don't destroy input                            
            points[:,0] -= self.xllcorner 
            points[:,1] -= self.yllcorner

        if is_list:
            points = points.tolist()
             
        return points

    def reconcile_zones(self, other):

        if other is None:
            # FIXME(Ole): Why would we do this?
            other = Geo_reference()
            #raise Exception('Expected georeference object, got None')
        
        if (self.zone == other.zone):
            pass        
        elif self.zone == DEFAULT_ZONE:
            self.zone = other.zone
        elif other.zone == DEFAULT_ZONE:
            other.zone = self.zone
        else:
            msg = ('Geospatial data must be in the same '
                   'ZONE to allow reconciliation. I got zone %d and %d'
                   % (self.zone, other.zone))
            raise ANUGAError(msg)

        # Should also reconcile hemisphere
        if (self.hemisphere == other.hemisphere):
            pass        
        elif self.hemisphere == DEFAULT_HEMISPHERE:
            self.hemisphere = other.hemisphere
        elif other.hemisphere == DEFAULT_HEMISPHERE:
            other.hemisphere = self.hemisphere
        else:
            msg = ('Geospatial data must be in the same '
                   'HEMISPHERE to allow reconciliation. I got hemisphere %d and %d'
                   % (self.hemisphere, other.hemisphere))
            raise ANUGAError(msg)        

    # FIXME (Ole): Do we need this back?    
    #def easting_northing2geo_reffed_point(self, x, y):
    #    return [x-self.xllcorner, y - self.xllcorner]

    #def easting_northing2geo_reffed_points(self, x, y):
    #    return [x-self.xllcorner, y - self.xllcorner]

    def get_origin(self):
        """Get origin of this geo_reference."""
      
        return (self.zone, self.xllcorner, self.yllcorner)

    def __repr__(self):
        return ('(zone=%i, easting=%f, northing=%f, hemisphere=%s)'
                % (self.zone, self.xllcorner, self.yllcorner, self.hemisphere))

    #def __cmp__(self, other):
    #    """Compare two geo_reference instances.#
    #
    #    self   this geo_reference instance
    #    other  another geo_reference instance to compare against#
    #
    #    Returns 0 if instances have the same attributes, else returns 1. 
    #
    #    Note: attributes are: zone, xllcorner, yllcorner.
    #    """

    #    # FIXME (DSG) add a tolerence
    #   if other is None:
    #        return 1
    #    cmp = 0
    #    if not (self.xllcorner == self.xllcorner):
    #        cmp = 1
    #    if not (self.yllcorner == self.yllcorner):
    #        cmp = 1
    #    if not (self.zone == self.zone):
    #        cmp = 1
    #    return cmp


def write_NetCDF_georeference(georef, outfile):
    """Write georeference info to a NetCDF file, usually a SWW file.

    georef   a georef instance or parameters to create a georef instance
    outfile  path to file to write

    Returns the normalised georef.
    """

    geo_ref = ensure_geo_reference(georef)
    geo_ref.write_NetCDF(outfile)
    return geo_ref


def ensure_geo_reference(origin):
    """Create a georef object from a tuple of attributes.

    origin  a georef instance or (zone, xllcorner, yllcorner)

    If origin is None, return None, so calling this function doesn't
    effect code logic.
    """

    if isinstance(origin, Geo_reference):
        geo_ref = origin
    elif origin is None:
        geo_ref = None
    else:
        if len(origin) == 1:
            geo_ref = Geo_reference(zone = origin)
        elif len(origin) == 2:
            geo_ref = Geo_reference(zone = -1, xllcorner=origin[0], yllcorner=origin[1])
        elif len(origin) == 3:
            geo_ref = Geo_reference(zone = origin[0], xllcorner=origin[1], yllcorner=origin[2])
        else:
            raise Exception(f'Invalid input {origin}, expected (zone), (xllcorner, yllcorner), or (zone, xllcorner, yllcorner).')


    return geo_ref


#-----------------------------------------------------------------------

if __name__ == "__main__":
    pass
