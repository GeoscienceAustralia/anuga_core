"""

"""


#FIXME: Ensure that all attributes of a georef are treated everywhere
#and unit test

import types, sys
from Numeric import array, Float, ArrayType, reshape, allclose
from anuga.utilities.numerical_tools import ensure_numeric
from anuga.utilities.anuga_exceptions import ANUGAError, TitleError, ParsingError, \
     ShapeError


DEFAULT_ZONE = -1
TITLE = '#geo reference' + "\n" #this title is referred to in the .xya format

DEFAULT_PROJECTION = 'UTM'
DEFAULT_DATUM = 'wgs84'
DEFAULT_UNITS = 'm'
DEFAULT_FALSE_EASTING = 500000
DEFAULT_FALSE_NORTHING = 10000000 #Default for southern hemisphere 

class Geo_reference:
    """
    """

    def __init__(self,
                 zone = DEFAULT_ZONE,
                 xllcorner = 0.0,
                 yllcorner = 0.0,
                 datum = DEFAULT_DATUM,
                 projection = DEFAULT_PROJECTION,
                 units = DEFAULT_UNITS,
                 false_easting = DEFAULT_FALSE_EASTING,
                 false_northing = DEFAULT_FALSE_NORTHING, 
                 NetCDFObject=None,
                 ASCIIFile=None,
                 read_title=None):
        """
        input:
        NetCDFObject - a handle to the netCDF file to be written to 
        ASCIIFile - a handle to the text file
        read_title - the title of the georeference text, if it was read in.
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
        self.false_easting = false_easting
        self.false_northing = false_northing        
        self.datum = datum
        self.projection = projection
        self.zone = zone        
        self.units = units
        self.xllcorner = xllcorner
        self.yllcorner = yllcorner        
            
        if NetCDFObject is not None:
            self.read_NetCDF(NetCDFObject)
  
        if ASCIIFile is not None:
            self.read_ASCII(ASCIIFile, read_title=read_title)

    def get_xllcorner(self):
        return self.xllcorner
    
    def get_yllcorner(self):
        return self.yllcorner
    
    def get_zone(self):
        return self.zone
        
    def write_NetCDF(self, outfile):
        outfile.xllcorner = self.xllcorner 
        outfile.yllcorner = self.yllcorner
        outfile.zone = self.zone

        outfile.false_easting = self.false_easting
        outfile.false_northing = self.false_northing

        outfile.datum = self.datum        
        outfile.projection = self.projection
        outfile.units = self.units

    def read_NetCDF(self, infile):
        self.xllcorner = infile.xllcorner[0]
        self.yllcorner = infile.yllcorner[0] 
        self.zone = infile.zone[0]

        
        # Fix some assertion failures
        if type(self.zone) == ArrayType and self.zone.shape == ():
            self.zone = self.zone[0]
        if type(self.xllcorner) == ArrayType and self.xllcorner.shape == ():
            self.xllcorner = self.xllcorner[0]
        if type(self.yllcorner) == ArrayType and self.yllcorner.shape == ():
            self.yllcorner = self.yllcorner[0]

        assert (type(self.xllcorner) == types.FloatType or\
                type(self.xllcorner) == types.IntType)
        assert (type(self.yllcorner) == types.FloatType or\
                type(self.yllcorner) == types.IntType)
        assert (type(self.zone) == types.IntType)
        
        try:
            self.false_easting = infile.false_easting[0]
            self.false_northing = infile.false_northing[0]
        
            self.datum = infile.datum        
            self.projection = infile.projection
            self.units = infile.units
        except:
            pass
        if (self.false_easting != DEFAULT_FALSE_EASTING):
            print "WARNING: False easting of %f specified." %self.false_easting
            print "Default false easting is %f." %DEFAULT_FALSE_EASTING
            print "ANUGA does not correct for differences in False Eastings."
            
        if (self.false_northing != DEFAULT_FALSE_NORTHING):
            print "WARNING: False northing of %f specified." \
                  %self.false_northing
            print "Default false northing is %f." %DEFAULT_FALSE_NORTHING
            print "ANUGA does not correct for differences in False Northings."
            
        if (self.datum.upper() != DEFAULT_DATUM.upper()):
            print "WARNING: Datum of %s specified." \
                  %self.datum
            print "Default Datum is %s." %DEFAULT_DATUM
            print "ANUGA does not correct for differences in datums."
            
        if (self.projection.upper() != DEFAULT_PROJECTION.upper()):
            print "WARNING: Projection of %s specified." \
                  %self.projection
            print "Default Projection is %s." %DEFAULT_PROJECTION
            print "ANUGA does not correct for differences in Projection."
            
        if (self.units.upper() != DEFAULT_UNITS.upper()):
            print "WARNING: Units of %s specified." \
                  %self.units
            print "Default units is %s." %DEFAULT_UNITS
            print "ANUGA does not correct for differences in units."
        
        
    def write_ASCII(self, fd):
        fd.write(TITLE)
        fd.write(str(self.zone) + "\n") 
        fd.write(str(self.xllcorner) + "\n") 
        fd.write(str(self.yllcorner) + "\n")

    def read_ASCII(self, fd,read_title=None):
        try:
            if read_title == None:
                read_title = fd.readline() # remove the title line
            if read_title[0:2].upper() != TITLE[0:2].upper():
                msg = 'File error.  Expecting line: %s.  Got this line: %s' \
                      %(TITLE, read_title)
                raise TitleError, msg 
            self.zone = int(fd.readline())
            self.xllcorner = float(fd.readline())
            self.yllcorner = float(fd.readline())
        except SyntaxError:
                msg = 'File error.  Got syntax error while parsing geo reference' 
                raise ParsingError, msg
            
        # Fix some assertion failures
        if(type(self.zone) == ArrayType and self.zone.shape == ()):
            self.zone = self.zone[0]
        if type(self.xllcorner) == ArrayType and self.xllcorner.shape == ():
            self.xllcorner = self.xllcorner[0]
        if type(self.yllcorner) == ArrayType and self.yllcorner.shape == ():
            self.yllcorner = self.yllcorner[0]
    
        assert (type(self.xllcorner) == types.FloatType)
        assert (type(self.yllcorner) == types.FloatType)
        assert (type(self.zone) == types.IntType)
        
        
    def change_points_geo_ref(self, points,
                              points_geo_ref=None):
        """
        Change the geo reference of a list or Numeric array of points to
        be this reference.(The reference used for this object)
        If the points do not have a geo ref, assume 'absolute' values
        """

        is_list = False
        if type(points) == types.ListType:
            is_list = True

        points = ensure_numeric(points, Float)
	
        if len(points.shape) == 1:
            #One point has been passed
            msg = 'Single point must have two elements'
            assert len(points) == 2, msg
            points = reshape(points, (1,2))

        msg = 'Input must be an N x 2 array or list of (x,y) values. '
        msg += 'I got an %d x %d array' %points.shape    
        assert points.shape[1] == 2, msg                

            
        if points_geo_ref is not self:
            #add point geo ref to points
            if not points_geo_ref is None:
                points[:,0] += points_geo_ref.xllcorner 
                points[:,1] += points_geo_ref.yllcorner 
        
            #subtract primary geo ref from points
            points[:,0] -= self.xllcorner 
            points[:,1] -= self.yllcorner

            
        if is_list:
            points = points.tolist()
            
        return points

    
    def is_absolute(self):
        """Return True if xllcorner==yllcorner==0 indicating that points
        in question are absolute.
        """

        return allclose([self.xllcorner, self.yllcorner], 0) 

        
    
    def get_absolute(self, points):
        """
        Given a set of points geo referenced to this instance,
        return the points as absolute values.
        """

        #if self.is_absolute():
        #    return points
        

        is_list = False
        if type(points) == types.ListType:
            is_list = True

        points = ensure_numeric(points, Float)
        if len(points.shape) == 1:
            #One point has been passed
            msg = 'Single point must have two elements'
            if not len(points) == 2:
                raise ShapeError, msg    
                #points = reshape(points, (1,2))


        msg = 'Input must be an N x 2 array or list of (x,y) values. '
        msg += 'I got an %d x %d array' %points.shape    
        if not points.shape[1] == 2:
            raise ShapeError, msg    
            
        
        # Add geo ref to points
        if not self.is_absolute():
            points[:,0] += self.xllcorner 
            points[:,1] += self.yllcorner

        
        if is_list:
            points = points.tolist()
             
        return points


    def reconcile_zones(self, other):
        if other is None:
            other = Geo_reference()
        if self.zone == other.zone or\
               self.zone == DEFAULT_ZONE and other.zone == DEFAULT_ZONE:
            pass
        elif self.zone == DEFAULT_ZONE:
            self.zone = other.zone
        elif other.zone == DEFAULT_ZONE:
            other.zone = self.zone            
        else:    
            msg = 'Both geospatial_data objects must be in the same '+\
                  'ZONE to allow reconciliation. I got zone %d and %d'\
                  %(self.zone, other.zone)
            raise ANUGAError, msg
    
    #def easting_northing2geo_reffed_point(self, x, y):
    #    return [x-self.xllcorner, y - self.xllcorner]

    #def easting_northing2geo_reffed_points(self, x, y):
    #    return [x-self.xllcorner, y - self.xllcorner]

    def get_origin(self):
        return (self.zone, self.xllcorner, self.yllcorner)
    
    def __repr__(self):
        return "(zone=%i easting=%f, northing=%f)" %(self.zone, self.xllcorner, self.yllcorner)
    
    def __cmp__(self,other):
        # FIXME (DSG) add a tolerence
        if other is None:
            return 1
        cmp = 0
        if not (self.xllcorner == self.xllcorner):
            cmp = 1
        if not (self.yllcorner == self.yllcorner):
            cmp = 1
        if not (self.zone == self.zone):
            cmp = 1
        return cmp

def write_NetCDF_georeference(origin, outfile):
    """
    Write georeferrence info to a netcdf file, usually sww.

    The origin can be a georef instance or parrameters for a geo_ref instance

    outfile is the name of the file to be written to.
    """
    geo_ref = ensure_geo_reference(origin)
    geo_ref.write_NetCDF(outfile)
    return geo_ref

def ensure_geo_reference(origin):
    """
    Given a list/tuple of zone, xllcorner and yllcorner of a geo-ref object,
    return a geo ref object.

    If the origin is None, return None, so calling this function doesn't
    effect code logic
    """
    if isinstance(origin, Geo_reference): 
        geo_ref = origin
    elif origin is None:
        geo_ref = None
    else:
        geo_ref = apply(Geo_reference,origin)       
    return geo_ref

    
#-----------------------------------------------------------------------

if __name__ == "__main__":
    pass
