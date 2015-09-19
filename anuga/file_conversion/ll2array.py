""" Load a lat long grid file, decimate it, and resave it.
"""

import numpy as num
import os

from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_float

def ll2array(filename, verbose=False):
        """Read latitude and longitude grid file from the following latlong format (.ll)
    
           ncols          11
	   nrows          12
	   latcorner      45
	   longcorner     90
	   lat increment  0.1
	   long increment 0.1
	   NODATA_value  -9999
	   """


        msg = 'Filename must be a text string'
        assert isinstance(filename, basestring), msg
        
        msg = 'Extension should be .ll or?'
        assert os.path.splitext(filename)[1] in ['.ll'], msg
     
        root = filename[:-4]
    
        #Read ll data
        datafile = open(filename)
    
        if verbose: log.critical('Reading data from %s' % (filename))
    
        lines = datafile.readlines()
        datafile.close()
    
        if verbose: log.critical('Got %d lines' % len(lines))
    

        ncols = int(lines[0].split()[1].strip())
        nrows = int(lines[1].split()[1].strip())
	latcorner = float(lines[2].split()[1].strip())
	longcorner = float(lines[3].split()[1].strip())

        latinc = float(lines[4].split()[2].strip())
        longinc = float(lines[5].split()[2].strip())
        
        # Checks suggested by Joaquim Luis
        # Our internal representation of xllcorner
        # and yllcorner is non-standard.
        latref = lines[2].split()
        if latref[0].strip() == 'latcorner':
            latcorner = float(latref[1].strip())
        else:
            msg = 'Unknown keyword: %s' % latref[0].strip()
            raise Exception, msg
    
        longref = lines[3].split()
        if longref[0].strip() == 'longcorner':
            longcorner = float(longref[1].strip())
        else:
            msg = 'Unknown keyword: %s' % longref[0].strip()
            raise Exception, msg
    
        NODATA_value = int(float(lines[6].split()[1].strip()))
    
        assert len(lines) == nrows + 7
    
  
        #Store data
        import numpy
    
        datafile = open(filename)
        Z = numpy.loadtxt(datafile, skiprows=7)
        datafile.close()
        
        #print Z.shape
        
        # For raster data we need to a flip and transpose
        Z = numpy.flipud(Z)

        # Transpose z to have y coordinates along the first axis and x coordinates
        # along the second axis
        Z = Z.transpose()
    
        lats = numpy.linspace(latcorner, latcorner+(nrows-1)*latinc, nrows)
	longs = numpy.linspace(longcorner, longcorner+(ncols-1)*longinc, ncols)
		
        # get rid of NODATA
        import numpy
        Z = numpy.where(Z == NODATA_value , numpy.nan, Z)
        
        return lats, longs, Z
        

