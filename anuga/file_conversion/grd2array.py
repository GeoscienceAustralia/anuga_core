""" Load a DEM file, decimate it, and resave it.
"""


import numpy as num
import os

from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_float

def grd2array(filename, verbose=False):
        """Read Digital Elevation model from the following ASCII format (.asc or .grd)
    
        Example:
        ncols         3121
        nrows         1800
        xllcorner     722000
        yllcorner     5893000
        cellsize      25
        NODATA_value  -9999
        138.3698 137.4194 136.5062 135.5558 ..........
    
    
        An accompanying file with same basename but extension .prj must exist
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



        msg = 'Filename must be a text string'
        assert isinstance(filename, str), msg
        
        msg = 'Extension should be .grd or asc'
        assert os.path.splitext(filename)[1] in ['.grd', '.asc'], msg
     

            
    
        root = filename[:-4]
    
    
        #Read DEM data
        datafile = open(filename)
    
        if verbose: log.critical('Reading data from %s' % (filename))
    
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
            raise Exception(msg)
    
        yref = lines[3].split()
        if yref[0].strip() == 'yllcorner':
            yllcorner = float(yref[1].strip()) # + 0.5*cellsize # Correct offset
        elif yref[0].strip() == 'yllcenter':
            yllcorner = float(yref[1].strip())
        else:
            msg = 'Unknown keyword: %s' % yref[0].strip()
            raise Exception(msg)
    
        NODATA_value = int(float(lines[5].split()[1].strip()))
    
        assert len(lines) == nrows + 6
    
  
        #Store data
        import numpy
    
        datafile = open(filename)
        Z = numpy.loadtxt(datafile, skiprows=6)
        datafile.close()
        
        #print Z.shape
        #print Z
        
        # For raster data we need to a flip and transpose
        Z = numpy.flipud(Z)

        # Transpose z to have y coordinates along the first axis and x coordinates
        # along the second axis
        Z = Z.transpose()
    
        x = num.linspace(xllcorner, xllcorner+cellsize*(ncols-1), ncols)
        y = num.linspace(yllcorner, yllcorner+cellsize*(nrows-1), nrows)
        
        # get rid of NODATA
        import numpy
        Z = numpy.where(Z == NODATA_value , numpy.nan, Z)
        
        return x,y, Z
        

