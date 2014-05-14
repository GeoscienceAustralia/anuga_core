#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$17/04/2012 11:32:04 AM$"


# external modules
import numpy as np
import anuga.utilities.log as log

import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


if __name__ == "__main__":
    """
    Example:
    ncols         3121
    nrows         1800
    xllcorner     722000
    yllcorner     5893000
    cellsize      25
    NODATA_value  -9999
    138.3698 137.4194 136.5062 135.5558 ..........
    """


    name_in = "V.asc"
    verbose = True


    #Read DEM data
    datafile = open(name_in)

    if verbose: log.critical('Reading DEM from %s' % (name_in))

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

    NODATA_value = float(lines[5].split()[1].strip())

    assert len(lines) == nrows + 6


    print 'rows ', nrows
    print 'cols ', ncols
    print 'cell size ',cellsize

    print 'xllcorner ',xllcorner
    print 'yllcorner ',yllcorner

    Z = np.zeros((nrows,ncols),dtype=np.float)

    n = len(lines[6:])
    for i, line in enumerate(lines[6:]):
        fields = line.split()
        if verbose and i % ((n+10)/10) == 0:
            log.critical('Processing row %d of %d' % (i, nrows))

        if len(fields) != ncols:
            msg = 'Wrong number of columns in file "%s" line %d\n' % (name_in, i)
            msg += 'I got %d elements, but there should have been %d\n' % (len(fields), ncols)
            raise Exception, msg

        Z[i, :] = np.array([float(x) for x in fields])


    print Z.shape


    
    ZZ =  np.where(Z > -9999.0 , Z, np.nan)
    

    z_min = int(np.nanmin(ZZ))
    z_max = int(np.nanmax(ZZ))+1

    print z_min,z_max

    xlen = ncols*cellsize
    ylen = nrows*cellsize
    x = np.arange(0.0, xlen, cellsize)
    y = np.arange(ylen, 0.0, -cellsize)
    X, Y = np.meshgrid(x, y)

    print X.shape
    print Y.shape

    plt.figure()


    levels = z_min + np.arange(11,dtype=np.float)/10.0*(z_max-z_min)


    print levels
    
    #CS = plt.contour(X, Y, Z)

    CS = plt.contour(X, Y, Z, levels)
                        #colors = ('r', 'g', 'b'),
                        #extend='both')

    #CS.cmap.set_under('yellow')
    #CS.cmap.set_over('yellow')

    cbar = plt.colorbar(CS)
    cbar.ax.set_ylabel('velocity')
    # Add the contour line levels to the colorbar
    cbar.add_lines(CS)

    plt.clabel(CS, inline=5, fontsize=10)
    plt.title('Simplest default with labels')





    plt.show()
    

    


    





