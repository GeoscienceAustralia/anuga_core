import unittest
import copy
import os
import numpy as num


            
from anuga.file_conversion.ll2array import ll2array
from anuga.coordinate_transforms.lat_long_UTM_conversion import UTMtoLL

#Aux for fit_interpolate.fit example
def linear_function(point):
    point = num.array(point)
    return point[:,0]+3*point[:,1]
    #return point[:,1]
    
    
def axes2points(x, y):
    """Generate all combinations of grid point coordinates from x and y axes

    Args:
        * x: x coordinates (array)
        * y: y coordinates (array)

    Returns:
        * P: Nx2 array consisting of coordinates for all
             grid points defined by x and y axes. The x coordinate
             will vary the fastest to match the way 2D numpy
             arrays are laid out by default ('C' order). That way,
             the x and y coordinates will match a corresponding
             2D array A when flattened (A.flat[:] or A.reshape(-1))

    Note:
        Example

        x = [1, 2, 3]
        y = [10, 20]

        P = [[1, 10],
             [2, 10],
             [3, 10],
             [1, 20],
             [2, 20],
             [3, 20]]
    """
    import numpy 
    
    # Reverse y coordinates to have them start at bottom of array
    y = numpy.flipud(y)

    # Repeat x coordinates for each y (fastest varying)
    X = numpy.kron(numpy.ones(len(y)), x)

    # Repeat y coordinates for each x (slowest varying)
    Y = numpy.kron(y, numpy.ones(len(x)))

    # Check
    N = len(X)
    assert len(Y) == N

    # Create Nx2 array of x and y coordinates
    X = numpy.reshape(X, (N, 1))
    Y = numpy.reshape(Y, (N, 1))
    P = numpy.concatenate((X, Y), axis=1)

    # Return
    return P

    
class Test_ll2array(unittest.TestCase):
 
     def test_ll2array(self):

        """ Format of ll file 
        ncols          11
        nrows          12
        latcorner      45
        longcorner     90
        lat spacing    0.1
	long spacing   0.1
        NODATA_value  -9999
        """
        
        x0 = 240000.0
        y0 = 7620000.0
	zone = 56
        
        ncols = 11  # Nx
        nrows = 12  # Ny
        xllcorner = x0
        yllcorner = y0
        latcorner, longcorner = UTMtoLL(y0, x0, zone)
	cellsize = 6000.0
        NODATA_value =  -9999

	# create lat long grid file from utm grids
        x = num.linspace(xllcorner, xllcorner+(ncols-1)*cellsize, ncols)
        y = num.linspace(yllcorner, yllcorner+(nrows-1)*cellsize, nrows)
	points = axes2points(x, y)
	yll = num.array([list(UTMtoLL(yi, x[0], zone)) for yi in y])
	lato = yll[:,0]
	xll = num.array([list(UTMtoLL(y[0], xi, zone)) for xi in x])
	longo = xll[:,1]
	latinc = (max(lato)-min(lato))/(nrows-1)
	longinc = (max(longo)-min(longo))/(ncols-1)
	
#	latinc = 0.0541605208492
#	longinc = 0.0578860943416

        #Create .ll file
        txt_file = 'test_ll.ll'
        datafile = open(txt_file,"w")
        datafile.write('ncols '+str(ncols)+"\n")
        datafile.write('nrows '+str(nrows)+"\n")
        datafile.write('latcorner '+str(latcorner)+"\n")
        datafile.write('longcorner '+str(longcorner)+"\n")
        datafile.write('lat spacing '+str(latinc)+"\n")
        datafile.write('long spacing '+str(longinc)+"\n")
        datafile.write('NODATA_value '+str(NODATA_value)+"\n")
        
 	datavalues = linear_function(points)
        datavalues = datavalues.reshape(nrows,ncols)
#        print datavalues, 'datavalues'
        #print datavalues.shape
        for row in datavalues:
            #print row
            datafile.write(" ".join(str(elem) for elem in row) + "\n")         
        datafile.close()       
        
        lats, longs, Z = ll2array(txt_file)

	lats_ex = num.linspace(latcorner, latcorner+(nrows-1)*latinc, nrows)
	longs_ex = num.linspace(longcorner, longcorner+(ncols-1)*longinc, ncols)

        Z_ex = [[    23100000., 23118000., 23136000., 23154000., 23172000., 23190000.,
                       23208000., 23226000., 23244000., 23262000., 23280000., 23298000.],
                     [ 23106000., 23124000., 23142000., 23160000., 23178000., 23196000.,
                       23214000., 23232000., 23250000., 23268000., 23286000., 23304000.],
                     [ 23112000., 23130000., 23148000., 23166000., 23184000., 23202000.,
                       23220000., 23238000., 23256000., 23274000., 23292000., 23310000.],
                     [ 23118000., 23136000., 23154000., 23172000., 23190000., 23208000.,
                       23226000., 23244000., 23262000., 23280000., 23298000., 23316000.],
                     [ 23124000., 23142000., 23160000., 23178000., 23196000., 23214000.,
                       23232000., 23250000., 23268000., 23286000., 23304000., 23322000.],
                     [ 23130000., 23148000., 23166000., 23184000., 23202000., 23220000.,
                       23238000., 23256000., 23274000., 23292000., 23310000., 23328000.],
                     [ 23136000., 23154000., 23172000., 23190000., 23208000., 23226000.,
                       23244000., 23262000., 23280000., 23298000., 23316000., 23334000.],
                     [ 23142000., 23160000., 23178000., 23196000., 23214000., 23232000.,
                       23250000., 23268000., 23286000., 23304000., 23322000., 23340000.],
                     [ 23148000., 23166000., 23184000., 23202000., 23220000., 23238000.,
                       23256000., 23274000., 23292000., 23310000., 23328000., 23346000.],
                     [ 23154000., 23172000., 23190000., 23208000., 23226000., 23244000.,
                       23262000., 23280000., 23298000., 23316000., 23334000., 23352000.],
                     [ 23160000., 23178000., 23196000., 23214000., 23232000., 23250000.,
                       23268000., 23286000., 23304000., 23322000., 23340000., 23358000.]]
        


        #print quantity.vertex_values
        
        assert num.allclose(Z, Z_ex)
        assert num.allclose(lats, lats_ex)
        assert num.allclose(longs,longs_ex)
        
        os.remove(txt_file)
      
 

################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_ll2array, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
        
