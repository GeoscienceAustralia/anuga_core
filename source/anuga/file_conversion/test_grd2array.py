import unittest
import copy
import os
import numpy as num


            
from grd2array import grd2array
                

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

    
class Test_grd2array(unittest.TestCase):
 
 
     def test_grd2array(self):



        """ Format of asc file 
        ncols         11
        nrows         12
        xllcorner     240000
        yllcorner     7620000
        cellsize      6000
        NODATA_value  -9999
        """
        
        x0 = 0.0
        y0 = 0.0
        
        ncols = 11  # Nx
        nrows = 12  # Ny
        xllcorner = x0
        yllcorner = y0
        cellsize  = 1.0
        NODATA_value =  -9999


        
        #Create .asc file
        #txt_file = tempfile.mktemp(".asc")from anuga.config import netcdf_float
        txt_file = 'test_asc.asc'
        datafile = open(txt_file,"w")
        datafile.write('ncols '+str(ncols)+"\n")
        datafile.write('nrows '+str(nrows)+"\n")
        datafile.write('xllcorner '+str(xllcorner)+"\n")
        datafile.write('yllcorner '+str(yllcorner)+"\n")
        datafile.write('cellsize '+str(cellsize)+"\n")
        datafile.write('NODATA_value '+str(NODATA_value)+"\n")
        
        x_ex = num.linspace(xllcorner, xllcorner+(ncols-1)*cellsize, ncols)
        y_ex = num.linspace(yllcorner, yllcorner+(nrows-1)*cellsize, nrows)
        points = axes2points(x_ex, y_ex)
        
        #print points
        #print x.shape, x
        #print y.shape, y
        
        datavalues = linear_function(points)
        #print datavalues 
        
        datavalues = datavalues.reshape(nrows,ncols)

        #print datavalues
        #print datavalues.shape
        for row in datavalues:
            #print row
            datafile.write(" ".join(str(elem) for elem in row) + "\n")         
        datafile.close()

        #print quantity.vertex_values
        #print quantity.centroid_values 

        x,y,Z = grd2array(txt_file)

        #print x
        #print y
        #print Z

         
        
        answer = [[  0.,  3.,  6.,  9., 12., 15., 18., 21., 24., 27., 30., 33.],
                 [  1.,  4.,  7., 10., 13., 16., 19., 22., 25., 28., 31., 34.],
                 [  2.,  5.,  8., 11., 14., 17., 20., 23., 26., 29., 32., 35.],
                 [  3.,  6.,  9., 12., 15., 18., 21., 24., 27., 30., 33., 36.],
                 [  4.,  7., 10., 13., 16., 19., 22., 25., 28., 31., 34., 37.],
                 [  5.,  8., 11., 14., 17., 20., 23., 26., 29., 32., 35., 38.],
                 [  6.,  9., 12., 15., 18., 21., 24., 27., 30., 33., 36., 39.],
                 [  7., 10., 13., 16., 19., 22., 25., 28., 31., 34., 37., 40.],
                 [  8., 11., 14., 17., 20., 23., 26., 29., 32., 35., 38., 41.],
                 [  9., 12., 15., 18., 21., 24., 27., 30., 33., 36., 39., 42.],
                 [ 10., 13., 16., 19., 22., 25., 28., 31., 34., 37., 40., 43.]]

        #print quantity.vertex_values
        
        assert num.allclose(Z, answer)
        assert num.allclose(x,x_ex)
        assert num.allclose(y,y_ex)
        
        
 

#################################################################################

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_grd2array, 'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
        
