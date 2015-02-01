#!/usr/bin/env python
"""
    Generic boundary conditions for a domain.
    
    A boundary represents the edge of the model, where inflow, outflow, and
    reflection can take place.
    
    The boundaries in this model can be applied universally across all
    domain models, without being tied to a particular implementation.
"""

import unittest
from math import sqrt, pi

from anuga.abstract_2d_finite_volumes.generic_boundary_conditions import *
from anuga.abstract_2d_finite_volumes.generic_domain import Generic_Domain
from anuga.abstract_2d_finite_volumes.quantity import Quantity
from anuga.config import epsilon

import numpy as num


class Test_Generic_Boundary_Conditions(unittest.TestCase):
    def setUp(self):
        pass
        #print "  Setting up"

    def tearDown(self):
        pass
        #print "  Tearing down"


    def test_generic(self):
        b = Boundary()

        try:
            b.evaluate()
        except:
            pass
        else:
            raise Exception('Should have raised exception')


    def test_dirichlet_empty(self):

        try:
            Bd = Dirichlet_boundary()
        except:
            pass
        else:
            raise Exception('Should have raised exception')

    def test_dirichlet(self):
        x = [3.14,0,0.1]
        Bd = Dirichlet_boundary(x)

        q = Bd.evaluate()
        assert num.allclose(q, x)


    def test_time(self):

        

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Generic_Domain(points, elements)
        domain.check_integrity()

        domain.conserved_quantities = ['stage', 'ymomentum']
        domain.evolved_quantities = ['stage', 'ymomentum']        
        domain.quantities['stage'] =\
                                   Quantity(domain, [[1,2,3], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])

        domain.quantities['ymomentum'] =\
                                   Quantity(domain, [[2,3,4], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])


        domain.check_integrity()

        #Test time bdry, you need to provide a domain and function
        try:
            T = Time_boundary(domain)
        except:
            pass
        else:
            raise Exception('Should have raised exception')

        #Test time bdry, you need to provide a function
        try:
            T = Time_boundary()
        except:
            pass
        else:
            raise Exception('Should have raised exception')


        def function(t):
            return [1.0, 0.0]
        
        T = Time_boundary(domain, function)

        from anuga.config import default_boundary_tag
        domain.set_boundary( {default_boundary_tag: T} )


        #FIXME: should not necessarily be true always.
        #E.g. with None as a boundary object.
        assert len(domain.boundary) == len(domain.boundary_objects)

        q = T.evaluate(0, 2)  #Vol=0, edge=2

        assert num.allclose(q, [1.0, 0.0])


    def test_time_space_boundary(self):


        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Generic_Domain(points, elements)
        domain.check_integrity()

        domain.conserved_quantities = ['stage', 'ymomentum']
        domain.evolved_quantities = ['stage', 'ymomentum']        
        domain.quantities['stage'] =\
                                   Quantity(domain, [[1,2,3], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])

        domain.quantities['ymomentum'] =\
                                   Quantity(domain, [[2,3,4], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])


        domain.check_integrity()

        #Test time space bdry, you need to provide a domain and function
        try:
            T = Time_space_boundary(domain)
        except:
            pass
        else:
            raise Exception('Should have raised exception')

        #Test time bdry, you need to provide a function
        try:
            T = Time_space_boundary()
        except:
            pass
        else:
            raise Exception('Should have raised exception')


        def function(t,x,y):
            return [x,y]
        
        T = Time_space_boundary(domain, function)

        from anuga.config import default_boundary_tag
        domain.set_boundary( {default_boundary_tag: T} )


        #FIXME: should not necessarily be true always.
        #E.g. with None as a boundary object.
        assert len(domain.boundary) == len(domain.boundary_objects)

        q = T.evaluate(0, 2)  #Vol=0, edge=2
        assert num.allclose(q, domain.get_edge_midpoint_coordinate(0,2))


        q = T.evaluate(1, 1)  #Vol=1, edge=1
        assert num.allclose(q, domain.get_edge_midpoint_coordinate(1,1))        





    def test_transmissive(self):


        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Generic_Domain(points, elements)
        domain.check_integrity()

        domain.conserved_quantities = ['stage', 'ymomentum']
        domain.evolved_quantities = ['stage', 'ymomentum']        
        domain.quantities['stage'] =\
                                   Quantity(domain, [[1,2,3], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])

        domain.quantities['ymomentum'] =\
                                   Quantity(domain, [[2,3,4], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])


        domain.check_integrity()

        #Test transmissve bdry
        try:
            T = Transmissive_boundary()
        except:
            pass
        else:
            raise Exception('Should have raised exception')

        T = Transmissive_boundary(domain)

        from anuga.config import default_boundary_tag
        domain.set_boundary( {default_boundary_tag: T} )


        #FIXME: should not necessarily be true always.
        #E.g. with None as a boundary object.
        assert len(domain.boundary) == len(domain.boundary_objects)

        q = T.evaluate(0, 2)  #Vol=0, edge=2

        assert num.allclose(q, [1.5, 2.5])


        # Now set the centroid_transmissive_bc flag to true
        domain.set_centroid_transmissive_bc(True)

        q = T.evaluate(0, 2)  #Vol=0, edge=2

        assert num.allclose(q, [2.0 ,3.0]) # centroid value



        

    def NOtest_fileboundary_time_only(self):
        """Test that boundary values can be read from file and interpolated
        This is using the .tms file format
        
        See also test_util for comprenhensive testing of the underlying 
        file_function and also tests in test_datamanager which tests 
        file_function using the sts format
        """
        #FIXME (Ole): This test was disabled 18 August 2008 as no
        # need for this was found. Rather I implemented an Exception
        # to catch possible errors in the model setup
        

        import time, os
        from math import sin, pi
        from anuga.config import time_format

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0, 0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0, 0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Generic_Domain(points, elements)
        domain.conserved_quantities = ['stage', 'ymomentum']
        domain.quantities['stage'] =\
                                   Quantity(domain, [[1,2,3], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])

        domain.quantities['ymomentum'] =\
                                   Quantity(domain, [[2,3,4], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])

        domain.check_integrity()


        #Write file
        filename = 'boundarytest' + str(time.time())
        fid = open(filename + '.txt', 'w')
        start = time.mktime(time.strptime('2000', '%Y'))
        dt = 5*60  #Five minute intervals
        for i in range(10):
            t = start + i*dt
            t_string = time.strftime(time_format, time.gmtime(t))

            fid.write('%s,%f %f\n' %(t_string, 1.0*i, sin(i*2*pi/10)))
        fid.close()


        #Convert ASCII file to NetCDF (Which is what we really like!)
        
        from anuga.shallow_water.data_manager import timefile2netcdf
        
        timefile2netcdf(filename, quantity_names = ['stage', 'ymomentum'])
        


        F = File_boundary(filename + '.tms', domain)

        
        os.remove(filename + '.txt')
        os.remove(filename + '.tms')        


        

        #Check that midpoint coordinates at boundary are correctly computed
        assert num.allclose( F.midpoint_coordinates,
                             [[1.0, 0.0], [0.0, 1.0], [3.0, 0.0],
                              [3.0, 1.0], [1.0, 3.0], [0.0, 3.0]])

        #assert allclose(F.midpoint_coordinates[(3,2)], [0.0, 3.0])
        #assert allclose(F.midpoint_coordinates[(3,1)], [1.0, 3.0])
        #assert allclose(F.midpoint_coordinates[(0,2)], [0.0, 1.0])
        #assert allclose(F.midpoint_coordinates[(0,0)], [1.0, 0.0])
        #assert allclose(F.midpoint_coordinates[(2,0)], [3.0, 0.0])
        #assert allclose(F.midpoint_coordinates[(2,1)], [3.0, 1.0])


        #Check time interpolation
        from anuga.config import default_boundary_tag
        domain.set_boundary( {default_boundary_tag: F} )

        domain.time = 5*30/2  #A quarter way through first step
        q = F.evaluate()
        assert num.allclose(q, [1.0/4, sin(2*pi/10)/4])


        domain.time = 2.5*5*60  #Half way between steps 2 and 3
        q = F.evaluate()
        assert num.allclose(q, [2.5, (sin(2*2*pi/10) + sin(3*2*pi/10))/2])



    def test_fileboundary_exception(self):
        """Test that boundary object complains if number of
        conserved quantities are wrong
        """


        import time, os
        from math import sin, pi
        from anuga.config import time_format

        a = [0.0, 0.0]
        b = [0.0, 2.0]
        c = [2.0,0.0]
        d = [0.0, 4.0]
        e = [2.0, 2.0]
        f = [4.0,0.0]

        points = [a, b, c, d, e, f]

        #bac, bce, ecf, dbe
        elements = [ [1,0,2], [1,2,4], [4,2,5], [3,1,4] ]

        domain = Generic_Domain(points, elements)
        domain.conserved_quantities = ['stage', 'xmomentum', 'ymomentum']
        domain.evolved_quantities = ['stage', 'xmomentum', 'ymomentum']
        domain.quantities['stage'] =\
                                   Quantity(domain, [[1,2,3], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])

        domain.quantities['xmomentum'] =\
                                   Quantity(domain, [[2,3,4], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])
        domain.quantities['ymomentum'] =\
                                   Quantity(domain, [[2,3,4], [5,5,5],
                                                     [0,0,9], [-6, 3, 3]])

        domain.check_integrity()

        #Write file (with only two values)
        filename = 'boundarytest' + str(time.time())
        fid = open(filename + '.txt', 'w')
        start = time.mktime(time.strptime('2000', '%Y'))
        dt = 5*60  #Five minute intervals
        for i in range(10):
            t = start + i*dt
            t_string = time.strftime(time_format, time.gmtime(t))

            fid.write('%s,%f %f\n' %(t_string, 1.0*i, sin(i*2*pi/10)))
        fid.close()


        #Convert ASCII file to NetCDF (Which is what we really like!)
        from anuga.file_conversion.file_conversion import timefile2netcdf
        
        timefile2netcdf(filename+'.txt', quantity_names = ['stage', 'xmomentum'])

        
        try:
            F = File_boundary(filename + '.tms',
                              domain)
        except:
            pass
        else:
            raise Exception('Should have raised an exception')
        
        os.remove(filename + '.txt')
        os.remove(filename + '.tms')        


#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_Generic_Boundary_Conditions, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
