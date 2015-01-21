#!/usr/bin/env python

import unittest
from math import sqrt

from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.abstract_2d_finite_volumes.generic_domain import Generic_Domain
from anuga.abstract_2d_finite_volumes.tag_region import *
#from anuga.config import epsilon

import numpy as num


"""
This is what the mesh in these tests look like;

3---7
|5 /|
| /6|
2---6
|3 /|
| /2|
1---5
|1 /|
| /0|
0---4
"""

def add_x_y(x, y):
    return x+y

def give_me_23(x, y):
    return 23.0

class Test_tag_region(unittest.TestCase):
    def setUp(self):
        pass


    def tearDown(self):
        pass


    def test_region_tags(self):
        """get values based on triangle lists."""

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.build_tagged_elements_dictionary({'bottom': [0,1],
                                                 'top': [4,5],
                                                 'all': [0,1,2,3,4,5]})

        #Set friction
        manning = 0.07
        domain.set_quantity('friction', manning)

        a = Set_tag_region('bottom', 'friction', 0.09)
        b = Set_tag_region('top', 'friction', 1.0)
        domain.set_tag_region([a, b])

        expected = [[ 0.09,  0.09,  0.09],
                    [ 0.09,  0.09,  0.09],
                    [ 0.07,  0.07,  0.07],
                    [ 0.07,  0.07,  0.07],
                    [ 1.0,   1.0,   1.0],
                    [ 1.0,   1.0,   1.0]]
        msg = ("\ndomain.quantities['friction']=%s\nexpected value=%s"
               % (str(domain.quantities['friction'].get_values()),
                  str(expected)))
        assert num.allclose(domain.quantities['friction'].get_values(),
                            expected), msg

        #c = Add_Value_To_region('all', 'friction', 10.0)
        domain.set_tag_region(Add_value_to_region('all', 'friction', 10.0))
        #print domain.quantities['friction'].get_values()
        assert num.allclose(domain.quantities['friction'].get_values(),
                            [[ 10.09, 10.09, 10.09],
                             [ 10.09, 10.09, 10.09],
                             [ 10.07, 10.07, 10.07],
                             [ 10.07, 10.07, 10.07],
                             [ 11.0,  11.0,  11.0],
                             [ 11.0,  11.0,  11.0]])

        # trying a function
        domain.set_tag_region(Set_tag_region('top', 'friction', add_x_y))
        #print domain.quantities['friction'].get_values()
        assert num.allclose(domain.quantities['friction'].get_values(),
                            [[ 10.09, 10.09, 10.09],
                             [ 10.09, 10.09, 10.09],
                             [ 10.07, 10.07, 10.07],
                             [ 10.07, 10.07, 10.07],
                             [ 5./3,  2.0,  2./3],
                             [ 1.0,  2./3,  2.0]])

        domain.set_quantity('elevation', 10.0)
        domain.set_quantity('stage', 10.0)
        domain.set_tag_region(Add_value_to_region('top', 'stage', 1.0,initial_quantity='elevation'))
        #print domain.quantities['stage'].get_values()
        assert num.allclose(domain.quantities['stage'].get_values(),
                            [[ 10., 10., 10.],
                             [ 10., 10., 10.],
                             [ 10., 10., 10.],
                             [ 10., 10., 10.],
                             [ 11.0,  11.0,  11.0],
                             [ 11.0,  11.0,  11.0]])

        
        domain.set_quantity('elevation', 10.0)
        domain.set_quantity('stage', give_me_23)
        #this works as well, (is cleaner, but doesn't work for regions)
        #domain.set_quantity('stage',
        #                    domain.quantities['stage'].vertex_values+ \
        #                    domain.quantities['elevation'].vertex_values)
        domain.set_tag_region(Add_quantities('top', 'elevation','stage'))
        #print domain.quantities['stage'].get_values()
        assert num.allclose(domain.quantities['elevation'].get_values(),
                            [[ 10., 10., 10.],
                             [ 10., 10., 10.],
                             [ 10., 10., 10.],
                             [ 10., 10., 10.],
                             [ 33.,  33.0,  33.],
                             [ 33.0,  33.,  33.]])
        
    def test_unique_vertices(self):
        """get values based on triangle lists."""

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.build_tagged_elements_dictionary({'bottom':[0,1],
                                                 'top':[4,5],
                                                 'all':[0,1,2,3,4,5]})

        #Set friction
        manning = 0.07
        domain.set_quantity('friction', manning)

        a = Set_tag_region('bottom', 'friction', 0.09, location = 'unique vertices')
        domain.set_tag_region(a)
        assert num.allclose(domain.quantities['friction'].get_values(),
                            [[ 0.09,  0.09,  0.09],
                             [ 0.09,  0.09,  0.09],
                             [ 0.09,  0.07,  0.09],
                             [ 0.07,  0.09,  0.07],
                             [ 0.07,  0.07,  0.07],
                             [ 0.07,  0.07,  0.07]])


    def test_unique_verticesII(self):
        """
        get values based on triangle lists.
        """

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.build_tagged_elements_dictionary({'bottom':[0,1],
                                                 'top':[4,5],
                                                 'all':[0,1,2,3,4,5]})

        #Set friction
        manning = 0.07
        domain.set_quantity('friction', manning)

        domain.set_tag_region(Add_value_to_region('bottom', 'friction', 1.0,initial_quantity='friction', location = 'unique vertices'))

        #print domain.quantities['friction'].get_values()
        assert num.allclose(domain.quantities['friction'].get_values(),\
                            [[ 1.07,  1.07,  1.07],
                             [ 1.07,  1.07,  1.07],
                             [ 1.07,  0.07,  1.07],
                             [ 0.07,  1.07,  0.07],
                             [ 0.07,  0.07,  0.07],
                             [ 0.07,  0.07,  0.07]])
                         
    def test_unique_vertices_average_loc_vert(self):
        """Get values based on triangle lists."""

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)

        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.build_tagged_elements_dictionary({'bottom': [0, 1],
                                                 'top': [4, 5],
                                                 'not_bottom': [2, 3, 4, 5]})

        #Set friction
        domain.set_quantity('friction', add_x_y)
        av_bottom = 2.0 / 3.0
        add = 60.0
        calc_frict = av_bottom + add
        domain.set_tag_region(Add_value_to_region('bottom', 'friction', add,
                                              initial_quantity='friction',
                                              location='vertices',
                                              average=True))

        frict_points = domain.quantities['friction'].get_values()

        expected = [calc_frict, calc_frict, calc_frict]
        msg = ('frict_points[0]=%s\nexpected=%s'
               % (str(frict_points[0]), str(expected)))
        assert num.allclose(frict_points[0], expected), msg

        msg = ('frict_points[1]=%s\nexpected=%s'
               % (str(frict_points[1]), str(expected)))
        assert num.allclose(frict_points[1], expected), msg
  
    def test_unique_vertices_average_loc_unique_vert(self):
        """
        get values based on triangle lists.
        """

        #Create basic mesh
        points, vertices, boundary = rectangular(1, 3)
        #Create shallow water domain
        domain = Domain(points, vertices, boundary)
        domain.build_tagged_elements_dictionary({'bottom':[0,1],
                                                 'top':[4,5],
                                                 'not_bottom':[2,3,4,5]})

        #Set friction
        domain.set_quantity('friction', add_x_y)
        av_bottom = 2.0/3.0
        add = 60.0
        calc_frict = av_bottom + add
        domain.set_tag_region(Add_value_to_region('bottom', 'friction', add,
                          initial_quantity='friction',
                           location='unique vertices',
                           average=True
                          ))

        #print domain.quantities['friction'].get_values()
        frict_points = domain.quantities['friction'].get_values()
        assert num.allclose(frict_points[0],\
                            [ calc_frict, calc_frict, calc_frict])
        assert num.allclose(frict_points[1],\
                            [ calc_frict, calc_frict, calc_frict])
        assert num.allclose(frict_points[2],\
                            [ calc_frict, 1.0 + 2.0/3.0, calc_frict])
        assert num.allclose(frict_points[3],\
                            [ 2.0/3.0,calc_frict, 1.0 + 2.0/3.0])
                                                
                         
#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_tag_region, 'test')    
    runner = unittest.TextTestRunner()
    runner.run(suite)
