#!/usr/bin/env python

from anuga.culvert_flows.culvert_polygons import *

import unittest
import os.path

from anuga.geometry.polygon import inside_polygon, polygon_area


class Test_poly(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_1(self):
        
        end_point0=[307138.813,6193474]
        end_point1=[307150.563,6193469]
        width=3 
        height=3
        number_of_barrels=1

        P = create_culvert_polygons(end_point0,
                                    end_point1,
                                    width=width,   
                                    height=height,
                                    number_of_barrels=number_of_barrels)
        
        # Compute area and check that it is greater than 0
        for key1 in ['exchange_polygon0',
                     'exchange_polygon1']:
            polygon = P[key1]
            area = polygon_area(polygon)
            
            msg = 'Polygon %s ' %(polygon)
            msg += ' has area = %f' % area
            assert area > 0.0, msg

            for key2 in ['enquiry_point0', 'enquiry_point1']:
                point = P[key2]
                assert not inside_polygon(point, polygon)
        

    def test_2(self):
        #end_point0=[307138.813,6193474]
        #end_point1=[307150.563,6193469]
        end_point0=[10., 5.]
        end_point1=[10., 10.]     
        width = 1
        height = 3.5 
        number_of_barrels=1

        P = create_culvert_polygons(end_point0,
                                    end_point1,
                                    width=width,   
                                    height=height,
                                    number_of_barrels=number_of_barrels)
        
        # Compute area and check that it is greater than 0
        for key1 in ['exchange_polygon0',
                    'exchange_polygon1']:
            polygon = P[key1]
            area = polygon_area(polygon)
            
            msg = 'Polygon %s ' % (polygon)
            msg += ' has area = %f' % area
            assert area > 0.0, msg

            for key2 in ['enquiry_point0', 'enquiry_point1']:
                point = P[key2]
                assert not inside_polygon(point, polygon)                        

    
               
#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_poly, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)
        
