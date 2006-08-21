#!/usr/bin/env python


import unittest
import os
from Numeric import zeros, array, allclose, concatenate
from math import sqrt, pi
import tempfile

from anuga.geospatial_data.geospatial_data import *
from anuga.coordinate_transforms.geo_reference import Geo_reference, TitleError
from anuga.coordinate_transforms.redfearn import degminsec2decimal_degrees


class Test_Geospatial_data(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass


    def test_0(self):
        """Basic points
        """
        from anuga.coordinate_transforms.geo_reference import Geo_reference
        
        points = [[1.0, 2.1], [3.0, 5.3]]
        G = Geospatial_data(points)

        assert allclose(G.data_points, [[1.0, 2.1], [3.0, 5.3]])

        #Check defaults
        assert G.attributes is None
        
        assert G.geo_reference.zone == Geo_reference().zone
        assert G.geo_reference.xllcorner == Geo_reference().xllcorner
        assert G.geo_reference.yllcorner == Geo_reference().yllcorner
        

    def test_1(self):
        points = [[1.0, 2.1], [3.0, 5.3]]
        attributes = [2, 4]
        G = Geospatial_data(points, attributes)        

        assert G.attributes.keys()[0] == 'attribute'
        assert allclose(G.attributes.values()[0], [2, 4])
        

    def test_2(self):
        from anuga.coordinate_transforms.geo_reference import Geo_reference
        points = [[1.0, 2.1], [3.0, 5.3]]
        attributes = [2, 4]
        G = Geospatial_data(points, attributes,
                            geo_reference=Geo_reference(56, 100, 200))

        assert G.geo_reference.zone == 56
        assert G.geo_reference.xllcorner == 100
        assert G.geo_reference.yllcorner == 200


    def test_get_attributes_1(self):
        from anuga.coordinate_transforms.geo_reference import Geo_reference
        points = [[1.0, 2.1], [3.0, 5.3]]
        attributes = [2, 4]
        G = Geospatial_data(points, attributes,
                            geo_reference=Geo_reference(56, 100, 200))


        P = G.get_data_points(absolute=False)
        assert allclose(P, [[1.0, 2.1], [3.0, 5.3]])        

        P = G.get_data_points(absolute=True)
        assert allclose(P, [[101.0, 202.1], [103.0, 205.3]])        

        V = G.get_attributes() #Simply get them
        assert allclose(V, [2, 4])

        V = G.get_attributes('attribute') #Get by name
        assert allclose(V, [2, 4])

    def test_get_attributes_2(self):
        """Multiple attributes
        """
        
        from anuga.coordinate_transforms.geo_reference import Geo_reference
        points = [[1.0, 2.1], [3.0, 5.3]]
        attributes = {'a0': [0, 0], 'a1': [2, 4], 'a2': [79.4, -7]}
        G = Geospatial_data(points, attributes,
                            geo_reference=Geo_reference(56, 100, 200),
                            default_attribute_name='a1')


        P = G.get_data_points(absolute=False)
        assert allclose(P, [[1.0, 2.1], [3.0, 5.3]])        
        
        V = G.get_attributes() #Get default attribute
        assert allclose(V, [2, 4])

        V = G.get_attributes('a0') #Get by name
        assert allclose(V, [0, 0])

        V = G.get_attributes('a1') #Get by name
        assert allclose(V, [2, 4])

        V = G.get_attributes('a2') #Get by name
        assert allclose(V, [79.4, -7])

        try:
            V = G.get_attributes('hdnoatedu') #Invalid
        except AssertionError:
            pass
        else:
            raise 'Should have raised exception' 

    def test_get_data_points(self):
        points_ab = [[12.5,34.7],[-4.5,-60.0]]
        x_p = -10
        y_p = -40
        geo_ref = Geo_reference(56, x_p, y_p)
        points_rel = geo_ref.change_points_geo_ref(points_ab)
        
        spatial = Geospatial_data(points_rel, geo_reference=geo_ref)

        results = spatial.get_data_points(absolute=False)
        
        assert allclose(results, points_rel)
        
        x_p = -1770
        y_p = 4.01
        geo_ref = Geo_reference(56, x_p, y_p)
        points_rel = geo_ref.change_points_geo_ref(points_ab)
        results = spatial.get_data_points \
                  ( geo_reference=geo_ref)
        
        assert allclose(results, points_rel)

        
    def test_set_geo_reference(self):
        points_ab = [[12.5,34.7],[-4.5,-60.0]]
        x_p = -10
        y_p = -40
        geo_ref = Geo_reference(56, x_p, y_p)
        points_rel = geo_ref.change_points_geo_ref(points_ab)
        spatial = Geospatial_data(points_rel)
        # since the geo_ref wasn't set
        assert not allclose( points_ab, spatial.get_data_points(absolute=True))
        
        spatial = Geospatial_data(points_rel, geo_reference=geo_ref)
        assert allclose( points_ab, spatial.get_data_points(absolute=True))
        
        x_p = 10
        y_p = 400
        new_geo_ref = Geo_reference(56, x_p, y_p)
        spatial.set_geo_reference(new_geo_ref)
        assert allclose( points_ab, spatial.get_data_points(absolute=True))
        
        
        
    def test_conversions_to_points_dict(self):
        """test conversions to points_dict
        """
        
        from anuga.coordinate_transforms.geo_reference import Geo_reference
        points = [[1.0, 2.1], [3.0, 5.3]]
        attributes = {'a0': [0, 0], 'a1': [2, 4], 'a2': [79.4, -7]}
        G = Geospatial_data(points, attributes,
                            geo_reference=Geo_reference(56, 100, 200),
                            default_attribute_name='a1')


        points_dict = geospatial_data2points_dictionary(G)

        assert points_dict.has_key('pointlist')
        assert points_dict.has_key('attributelist')        
        assert points_dict.has_key('geo_reference')

        assert allclose( points_dict['pointlist'], points )

        A = points_dict['attributelist']
        assert A.has_key('a0')
        assert A.has_key('a1')
        assert A.has_key('a2')        

        assert allclose( A['a0'], [0, 0] )
        assert allclose( A['a1'], [2, 4] )        
        assert allclose( A['a2'], [79.4, -7] )


        geo = points_dict['geo_reference']
        assert geo is G.geo_reference


    def test_conversions_from_points_dict(self):
        """test conversions from points_dict
        """

        from anuga.coordinate_transforms.geo_reference import Geo_reference
        
        points = [[1.0, 2.1], [3.0, 5.3]]
        attributes = {'a0': [0, 0], 'a1': [2, 4], 'a2': [79.4, -7]}

        points_dict = {}
        points_dict['pointlist'] = points
        points_dict['attributelist'] = attributes
        points_dict['geo_reference'] = Geo_reference(56, 100, 200)
        

        G = points_dictionary2geospatial_data(points_dict)

        P = G.get_data_points(absolute=False)
        assert allclose(P, [[1.0, 2.1], [3.0, 5.3]])        
        
        #V = G.get_attribute_values() #Get default attribute
        #assert allclose(V, [2, 4])

        V = G.get_attributes('a0') #Get by name
        assert allclose(V, [0, 0])

        V = G.get_attributes('a1') #Get by name
        assert allclose(V, [2, 4])

        V = G.get_attributes('a2') #Get by name
        assert allclose(V, [79.4, -7])

    def test_add(self):
        """ test the addition of two geospatical objects
            no geo_reference see next test
        """
        points = [[1.0, 2.1], [3.0, 5.3]]
        attributes = {'depth':[2, 4], 'elevation':[6.1, 5]}
        attributes1 = {'depth':[2, 4], 'elevation':[2.5, 1]}
        G1 = Geospatial_data(points, attributes)        
        G2 = Geospatial_data(points, attributes1) 
        
#        g3 = geospatial_data2points_dictionary(G1)
#        print 'g3=', g3
        
        G = G1 + G2

        assert G.attributes.has_key('depth')
        assert G.attributes.has_key('elevation')
        assert allclose(G.attributes['depth'], [2, 4, 2, 4])
        assert allclose(G.attributes['elevation'], [6.1, 5, 2.5, 1])
        assert allclose(G.get_data_points(), [[1.0, 2.1], [3.0, 5.3],
                                              [1.0, 2.1], [3.0, 5.3]])
        
    def test_add_with_geo (self):
        """
        Difference in Geo_reference resolved
        """
        points1 = [[1.0, 2.1], [3.0, 5.3]]
        points2 = [[5.0, 6.1], [6.0, 3.3]]
        attributes1 = [2, 4]
        attributes2 = [5, 76]
        geo_ref1= Geo_reference(55, 1.0, 2.0)
        geo_ref2 = Geo_reference(zone=55,
                                 xllcorner=0.1,
                                 yllcorner=3.0,
                                 datum='wgs84',
                                 projection='UTM',
                                 units='m')
                                
        G1 = Geospatial_data(points1, attributes1, geo_ref1)
        G2 = Geospatial_data(points2, attributes2, geo_ref2)

        #Check that absolute values are as expected
        P1 = G1.get_data_points(absolute=True)
        assert allclose(P1, [[2.0, 4.1], [4.0, 7.3]])

        P2 = G2.get_data_points(absolute=True)
        assert allclose(P2, [[5.1, 9.1], [6.1, 6.3]])        
        
        G = G1 + G2
        
        assert allclose(G.get_geo_reference().get_xllcorner(), 0.1)
        assert allclose(G.get_geo_reference().get_yllcorner(), 2.0)

        P = G.get_data_points(absolute=True)

        P_relative = G.get_data_points(absolute=False)
        
        assert allclose(P_relative, P - [0.1, 2.0])

        assert allclose(P, concatenate( (P1,P2) ))
        assert allclose(P, [[2.0, 4.1], [4.0, 7.3],
                            [5.1, 9.1], [6.1, 6.3]])
        


        

    def test_add_with_geo_absolute (self):
        """
        Difference in Geo_reference resolved
        """
        points1 = array([[2.0, 4.1], [4.0, 7.3]])
        points2 = array([[5.1, 9.1], [6.1, 6.3]])        
        attributes1 = [2, 4]
        attributes2 = [5, 76]
        geo_ref1= Geo_reference(55, 1.0, 2.0)
        geo_ref2 = Geo_reference(55, 2.0, 3.0)

        
                                
        G1 = Geospatial_data(points1 - [geo_ref1.get_xllcorner(), geo_ref1.get_yllcorner()],
                             attributes1, geo_ref1)
        
        G2 = Geospatial_data(points2 - [geo_ref2.get_xllcorner(), geo_ref2.get_yllcorner()],
                             attributes2, geo_ref2)

        #Check that absolute values are as expected
        P1 = G1.get_data_points(absolute=True)
        assert allclose(P1, points1)

        P1 = G1.get_data_points(absolute=False)
        assert allclose(P1, points1 - [geo_ref1.get_xllcorner(), geo_ref1.get_yllcorner()])        

        P2 = G2.get_data_points(absolute=True)
        assert allclose(P2, points2)

        P2 = G2.get_data_points(absolute=False)
        assert allclose(P2, points2 - [geo_ref2.get_xllcorner(), geo_ref2.get_yllcorner()])                
        
        G = G1 + G2
        
        assert allclose(G.get_geo_reference().get_xllcorner(), 1.0)
        assert allclose(G.get_geo_reference().get_yllcorner(), 2.0)

        P = G.get_data_points(absolute=True)

        P_relative = G.get_data_points(absolute=False)
        
        assert allclose(P_relative, [[1.0, 2.1], [3.0, 5.3], [4.1, 7.1], [5.1, 4.3]])

        assert allclose(P, concatenate( (points1,points2) ))
                           
        


    def test_create_from_xya_file(self):
        """Check that object can be created from a points file (.pts and .xya)
        """

        points = [[1.0, 2.1], [3.0, 5.3], [5.0, 6.1], [6.0, 3.3]]
        attributes = [2, 4, 5, 76]
        '''
        # Use old pointsdict format
        pointsdict = {'pointlist': points,
                      'attributelist': {'att1': attributes,
                                        'att2': array(attributes) + 1}} 
        '''
        att_dict = {'att1': attributes,
                    'att2': array(attributes) +1}
                    
        # Create points as an xya file 
        FN = 'test_points.xya'
        G1 = Geospatial_data(points, att_dict)
        G1.export_points_file(FN)
#        G1.export_points_file(ofile)

        #Create object from file
        G = Geospatial_data(file_name = FN)
        
        assert allclose(G.get_data_points(), points)
        assert allclose(G.get_attributes('att1'), attributes)
        assert allclose(G.get_attributes('att2'), array(attributes) + 1)
        
        os.remove(FN)

    def test_create_from_xya_file1(self):
        """
        Check that object can be created from an Absolute xya file
        """

        points = [[1.0, 2.1], [3.0, 5.3], [5.0, 6.1], [6.0, 3.3]]
        attributes = [2, 4, 5, 76]

        att_dict = {'att1': attributes,
                    'att2': array(attributes) +1}

        geo_ref = Geo_reference(56, 10, 5)
                    
        # Create points as an xya file 
        FN = 'test_points.xya'
        G1 = Geospatial_data(points, att_dict, geo_ref)

        G1.export_points_file(FN, absolute=True)

        #Create object from file
        G = Geospatial_data(file_name = FN)
        
        assert allclose(G.get_data_points(absolute=True), 
                        G1.get_data_points(absolute=True))
        assert allclose(G.get_attributes('att1'), attributes)
        assert allclose(G.get_attributes('att2'), array(attributes) + 1)
        
        os.remove(FN)
        
    def test_loadxya(self):
        """
        comma delimited
        """
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("elevation  , speed \n\
1.0, 0.0, 10.0, 0.0\n\
0.0, 1.0, 0.0, 10.0\n\
1.0, 0.0, 10.4, 40.0\n")
        file.close()
        results = Geospatial_data(fileName, delimiter=',')
        os.remove(fileName)
#        print 'data', results.get_data_points()
        assert allclose(results.get_data_points(), [[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        assert allclose(results.get_attributes(attribute_name='elevation'), [10.0, 0.0, 10.4])
        assert allclose(results.get_attributes(attribute_name='speed'), [0.0, 10.0, 40.0])

    def test_loadxya2(self):
        """
        space delimited
        """
        import os
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation   speed \n\
1.0 0.0 10.0 0.0\n\
0.0 1.0 0.0 10.0\n\
1.0 0.0 10.4 40.0\n")
        file.close()

        results = Geospatial_data(fileName, delimiter=' ')

        os.remove(fileName)

        assert allclose(results.get_data_points(), [[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        assert allclose(results.get_attributes(attribute_name='elevation'), [10.0, 0.0, 10.4])
        assert allclose(results.get_attributes(attribute_name='speed'), [0.0, 10.0, 40.0])
     
    def test_loadxya3(self):
        """
        space delimited
        """
        import os
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation   speed \n\
1.0 0.0 10.0 0.0\n\
0.0 1.0 0.0 10.0\n\
1.0 0.0 10.4 40.0\n\
#geocrap\n\
56\n\
56.6\n\
3\n")
        file.close()

        results = Geospatial_data(fileName, delimiter=' ')

        os.remove(fileName)
        assert allclose(results.get_data_points(absolute=False), [[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        assert allclose(results.get_attributes(attribute_name='elevation'), [10.0, 0.0, 10.4])
        assert allclose(results.get_attributes(attribute_name='speed'), [0.0, 10.0, 40.0])

    def BADtest_loadxya4(self):
        """
        comma delimited
        """
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("elevation  , speed \n\
1.0, 0.0, splat, 0.0\n\
0.0, 1.0, 0.0, 10.0\n\
1.0, 0.0, 10.4, 40.0\n")
        file.close()
        results = Geospatial_data(fileName, delimiter=',')
        os.remove(fileName)
#        print 'data', results.get_data_points()
        assert allclose(results.get_data_points(), [[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        assert allclose(results.get_attributes(attribute_name='elevation'), ["splat", 0.0, 10.4])
        assert allclose(results.get_attributes(attribute_name='speed'), [0.0, 10.0, 40.0])
        
    def test_read_write_points_file_bad2(self):
        att_dict = {}
        pointlist = array([[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        att_dict['elevation'] = array([10.0, 0.0, 10.4])
        att_dict['brightness'] = array([10.0, 0.0, 10.4])
        geo_reference=Geo_reference(56,1.9,1.9)
        
        G = Geospatial_data(pointlist, att_dict, geo_reference)
        
        try:
            G.export_points_file("_???/yeah.xya")
            
        except IOError:
            pass
        else:
            msg = 'bad points file extension did not raise error!'
            raise msg
#            self.failUnless(0 == 1,
#                        'bad points file extension did not raise error!')
                   
    def test_loadxy_bad(self):
        import os
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation   \n\
1.0 0.0 10.0 0.0\n\
0.0 1.0 0.0 10.0\n\
1.0 0.0 10.4 40.0\n")
        file.close()
        #print fileName
        try:
            results = Geospatial_data(fileName, delimiter=' ')
        except IOError:
            pass
        else:
            msg = 'bad xya file did not raise error!'
            raise msg
#            self.failUnless(0 == 1,
#                        'bad xya file did not raise error!')
        os.remove(fileName)
       
    def test_loadxy_bad2(self):
        import os
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("elevation\n\
1.0 0.0 10.0 \n\
0.0 1.0\n\
1.0 \n")
        file.close()
        #print fileName
        try:
            results = Geospatial_data(fileName, delimiter=' ')
        except IOError:
            pass
        else:
            msg = 'bad xya file did not raise error!'
            raise msg
        os.remove(fileName)
   
    def test_loadxy_bad3(self):
        """
        specifying wrong delimiter
        """
        import os
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation  , speed \n\
1.0, 0.0, 10.0, 0.0\n\
0.0, 1.0, 0.0, 10.0\n\
1.0, 0.0, 10.4, 40.0\n")
        file.close()
        try:
            results = Geospatial_data(fileName, delimiter=' ')
        except IOError:
            pass
        else:
            msg = 'bad xya file did not raise error!'
            raise msg
        os.remove(fileName)
     
    def test_loadxy_bad4(self):
        """
         specifying wrong delimiter
        """
        import os
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation   speed \n\
1.0 0.0 10.0 0.0\n\
0.0 1.0 0.0 10.0\n\
1.0 0.0 10.4 40.0\n\
#geocrap\n\
56\n\
56.6\n\
3\n"
)
        file.close()
        try:
            results = Geospatial_data(fileName, delimiter=',')
        except IOError:
            pass
        else:
            msg = 'bad xya file did not raise error!'
            raise msg

        os.remove(fileName)

    def test_loadxy_bad5(self):
        """
        specifying wrong delimiter
        """
        import os
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation   speed \n\
1.0 0.0 10.0 0.0\n\
0.0 1.0 0.0 10.0\n\
1.0 0.0 10.4 40.0\n\
#geocrap\n\
crap")
        file.close()
        try:
#            dict = import_points_file(fileName,delimiter=' ')
#            results = Geospatial_data()
            results = Geospatial_data(fileName, delimiter=' ', verbose=True)
#            results.import_points_file(fileName, delimiter=' ')
        except IOError:
            pass
        else:
            msg = 'bad xya file did not raise error!'
            raise msg

#            self.failUnless(0 ==1,
#                        'bad xya file did not raise error!')    
        os.remove(fileName)              

    def test_loadxy_bad_no_file_xya(self):
        import os
       
        fileName = tempfile.mktemp(".xya")
        try:
            results = Geospatial_data(fileName, delimiter=' ')
        except IOError:
            pass
        else:
            msg = 'imaginary file did not raise error!'
            raise msg

#        except IOError:
#            pass
#        else:
#            self.failUnless(0 == 1,
#                        'imaginary file did not raise error!')

                        
  ###################### .XYA ##############################
        
    def test_export_xya_file(self):
#        dict = {}
        att_dict = {}
        pointlist = array([[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        att_dict['elevation'] = array([10.0, 0.0, 10.4])
        att_dict['brightness'] = array([10.0, 0.0, 10.4])
#        dict['attributelist'] = att_dict
        geo_reference=Geo_reference(56,1.9,1.9)
        
        
        fileName = tempfile.mktemp(".xya")
        G = Geospatial_data(pointlist, att_dict, geo_reference)
        G.export_points_file(fileName, False)

#        dict2 = import_points_file(fileName)
        results = Geospatial_data(file_name = fileName)
        #print "fileName",fileName 
        os.remove(fileName)
        
        assert allclose(results.get_data_points(absolute=False),[[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        assert allclose(results.get_attributes(attribute_name='elevation'), [10.0, 0.0, 10.4])
        answer = [10.0, 0.0, 10.4]
        assert allclose(results.get_attributes(attribute_name='brightness'), answer)
        #print "dict2['geo_reference']",dict2['geo_reference'] 
        self.failUnless(results.get_geo_reference() == geo_reference,
                         'test_writepts failed. Test geo_reference')

    def test_export_xya_file2(self):
        """test absolute xya file
        """
        att_dict = {}
        pointlist = array([[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        att_dict['elevation'] = array([10.0, 0.0, 10.4])
        att_dict['brightness'] = array([10.0, 0.0, 10.4])
        
        fileName = tempfile.mktemp(".xya")
        G = Geospatial_data(pointlist, att_dict)
        G.export_points_file(fileName)
        results = Geospatial_data(file_name = fileName)
#        dict2 = import_points_file(fileName)
        os.remove(fileName)
        
        assert allclose(results.get_data_points(False),[[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        assert allclose(results.get_attributes('elevation'), [10.0, 0.0, 10.4])
        answer = [10.0, 0.0, 10.4]
        assert allclose(results.get_attributes('brightness'), answer)

    def test_export_xya_file3(self):
        """test absolute xya file with geo_ref
        """
        att_dict = {}
        pointlist = array([[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        att_dict['elevation'] = array([10.0, 0.0, 10.4])
        att_dict['brightness'] = array([10.0, 0.0, 10.4])
        geo_reference=Geo_reference(56,1.9,1.9)
        
        
        fileName = tempfile.mktemp(".xya")
        G = Geospatial_data(pointlist, att_dict, geo_reference)
        
        G.export_points_file(fileName, absolute=True)
        
        results = Geospatial_data(file_name = fileName)
        os.remove(fileName)

        assert allclose(results.get_data_points(),
                        [[2.9, 1.9],[1.9, 2.9],[2.9, 1.9]])
        assert allclose(results.get_attributes(attribute_name='elevation'),
                         [10.0, 0.0, 10.4])
        answer = [10.0, 0.0, 10.4]
        assert allclose(results.get_attributes(attribute_name='brightness'), answer)
        self.failUnless(results.get_geo_reference() == geo_reference,
                         'test_writepts failed. Test geo_reference')                         
                        
                        
                        
    def test_new_export_pts_file(self):
        att_dict = {}
        pointlist = array([[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        att_dict['elevation'] = array([10.1, 0.0, 10.4])
        att_dict['brightness'] = array([10.0, 1.0, 10.4])
        
        fileName = tempfile.mktemp(".pts")
        
        G = Geospatial_data(pointlist, att_dict)
        
        G.export_points_file(fileName)

        results = Geospatial_data(file_name = fileName)

        os.remove(fileName)
        
        assert allclose(results.get_data_points(),[[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        assert allclose(results.get_attributes(attribute_name='elevation'), [10.1, 0.0, 10.4])
        answer = [10.0, 1.0, 10.4]
        assert allclose(results.get_attributes(attribute_name='brightness'), answer)

    def test_new_export_absolute_pts_file(self):
        att_dict = {}
        pointlist = array([[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        att_dict['elevation'] = array([10.1, 0.0, 10.4])
        att_dict['brightness'] = array([10.0, 1.0, 10.4])
        geo_ref = Geo_reference(50, 25, 55)
        
        fileName = tempfile.mktemp(".pts")
        
        G = Geospatial_data(pointlist, att_dict, geo_ref)
        
        G.export_points_file(fileName, absolute=True)

        results = Geospatial_data(file_name = fileName)

        os.remove(fileName)
        
        assert allclose(results.get_data_points(), G.get_data_points(True))
        assert allclose(results.get_attributes(attribute_name='elevation'), [10.1, 0.0, 10.4])
        answer = [10.0, 1.0, 10.4]
        assert allclose(results.get_attributes(attribute_name='brightness'), answer)

    def test_loadpts(self):
        
        from Scientific.IO.NetCDF import NetCDFFile

        fileName = tempfile.mktemp(".pts")
        # NetCDF file definition
        outfile = NetCDFFile(fileName, 'w')
        
        # dimension definitions
        outfile.createDimension('number_of_points', 3)    
        outfile.createDimension('number_of_dimensions', 2) #This is 2d data
    
        # variable definitions
        outfile.createVariable('points', Float, ('number_of_points',
                                                 'number_of_dimensions'))
        outfile.createVariable('elevation', Float, ('number_of_points',))
    
        # Get handles to the variables
        points = outfile.variables['points']
        elevation = outfile.variables['elevation']
 
        points[0, :] = [1.0,0.0]
        elevation[0] = 10.0 
        points[1, :] = [0.0,1.0]
        elevation[1] = 0.0  
        points[2, :] = [1.0,0.0]
        elevation[2] = 10.4    

        outfile.close()
        
        results = Geospatial_data(file_name = fileName)
        os.remove(fileName)
        answer =  [[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]]
        assert allclose(results.get_data_points(), [[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        assert allclose(results.get_attributes(attribute_name='elevation'), [10.0, 0.0, 10.4])
        
    def test_writepts(self):
        att_dict = {}
        pointlist = array([[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        att_dict['elevation'] = array([10.0, 0.0, 10.4])
        att_dict['brightness'] = array([10.0, 0.0, 10.4])
        geo_reference=Geo_reference(56,1.9,1.9)
        
        fileName = tempfile.mktemp(".pts")
        
        G = Geospatial_data(pointlist, att_dict, geo_reference)
        
        G.export_points_file(fileName, False)
        
        results = Geospatial_data(file_name = fileName)

        os.remove(fileName)

        assert allclose(results.get_data_points(False),[[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        assert allclose(results.get_attributes('elevation'), [10.0, 0.0, 10.4])
        answer = [10.0, 0.0, 10.4]
        assert allclose(results.get_attributes('brightness'), answer)

        
        self.failUnless(geo_reference == geo_reference,
                         'test_writepts failed. Test geo_reference')
        
 ########################## BAD .PTS ##########################          

    def test_load_bad_no_file_pts(self):
        import os
        import tempfile
       
        fileName = tempfile.mktemp(".pts")
        #print fileName
        try:
            results = Geospatial_data(file_name = fileName)
#            dict = import_points_file(fileName)
        except IOError:
            pass
        else:
            msg = 'imaginary file did not raise error!'
            raise msg
#            self.failUnless(0 == 1,
#                        'imaginary file did not raise error!')


    def test_create_from_pts_file(self):
        
        from Scientific.IO.NetCDF import NetCDFFile

#        fileName = tempfile.mktemp(".pts")
        FN = 'test_points.pts'
        # NetCDF file definition
        outfile = NetCDFFile(FN, 'w')
        
        # dimension definitions
        outfile.createDimension('number_of_points', 3)    
        outfile.createDimension('number_of_dimensions', 2) #This is 2d data
    
        # variable definitions
        outfile.createVariable('points', Float, ('number_of_points',
                                                 'number_of_dimensions'))
        outfile.createVariable('elevation', Float, ('number_of_points',))
    
        # Get handles to the variables
        points = outfile.variables['points']
        elevation = outfile.variables['elevation']
 
        points[0, :] = [1.0,0.0]
        elevation[0] = 10.0 
        points[1, :] = [0.0,1.0]
        elevation[1] = 0.0  
        points[2, :] = [1.0,0.0]
        elevation[2] = 10.4    

        outfile.close()

        G = Geospatial_data(file_name = FN)

        assert allclose(G.get_geo_reference().get_xllcorner(), 0.0)
        assert allclose(G.get_geo_reference().get_yllcorner(), 0.0)

        assert allclose(G.get_data_points(), [[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        assert allclose(G.get_attributes(), [10.0, 0.0, 10.4])
        os.remove(FN)

    def test_create_from_pts_file_with_geo(self):
        """This test reveals if Geospatial data is correctly instantiated from a pts file.
        """
        
        from Scientific.IO.NetCDF import NetCDFFile

        FN = 'test_points.pts'
        # NetCDF file definition
        outfile = NetCDFFile(FN, 'w')

        # Make up an arbitrary georef
        xll = 0.1
        yll = 20
        geo_reference=Geo_reference(56, xll, yll)
        geo_reference.write_NetCDF(outfile)

        # dimension definitions
        outfile.createDimension('number_of_points', 3)    
        outfile.createDimension('number_of_dimensions', 2) #This is 2d data
    
        # variable definitions
        outfile.createVariable('points', Float, ('number_of_points',
                                                 'number_of_dimensions'))
        outfile.createVariable('elevation', Float, ('number_of_points',))
    
        # Get handles to the variables
        points = outfile.variables['points']
        elevation = outfile.variables['elevation']

        points[0, :] = [1.0,0.0]
        elevation[0] = 10.0 
        points[1, :] = [0.0,1.0]
        elevation[1] = 0.0  
        points[2, :] = [1.0,0.0]
        elevation[2] = 10.4    

        outfile.close()

        G = Geospatial_data(file_name = FN)

        assert allclose(G.get_geo_reference().get_xllcorner(), xll)
        assert allclose(G.get_geo_reference().get_yllcorner(), yll)

        assert allclose(G.get_data_points(), [[1.0+xll, 0.0+yll],
                                              [0.0+xll, 1.0+yll],
                                              [1.0+xll, 0.0+yll]])
        
        assert allclose(G.get_attributes(), [10.0, 0.0, 10.4])
        os.remove(FN)

        
    def test_add_(self):
        '''adds an xya and pts files, reads the files and adds them
           checking results are correct
        '''

        # create files
        att_dict1 = {}
        pointlist1 = array([[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        att_dict1['elevation'] = array([-10.0, 0.0, 10.4])
        att_dict1['brightness'] = array([10.0, 0.0, 10.4])
        geo_reference1 = Geo_reference(56, 2.0, 1.0)
        
        att_dict2 = {}
        pointlist2 = array([[2.0, 1.0],[1.0, 2.0],[2.0, 1.0]])
        att_dict2['elevation'] = array([1.0, 15.0, 1.4])
        att_dict2['brightness'] = array([14.0, 1.0, -12.4])
        geo_reference2 = Geo_reference(56, 1.0, 2.0) 

        G1 = Geospatial_data(pointlist1, att_dict1, geo_reference1)
        G2 = Geospatial_data(pointlist2, att_dict2, geo_reference2)
        
        fileName1 = tempfile.mktemp(".xya")
        fileName2 = tempfile.mktemp(".pts")

        #makes files
        G1.export_points_file(fileName1)
        G2.export_points_file(fileName2)
        
        # add files
        
        G3 = Geospatial_data(file_name = fileName1)
        G4 = Geospatial_data(file_name = fileName2)
        
        G = G3 + G4

        
        #read results
#        print'res', G.get_data_points()
#        print'res1', G.get_data_points(False)
        assert allclose(G.get_data_points(),
                        [[ 3.0, 1.0], [ 2.0, 2.0],
                         [ 3.0, 1.0], [ 3.0, 3.0],
                         [ 2.0, 4.0], [ 3.0, 3.0]])
                         
        assert allclose(G.get_attributes(attribute_name='elevation'),
                        [-10.0, 0.0, 10.4, 1.0, 15.0, 1.4])
        
        answer = [10.0, 0.0, 10.4, 14.0, 1.0, -12.4]
        assert allclose(G.get_attributes(attribute_name='brightness'), answer)
        
        self.failUnless(G.get_geo_reference() == geo_reference1,
                         'test_writepts failed. Test geo_reference')
                         
        os.remove(fileName1)
        os.remove(fileName2)
        
    def test_ensure_absolute(self):
        points = [[2.0, 0.0],[1.0, 1.0],
                         [2.0, 0.0],[2.0, 2.0],
                         [1.0, 3.0],[2.0, 2.0]]
        new_points = ensure_absolute(points)
        
        assert allclose(new_points, points)

        points = array([[2.0, 0.0],[1.0, 1.0],
                         [2.0, 0.0],[2.0, 2.0],
                         [1.0, 3.0],[2.0, 2.0]])
        new_points = ensure_absolute(points)
        
        assert allclose(new_points, points)
        
        ab_points = array([[2.0, 0.0],[1.0, 1.0],
                         [2.0, 0.0],[2.0, 2.0],
                         [1.0, 3.0],[2.0, 2.0]])
        
        mesh_origin = (56, 290000, 618000) #zone, easting, northing

        data_points = zeros((ab_points.shape), Float)
        #Shift datapoints according to new origins
        for k in range(len(ab_points)):
            data_points[k][0] = ab_points[k][0] - mesh_origin[1]
            data_points[k][1] = ab_points[k][1] - mesh_origin[2]
        #print "data_points",data_points     
        new_points = ensure_absolute(data_points,
                                             geo_reference=mesh_origin)
        #print "new_points",new_points
        #print "ab_points",ab_points
           
        assert allclose(new_points, ab_points)

        geo = Geo_reference(56,67,-56)

        data_points = geo.change_points_geo_ref(ab_points)   
        new_points = ensure_absolute(data_points,
                                             geo_reference=geo)
        #print "new_points",new_points
        #print "ab_points",ab_points
           
        assert allclose(new_points, ab_points)


        geo_reference = Geo_reference(56, 100, 200)
        ab_points = [[1.0, 2.1], [3.0, 5.3]]
        points = geo_reference.change_points_geo_ref(ab_points)
        attributes = [2, 4]
        #print "geo in points", points
        G = Geospatial_data(points, attributes,
                            geo_reference=geo_reference)
          
        new_points = ensure_absolute(G)
        #print "new_points",new_points
        #print "ab_points",ab_points
           
        assert allclose(new_points, ab_points)

       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation   speed \n\
1.0 0.0 10.0 0.0\n\
0.0 1.0 0.0 10.0\n\
1.0 0.0 10.4 40.0\n\
#geocrap\n\
56\n\
10\n\
20\n")
        file.close()
        
        ab_points = ensure_absolute(fileName)
        actual =  [[11, 20.0],[10.0, 21.0],[11.0, 20.0]]
        assert allclose(ab_points, actual)
        os.remove(fileName)

        
    def test_ensure_geospatial(self):
        points = [[2.0, 0.0],[1.0, 1.0],
                         [2.0, 0.0],[2.0, 2.0],
                         [1.0, 3.0],[2.0, 2.0]]
        new_points = ensure_geospatial(points)
        
        assert allclose(new_points.get_data_points(absolute = True), points)

        points = array([[2.0, 0.0],[1.0, 1.0],
                         [2.0, 0.0],[2.0, 2.0],
                         [1.0, 3.0],[2.0, 2.0]])
        new_points = ensure_geospatial(points)
        
        assert allclose(new_points.get_data_points(absolute = True), points)
        
        ab_points = array([[2.0, 0.0],[1.0, 1.0],
                         [2.0, 0.0],[2.0, 2.0],
                         [1.0, 3.0],[2.0, 2.0]])
        
        mesh_origin = (56, 290000, 618000) #zone, easting, northing

        data_points = zeros((ab_points.shape), Float)
        #Shift datapoints according to new origins
        for k in range(len(ab_points)):
            data_points[k][0] = ab_points[k][0] - mesh_origin[1]
            data_points[k][1] = ab_points[k][1] - mesh_origin[2]
        #print "data_points",data_points     
        new_geospatial = ensure_geospatial(data_points,
                                             geo_reference=mesh_origin)
        new_points = new_geospatial.get_data_points(absolute=True)
        #print "new_points",new_points
        #print "ab_points",ab_points
           
        assert allclose(new_points, ab_points)

        geo = Geo_reference(56,67,-56)

        data_points = geo.change_points_geo_ref(ab_points)   
        new_geospatial = ensure_geospatial(data_points,
                                             geo_reference=geo)
        new_points = new_geospatial.get_data_points(absolute=True)
        #print "new_points",new_points
        #print "ab_points",ab_points
           
        assert allclose(new_points, ab_points)


        geo_reference = Geo_reference(56, 100, 200)
        ab_points = [[1.0, 2.1], [3.0, 5.3]]
        points = geo_reference.change_points_geo_ref(ab_points)
        attributes = [2, 4]
        #print "geo in points", points
        G = Geospatial_data(points, attributes,
                            geo_reference=geo_reference)
          
        new_geospatial  = ensure_geospatial(G)
        new_points = new_geospatial.get_data_points(absolute=True)
        #print "new_points",new_points
        #print "ab_points",ab_points
           
        assert allclose(new_points, ab_points)
        
    def test_isinstance(self):

        import os
       
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation   speed \n\
1.0 0.0 10.0 0.0\n\
0.0 1.0 0.0 10.0\n\
1.0 0.0 10.4 40.0\n\
#geocrap\n\
56\n\
56.6\n\
3\n")
        file.close()

        results = Geospatial_data(fileName)
        assert allclose(results.get_data_points(absolute=False), \
                        [[1.0, 0.0],[0.0, 1.0],[1.0, 0.0]])
        assert allclose(results.get_attributes(attribute_name='elevation'), \
                        [10.0, 0.0, 10.4])
        assert allclose(results.get_attributes(attribute_name='speed'), \
                        [0.0, 10.0, 40.0])

        os.remove(fileName)
        
    def test_delimiter(self):
        
        try:
            G = Geospatial_data(delimiter=',')
#            results = Geospatial_data(file_name = fileName)
#            dict = import_points_file(fileName)
        except ValueError:
            pass
        else:
            msg = 'Instance with No fileName but has a delimiter\
                  did not raise error!'
            raise msg

    def test_no_constructors(self):
        
        try:
            G = Geospatial_data()
#            results = Geospatial_data(file_name = fileName)
#            dict = import_points_file(fileName)
        except ValueError:
            pass
        else:
            msg = 'Instance must have a filename or data points'
            raise msg        

    def test_check_geo_reference(self):
        """
        checks geo reference details are OK. eg can be called '#geo reference'
        if not throws a clear error message
        """
        import os
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation  \n\
1.0 0.0 10.0\n\
0.0 1.0 0.0\n\
1.0 0.0 10.4\n\
#ge oreference\n\
56\n\
1.1\n\
1.0\n")

        file.close()
        results = Geospatial_data(fileName)
        assert allclose(results.get_geo_reference().get_xllcorner(), 1.1)
        assert allclose(results.get_geo_reference().get_yllcorner(), 1.0)

        os.remove(fileName)
        
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation  \n\
1.0 0.0 10.0\n\
0.0 1.0 0.0\n\
1.0 0.0 10.4\n")

        file.close()
        results = Geospatial_data(fileName)
       
        os.remove(fileName)
        
    def test_check_geo_reference1(self):
        """
        checks geo reference details are OK. eg can be called '#geo reference'
        if not throws a clear error message
        """
        import os
        fileName = tempfile.mktemp(".xya")
        file = open(fileName,"w")
        file.write("  elevation  \n\
1.0 0.0 10.0\n\
0.0 1.0 0.0\n\
1.0 0.0 10.4\n\
#geo t t\n\
56\n\
1.1\n"
)
        file.close()

        try:
            results = Geospatial_data(fileName, delimiter = " ")
        except IOError:
            pass
        else:
            msg = 'Geo reference data format is incorrect'
            raise msg        


        os.remove(fileName)


        
    def test_lat_long(self):
        lat_gong = degminsec2decimal_degrees(-34,30,0.)
        lon_gong = degminsec2decimal_degrees(150,55,0.)
        
        lat_2 = degminsec2decimal_degrees(-34,00,0.)
        lon_2 = degminsec2decimal_degrees(150,00,0.)
        
        lats = [lat_gong, lat_2]
        longs = [lon_gong, lon_2]
        gsd = Geospatial_data(latitudes=lats, longitudes=longs)

        points = gsd.get_data_points(absolute=True)
        
        assert allclose(points[0][0], 308728.009)
        assert allclose(points[0][1], 6180432.601)
        assert allclose(points[1][0],  222908.705)
        assert allclose(points[1][1], 6233785.284)
        self.failUnless(gsd.get_geo_reference().get_zone() == 56,
                        'Bad zone error!')
        
        try:
            results = Geospatial_data(latitudes=lats)
        except ValueError:
            pass
        else:
            self.failUnless(0 ==1,  'Error not thrown error!')
        try:
            results = Geospatial_data(latitudes=lats)
        except ValueError:
            pass
        else:
            self.failUnless(0 ==1,  'Error not thrown error!')
        try:
            results = Geospatial_data(longitudes=lats)
        except ValueError:
            pass
        else:
            self.failUnless(0 ==1, 'Error not thrown error!')
        try:
            results = Geospatial_data(latitudes=lats, longitudes=longs,
                                      geo_reference="p")
        except ValueError:
            pass
        else:
            self.failUnless(0 ==1,  'Error not thrown error!')
            
        try:
            results = Geospatial_data(latitudes=lats, longitudes=longs,
                                      data_points=12)
        except ValueError:
            pass
        else:
            self.failUnless(0 ==1,  'Error not thrown error!')

    def test_lat_long2(self):
        lat_gong = degminsec2decimal_degrees(-34,30,0.)
        lon_gong = degminsec2decimal_degrees(150,55,0.)
        
        lat_2 = degminsec2decimal_degrees(-34,00,0.)
        lon_2 = degminsec2decimal_degrees(150,00,0.)
        
        points = [[lat_gong, lon_gong], [lat_2, lon_2]]
        gsd = Geospatial_data(data_points=points, points_are_lats_longs=True)

        points = gsd.get_data_points(absolute=True)
        
        assert allclose(points[0][0], 308728.009)
        assert allclose(points[0][1], 6180432.601)
        assert allclose(points[1][0],  222908.705)
        assert allclose(points[1][1], 6233785.284)
        self.failUnless(gsd.get_geo_reference().get_zone() == 56,
                        'Bad zone error!')

        try:
            results = Geospatial_data(points_are_lats_longs=True)
        except ValueError:
            pass
        else:
            self.failUnless(0 ==1,  'Error not thrown error!')

    def test_len(self):
        
        points = [[1.0, 2.1], [3.0, 5.3]]
        G = Geospatial_data(points)
        self.failUnless(2 ==len(G),  'Len error!')
        
        points = [[1.0, 2.1]]
        G = Geospatial_data(points)
        self.failUnless(1 ==len(G),  'Len error!')

        points = [[1.0, 2.1], [3.0, 5.3], [3.0, 5.3], [3.0, 5.3]]
        G = Geospatial_data(points)
        self.failUnless(4 ==len(G),  'Len error!')
         
if __name__ == "__main__":

    #suite = unittest.makeSuite(Test_Geospatial_data, 'test_ensure_geospatial')
    suite = unittest.makeSuite(Test_Geospatial_data, 'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

    
