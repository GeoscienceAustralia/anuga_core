#!/usr/bin/env python

#TEST
import sys
import unittest

import numpy as num

try:
    from anuga.alpha_shape.alpha_shape import *
except ImportError:  
    from alpha_shape import *

class TestCase(unittest.TestCase):

    def setUp(self):
        pass

         
    def tearDown(self):
        pass

    def test_delaunay(self):
        #print "test_delaunay"
        a = [0.0, 0.0]
        b = [1.0, 0.0]
        c = [2.0, 0.0]
        d = [2.0, 2.0]
        e = [1.0, 2.0]
        f = [0.0, 2.0]

        alpha = Alpha_Shape([a,b,c,d,e,f])
        result = alpha.get_delaunay()
        answer = [(0, 1, 5), (5, 1, 4), (4, 2, 3), (2, 4, 1)]
        assert num.allclose(answer, result) 

    def test_3_points_on_line(self):
        #print "test_delaunay"
        a = [0.0, 0.0]
        b = [1.0, 0.0]
        c = [2.0, 0.0]

        try:
            alpha = Alpha_Shape([a,b,c])
        except PointError:
            pass
        else:
            self.assertTrue(0==1, \
                        'point list with 2 points did not raise an error!')



    def test_alpha_1(self):
        #print "test_alpha" 
        a = [0.0, 0.0]
        b = [1.0, 0.0]
        c = [2.0, 0.0]
        d = [2.0, 2.0]
        e = [1.0, 2.0]
        f = [0.0, 2.0]

        alpha = Alpha_Shape([a,b,c,d,e,f])
        result = alpha.get_boundary()
        #print "result",result
        answer = [(5, 0), (0, 1), (4, 5), (2, 3), (3, 4), (1, 2)]
        assert num.allclose(answer, result) 


    def test_alpha_2(self):
        #print "test_alpha" 
        a = [0.0, 0.0]
        b = [2.0, 0.0]
        c = [4.0, 0.0]
        cd = [3.0,2.0]
        d = [4.0, 4.0]
        e = [2.0, 4.0]
        f = [0.0, 4.0]

        alpha = Alpha_Shape([a,b,c,cd,d,e,f])
        result = alpha.get_boundary()
        #print "result",result
        answer = [(0, 3), (3, 6), (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]
        assert num.allclose(answer, result) 


    def test_alpha_3(self):
        #print "test_alpha" 
        a = [0.0, 0.0]
        b = [1.0, 0.0]
        c = [2.0, 0.0]
        d = [2.0, 2.0]
        e = [1.0, 2.0]
        f = [0.0, 2.0]
        alf = 1.5

        alpha = Alpha_Shape([a,b,c,d,e,f],alf)
        result = alpha.get_boundary()
        #print "result",result
        answer = [(5, 0), (0, 1), (4, 5), (2, 3), (3, 4), (1, 2)]
        assert num.allclose(answer, result) 


    def test_boundary_1(self):
        a = [0.0, 0.0]
        b = [2.0, 0.0]
        c = [4.0, 0.0]
        cd = [3.0,2.0]
        d = [4.0, 4.0]
        e = [2.0, 4.0]
        f = [0.0, 4.0]

        alpha = Alpha_Shape([a,b,c,cd,d,e,f])
        alpha.set_boundary_type(0,1,0,0)
        result = alpha.get_boundary()
        #print "result",result
        answer = [(0, 3), (3, 6), (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)]
        assert num.allclose(answer, result) 
   


    def test_boundary_2(self):
        a = [0.0, 0.0]
        b = [2.0, 0.0]
        c = [4.0, 0.0]
        cd = [3.0,2.0]
        d = [4.0, 4.0]
        e = [2.0, 4.0]
        f = [0.0, 4.0]

        alpha = Alpha_Shape([a,b,c,cd,d,e,f])
        alpha.set_boundary_type(0,0,1,0)
        result = alpha.get_boundary()
        #print "result",result
        answer = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 0)]
        assert num.allclose(answer, result) 


    def test_boundary_3(self):
        a = [452.0, -300.0]
        b = [637.0, -322.0]
        c = [844.0, -474.0]
        d = [900.0, -585.0]
        e = [726.0, -492.0]
        f = [481.0, -354.0]
        g = [744.0, -388.0]

        alpha = Alpha_Shape([a,b,c,d,e,f,g])
        alpha.set_boundary_type(0,0,0,1)
        result = alpha.get_boundary()
        answer = [(1, 0), (0, 5), (3, 2), (4, 3), (2, 6), (6, 1), (5, 4)]
        assert num.allclose(answer, result)


        
    def test_not_enough_points(self):
        #print "test_not_enought_points"
        a = [0.0, 0.0]
        b = [1.0, 0.0]

        try:
            alpha = Alpha_Shape([a,b])
        except PointError:
            pass
        else:
            self.assertTrue(0==1,
                        'point list with 2 points did not raise an error!')

## FIXME?  Running from the Comand line isn't in vogue these days
#  The test was breaking when test_all at the inundation level was running
# was running it. issue - not running the test in this directory
    def FIXtest_alpha_stand_alone(self):
        #print "test_alpha_stand_alone"
        import os
        import tempfile

        fileName = tempfile.mktemp(".csv")
        file = open(fileName,"w")
        file.write("x,y\n\
0.0, 0.0\n\
1.0, 0.0\n\
2.0, 0.0\n\
2.0, 2.0\n\
1.0, 2.0\n\
0.0, 2.0\n")
        file.close()

        output_file_name = tempfile.mktemp(".bnd")
        #(ho, output_file_name) = tempfile.mkstemp(".bnd")
        command = sys.executable  + 'alpha_shape'+ os.sep +'alpha_shape.py ' \
                  + fileName  + ' ' + output_file_name
        #print "command", command
        os.system(command)
        os.remove(fileName)
        
        file = open(output_file_name,"r")
        lFile = file.read().split('\n')
        file.close()
        os.remove(output_file_name)
        self.assertTrue(lFile[1] == "5,0" and
                        lFile[2] == "0,1" and 
                        lFile[3] == "4,5" and 
                        lFile[4] == "2,3" and 
                        lFile[5] == "3,4" and 
                        lFile[6] == "1,2" and 
                        lFile[7] == ""
                        ,
                        'boundary file is wrong')


    def test_expand_pinch(self):
        a = [1.0, 1.0]
        b = [1.0, 5.0]
        c = [4.0, 3.0]
        d = [8.0, 5.0]
        e = [8.0, 1.0]

        alpha = Alpha_Shape([a,b,c,d,e])
        alpha.set_boundary_type(raw_boundary=True,
                          remove_holes=False,
                          smooth_indents=False,
                          expand_pinch=True)
        result = alpha.get_boundary()
        #print "result",result
        answer = [(1, 0), (4, 3), (0, 4), (3, 1)]
        assert num.allclose(answer, result)
  
    def test_sharp_indents(self):
        a = [3.0, 1.0]
        b = [5.0, 3.0]
        c = [4.0, 4.0]
        d = [3.0, 5.0]
        e = [1.0, 3.0]

        alpha = Alpha_Shape([a,b,c,d,e])
        alpha.set_boundary_type(raw_boundary=True,
                          remove_holes=False,
                          smooth_indents=True,
                          expand_pinch=False)
        result = alpha.get_boundary()
        #print "result",result
        answer = [(3, 4), (2, 3), (0, 1), (1, 2), (4, 0)]
        assert num.allclose(answer, result)

        
    def test_small_islands(self):
        """
        I couldn't find a small data set that could test this feature...
        """
        alpha = Alpha_Shape([(191.0,92.0), \
(166.0,79.0), \
(142.0,63.0), \
(124.0,44.0), \
(134.0,19.0), \
(154.0,-9.0), \
(182.0,-21.0), \
(200.0,-3.0), \
(227.0,13.0), \
(237.0,59.0), \
(208.0,86.0), \
(243.0,34.0), \
(209.0,35.0), \
(171.0,15.0), \
(159.0,41.0), \
(189.0,56.0), \
(235.0,-6.0), \
(221.0,-17.0), \
(204.0,-30.0), \
(184.0,-41.0), \
(160.0,-37.0), \
(144.0,-24.0), \
(128.0,-8.0), \
(113.0,9.0), \
(102.0,29.0), \
(95.0,46.0), \
(112.0,61.0), \
(127.0,72.0), \
(145.0,85.0), \
(164.0,100.0), \
(187.0,112.0), \
(210.0,107.0), \
(227.0,94.0), \
(244.0,74.0), \
(261.0,49.0), \
(266.0,20.0), \
(256.0,1.0), \
(244.0,-28.0), \
(228.0,-40.0), \
(208.0,-50.0), \
(192.0,-56.0), \
(169.0,-64.0), \
(155.0,-58.0), \
(141.0,-46.0), \
(129.0,-35.0), \
(113.0,-17.0), \
(103.0,-3.0), \
(90.0,9.0), \
(78.0,26.0), \
(70.0,46.0), \
(86.0,62.0), \
(101.0,79.0), \
(118.0,89.0), \
(135.0,102.0), \
(152.0,115.0), \
(177.0,125.0), \
(192.0,130.0), \
(209.0,125.0), \
(227.0,116.0), \
(234.0,104.0), \
(249.0,92.0), \
(258.0,74.0), \
(264.0,60.0), \
(280.0,42.0), \
(279.0,26.0), \
(274.0,3.0), \
(270.0,-19.0), \
(259.0,-34.0), \
(244.0,-45.0), \
(229.0,-60.0), \
(186.0,4.0), \
(203.0,20.0), \
(215.0,75.0), \
(228.0,49.0), \
(151.0,21.0), \
(163.0,60.0), \
(178.0,70.0), \
(196.0,45.0), \
(261.0,37.0), \
(113.0,36.0), \
(241.0,12.0), \
(196.0,75.0), \
(210.0,58.0), \
(163.248648097,111.633266713), \
(150.844951796,94.8282588211), \
(130.438870783,84.0250394618), \
(112.033385949,72.8217008669), \
(98.8294511765,61.2182430364), \
(184.855086816,-50.4150236771), \
(171.251032808,-52.0155006192), \
(156.446621093,-49.2146659705), \
(148.044117147,-42.8127582019), \
(264.878933922,-2.80083464873), \
(260.077503096,-12.4036963015), \
(255.27607227,-28.0083464873), \
(163.248648097,-22.0065579543), \
(128.03815537,6.00178853298), \
(221.265937249,0.0), \
(202.060213944,-13.6040540081), \
(166.449601981,-2.80083464873), \
(149.244474854,6.80202700405), \
(182.454371403,29.2087041938), \
(147.243878676,30.4090619004), \
(126.437678428,28.0083464873), \
(133.639824668,50.0149044415), \
(176.852702105,43.612996673), \
(123.636843779,-24.4072733675), \
(73.6219393379,35.2104927268), \
(212.463314068,-39.2116850822), \
(254.875953034,66.4197930983), \
(231.669037373,71.6213431603), \
(85.6255164039,23.6070348964), \
(98.4293319409,16.4048886568), \
(223.266533427,-52.0155006192), \
(208.062002477,-57.2170506811), \
(243.672614439,86.425754875), \
(214.06379101,113.633862891), \
(205.261167828,115.234339833), \
(194.057829233,120.836009131), \
(85.2253971684,50.8151429126), \
(84.0250394618,34.0101350202), \
(194.457948469,-47.6141890283), \
(232.86939508,-20.4060810121), \
(252.875356856,23.2069156609), \
(232.86939508,27.2081080162), \
(214.463910245,13.6040540081), \
(182.054252167,18.4054848345), \
(163.248648097,29.6088234294), \
(218.865221836,42.8127582019), \
(176.45258287,97.2289742343), \
(154.04590568,70.8211046892), \
(200.059617766,110.833028242)])
        alpha.set_boundary_type(raw_boundary=False ,
                          remove_holes=True,
                          smooth_indents=True,
                          expand_pinch=True)
        result = alpha.get_boundary()
        #print "result",result
        answer = [(45, 106), (44, 106), (46, 47), (44, 43), (48, 111), \
                  (47, 111), (43, 42), (42, 41), (41, 88), (40, 88), \
                  (107, 48), (49, 107), (49, 50), (50, 51), (52, 53), \
                  (51, 52), (53, 54), (55, 83), (54, 83), (40, 114), \
                  (67, 68), (68, 69), (67, 66), (65, 66), (64, 65), \
                  (69, 114), (57, 56), (56, 55), (58, 57), (64, 63), \
                  (61, 62), (62, 63), (59, 60), (58, 59), (60, 61), \
                  (46, 45)]
        assert num.allclose(answer, result)    

#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(TestCase,'test')
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite)
