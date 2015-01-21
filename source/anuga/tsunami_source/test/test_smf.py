import unittest
import os
import numpy as num
from anuga.tsunami_source.smf import slide_tsunami, slump_tsunami, Double_gaussian
from anuga import Domain

class Test_smf(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('test.msh')
        except:
            pass


    def test_Double_gaussian(self):
        a3D = 5.6
        wavelength = 10.2
        width = 6.6
        x0 = 0.0
        y0 = 0.0
        alpha = 0.0
        kappa = 3.0
        kappad = 0.8
        dx = 0.5
        zsmall = 0.01
        scale = 1.0

        dg = Double_gaussian(a3D, wavelength, width, \
                             x0, y0, alpha, \
                             kappa, kappad, zsmall, dx, scale)

        assert num.allclose(dg.a3D, a3D)
        assert num.allclose(dg.wavelength, wavelength)
        assert num.allclose(dg.width, width)
        assert num.allclose(dg.x0, x0)
        assert num.allclose(dg.y0, y0)
        assert num.allclose(dg.alpha, alpha)
        assert num.allclose(dg.kappa, kappa)
        assert num.allclose(dg.kappad, kappad)
        assert num.allclose(dg.dx, dx)

    
    def test_slide_tsunami(self):

        len = 600.0
        dep = 150.0
        th = 9.0
        thk = 15.0
        wid = 340.0
        kappa = 3.0
        kappad = 0.8
        x0 = 100000.

        slide = slide_tsunami(length=len, depth=dep, slope=th, x0=x0, \
                              width = wid, thickness=thk, kappa=kappa, kappad=kappad, \
                              verbose=False)

        assert num.allclose(slide.a3D, 0.07775819)
        assert num.allclose(slide.wavelength, 2938.66695708)
        assert num.allclose(slide.width, 340.0)
        assert num.allclose(slide.y0, 0.0)
        assert num.allclose(slide.alpha, 0.0)


    def test_slump_tsunami(self):

        length = 4500.0
        thk = 760.0
        wid = 4500.0
        dep = 1200.0
        rad = 3330
        dp = 0.23
        th = 12
        alpha = 0.0
        x0 = 0
        y0 = 0
        gamma = 1.85

        slump = slump_tsunami(length, dep, th, wid, thk, rad, dp, x0, y0, alpha, gamma, scale=1.0)

        assert num.allclose(slump.a3D, 9.82538623)
        assert num.allclose(slump.wavelength, 3660.37606554)
        assert num.allclose(slump.width, 4500.0)
        assert num.allclose(slump.x0, 0.0)
        assert num.allclose(slump.y0, 0.0)
        assert num.allclose(slump.alpha, 0.0)
        assert num.allclose(slump.kappa, 3.0)
        assert num.allclose(slump.kappad, 1.0)

    def test_slide_tsunami_domain(self):

        length = 600.0
        dep = 150.0
        th = 9.0
        thk = 15.0
        wid = 340.0
        kappa = 3.0
        kappad = 0.8
        x0 = 100000.
        y0 = x0
        
        from anuga.pmesh.mesh_interface import create_mesh_from_regions
        polygon = [[0,0],[200000,0],[200000,200000],[0,200000]]
        create_mesh_from_regions(polygon,
                                 {'e0': [0], 'e1': [1], 'e2': [2], 'e3': [3]},
                                 maximum_triangle_area=5000000000,
                                 filename='test.msh',
                                 verbose = False)

        domain = Domain('test.msh', use_cache = True, verbose = False)

        slide = slide_tsunami(length, dep, th, x0, y0, \
                              wid, thk, kappa, kappad, \
                              domain=domain,verbose=False)

        domain.set_quantity('stage', slide)
        stage = domain.get_quantity('stage')
        w = stage.get_values()

##        check = [[-0.0 -0.0 -0.0],
##                 [-.189709745 -517.877716 -0.0],
##                 [-0.0 -0.0 -2.7695931e-08],
##                 [-0.0 -2.7695931e-08 -1.897097e-01]
##                 [-0.0 -517.877716 -0.0],
##                 [-0.0 -0.0 -0.0],
##                 [-0.0 -0.0 -0.0],
##                 [-0.0 -0.0 -0.0]]

        assert num.allclose(num.min(w), -517.877771593)
        assert num.allclose(num.max(w), 0.0)
        assert num.allclose(slide.a3D, 518.38797486)


#-------------------------------------------------------------

if __name__ == "__main__":
    suite = unittest.makeSuite(Test_smf,'test')
    runner = unittest.TextTestRunner()
    runner.run(suite)

