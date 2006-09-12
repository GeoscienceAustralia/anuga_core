import unittest
from Numeric import allclose
from smf import slide_tsunami, slump_tsunami, Double_gaussian

class Test_smf(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
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

        dg = Double_gaussian(a3D=a3D, wavelength=wavelength, width=width, \
                             x0=x0, y0=y0, alpha=alpha, dx = dx, \
                             kappa=kappa, kappad = kappad, zsmall = zsmall)

        assert allclose(dg.a3D, a3D)
        assert allclose(dg.wavelength, wavelength)
        assert allclose(dg.width, width)
        assert allclose(dg.x0, x0)
        assert allclose(dg.y0, y0)
        assert allclose(dg.alpha, alpha)
        assert allclose(dg.kappa, kappa)
        assert allclose(dg.kappad, kappad)
        assert allclose(dg.dx, dx)


    def test_slide_tsunami(self):

        len = 600.0
        dep = 150.0
        th = 9.0
        thk = 15.0
        wid = 340.0

        slide = slide_tsunami(length=len, depth=dep, slope=th, \
                              width = wid, thickness=thk)

        assert allclose(slide.a3D, 0.07775819)
        assert allclose(slide.wavelength, 2938.66695708)
        assert allclose(slide.width, 340.0)
        assert allclose(slide.x0, 0.0)
        assert allclose(slide.y0, 0.0)
        assert allclose(slide.alpha, 0.0)
        assert allclose(slide.kappa, 3.0)
        assert allclose(slide.kappad, 0.8)


    def test_slump_tsunami(self):

        len = 4500.0
        thk = 760.0
        wid = 4500.0
        dep = 1200.0
        rad = 3330
        dp = 0.23
        th = 12
        alpha = 0.0
        x0 = 0
        y0 = 0

        slump = slump_tsunami(length=len, depth=dep, slope=th, thickness=thk,\
                  radius=rad, dphi=dp, x0=x0, y0=y0, alpha=alpha)

        assert allclose(slump.a3D, 9.82538623)
        assert allclose(slump.wavelength, 3660.37606554)
        assert allclose(slump.width, 4500.0)
        assert allclose(slump.x0, 0.0)
        assert allclose(slump.y0, 0.0)
        assert allclose(slump.alpha, 0.0)
        assert allclose(slump.kappa, 3.0)
        assert allclose(slump.kappad, 1.0)


#-------------------------------------------------------------
if __name__ == "__main__":
    suite = unittest.makeSuite(Test_smf,'test_Double_gaussian')
    runner = unittest.TextTestRunner()
    runner.run(suite)

