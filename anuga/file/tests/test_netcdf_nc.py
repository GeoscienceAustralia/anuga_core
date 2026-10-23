"""Tests for anuga/file/netcdf.py — Write_nc, write_elevation_nc, nc_lon_lat_header."""

import os
import tempfile
import unittest
import numpy as num

from anuga.file.netcdf import NetCDFFile, Write_nc, write_elevation_nc, nc_lon_lat_header
from anuga.config import netcdf_mode_r, netcdf_mode_w


class Test_NetCDFFile(unittest.TestCase):

    def test_netcdffile_write_mode(self):
        """NetCDFFile in 'w' mode creates a writable Dataset."""
        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
            fname = f.name
        try:
            fid = NetCDFFile(fname, netcdf_mode_w)
            fid.close()
            self.assertTrue(os.path.exists(fname))
        finally:
            os.remove(fname)

    def test_netcdffile_wl_mode(self):
        """NetCDFFile with legacy 'wl' mode creates a NETCDF3_64BIT file."""
        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
            fname = f.name
        try:
            fid = NetCDFFile(fname, 'wl')
            fid.close()
            self.assertTrue(os.path.exists(fname))
        finally:
            os.remove(fname)

    def test_netcdffile_read_mode(self):
        """NetCDFFile in 'r' mode can read a previously written file."""
        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
            fname = f.name
        try:
            fid = NetCDFFile(fname, netcdf_mode_w)
            fid.close()
            fid = NetCDFFile(fname, netcdf_mode_r)
            fid.close()
        finally:
            os.remove(fname)


class Test_nc_lon_lat_header(unittest.TestCase):

    def test_writes_lon_lat_dimensions(self):
        """nc_lon_lat_header creates LON/LAT dimensions and variables."""
        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
            fname = f.name
        try:
            lon = num.array([150.0, 151.0, 152.0])
            lat = num.array([-34.0, -33.0])
            fid = NetCDFFile(fname, netcdf_mode_w)
            nc_lon_lat_header(fid, lon, lat)
            self.assertIn('LON', fid.variables)
            self.assertIn('LAT', fid.variables)
            self.assertEqual(len(fid.variables['LON']), 3)
            self.assertEqual(len(fid.variables['LAT']), 2)
            num.testing.assert_allclose(fid.variables['LON'][:], lon)
            num.testing.assert_allclose(fid.variables['LAT'][:], lat)
            fid.close()
        finally:
            os.remove(fname)


class Test_write_elevation_nc(unittest.TestCase):

    def test_writes_elevation_variable(self):
        """write_elevation_nc creates an ELEVATION variable with correct shape."""
        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
            fname = f.name
        try:
            lon = num.array([150.0, 151.0, 152.0])
            lat = num.array([-34.0, -33.0])
            depth_vector = num.array([10.0, 20.0, 30.0, 15.0, 25.0, 35.0])
            write_elevation_nc(fname, lon, lat, depth_vector)

            fid = NetCDFFile(fname, netcdf_mode_r)
            self.assertIn('ELEVATION', fid.variables)
            elev = fid.variables['ELEVATION'][:]
            self.assertEqual(elev.shape, (len(lat), len(lon)))
            fid.close()
        finally:
            os.remove(fname)


class Test_Write_nc(unittest.TestCase):

    def _make_write_nc(self, fname, quantity_name='HA'):
        lon = num.array([150.0, 151.0])
        lat = num.array([-34.0, -33.0])
        return Write_nc(quantity_name, fname, time_step_count=3,
                        time_step=0.5, lon=lon, lat=lat)

    def test_init_HA(self):
        """Write_nc initialises for HA quantity without error."""
        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
            fname = f.name
        try:
            wnc = self._make_write_nc(fname, 'HA')
            self.assertEqual(wnc.quantity_name, 'HA')
            self.assertAlmostEqual(wnc.quantity_multiplier, 100.0)
            wnc.close()
        finally:
            os.remove(fname)

    def test_init_UA(self):
        """Write_nc initialises for UA quantity."""
        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
            fname = f.name
        try:
            wnc = self._make_write_nc(fname, 'UA')
            self.assertAlmostEqual(wnc.quantity_multiplier, 100.0)
            wnc.close()
        finally:
            os.remove(fname)

    def test_init_VA(self):
        """Write_nc initialises for VA quantity (negative multiplier)."""
        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
            fname = f.name
        try:
            wnc = self._make_write_nc(fname, 'VA')
            self.assertAlmostEqual(wnc.quantity_multiplier, -100.0)
            wnc.close()
        finally:
            os.remove(fname)

    def test_store_timestep(self):
        """store_timestep writes a 2-D slice correctly scaled."""
        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
            fname = f.name
        try:
            wnc = self._make_write_nc(fname, 'HA')
            # lat=2, lon=2 → 2×2 slice
            slice0 = num.array([[1.0, 2.0], [3.0, 4.0]])
            slice1 = num.array([[0.5, 1.5], [2.5, 3.5]])
            wnc.store_timestep(slice0)
            wnc.store_timestep(slice1)
            wnc.close()

            fid = NetCDFFile(fname, netcdf_mode_r)
            stored = fid.variables['HA'][:]
            # multiplier for HA is 100
            num.testing.assert_allclose(stored[0], slice0 * 100.0, rtol=1e-5)
            num.testing.assert_allclose(stored[1], slice1 * 100.0, rtol=1e-5)
            fid.close()
        finally:
            os.remove(fname)

    def test_time_increments(self):
        """store_timestep assigns correct time values."""
        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
            fname = f.name
        try:
            wnc = self._make_write_nc(fname, 'HA')
            s = num.ones((2, 2))
            wnc.store_timestep(s)
            wnc.store_timestep(s)
            wnc.close()

            fid = NetCDFFile(fname, netcdf_mode_r)
            times = fid.variables['TIME'][:]
            self.assertAlmostEqual(times[0], 0.0)
            self.assertAlmostEqual(times[1], 0.5)
            fid.close()
        finally:
            os.remove(fname)


if __name__ == '__main__':
    unittest.main()
