"""Tests for anuga/file/sts.py — Write_sts and create_sts_boundary."""

import os
import tempfile
import unittest
import numpy as num

from anuga.file.netcdf import NetCDFFile
from anuga.file.sts import Write_sts, create_sts_boundary
from anuga.config import netcdf_mode_w, netcdf_mode_r


def _make_sts_file(fname, n_points=4, n_times=3):
    """Write a minimal STS file and return it for read-back tests."""
    times = num.linspace(0.0, 2.0, n_times)
    points_utm = num.array([[100.0, 200.0],
                             [110.0, 200.0],
                             [110.0, 210.0],
                             [100.0, 210.0]], dtype=float)[:n_points]
    elevation = num.zeros(n_points)
    zone = 55

    fid = NetCDFFile(fname, netcdf_mode_w)
    ws = Write_sts()
    ws.store_header(fid, times, n_points)
    ws.store_points(fid, points_utm, elevation, zone=zone)
    for i in range(n_times):
        ws.store_quantities(fid,
                            slice_index=i,
                            stage=num.ones(n_points) * float(i),
                            xmomentum=num.zeros(n_points),
                            ymomentum=num.zeros(n_points))
    fid.close()


class Test_Write_sts_header(unittest.TestCase):

    def test_store_header_with_time_array(self):
        """store_header with a time array stores relative times."""
        with tempfile.NamedTemporaryFile(suffix='.sts', delete=False) as f:
            fname = f.name
        try:
            times = num.array([10.0, 11.0, 12.0])
            fid = NetCDFFile(fname, netcdf_mode_w)
            ws = Write_sts()
            ws.store_header(fid, times, number_of_points=3)
            # starttime should be first time value
            self.assertAlmostEqual(fid.starttime, 10.0)
            # relative times stored
            stored_times = fid.variables['time'][:]
            num.testing.assert_allclose(stored_times, [0.0, 1.0, 2.0])
            fid.close()
        finally:
            os.remove(fname)

    def test_store_header_empty_times(self):
        """store_header with an empty time list sets starttime=0."""
        with tempfile.NamedTemporaryFile(suffix='.sts', delete=False) as f:
            fname = f.name
        try:
            fid = NetCDFFile(fname, netcdf_mode_w)
            ws = Write_sts()
            ws.store_header(fid, [], number_of_points=2)
            self.assertEqual(fid.starttime, 0)
            fid.close()
        finally:
            os.remove(fname)

    def test_store_header_scalar_starttime(self):
        """store_header with a scalar stores it as starttime, zero timesteps."""
        with tempfile.NamedTemporaryFile(suffix='.sts', delete=False) as f:
            fname = f.name
        try:
            fid = NetCDFFile(fname, netcdf_mode_w)
            ws = Write_sts()
            ws.store_header(fid, 5.0, number_of_points=2)
            self.assertAlmostEqual(fid.starttime, 5.0)
            # No time steps stored
            self.assertEqual(len(fid.variables['time']), 0)
            fid.close()
        finally:
            os.remove(fname)


class Test_Write_sts_store_points(unittest.TestCase):

    def test_store_points_with_zone(self):
        """store_points writes x, y, elevation and geo_reference."""
        with tempfile.NamedTemporaryFile(suffix='.sts', delete=False) as f:
            fname = f.name
        try:
            n = 3
            times = num.array([0.0, 1.0])
            points = num.array([[300000.0, 6000000.0],
                                 [300100.0, 6000000.0],
                                 [300000.0, 6000100.0]])
            elevation = num.array([0.0, -1.0, -2.0])

            fid = NetCDFFile(fname, netcdf_mode_w)
            ws = Write_sts()
            ws.store_header(fid, times, n)
            ws.store_points(fid, points, elevation, zone=55)

            # Check elevation range attributes updated
            self.assertAlmostEqual(fid.variables['elevation_range'][0],
                                   min(elevation))
            self.assertAlmostEqual(fid.variables['elevation_range'][1],
                                   max(elevation))
            fid.close()
        finally:
            os.remove(fname)

    def test_store_points_with_existing_georef(self):
        """store_points with points_georeference uses it directly."""
        from anuga.coordinate_transforms.geo_reference import Geo_reference
        with tempfile.NamedTemporaryFile(suffix='.sts', delete=False) as f:
            fname = f.name
        try:
            n = 2
            times = num.array([0.0, 1.0])
            geo_ref = Geo_reference(zone=55, xllcorner=100.0, yllcorner=200.0)
            points = num.array([[0.0, 0.0], [10.0, 0.0]])
            elevation = num.array([0.0, 0.0])

            fid = NetCDFFile(fname, netcdf_mode_w)
            ws = Write_sts()
            ws.store_header(fid, times, n)
            ws.store_points(fid, points, elevation,
                            points_georeference=geo_ref)
            fid.close()
        finally:
            os.remove(fname)


class Test_Write_sts_store_quantities(unittest.TestCase):

    def test_store_quantities_updates_range(self):
        """store_quantities updates stage_range correctly."""
        with tempfile.NamedTemporaryFile(suffix='.sts', delete=False) as f:
            fname = f.name
        try:
            n = 3
            times = num.array([0.0, 1.0, 2.0])
            fid = NetCDFFile(fname, netcdf_mode_w)
            ws = Write_sts()
            ws.store_header(fid, times, n)

            ws.store_quantities(fid, slice_index=0,
                                stage=num.array([1.0, 2.0, 3.0]),
                                xmomentum=num.zeros(n),
                                ymomentum=num.zeros(n))
            ws.store_quantities(fid, slice_index=1,
                                stage=num.array([0.5, 1.5, 4.0]),
                                xmomentum=num.zeros(n),
                                ymomentum=num.zeros(n))

            stage_range = fid.variables['stage_range'][:]
            self.assertAlmostEqual(stage_range[0], 0.5)   # global min
            self.assertAlmostEqual(stage_range[1], 4.0)   # global max
            fid.close()
        finally:
            os.remove(fname)

    def test_store_quantities_with_time(self):
        """store_quantities with time= argument advances the time dimension."""
        with tempfile.NamedTemporaryFile(suffix='.sts', delete=False) as f:
            fname = f.name
        try:
            n = 2
            fid = NetCDFFile(fname, netcdf_mode_w)
            ws = Write_sts()
            ws.store_header(fid, 0.0, n)   # scalar starttime → 0 timesteps initially
            ws.store_quantities(fid, time=1.0,
                                stage=num.array([2.0, 3.0]),
                                xmomentum=num.zeros(n),
                                ymomentum=num.zeros(n))
            self.assertEqual(len(fid.variables['time']), 1)
            self.assertAlmostEqual(fid.variables['time'][0], 1.0)
            fid.close()
        finally:
            os.remove(fname)

    def test_store_quantities_missing_key_raises(self):
        """store_quantities raises if a required quantity is missing."""
        with tempfile.NamedTemporaryFile(suffix='.sts', delete=False) as f:
            fname = f.name
        try:
            n = 2
            fid = NetCDFFile(fname, netcdf_mode_w)
            ws = Write_sts()
            ws.store_header(fid, num.array([0.0, 1.0]), n)
            with self.assertRaises(Exception):
                ws.store_quantities(fid, slice_index=0,
                                    stage=num.zeros(n))  # missing xmomentum, ymomentum
            fid.close()
        finally:
            os.remove(fname)


class Test_create_sts_boundary(unittest.TestCase):

    def test_returns_list_of_points(self):
        """create_sts_boundary returns a list of [x, y] pairs from the .sts file."""
        with tempfile.NamedTemporaryFile(suffix='.sts', delete=False) as f:
            fname = f.name
        try:
            _make_sts_file(fname, n_points=4, n_times=3)
            boundary = create_sts_boundary(fname)
            self.assertIsInstance(boundary, list)
            self.assertEqual(len(boundary), 4)
            # Each point should be a two-element list
            for pt in boundary:
                self.assertEqual(len(pt), 2)
        finally:
            os.remove(fname)

    def test_accepts_name_without_extension(self):
        """create_sts_boundary appends .sts if not present."""
        with tempfile.NamedTemporaryFile(suffix='.sts', delete=False) as f:
            fname = f.name
        try:
            _make_sts_file(fname, n_points=3, n_times=2)
            # Strip the .sts extension
            name_no_ext = fname[:-4]
            boundary = create_sts_boundary(name_no_ext)
            self.assertEqual(len(boundary), 3)
        finally:
            os.remove(fname)

    def test_missing_file_raises_oserror(self):
        """create_sts_boundary raises OSError for a non-existent file."""
        with self.assertRaises(OSError):
            create_sts_boundary('/nonexistent/path/missing.sts')


if __name__ == '__main__':
    unittest.main()
