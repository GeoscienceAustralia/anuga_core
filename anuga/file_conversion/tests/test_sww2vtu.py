"""Tests for sww2vtu — SWW to VTU/PVD conversion."""

import os
import tempfile
import unittest
import xml.etree.ElementTree as ET

import numpy as np

import anuga
from anuga.file_conversion.sww2vtu import sww2vtu, _b64_encode, _write_vtu, _write_pvd


def _make_sww(tmpdir, name='test.sww', n_steps=3):
    """Create a minimal SWW file using ANUGA's domain output."""
    sww_path = os.path.join(tmpdir, name)
    domain = anuga.rectangular_cross_domain(3, 3, len1=3.0, len2=3.0)
    domain.set_flow_algorithm('DE0')
    domain.set_name(os.path.splitext(name)[0])
    domain.set_datadir(tmpdir)
    domain.set_quantity('elevation', lambda x, y: 0.1 * x)
    domain.set_quantity('stage', 1.0)
    domain.set_quantity('friction', 0.03)
    domain.set_boundary({'left': anuga.Reflective_boundary(domain),
                         'right': anuga.Reflective_boundary(domain),
                         'top': anuga.Reflective_boundary(domain),
                         'bottom': anuga.Reflective_boundary(domain)})
    for _ in domain.evolve(yieldstep=0.5, finaltime=(n_steps - 1) * 0.5):
        pass
    return sww_path


class Test_b64_encode(unittest.TestCase):

    def test_round_trip(self):
        import base64
        import struct
        arr = np.array([1.0, 2.0, 3.0], dtype=np.float32)
        encoded = _b64_encode(arr)
        raw = base64.b64decode(encoded)
        nbytes = struct.unpack('<I', raw[:4])[0]
        self.assertEqual(nbytes, arr.nbytes)
        decoded = np.frombuffer(raw[4:], dtype=np.float32)
        np.testing.assert_array_equal(decoded, arr)


class Test_write_vtu(unittest.TestCase):

    def test_valid_xml(self):
        tmpdir = tempfile.mkdtemp()
        path = os.path.join(tmpdir, 'out.vtu')
        x = np.array([0.0, 1.0, 0.5])
        y = np.array([0.0, 0.0, 1.0])
        z = np.zeros(3)
        tris = np.array([[0, 1, 2]], dtype=np.int32)
        point_data = {'stage': np.array([1.0, 1.1, 1.05], dtype=np.float32)}
        _write_vtu(path, x, y, z, tris, point_data)
        self.assertTrue(os.path.isfile(path))
        tree = ET.parse(path)
        root = tree.getroot()
        self.assertEqual(root.tag, 'VTKFile')

    def test_multiple_quantities(self):
        tmpdir = tempfile.mkdtemp()
        path = os.path.join(tmpdir, 'out.vtu')
        x = np.array([0.0, 1.0, 0.5])
        y = np.array([0.0, 0.0, 1.0])
        z = np.zeros(3)
        tris = np.array([[0, 1, 2]], dtype=np.int32)
        point_data = {
            'stage': np.ones(3, dtype=np.float32),
            'elevation': np.zeros(3, dtype=np.float32),
        }
        _write_vtu(path, x, y, z, tris, point_data)
        content = open(path).read()
        self.assertIn('stage', content)
        self.assertIn('elevation', content)


class Test_write_pvd(unittest.TestCase):

    def test_valid_xml(self):
        tmpdir = tempfile.mkdtemp()
        path = os.path.join(tmpdir, 'out.pvd')
        _write_pvd(path, ['out_0000.vtu', 'out_0001.vtu'], [0.0, 1.0])
        self.assertTrue(os.path.isfile(path))
        tree = ET.parse(path)
        root = tree.getroot()
        self.assertEqual(root.tag, 'VTKFile')
        datasets = root.findall('.//DataSet')
        self.assertEqual(len(datasets), 2)


class Test_sww2vtu(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.sww_path = _make_sww(self.tmpdir, n_steps=3)

    def test_produces_pvd(self):
        outdir = os.path.join(self.tmpdir, 'vtu_out')
        pvd = sww2vtu(self.sww_path, output_dir=outdir)
        self.assertTrue(os.path.isfile(pvd))
        self.assertTrue(pvd.endswith('.pvd'))

    def test_produces_vtu_files(self):
        outdir = os.path.join(self.tmpdir, 'vtu_out2')
        sww2vtu(self.sww_path, output_dir=outdir)
        vtu_files = [f for f in os.listdir(outdir) if f.endswith('.vtu')]
        self.assertGreater(len(vtu_files), 0)

    def test_vtu_is_valid_xml(self):
        outdir = os.path.join(self.tmpdir, 'vtu_out3')
        sww2vtu(self.sww_path, output_dir=outdir)
        vtu_files = sorted(f for f in os.listdir(outdir) if f.endswith('.vtu'))
        tree = ET.parse(os.path.join(outdir, vtu_files[0]))
        root = tree.getroot()
        self.assertEqual(root.tag, 'VTKFile')

    def test_verbose(self):
        outdir = os.path.join(self.tmpdir, 'vtu_out4')
        pvd = sww2vtu(self.sww_path, output_dir=outdir, verbose=True)
        self.assertTrue(os.path.isfile(pvd))

    def test_z_scale(self):
        outdir = os.path.join(self.tmpdir, 'vtu_out5')
        pvd = sww2vtu(self.sww_path, output_dir=outdir, z_scale=10.0)
        self.assertTrue(os.path.isfile(pvd))

    def test_absolute_coords(self):
        outdir = os.path.join(self.tmpdir, 'vtu_out6')
        pvd = sww2vtu(self.sww_path, output_dir=outdir, absolute_coords=True)
        self.assertTrue(os.path.isfile(pvd))

    def test_default_output_dir(self):
        """When output_dir is None the VTU files land next to the SWW."""
        pvd = sww2vtu(self.sww_path)
        self.assertTrue(os.path.isfile(pvd))
        self.assertEqual(os.path.dirname(pvd),
                         os.path.dirname(os.path.abspath(self.sww_path)))


if __name__ == '__main__':
    unittest.main()
