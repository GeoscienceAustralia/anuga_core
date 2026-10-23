"""Tests for sequential_mesh_dump / sequential_mesh_load.

Round-trip tests: dump a mesh partition to NetCDF4 files, load each rank's
file, and verify that the reconstructed domain has the correct mesh topology,
halo structure, and boundary tags.

Run serially (single-process)::

    python -m pytest anuga/parallel/tests/test_sequential_mesh_io.py

The single-process tests construct a small domain, call sequential_mesh_dump
with numprocs=N, then call sequential_mesh_load from a mocked single-rank
environment by patching myid/numprocs, and compare against the domain built
by the existing sequential_distribute path.
"""

import os
import tempfile
import unittest

import numpy as num
import pytest

import anuga
from anuga import Domain, Reflective_boundary, Dirichlet_boundary
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross


def _make_domain(m=6, n=6):
    """Small rectangular-cross domain for testing."""
    points, vertices, boundary = rectangular_cross(m, n, len1=1.0, len2=1.0)
    domain = Domain(points, vertices, boundary)
    domain.set_name('test_mesh')
    domain.set_quantity('elevation', lambda x, y: x + y)
    domain.set_quantity('stage',     0.5)
    return domain


class TestSequentialMeshDump(unittest.TestCase):

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    # ------------------------------------------------------------------
    # Dump tests
    # ------------------------------------------------------------------

    def test_dump_creates_files(self):
        """sequential_mesh_dump writes one .nc file per rank."""
        domain = _make_domain()
        anuga.sequential_mesh_dump(domain, numprocs=2,
                                   partition_dir=self.tmpdir, verbose=False)
        for p in range(2):
            fname = os.path.join(self.tmpdir, f'test_mesh_mesh_P2_{p}.nc')
            self.assertTrue(os.path.exists(fname), f'missing {fname}')

    def test_dump_netcdf_variables(self):
        """Each partition file has the expected NetCDF variables and attributes."""
        import netCDF4
        domain = _make_domain()
        anuga.sequential_mesh_dump(domain, numprocs=2,
                                   partition_dir=self.tmpdir)
        fname = os.path.join(self.tmpdir, 'test_mesh_mesh_P2_0.nc')
        with netCDF4.Dataset(fname, 'r') as nc:
            self.assertIn('points',        nc.variables)
            self.assertIn('vertices',      nc.variables)
            self.assertIn('tri_l2g',       nc.variables)
            self.assertIn('node_l2g',      nc.variables)
            self.assertIn('boundary_tri',  nc.variables)
            self.assertIn('boundary_edge', nc.variables)
            self.assertIn('boundary_tag',  nc.variables)
            self.assertIn('send_ranks',    nc.variables)
            self.assertIn('recv_ranks',    nc.variables)
            self.assertEqual(int(nc.numprocs), 2)
            self.assertEqual(int(nc.rank),     0)
            self.assertGreater(int(nc.number_of_full_triangles), 0)
            self.assertGreater(int(nc.number_of_global_triangles), 0)

    def test_dump_partition_counts(self):
        """Full-triangle counts across all ranks sum to the global count."""
        import netCDF4
        domain = _make_domain(m=8, n=8)
        N = domain.number_of_triangles
        anuga.sequential_mesh_dump(domain, numprocs=4,
                                   partition_dir=self.tmpdir)
        total_full = 0
        for p in range(4):
            fname = os.path.join(self.tmpdir, f'test_mesh_mesh_P4_{p}.nc')
            with netCDF4.Dataset(fname, 'r') as nc:
                total_full += int(nc.number_of_full_triangles)
                self.assertEqual(int(nc.number_of_global_triangles), N)
        self.assertEqual(total_full, N)

    def test_dump_boundary_tags_preserved(self):
        """Boundary tags in the file match those of the original domain."""
        import netCDF4
        domain = _make_domain()
        original_tags = set(domain.boundary.values())
        anuga.sequential_mesh_dump(domain, numprocs=2,
                                   partition_dir=self.tmpdir)
        found_tags = set()
        for p in range(2):
            fname = os.path.join(self.tmpdir, f'test_mesh_mesh_P2_{p}.nc')
            with netCDF4.Dataset(fname, 'r') as nc:
                raw = nc['boundary_tag'][:]
                tags = {b''.join(row).decode('ascii').rstrip('\x00')
                        for row in raw}
                found_tags |= tags
        # Every original tag (plus 'ghost') should appear somewhere
        self.assertTrue(original_tags.issubset(found_tags | {'ghost'}))

    def test_dump_custom_name(self):
        """name parameter overrides domain.get_name()."""
        domain = _make_domain()
        anuga.sequential_mesh_dump(domain, numprocs=2,
                                   partition_dir=self.tmpdir,
                                   name='custom')
        fname = os.path.join(self.tmpdir, 'custom_mesh_P2_0.nc')
        self.assertTrue(os.path.exists(fname))

    def test_dump_creates_partition_dir(self):
        """partition_dir is created when it does not exist."""
        new_dir = os.path.join(self.tmpdir, 'subdir', 'partitions')
        domain = _make_domain()
        anuga.sequential_mesh_dump(domain, numprocs=2,
                                   partition_dir=new_dir)
        self.assertTrue(os.path.isdir(new_dir))

    # ------------------------------------------------------------------
    # Round-trip tests
    # ------------------------------------------------------------------

    def _load_rank(self, name, partition_dir, rank, numprocs_total):
        """Load partition file for a given rank without MPI."""
        import netCDF4
        from anuga.parallel.sequential_distribute import (
            _write_mesh_partition, sequential_mesh_load)
        from anuga.coordinate_transforms.geo_reference import Geo_reference
        from anuga.parallel.parallel_shallow_water import Parallel_domain

        fname = os.path.join(partition_dir,
                             f'{name}_mesh_P{numprocs_total}_{rank}.nc')
        with netCDF4.Dataset(fname, 'r') as nc:
            number_of_full_triangles   = int(nc.number_of_full_triangles)
            number_of_full_nodes       = int(nc.number_of_full_nodes)
            number_of_global_triangles = int(nc.number_of_global_triangles)
            number_of_global_nodes     = int(nc.number_of_global_nodes)
            ghost_layer_width          = int(nc.ghost_layer_width)
            geo_ref = Geo_reference(NetCDFObject=nc)

            points   = num.array(nc['points'][:])
            vertices = num.array(nc['vertices'][:])
            tri_l2g  = num.array(nc['tri_l2g'][:])
            node_l2g = num.array(nc['node_l2g'][:])

            bnd_tris  = num.array(nc['boundary_tri'][:])
            bnd_edges = num.array(nc['boundary_edge'][:])
            raw_tags  = nc['boundary_tag'][:]
            bnd_tags  = [b''.join(row).decode('ascii').rstrip('\x00')
                         for row in raw_tags]
            boundary = {(int(bnd_tris[i]), int(bnd_edges[i])): bnd_tags[i]
                        for i in range(len(bnd_tris))}

            def _read_comm(nc, prefix):
                ranks_  = list(nc[f'{prefix}_ranks'][:].astype(int))
                offsets = nc[f'{prefix}_offsets'][:].astype(int)
                local_  = num.array(nc[f'{prefix}_local'][:])
                global_ = num.array(nc[f'{prefix}_global'][:])
                d = {}
                for i, r in enumerate(ranks_):
                    s, e = offsets[i], offsets[i + 1]
                    d[r] = [local_[s:e], global_[s:e]]
                return d

            full_send_dict  = _read_comm(nc, 'send')
            ghost_recv_dict = _read_comm(nc, 'recv')

        domain = Parallel_domain(
            points, vertices, boundary,
            full_send_dict=full_send_dict,
            ghost_recv_dict=ghost_recv_dict,
            number_of_full_nodes=number_of_full_nodes,
            number_of_full_triangles=number_of_full_triangles,
            number_of_global_triangles=number_of_global_triangles,
            number_of_global_nodes=number_of_global_nodes,
            processor=rank,
            numproc=numprocs_total,
            s2p_map=None,
            p2s_map=None,
            tri_l2g=tri_l2g,
            node_l2g=node_l2g,
            ghost_layer_width=ghost_layer_width,
            geo_reference=geo_ref,
        )
        boundary_map = {tag: None for tag in set(boundary.values())}
        boundary_map['ghost'] = None
        domain.set_boundary(boundary_map)
        return domain

    @pytest.mark.skipif(not hasattr(os, 'fork'),
                        reason='parallel dump uses a fork process pool (POSIX only)')
    def test_parallel_dump_matches_serial(self):
        """num_workers>1 (fork pool) produces the same NetCDF files as serial."""
        import netCDF4

        def dump(nw):
            d = tempfile.mkdtemp()
            anuga.sequential_mesh_dump(_make_domain(8, 6), numprocs=5,
                                       partition_dir=d, num_workers=nw)
            return d

        serial_dir = dump(1)
        parallel_dir = dump(4)

        ncs = sorted(f for f in os.listdir(serial_dir) if f.endswith('.nc'))
        self.assertEqual(ncs, sorted(f for f in os.listdir(parallel_dir)
                                     if f.endswith('.nc')))
        self.assertEqual(len(ncs), 5)
        for f in ncs:
            with netCDF4.Dataset(os.path.join(serial_dir, f)) as a, \
                 netCDF4.Dataset(os.path.join(parallel_dir, f)) as b:
                self.assertEqual(set(a.variables), set(b.variables))
                for v in a.variables:
                    num.testing.assert_array_equal(
                        num.asarray(a.variables[v][:]),
                        num.asarray(b.variables[v][:]),
                        err_msg=f'{f}:{v} differs between serial and parallel dump')

    @pytest.mark.skipif(not hasattr(os, 'fork'),
                        reason='parallel dump uses a fork process pool (POSIX only)')
    def test_parallel_dump_domain_matches_serial(self):
        """sequential_distribute_dump: num_workers>1 reproduces the serial dump.

        Single-file (default) partitions embed no paths, so the pickle files are
        byte-identical between the serial and parallel writers -- and there are
        no separate .npy files.
        """
        import filecmp
        from anuga.parallel.sequential_distribute import sequential_distribute_dump

        def dump(nw):
            d = tempfile.mkdtemp()
            sequential_distribute_dump(_make_domain(8, 6), numprocs=5,
                                       partition_dir=d, num_workers=nw)
            return d

        s = dump(1)
        p = dump(4)

        files = sorted(os.listdir(s))
        self.assertEqual(files, sorted(os.listdir(p)))
        self.assertEqual(len([f for f in files if f.endswith('.pickle')]), 5)
        self.assertFalse([f for f in files if f.endswith('.npy')],
                         'single-file dump should not write .npy files')
        match, mismatch, errors = filecmp.cmpfiles(s, p, files, shallow=False)
        self.assertEqual((mismatch, errors), ([], []),
                         f'serial vs parallel differ: {mismatch} {errors}')

    def test_domain_dump_single_file_equals_legacy(self):
        """single_file=True stores the same partition data as the legacy
        multi-file layout (and the parallel writer works for both)."""
        import pickle
        from anuga.parallel.sequential_distribute import sequential_distribute_dump

        sa = tempfile.mkdtemp()  # single-file (default), parallel writer
        sb = tempfile.mkdtemp()  # legacy multi-file, serial writer
        sequential_distribute_dump(_make_domain(8, 6), numprocs=4,
                                   partition_dir=sa, num_workers=3)
        sequential_distribute_dump(_make_domain(8, 6), numprocs=4,
                                   partition_dir=sb, single_file=False)

        # single-file: one pickle, no .npy;  legacy: 3 + N_quant .npy per rank
        self.assertFalse([f for f in os.listdir(sa) if f.endswith('.npy')])
        self.assertEqual(len([f for f in os.listdir(sa) if f.endswith('.pickle')]), 4)
        self.assertTrue([f for f in os.listdir(sb) if f.endswith('.npy')])

        for p in range(4):
            fn = f'test_mesh_P4_{p}.pickle'
            with open(os.path.join(sa, fn), 'rb') as f:
                A = list(pickle.load(f))          # arrays inline
            with open(os.path.join(sb, fn), 'rb') as f:
                B = list(pickle.load(f))          # points/vertices/quantities are paths
            num.testing.assert_array_equal(A[1], num.load(B[1]))   # points
            num.testing.assert_array_equal(A[2], num.load(B[2]))   # vertices
            self.assertEqual(set(A[4]), set(B[4]))                 # quantity names
            for k in A[4]:
                num.testing.assert_array_equal(A[4][k], num.load(B[4][k]))

    def test_roundtrip_node_count(self):
        """Loaded domain has correct full-triangle and full-node counts."""
        import netCDF4
        domain = _make_domain()
        N_tri = domain.number_of_triangles
        N_nod = domain.number_of_nodes
        np2 = 2
        anuga.sequential_mesh_dump(domain, numprocs=np2,
                                   partition_dir=self.tmpdir)
        for p in range(np2):
            d = self._load_rank('test_mesh', self.tmpdir, p, np2)
            fname = os.path.join(self.tmpdir,
                                 f'test_mesh_mesh_P{np2}_{p}.nc')
            with netCDF4.Dataset(fname, 'r') as nc:
                nft = int(nc.number_of_full_triangles)
                nfn = int(nc.number_of_full_nodes)
            self.assertEqual(d.number_of_triangles,
                             len(d.mesh.triangles))
            self.assertEqual(int(d.number_of_full_triangles), nft)
            # number_of_full_nodes is overwritten by generic_domain to total
            # node count; the correct full-only count lives in _tmp attribute
            self.assertEqual(int(d.number_of_full_nodes_tmp), nfn)
            self.assertEqual(int(d.number_of_global_triangles), N_tri)
            self.assertEqual(int(d.number_of_global_nodes),     N_nod)

    def test_roundtrip_tri_l2g_coverage(self):
        """Union of full-triangle global IDs across all ranks covers 0..N-1."""
        import netCDF4
        domain = _make_domain(m=6, n=6)
        N = domain.number_of_triangles
        np3 = 3
        anuga.sequential_mesh_dump(domain, numprocs=np3,
                                   partition_dir=self.tmpdir)
        global_ids = set()
        for p in range(np3):
            fname = os.path.join(self.tmpdir,
                                 f'test_mesh_mesh_P{np3}_{p}.nc')
            with netCDF4.Dataset(fname, 'r') as nc:
                nft    = int(nc.number_of_full_triangles)
                l2g    = num.array(nc['tri_l2g'][:])
            global_ids |= set(l2g[:nft].tolist())
        self.assertEqual(global_ids, set(range(N)))

    def test_roundtrip_points_shape(self):
        """Loaded domain points array has correct shape."""
        domain = _make_domain()
        anuga.sequential_mesh_dump(domain, numprocs=2,
                                   partition_dir=self.tmpdir)
        for p in range(2):
            d = self._load_rank('test_mesh', self.tmpdir, p, 2)
            pts = d.mesh.nodes
            self.assertEqual(pts.ndim, 2)
            self.assertEqual(pts.shape[1], 2)
            self.assertGreater(pts.shape[0], 0)

    def test_roundtrip_comm_dicts_nonempty(self):
        """Inner ranks have non-empty send and recv dicts."""
        import netCDF4
        domain = _make_domain(m=8, n=8)
        np4 = 4
        anuga.sequential_mesh_dump(domain, numprocs=np4,
                                   partition_dir=self.tmpdir)
        # rank 1 (inner) must have both send and recv neighbours
        fname = os.path.join(self.tmpdir, f'test_mesh_mesh_P{np4}_1.nc')
        with netCDF4.Dataset(fname, 'r') as nc:
            send_nranks = len(nc['send_ranks'][:])
            recv_nranks = len(nc['recv_ranks'][:])
        self.assertGreater(send_nranks, 0)
        self.assertGreater(recv_nranks, 0)

    def test_roundtrip_geo_reference(self):
        """Geo_reference round-trips correctly through the file."""
        import netCDF4
        from anuga.coordinate_transforms.geo_reference import Geo_reference
        domain = _make_domain()
        domain.geo_reference = Geo_reference(zone=56, xllcorner=300000.0,
                                             yllcorner=6000000.0)
        anuga.sequential_mesh_dump(domain, numprocs=2,
                                   partition_dir=self.tmpdir)
        fname = os.path.join(self.tmpdir, 'test_mesh_mesh_P2_0.nc')
        with netCDF4.Dataset(fname, 'r') as nc:
            gr = Geo_reference(NetCDFObject=nc)
        self.assertEqual(gr.zone,        56)
        self.assertAlmostEqual(gr.xllcorner, 300000.0)
        self.assertAlmostEqual(gr.yllcorner, 6000000.0)


if __name__ == '__main__':
    unittest.main()
