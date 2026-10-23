"""Tests for uniform mesh refinement: uniform_refine_domain,
sequential_mesh_refine, and create_parallel_mesh.

Serial tests only — no MPI required.
"""

import os
import tempfile
import unittest

import numpy as num
import netCDF4

import anuga
from anuga import Domain
from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross


def _make_domain(m=4, n=4):
    """Small rectangular-cross domain for testing."""
    points, vertices, boundary = rectangular_cross(m, n, len1=1.0, len2=1.0)
    domain = Domain(points, vertices, boundary)
    domain.set_name('test_refine')
    return domain


class TestUniformRefineDomain(unittest.TestCase):
    """Tests for uniform_refine_domain (sequential, single-process)."""

    def test_triangle_count(self):
        domain = _make_domain(3, 3)
        N0 = domain.number_of_triangles
        refined = anuga.uniform_refine_domain(domain)
        self.assertEqual(refined.number_of_triangles, 4 * N0)

    def test_node_count(self):
        domain = _make_domain(3, 3)
        M0 = domain.number_of_nodes
        N0 = domain.number_of_triangles
        refined = anuga.uniform_refine_domain(domain)
        M1 = refined.number_of_nodes
        # New nodes = midpoints of unique edges; for a closed triangulation
        # with E edges: E = (3N0 + B) / 2 where B = boundary edges.
        self.assertGreater(M1, M0)
        self.assertLess(M1, M0 + 3 * N0)  # at most 3 new per triangle

    def test_boundary_tag_count_doubles(self):
        """Each original boundary edge generates two child boundary edges."""
        domain = _make_domain(3, 3)
        n_bnd = len(domain.boundary)
        refined = anuga.uniform_refine_domain(domain)
        self.assertEqual(len(refined.boundary), 2 * n_bnd)

    def test_boundary_tags_preserved(self):
        domain = _make_domain(3, 3)
        orig_tags = set(domain.boundary.values())
        refined = anuga.uniform_refine_domain(domain)
        refined_tags = set(refined.boundary.values())
        self.assertEqual(refined_tags, orig_tags)

    def test_original_nodes_unchanged(self):
        """Original node coordinates must not move after refinement."""
        domain = _make_domain(3, 3)
        orig_nodes = domain.mesh.nodes.copy()
        refined = anuga.uniform_refine_domain(domain)
        M0 = len(orig_nodes)
        self.assertTrue(
            num.allclose(refined.mesh.nodes[:M0], orig_nodes),
            'Original node coordinates changed after refinement')

    def test_new_nodes_are_midpoints(self):
        """Every new node must be the midpoint of some original edge."""
        domain = _make_domain(2, 2)
        orig_nodes = domain.mesh.nodes
        orig_tris  = domain.triangles
        M0 = len(orig_nodes)

        # Build set of original edge midpoints
        orig_midpoints = set()
        for v0, v1, v2 in orig_tris:
            for a, b in ((v0, v1), (v1, v2), (v0, v2)):
                mp = tuple(((orig_nodes[a] + orig_nodes[b]) * 0.5).tolist())
                orig_midpoints.add(mp)

        refined = anuga.uniform_refine_domain(domain)
        for mid_idx in range(M0, len(refined.mesh.nodes)):
            mid_pt = tuple(refined.mesh.nodes[mid_idx].tolist())
            self.assertIn(mid_pt, orig_midpoints,
                          f'Node {mid_idx} at {mid_pt} is not a midpoint of any original edge')

    def test_double_refinement(self):
        domain = _make_domain(2, 2)
        N0 = domain.number_of_triangles
        r1 = anuga.uniform_refine_domain(domain)
        r2 = anuga.uniform_refine_domain(r1)
        self.assertEqual(r2.number_of_triangles, 16 * N0)
        self.assertEqual(set(r2.boundary.values()), set(domain.boundary.values()))


class TestSequentialMeshRefine(unittest.TestCase):
    """Tests for sequential_mesh_refine (offline partition refinement)."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()
        self.domain = _make_domain(4, 4)
        self.N0 = self.domain.number_of_triangles
        anuga.sequential_mesh_dump(self.domain, numprocs=2,
                                    partition_dir=self.tmpdir, name='test')

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_files_created(self):
        anuga.sequential_mesh_refine('test', numprocs=2, levels=1,
                                      partition_dir=self.tmpdir)
        for p in range(2):
            fname = os.path.join(self.tmpdir, f'test_mesh_P2_{p}.nc')
            self.assertTrue(os.path.exists(fname))

    def test_global_triangle_count(self):
        anuga.sequential_mesh_refine('test', numprocs=2, levels=1,
                                      partition_dir=self.tmpdir)
        for p in range(2):
            fname = os.path.join(self.tmpdir, f'test_mesh_P2_{p}.nc')
            with netCDF4.Dataset(fname, 'r') as nc:
                self.assertEqual(int(nc.number_of_global_triangles),
                                  4 * self.N0)

    def test_full_triangle_sum(self):
        """Full-triangle counts across all ranks must sum to 4×N0."""
        anuga.sequential_mesh_refine('test', numprocs=2, levels=1,
                                      partition_dir=self.tmpdir)
        total = 0
        for p in range(2):
            fname = os.path.join(self.tmpdir, f'test_mesh_P2_{p}.nc')
            with netCDF4.Dataset(fname, 'r') as nc:
                total += int(nc.number_of_full_triangles)
        self.assertEqual(total, 4 * self.N0)

    def test_two_levels(self):
        """Two refinement levels produce 16×N0 triangles."""
        anuga.sequential_mesh_refine('test', numprocs=2, levels=2,
                                      partition_dir=self.tmpdir)
        total = 0
        for p in range(2):
            fname = os.path.join(self.tmpdir, f'test_mesh_P2_{p}.nc')
            with netCDF4.Dataset(fname, 'r') as nc:
                total += int(nc.number_of_full_triangles)
                self.assertEqual(int(nc.number_of_global_triangles), 16 * self.N0)
        self.assertEqual(total, 16 * self.N0)

    def test_vertex_indices_in_bounds(self):
        """All vertex indices must be valid local node indices."""
        anuga.sequential_mesh_refine('test', numprocs=2, levels=1,
                                      partition_dir=self.tmpdir)
        for p in range(2):
            fname = os.path.join(self.tmpdir, f'test_mesh_P2_{p}.nc')
            with netCDF4.Dataset(fname, 'r') as nc:
                M_local = len(nc['points'][:])
                vertices = num.array(nc['vertices'][:])
            self.assertTrue(num.all(vertices >= 0))
            self.assertTrue(num.all(vertices < M_local),
                             f'rank {p}: vertex index {vertices.max()} >= M_local {M_local}')

    def test_node_l2g_in_bounds(self):
        """All node_l2g values must be in [0, M_global)."""
        anuga.sequential_mesh_refine('test', numprocs=2, levels=1,
                                      partition_dir=self.tmpdir)
        for p in range(2):
            fname = os.path.join(self.tmpdir, f'test_mesh_P2_{p}.nc')
            with netCDF4.Dataset(fname, 'r') as nc:
                M_global = int(nc.number_of_global_nodes)
                node_l2g = num.array(nc['node_l2g'][:])
            self.assertTrue(num.all(node_l2g >= 0))
            self.assertTrue(num.all(node_l2g < M_global))

    def test_boundary_tags_preserved(self):
        orig_tags = set(self.domain.boundary.values())
        anuga.sequential_mesh_refine('test', numprocs=2, levels=1,
                                      partition_dir=self.tmpdir)
        for p in range(2):
            fname = os.path.join(self.tmpdir, f'test_mesh_P2_{p}.nc')
            with netCDF4.Dataset(fname, 'r') as nc:
                raw_tags = nc['boundary_tag'][:]
                tags = {b''.join(r).decode('ascii').rstrip('\x00')
                        for r in raw_tags}
            self.assertTrue(tags <= orig_tags | {'ghost'},
                             f'rank {p}: unexpected tags {tags - orig_tags}')

    def test_midpoint_coordinates_correct(self):
        """For a single-rank partition, new nodes are midpoints of original edges."""
        domain2 = _make_domain(2, 2)
        tmpdir2 = tempfile.mkdtemp()
        try:
            anuga.sequential_mesh_dump(domain2, numprocs=1,
                                        partition_dir=tmpdir2, name='t1')
            with netCDF4.Dataset(os.path.join(tmpdir2, 't1_mesh_P1_0.nc'), 'r') as nc:
                pts_orig = num.array(nc['points'][:])
                verts_orig = num.array(nc['vertices'][:])
            M0 = len(pts_orig)

            orig_midpoints = set()
            for v0, v1, v2 in verts_orig:
                for a, b in ((v0, v1), (v1, v2), (v0, v2)):
                    mp = tuple(((pts_orig[a] + pts_orig[b]) * 0.5).tolist())
                    orig_midpoints.add(mp)

            anuga.sequential_mesh_refine('t1', numprocs=1, levels=1,
                                          partition_dir=tmpdir2)
            with netCDF4.Dataset(os.path.join(tmpdir2, 't1_mesh_P1_0.nc'), 'r') as nc:
                pts_ref = num.array(nc['points'][:])

            self.assertTrue(num.allclose(pts_ref[:M0], pts_orig))
            for mid_idx in range(M0, len(pts_ref)):
                mid_pt = tuple(pts_ref[mid_idx].tolist())
                self.assertIn(mid_pt, orig_midpoints,
                               f'Node {mid_idx} is not a midpoint of any original edge')
        finally:
            import shutil
            shutil.rmtree(tmpdir2, ignore_errors=True)

    def test_output_dir_separate(self):
        """output_dir places refined files separately from input."""
        outdir = tempfile.mkdtemp()
        try:
            anuga.sequential_mesh_refine('test', numprocs=2, levels=1,
                                          partition_dir=self.tmpdir,
                                          output_dir=outdir)
            for p in range(2):
                self.assertTrue(
                    os.path.exists(os.path.join(outdir, f'test_mesh_P2_{p}.nc')))
        finally:
            import shutil
            shutil.rmtree(outdir, ignore_errors=True)


class TestCreateParallelMesh(unittest.TestCase):
    """Tests for create_parallel_mesh."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_no_refinement(self):
        domain = _make_domain(3, 3)
        N0 = domain.number_of_triangles
        name = anuga.create_parallel_mesh(domain, numprocs=2,
                                           refinement_levels=0,
                                           partition_dir=self.tmpdir)
        self.assertEqual(name, 'test_refine')
        for p in range(2):
            fname = os.path.join(self.tmpdir, f'{name}_mesh_P2_{p}.nc')
            self.assertTrue(os.path.exists(fname))
            with netCDF4.Dataset(fname, 'r') as nc:
                self.assertEqual(int(nc.number_of_global_triangles), N0)

    def test_one_level_refinement(self):
        domain = _make_domain(3, 3)
        N0 = domain.number_of_triangles
        name = anuga.create_parallel_mesh(domain, numprocs=2,
                                           refinement_levels=1,
                                           partition_dir=self.tmpdir)
        total = 0
        for p in range(2):
            fname = os.path.join(self.tmpdir, f'{name}_mesh_P2_{p}.nc')
            with netCDF4.Dataset(fname, 'r') as nc:
                total += int(nc.number_of_full_triangles)
                self.assertEqual(int(nc.number_of_global_triangles), 4 * N0)
        self.assertEqual(total, 4 * N0)

    def test_two_level_refinement(self):
        domain = _make_domain(3, 3)
        N0 = domain.number_of_triangles
        anuga.create_parallel_mesh(domain, numprocs=3,
                                    refinement_levels=2,
                                    name='mesh2',
                                    partition_dir=self.tmpdir)
        total = 0
        for p in range(3):
            fname = os.path.join(self.tmpdir, f'mesh2_mesh_P3_{p}.nc')
            self.assertTrue(os.path.exists(fname))
            with netCDF4.Dataset(fname, 'r') as nc:
                total += int(nc.number_of_full_triangles)
        self.assertEqual(total, 16 * N0)

    def test_returns_name(self):
        domain = _make_domain(2, 2)
        name = anuga.create_parallel_mesh(domain, numprocs=2,
                                           name='mymesh',
                                           partition_dir=self.tmpdir)
        self.assertEqual(name, 'mymesh')

    def test_no_intermediate_files_in_output_dir(self):
        """With refinement, only final-level files appear in partition_dir."""
        domain = _make_domain(2, 2)
        anuga.create_parallel_mesh(domain, numprocs=2,
                                    refinement_levels=2,
                                    name='mesh',
                                    partition_dir=self.tmpdir)
        nc_files = [f for f in os.listdir(self.tmpdir) if f.endswith('.nc')]
        self.assertEqual(len(nc_files), 2)  # exactly one .nc per rank

    def test_accepts_basic_mesh(self):
        """create_parallel_mesh works when passed a Basic_mesh instead of Domain."""
        from anuga import Basic_mesh
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross

        points, vertices, boundary = rectangular_cross(3, 3, len1=1.0, len2=1.0)
        mesh = Basic_mesh(points, vertices, boundary)
        N0 = mesh.number_of_triangles

        name = anuga.create_parallel_mesh(mesh, numprocs=2,
                                           refinement_levels=1,
                                           name='bm',
                                           partition_dir=self.tmpdir)
        self.assertEqual(name, 'bm')
        total = 0
        for p in range(2):
            fname = os.path.join(self.tmpdir, f'bm_mesh_P2_{p}.nc')
            self.assertTrue(os.path.exists(fname))
            with netCDF4.Dataset(fname, 'r') as nc:
                total += int(nc.number_of_full_triangles)
                self.assertEqual(int(nc.number_of_global_triangles), 4 * N0)
        self.assertEqual(total, 4 * N0)

    def test_basic_mesh_default_name(self):
        """Basic_mesh with no name defaults to 'mesh' for partition files."""
        from anuga import Basic_mesh
        from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross

        points, vertices, boundary = rectangular_cross(2, 2, len1=1.0, len2=1.0)
        mesh = Basic_mesh(points, vertices, boundary)

        name = anuga.create_parallel_mesh(mesh, numprocs=2,
                                           partition_dir=self.tmpdir)
        self.assertEqual(name, 'mesh')
        for p in range(2):
            fname = os.path.join(self.tmpdir, f'mesh_mesh_P2_{p}.nc')
            self.assertTrue(os.path.exists(fname))

    def test_parallel_workers_match_serial(self):
        """num_workers > 1 produces identical output to the serial path."""
        domain_s = _make_domain(4, 4)
        domain_p = _make_domain(4, 4)
        N0 = domain_s.number_of_triangles

        tmpdir_s = tempfile.mkdtemp()
        tmpdir_p = tempfile.mkdtemp()
        try:
            anuga.create_parallel_mesh(domain_s, numprocs=3,
                                        refinement_levels=1, name='ms',
                                        partition_dir=tmpdir_s,
                                        num_workers=1)
            anuga.create_parallel_mesh(domain_p, numprocs=3,
                                        refinement_levels=1, name='ms',
                                        partition_dir=tmpdir_p,
                                        num_workers=3)
            for p in range(3):
                fname_s = os.path.join(tmpdir_s, f'ms_mesh_P3_{p}.nc')
                fname_p = os.path.join(tmpdir_p, f'ms_mesh_P3_{p}.nc')
                with netCDF4.Dataset(fname_s, 'r') as ns, \
                     netCDF4.Dataset(fname_p, 'r') as np_:
                    self.assertEqual(int(ns.number_of_full_triangles),
                                     int(np_.number_of_full_triangles))
                    self.assertEqual(int(ns.number_of_global_triangles),
                                     4 * N0)
                    self.assertEqual(int(np_.number_of_global_triangles),
                                     4 * N0)
                    pts_s = num.array(ns['points'][:])
                    pts_p = num.array(np_['points'][:])
                    self.assertEqual(pts_s.shape, pts_p.shape)
        finally:
            import shutil
            shutil.rmtree(tmpdir_s, ignore_errors=True)
            shutil.rmtree(tmpdir_p, ignore_errors=True)


if __name__ == '__main__':
    unittest.main()
