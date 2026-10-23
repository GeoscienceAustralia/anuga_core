"""Tests for reorder_domain() and domain.reorder().

Verifies that reordering triangles for cache locality:
  - produces a valid permuted mesh (neighbours consistent)
  - preserves centroid_values after reorder
  - gives numerically identical simulation results
"""

import tempfile
import unittest

import numpy as np
import pytest

import anuga
from anuga import rectangular_cross_domain, Dirichlet_boundary, Reflective_boundary
from anuga.parallel.partitioning import reorder_domain


def topography(x, y):
    return -x / 2


def run_simulation(use_reorder=False, method='hilbert'):
    domain = rectangular_cross_domain(20, 20, 1.0, 1.0)
    domain.set_quantity('elevation', topography)
    domain.set_quantity('friction', 0.0)
    domain.set_quantity('stage', expression='elevation')
    domain.set_flow_algorithm('DE1')
    domain.set_quantities_to_be_stored(None)
    domain.set_name('test_reorder')
    domain.set_datadir(tempfile.mkdtemp())

    if use_reorder:
        reorder_domain(domain, method=method)

    Br = Reflective_boundary(domain)
    Bd = Dirichlet_boundary([-0.2, 0., 0.])
    domain.set_boundary({'left': Br, 'right': Bd, 'top': Br, 'bottom': Br})

    for t in domain.evolve(yieldstep=0.25, finaltime=1.0):
        pass

    # Return stage and elevation sorted by centroid (x, y) for order-independent comparison
    coords = domain.centroid_coordinates.copy()
    stage = domain.quantities['stage'].centroid_values.copy()
    elev  = domain.quantities['elevation'].centroid_values.copy()
    idx = np.lexsort((coords[:, 1], coords[:, 0]))
    return stage[idx], elev[idx], coords[idx]


class TestReorderDomain(unittest.TestCase):

    def test_hilbert_reorder_preserves_simulation(self):
        """Hilbert reorder gives same final stage as no-reorder."""
        stage_ref, elev_ref, coords_ref = run_simulation(use_reorder=False)
        stage_hl,  elev_hl,  coords_hl  = run_simulation(use_reorder=True, method='hilbert')

        np.testing.assert_allclose(coords_ref, coords_hl, atol=1e-10,
                                   err_msg='centroid coords differ after hilbert reorder')
        np.testing.assert_allclose(elev_ref, elev_hl, atol=1e-10,
                                   err_msg='elevation differs after hilbert reorder')
        np.testing.assert_allclose(stage_ref, stage_hl, atol=1e-10,
                                   err_msg='stage differs after hilbert reorder')

    def test_morton_reorder_preserves_simulation(self):
        """Morton reorder gives same final stage as no-reorder."""
        stage_ref, elev_ref, coords_ref = run_simulation(use_reorder=False)
        stage_mt,  elev_mt,  coords_mt  = run_simulation(use_reorder=True, method='morton')

        np.testing.assert_allclose(coords_ref, coords_mt, atol=1e-10,
                                   err_msg='centroid coords differ after morton reorder')
        np.testing.assert_allclose(stage_ref, stage_mt, atol=1e-10,
                                   err_msg='stage differs after morton reorder')

    def test_reorder_mesh_neighbours_consistent(self):
        """After reorder every neighbour pair shares exactly 2 nodes."""
        domain = rectangular_cross_domain(15, 15)
        reorder_domain(domain, method='hilbert')

        N    = domain.number_of_triangles
        tris = domain.triangles
        nbrs = domain.mesh.neighbours

        for i in range(N):
            for j in range(3):
                k = nbrs[i, j]
                if k < 0:
                    continue
                shared = set(tris[i]) & set(tris[k])
                self.assertEqual(len(shared), 2,
                    f'Triangle {i} edge {j}: neighbour {k} shares {len(shared)} nodes '
                    f'(expected 2) — neighbours stale after reorder')

    def test_reorder_tri_full_flag_sequential(self):
        """tri_full_flag stays all-ones on a sequential domain."""
        domain = rectangular_cross_domain(10, 10)
        reorder_domain(domain, method='morton')
        N = domain.number_of_triangles
        np.testing.assert_array_equal(domain.tri_full_flag, np.ones(N, int),
                                      err_msg='tri_full_flag corrupted by reorder')

    def test_reorder_centroid_coordinates_permuted(self):
        """domain.centroid_coordinates are reordered consistently with quantities."""
        domain = rectangular_cross_domain(10, 10)
        # Store original centroid coords and elevation centroid values
        coords_before = domain.centroid_coordinates.copy()
        elev_before   = domain.quantities['elevation'].centroid_values.copy()

        from anuga.parallel.partitioning import hilbert_partition
        epart_order, _ = hilbert_partition(domain, 1)
        domain.reorder(epart_order)

        coords_after = domain.centroid_coordinates
        elev_after   = domain.quantities['elevation'].centroid_values

        # After reorder: new[i] == old[epart_order[i]]
        np.testing.assert_allclose(coords_after, coords_before[epart_order], atol=1e-12,
                                   err_msg='centroid_coordinates not reordered consistently')
        np.testing.assert_allclose(elev_after, elev_before[epart_order], atol=1e-12,
                                   err_msg='elevation centroid_values not reordered consistently')

    def test_metis_reorder_preserves_simulation(self):
        """Metis reorder gives same final stage as no-reorder."""
        pytest.importorskip('pymetis')
        stage_ref, _, coords_ref = run_simulation(use_reorder=False)
        stage_mt,  _, coords_mt  = run_simulation(use_reorder=True, method='metis')

        np.testing.assert_allclose(coords_ref, coords_mt, atol=1e-10)
        np.testing.assert_allclose(stage_ref, stage_mt, atol=1e-10,
                                   err_msg='stage differs after metis reorder')


if __name__ == '__main__':
    unittest.main()
