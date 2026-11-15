import subprocess
import sys

import pytest
import numpy as np

from anuga.parallel.partitioning import morton_order_from_points
from anuga.parallel.partitioning import hilbert_order_from_points

def test_morton_order_from_points_basic():
    """Test basic Morton ordering with simple points."""
    points = np.array([[0.0, 0.0], [1.0, 1.0], [0.5, 0.5]])
    order = morton_order_from_points(points)
    assert len(order) == 3
    assert np.all(np.isin(order, [0, 1, 2]))


def test_morton_order_from_points_sorted():
    """Test that returned order is valid permutation."""
    points = np.array([[0.1, 0.2], [0.5, 0.5], [0.9, 0.8]])
    order = morton_order_from_points(points)
    assert np.array_equal(np.sort(order), np.arange(len(points)))


def test_morton_order_from_points_single_point():
    """Test with single point."""
    points = np.array([[0.5, 0.5]])
    order = morton_order_from_points(points)
    assert len(order) == 1
    assert order[0] == 0



def test_morton_order_from_points_invalid_shape():
    """Test error handling for invalid point shapes."""
    with pytest.raises(ValueError, match="points must have shape"):
        morton_order_from_points(np.array([0.0, 1.0]))
    with pytest.raises(ValueError, match="points must have shape"):
        morton_order_from_points(np.array([[0.0, 1.0, 0.5]]))


def test_morton_order_from_points_degenerate():
    """Test with degenerate points (zero span in one direction)."""
    points = np.array([[0.5, 0.0], [0.5, 1.0], [0.5, 0.5]])
    order = morton_order_from_points(points)
    assert len(order) == 3
    assert np.array_equal(np.sort(order), np.arange(3))


def test_morton_order_from_points_large_array():
    """Test with larger array of points."""
    rng = np.random.default_rng(42)
    points = rng.random((100, 2))
    order = morton_order_from_points(points)
    assert len(order) == 100
    assert np.array_equal(np.sort(order), np.arange(100))


def test_morton_order_from_points_type_conversion():
    """Test that function converts input to float64."""
    points = np.array([[0, 0], [1, 1]], dtype=np.int32)
    order = morton_order_from_points(points)
    assert len(order) == 2
    assert np.array_equal(np.sort(order), np.arange(2))


def test_hilbert_order_from_points_basic():
    """Test basic Hilbert ordering with simple points."""
    points = np.array([[0.0, 0.0], [1.0, 1.0], [0.5, 0.5]])
    order = hilbert_order_from_points(points)
    assert len(order) == 3
    assert np.all(np.isin(order, [0, 1, 2]))


def test_hilbert_order_from_points_sorted():
    """Test that returned order is valid permutation."""
    points = np.array([[0.1, 0.2], [0.5, 0.5], [0.9, 0.8]])
    order = hilbert_order_from_points(points)
    assert np.array_equal(np.sort(order), np.arange(len(points)))


def test_hilbert_order_from_points_single_point():
    """Test with single point."""
    points = np.array([[0.5, 0.5]])
    order = hilbert_order_from_points(points)
    assert len(order) == 1
    assert order[0] == 0


def test_hilbert_order_from_points_invalid_shape():
    """Test error handling for invalid point shapes."""
    with pytest.raises(ValueError, match="points must have shape"):
        hilbert_order_from_points(np.array([0.0, 1.0]))
    with pytest.raises(ValueError, match="points must have shape"):
        hilbert_order_from_points(np.array([[0.0, 1.0, 0.5]]))


def test_hilbert_order_from_points_degenerate():
    """Test with degenerate points (zero span in one direction)."""
    points = np.array([[0.5, 0.0], [0.5, 1.0], [0.5, 0.5]])
    order = hilbert_order_from_points(points)
    assert len(order) == 3
    assert np.array_equal(np.sort(order), np.arange(3))


def test_hilbert_order_from_points_large_array():
    """Test with larger array of points."""
    rng = np.random.default_rng(42)
    points = rng.random((100, 2))
    order = hilbert_order_from_points(points)
    assert len(order) == 100
    assert np.array_equal(np.sort(order), np.arange(100))


def test_hilbert_order_from_points_type_conversion():
    """Test that function converts input to float64."""
    points = np.array([[0, 0], [1, 1]], dtype=np.int32)
    order = hilbert_order_from_points(points)
    assert len(order) == 2
    assert np.array_equal(np.sort(order), np.arange(2))


def test_hilbert_order_from_points_custom_p():
    """Test with custom precision parameter."""
    points = np.array([[0.0, 0.0], [1.0, 1.0], [0.5, 0.5]])
    order = hilbert_order_from_points(points, p=8)
    assert len(order) == 3
    assert np.array_equal(np.sort(order), np.arange(3))


def test_hilbert_order_from_points_all_same():
    """Test with all identical points."""
    points = np.array([[0.5, 0.5], [0.5, 0.5], [0.5, 0.5]])
    order = hilbert_order_from_points(points)
    assert len(order) == 3
    assert np.array_equal(np.sort(order), np.arange(3))


# ---------------------------------------------------------------------------
# pymetis is optional: pip does not install it on Windows
# ---------------------------------------------------------------------------

def _small_domain():
    import anuga
    domain = anuga.rectangular_cross_domain(4, 4)
    domain.set_name('test_no_pymetis')
    domain.set_store(False)
    domain.set_quantity('elevation', 0.0)
    domain.set_quantity('stage', 1.0)
    domain.set_boundary({'left': anuga.Reflective_boundary(domain),
                         'right': anuga.Reflective_boundary(domain),
                         'top': anuga.Reflective_boundary(domain),
                         'bottom': anuga.Reflective_boundary(domain)})
    return domain


def test_without_pymetis_metis_partition_says_what_to_do(monkeypatch):
    from anuga.parallel.partitioning import metis_partition
    domain = _small_domain()
    monkeypatch.setitem(sys.modules, 'pymetis', None)   # makes the import fail

    with pytest.raises(ImportError, match='conda-forge pymetis.*hilbert'):
        metis_partition(domain, 2)


def test_without_pymetis_one_part_and_other_schemes_still_work(monkeypatch):
    from anuga.parallel.partitioning import (metis_partition, hilbert_partition,
                                             morton_partition, rcm_partition,
                                             reorder_domain)
    domain = _small_domain()
    n_tri = domain.number_of_triangles
    monkeypatch.setitem(sys.modules, 'pymetis', None)

    order, counts = metis_partition(domain, 1)
    assert np.array_equal(order, np.arange(n_tri))
    assert list(counts) == [n_tri]

    for partition in (hilbert_partition, morton_partition, rcm_partition):
        order, counts = partition(domain, 2)
        assert np.array_equal(np.sort(order), np.arange(n_tri))
        assert sum(counts) == n_tri

    reorder_domain(domain, method='hilbert')


def test_without_pymetis_serial_run_works(monkeypatch):
    monkeypatch.setitem(sys.modules, 'pymetis', None)
    domain = _small_domain()
    for _ in domain.evolve(yieldstep=0.1, finaltime=0.1):
        pass
    assert domain.get_time() == pytest.approx(0.1)


def test_import_anuga_without_pymetis():
    code = "import sys; sys.modules['pymetis'] = None; import anuga; print(anuga.__version__)"
    result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
