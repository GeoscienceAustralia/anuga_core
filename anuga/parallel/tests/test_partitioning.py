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
