"""Tests for where_close.py"""
import unittest
import numpy as num
from anuga.utilities.where_close import where_close


class Test_where_close(unittest.TestCase):

    # ------------------------------------------------------------------
    # Float comparisons
    # ------------------------------------------------------------------

    def test_float_arrays_equal(self):
        x = [1.0, 2.0, 3.0]
        y = [1.0, 2.0, 3.0]
        result = where_close(x, y)
        num.testing.assert_array_equal(result, [1, 1, 1])

    def test_float_arrays_not_equal(self):
        x = [1.0, 2.0, 3.0]
        y = [1.1, 2.0, 4.0]
        result = where_close(x, y)
        num.testing.assert_array_equal(result, [0, 1, 0])

    def test_float_within_tolerance(self):
        # 2.0 and 2.000000000001 should be close
        x = [2.0]
        y = [2.000000000001]
        result = where_close(x, y)
        self.assertEqual(result[0], 1)

    def test_float_scalar_broadcast(self):
        x = [1.0, 2.0, 3.0]
        y = 2.0
        result = where_close(x, y)
        num.testing.assert_array_equal(result, [0, 1, 0])

    def test_float_negative_values(self):
        x = [-32.0]
        y = [-32.0]
        result = where_close(x, y)
        self.assertEqual(result[0], 1)

    def test_2d_float_array(self):
        x = num.array([[1.0, -1.0], [2.0, 3.0]])
        y = num.array([[1.0, -1.0], [2.0, 4.0]])
        result = where_close(x, y)
        self.assertEqual(result[0, 0], 1)
        self.assertEqual(result[0, 1], 1)
        self.assertEqual(result[1, 0], 1)
        self.assertEqual(result[1, 1], 0)

    # ------------------------------------------------------------------
    # Integer comparisons
    # ------------------------------------------------------------------

    def test_integer_arrays_equal(self):
        x = num.array([1, 5, 7])
        y = num.array([1, 5, 7])
        result = where_close(x, y)
        num.testing.assert_array_equal(result, [1, 1, 1])

    def test_integer_arrays_not_equal(self):
        x = num.array([1, 5, 7, -2, 10])
        y = num.array([1, -5, 17, -2, 0])
        result = where_close(x, y)
        num.testing.assert_array_equal(result, [1, 0, 0, 1, 0])

    def test_integer_strict_equality(self):
        # Integers use strict equality even if they're close
        x = num.array([5])
        y = num.array([6])
        result = where_close(x, y)
        self.assertEqual(result[0], 0)

    # ------------------------------------------------------------------
    # Mixed float/integer (should use float path)
    # ------------------------------------------------------------------

    def test_mixed_int_float(self):
        x = num.array([1, 2])
        y = num.array([1.0, 2.0])
        result = where_close(x, y)
        num.testing.assert_array_equal(result, [1, 1])

    # ------------------------------------------------------------------
    # Custom tolerances
    # ------------------------------------------------------------------

    def test_custom_atol(self):
        x = [1.0]
        y = [1.05]
        # Default rtol=1e-5, atol=1e-8 → not close
        result_default = where_close(x, y)
        self.assertEqual(result_default[0], 0)
        # With atol=0.1 → close
        result_custom = where_close(x, y, atol=0.1)
        self.assertEqual(result_custom[0], 1)

    def test_custom_rtol(self):
        x = [100.0]
        y = [101.0]
        # rtol=0.01 → not close (1% tolerance but diff = 1%)
        result = where_close(x, y, rtol=0.005, atol=0.0)
        self.assertEqual(result[0], 0)
        result2 = where_close(x, y, rtol=0.02, atol=0.0)
        self.assertEqual(result2[0], 1)

    # ------------------------------------------------------------------
    # Invalid type raises ValueError
    # ------------------------------------------------------------------

    def test_complex_input_raises(self):
        x = num.array([1+2j])
        y = num.array([1+2j])
        with self.assertRaises(ValueError):
            where_close(x, y)


if __name__ == '__main__':
    unittest.main()
