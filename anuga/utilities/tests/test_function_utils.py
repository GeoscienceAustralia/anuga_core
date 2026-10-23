import unittest
import numpy as num

from anuga.utilities.function_utils import determine_function_type


class Test_Function_Utils(unittest.TestCase):

    def test_determine_function_type_scalar(self):

        type = determine_function_type(1.0)

        assert type == 'scalar'

        type = determine_function_type(1)

        assert type == 'scalar'

    def test_determine_function_type_exception(self):

        type = determine_function_type([])

        assert type == 'array'

    def test_determine_function_type_time_only(self):

        def myfunc(t):
            return t*2

        type = determine_function_type(myfunc)

        assert type == 't'

    def test_determine_function_type_spatial_only(self):

        def myfunc(x, y):
            return x+y

        type = determine_function_type(myfunc)

        assert type == 'x,y'

    def test_determine_function_type_spatial_only(self):

        def myfunc(x, y, t):
            return x+y+t

        type = determine_function_type(myfunc)

        assert type == 'x,y,t'


# -------------------------------------------------------------

class Test_Function_Utils_extra(unittest.TestCase):
    """Tests for uncovered paths in function_utils."""

    def test_determine_function_type_none(self):
        """None returns None (line 26)."""
        result = determine_function_type(None)
        self.assertIsNone(result)

    def test_determine_function_type_ndarray(self):
        """ndarray returns 'array' (lines 81-82)."""
        import numpy as np
        result = determine_function_type(np.array([1.0, 2.0]))
        self.assertEqual(result, 'array')

    def test_determine_function_type_t_valueerror(self):
        """Function of t that raises ValueError returns 't' (lines 57-59)."""
        def f(t):
            if isinstance(t, float):
                raise ValueError('out of range')
            raise TypeError('not a float')
        result = determine_function_type(f)
        self.assertEqual(result, 't')

    def test_evaluate_temporal_function_scalar(self):
        """evaluate_temporal_function with non-callable returns value (line 115)."""
        from anuga.utilities.function_utils import evaluate_temporal_function
        result = evaluate_temporal_function(42.0, t=1.0)
        self.assertAlmostEqual(result, 42.0)

    def test_evaluate_temporal_function_callable(self):
        """evaluate_temporal_function with callable (line 89)."""
        from anuga.utilities.function_utils import evaluate_temporal_function
        result = evaluate_temporal_function(lambda t: t * 2.0, t=3.0)
        self.assertAlmostEqual(result, 6.0)

    def test_determine_function_type_all_raises(self):
        """Function raising TypeError for all call forms raises Exception (lines 55-56)."""
        def bad(*args, **kwargs):
            raise TypeError('not callable')
        with self.assertRaises(Exception):
            determine_function_type(bad)

    def test_determine_function_type_xy_valueerror(self):
        """f(x,y,t) TypeError but f(x,y) ValueError returns 'x,y' (lines 69-70)."""
        def f(*args):
            if len(args) == 3:
                raise TypeError('no t')
            raise ValueError('out of range')
        result = determine_function_type(f)
        self.assertEqual(result, 'x,y')

    def test_determine_function_type_xyt_valueerror(self):
        """f(x,y,t) raises ValueError returns 'x,y,t' (lines 73-74)."""
        def f(*args):
            raise ValueError('out of range')
        result = determine_function_type(f)
        self.assertEqual(result, 'x,y,t')

    def test_determine_function_type_t_modeltime_early(self):
        """f(t) raises Modeltime_too_early returns 't' (lines 60-62)."""
        from anuga.fit_interpolate.interpolate import Modeltime_too_early
        def f(*args):
            if len(args) != 1:
                raise TypeError('not f(t)')
            raise Modeltime_too_early('too early')
        result = determine_function_type(f)
        self.assertEqual(result, 't')

    def test_determine_function_type_t_modeltime_late(self):
        """f(t) raises Modeltime_too_late returns 't' (lines 63-65)."""
        from anuga.fit_interpolate.interpolate import Modeltime_too_late
        def f(*args):
            if len(args) != 1:
                raise TypeError('not f(t)')
            raise Modeltime_too_late('too late')
        result = determine_function_type(f)
        self.assertEqual(result, 't')

    def test_evaluate_temporal_modeltime_early_default_none_raises(self):
        """evaluate_temporal_function Modeltime_too_early raises when default_left=None (lines 93-94)."""
        from anuga.utilities.function_utils import evaluate_temporal_function
        from anuga.fit_interpolate.interpolate import Modeltime_too_early

        def f(t):
            raise Modeltime_too_early('too early')

        with self.assertRaises(Modeltime_too_early):
            evaluate_temporal_function(f, t=0.0)

    def test_evaluate_temporal_modeltime_early_default_callable(self):
        """evaluate_temporal_function uses callable default_left_value (line 102)."""
        import warnings
        from anuga.utilities.function_utils import evaluate_temporal_function
        from anuga.fit_interpolate.interpolate import Modeltime_too_early

        def f(t):
            raise Modeltime_too_early('too early')

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            result = evaluate_temporal_function(f, t=0.0,
                                                default_left_value=lambda t: 99.0)
        self.assertAlmostEqual(result, 99.0)

    def test_evaluate_temporal_modeltime_late_default_none_raises(self):
        """evaluate_temporal_function Modeltime_too_late raises when default_right=None (lines 106-107)."""
        from anuga.utilities.function_utils import evaluate_temporal_function
        from anuga.fit_interpolate.interpolate import Modeltime_too_late

        def f(t):
            raise Modeltime_too_late('too late')

        with self.assertRaises(Modeltime_too_late):
            evaluate_temporal_function(f, t=100.0)

    def test_evaluate_temporal_modeltime_late_default_scalar(self):
        """evaluate_temporal_function uses scalar default_right_value (line 115)."""
        import warnings
        from anuga.utilities.function_utils import evaluate_temporal_function
        from anuga.fit_interpolate.interpolate import Modeltime_too_late

        def f(t):
            raise Modeltime_too_late('too late')

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', UserWarning)
            result = evaluate_temporal_function(f, t=100.0, default_right_value=42.0)
        self.assertAlmostEqual(result, 42.0)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test_Function_Utils)
    runner = unittest.TextTestRunner()  # verbosity=2)
    runner.run(suite)
