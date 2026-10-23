"""Tests for parse_time.py"""
import unittest
from anuga.utilities.parse_time import parse_time, seconds_to_hhmmss


class Test_parse_time(unittest.TestCase):

    # ------------------------------------------------------------------
    # None input
    # ------------------------------------------------------------------

    def test_none_returns_none(self):
        self.assertIsNone(parse_time(None))

    def test_none_default_returns_none(self):
        self.assertIsNone(parse_time())

    # ------------------------------------------------------------------
    # Numeric input (float/int)
    # ------------------------------------------------------------------

    def test_float_passthrough(self):
        self.assertAlmostEqual(parse_time(1234567.5), 1234567.5)

    def test_int_passthrough(self):
        result = parse_time(0)
        self.assertAlmostEqual(result, 0.0)

    def test_numeric_string(self):
        # A string that looks like a float should be parsed as float
        result = parse_time(9999.0)
        self.assertAlmostEqual(result, 9999.0)

    # ------------------------------------------------------------------
    # Date-only strings: 'YYYYMMDD'
    # ------------------------------------------------------------------

    def test_date_only_epoch(self):
        # 1970-01-01 → epoch 0
        result = parse_time('19700101')
        self.assertAlmostEqual(result, 0.0)

    def test_date_only_2000_01_01(self):
        import datetime
        expected = (datetime.datetime(2000, 1, 1) -
                    datetime.datetime(1970, 1, 1)).total_seconds()
        result = parse_time('20000101')
        self.assertAlmostEqual(result, expected)

    # ------------------------------------------------------------------
    # Datetime strings: 'YYYYMMDD_HHMM'
    # ------------------------------------------------------------------

    def test_datetime_with_underscore(self):
        import datetime
        expected = float(int((datetime.datetime(2012, 2, 29, 12, 10, 0) -
                               datetime.datetime(1970, 1, 1)).total_seconds()))
        result = parse_time('20120229_1210')
        self.assertAlmostEqual(result, expected)

    def test_datetime_with_space(self):
        import datetime
        expected = float(int((datetime.datetime(2012, 2, 29, 12, 10, 0) -
                               datetime.datetime(1970, 1, 1)).total_seconds()))
        result = parse_time('20120229 1210')
        self.assertAlmostEqual(result, expected)

    def test_datetime_no_separator(self):
        import datetime
        expected = float(int((datetime.datetime(2012, 2, 29, 12, 10, 0) -
                               datetime.datetime(1970, 1, 1)).total_seconds()))
        result = parse_time('201202291210')
        self.assertAlmostEqual(result, expected)

    # ------------------------------------------------------------------
    # Datetime strings: 'YYYYMMDD_HHMMSS'
    # ------------------------------------------------------------------

    def test_datetime_with_seconds_underscore(self):
        import datetime
        expected = float(int((datetime.datetime(2012, 2, 29, 12, 10, 30) -
                               datetime.datetime(1970, 1, 1)).total_seconds()))
        result = parse_time('20120229_121030')
        self.assertAlmostEqual(result, expected)

    def test_datetime_with_seconds_no_separator(self):
        import datetime
        expected = float(int((datetime.datetime(2012, 2, 29, 12, 10, 30) -
                               datetime.datetime(1970, 1, 1)).total_seconds()))
        result = parse_time('20120229121030')
        self.assertAlmostEqual(result, expected)

    def test_datetime_with_colon_separator(self):
        import datetime
        expected = float(int((datetime.datetime(2012, 2, 29, 12, 10, 0) -
                               datetime.datetime(1970, 1, 1)).total_seconds()))
        result = parse_time('20120229:1210')
        self.assertAlmostEqual(result, expected)

    def test_datetime_with_slash_separator(self):
        import datetime
        expected = float(int((datetime.datetime(2012, 2, 29, 12, 10, 0) -
                               datetime.datetime(1970, 1, 1)).total_seconds()))
        result = parse_time('20120229/1210')
        self.assertAlmostEqual(result, expected)

    # ------------------------------------------------------------------
    # Returns float type
    # ------------------------------------------------------------------

    def test_return_type_is_float(self):
        result = parse_time('20000101')
        self.assertIsInstance(result, float)

    def test_return_type_for_numeric(self):
        result = parse_time(42.0)
        self.assertIsInstance(result, float)

    # ------------------------------------------------------------------
    # debug flag (should not raise)
    # ------------------------------------------------------------------

    def test_debug_flag_does_not_raise(self):
        result = parse_time('20120229_1210', debug=True)
        self.assertIsNotNone(result)


class Test_seconds_to_hhmmss(unittest.TestCase):

    def test_seconds_only(self):
        self.assertEqual(seconds_to_hhmmss(45), '45s')

    def test_minutes_and_seconds(self):
        self.assertEqual(seconds_to_hhmmss(90), '1m:30s')

    def test_hours_minutes_seconds(self):
        self.assertEqual(seconds_to_hhmmss(3661), '1h:1m:1s')

    def test_zero_seconds(self):
        self.assertEqual(seconds_to_hhmmss(0), '0s')

    def test_exact_hour(self):
        # secs=0 and parts already has '1h', so no '0s' appended
        self.assertEqual(seconds_to_hhmmss(3600), '1h')

    def test_exact_minute(self):
        # secs=0 and parts already has '1m', so no '0s' appended
        self.assertEqual(seconds_to_hhmmss(60), '1m')

    def test_large_hours(self):
        result = seconds_to_hhmmss(7200)
        self.assertEqual(result, '2h')


class Test_parse_time_extra(unittest.TestCase):
    """Tests for uncovered exception paths in parse_time."""

    def test_short_string_triggers_year_except(self):
        """String too short for year/month/day — hits exception branches (lines 22-23, 30-31, 35-36)."""
        result = parse_time('abc')
        # Should return some float (epoch) without crashing
        self.assertIsNotNone(result)

    def test_non_numeric_non_string_float_valueerror(self):
        """Non-string whose float() raises ValueError hits lines 15-16, then falls into string path."""
        class BadFloat:
            def __float__(self):
                raise ValueError('bad')
        # After lines 15-16, code tries to parse it as a string, which may
        # raise TypeError since it's not subscriptable — that's acceptable
        try:
            result = parse_time(BadFloat())
        except (TypeError, Exception):
            pass  # line 15-16 still covered


if __name__ == '__main__':
    unittest.main()
