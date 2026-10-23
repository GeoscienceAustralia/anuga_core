#!/usr/bin/env python
#
# Tests for anuga.utilities.tk_scaling — the HiDPI scale detection.
#
# These are pure logic tests: the DPI sources are monkeypatched, so no display,
# no X server and no xrandr/xrdb binaries are needed.

import unittest
from unittest import mock

from anuga.utilities import tk_scaling


class _FakeRoot:
    """Stand-in for a Tk root; only used when the X sources give nothing."""

    def winfo_screenwidth(self):
        return 1920

    def winfo_screenmmwidth(self):
        return 508          # exactly 96 DPI

    def winfo_fpixels(self, spec):
        return 96.0


class Test_detect_scale(unittest.TestCase):

    def setUp(self):
        self.root = _FakeRoot()

    def _detect(self, xft=None, panel=None, env=None):
        with mock.patch.object(tk_scaling, '_xft_dpi', return_value=xft), \
             mock.patch.object(tk_scaling, '_panel_dpi', return_value=panel), \
             mock.patch.dict('os.environ', env or {}, clear=False):
            if env is None:
                with mock.patch.dict('os.environ', {'ANUGA_GUI_SCALE': ''}):
                    return tk_scaling.detect_scale(self.root)
            return tk_scaling.detect_scale(self.root)

    def test_normal_dpi_is_a_no_op(self):
        """A plain 96 DPI display must not be scaled at all."""
        self.assertEqual(self._detect(xft=96.0, panel=96.0), 1.0)

    def test_configured_xft_dpi_is_honoured(self):
        self.assertEqual(self._detect(xft=192.0, panel=192.0), 2.0)

    def test_dense_panel_lifts_a_conservative_xft_setting(self):
        """287 DPI panel + Xft.dpi 192 -> between the two, not just 2.0.

        This is the laptop case that prompted the change: the desktop is
        configured for 2x, but on a ~287 DPI panel that still leaves the UI
        physically small.
        """
        scale = self._detect(xft=192.0, panel=286.9)
        self.assertGreater(scale, 2.0)
        self.assertLessEqual(scale, 3.0)

    def test_panel_never_shrinks_the_users_setting(self):
        """A low-density panel must not drag a deliberate 2x setting down."""
        self.assertEqual(self._detect(xft=192.0, panel=96.0), 2.0)

    def test_env_override_wins(self):
        scale = self._detect(xft=192.0, panel=286.9,
                             env={'ANUGA_GUI_SCALE': '1.25'})
        self.assertEqual(scale, 1.25)

    def test_env_override_can_disable_scaling(self):
        scale = self._detect(xft=192.0, panel=286.9,
                             env={'ANUGA_GUI_SCALE': '1'})
        self.assertEqual(scale, 1.0)

    def test_junk_override_is_ignored(self):
        """A malformed override must fall back to detection, not crash."""
        scale = self._detect(xft=192.0, panel=192.0,
                             env={'ANUGA_GUI_SCALE': 'huge'})
        self.assertEqual(scale, 2.0)

    def test_scale_is_clamped(self):
        """An absurdly dense reading cannot blow the UI up without bound."""
        self.assertLessEqual(self._detect(xft=1200.0, panel=1200.0), 3.0)


if __name__ == '__main__':
    unittest.main()
