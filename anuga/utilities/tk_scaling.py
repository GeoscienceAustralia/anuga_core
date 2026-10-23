"""HiDPI scaling helper for ANUGA's Tkinter GUIs.

On Wayland/XWayland sessions (and some misconfigured X servers) Tk reports a
bogus 96 DPI via ``winfo_fpixels('1i')`` -- and back-fills
``winfo_screenmmwidth`` to be consistent with it -- so estimating DPI from
either source also fails.  The result is that GUIs render physically tiny on
HiDPI panels.

``detect_dpi`` instead reads the user's *configured* DPI from the X resource
``Xft.dpi`` (the value GTK/Qt honour for HiDPI scaling, e.g. 192 -> 2x), and
``apply_hidpi_scaling`` uses it to enlarge Tk's named fonts so the whole UI
scales up.  On a normal 96 DPI display the scale resolves to 1.0 and nothing
changes, so this is safe to call unconditionally at startup.
"""

import os
import re
import subprocess

# Tk is imported defensively so this module stays importable without it.
#
# conda-forge's `tk` does not depend on `xorg-libx11`; it expects libX11.so.6
# from outside the prefix, so whether `_tkinter` imports depends on whether
# something else in the environment happened to pull libX11 in. That made the
# 4.0.0 conda-forge build fail on linux_64/py3.10 only, while the other four
# Python versions passed by luck of the solve (#266). Any headless or minimal
# Linux image hits the same thing.
#
# Nothing here needs a display: detect_dpi/detect_scale read Xft.dpi and xrandr
# and fall back to arithmetic, so they are worth keeping testable and callable
# where Tk will never exist. Only apply_hidpi_scaling touches Tk objects, and it
# already requires a live root passed in by the caller.
try:
    import tkinter as tk
    import tkinter.font as tkfont
    TclError = tk.TclError
except ImportError:  # no Tk, or Tk present but its C extension cannot load
    tk = None
    tkfont = None

    class TclError(Exception):
        """Placeholder so ``except TclError`` parses when Tk is absent.

        Unreachable in that case: the calls that raise it need a Tk root, which
        cannot exist without Tk.
        """


# Named fonts every Tk/ttk widget references by default; resizing these resizes
# the UI regardless of theme.
_NAMED_FONTS = ('TkDefaultFont', 'TkTextFont', 'TkFixedFont', 'TkMenuFont',
                'TkHeadingFont', 'TkCaptionFont', 'TkSmallCaptionFont',
                'TkIconFont', 'TkTooltipFont')


def _xft_dpi():
    """The user's configured Xft.dpi from ``xrdb -query``, or None."""
    try:
        out = subprocess.run(['xrdb', '-query'], capture_output=True,
                             text=True, timeout=2).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    m = re.search(r'^Xft\.dpi:\s*([0-9.]+)', out, re.M)
    return float(m.group(1)) if m else None


def detect_dpi(root):
    """Best-effort *true* screen DPI, working around XWayland's bogus 96.

    Order: the user's configured Xft.dpi (preferred -- it reflects intended
    scaling, not raw panel density), then the QT_FONT_DPI env hint, then Tk's
    own (often bogus) numbers as a last resort.
    """
    dpi = _xft_dpi()
    if dpi and dpi > 0:
        return dpi
    try:
        v = float(os.environ.get('QT_FONT_DPI', ''))
        if v > 0:
            return v
    except ValueError:
        pass
    try:
        px = root.winfo_screenwidth()
        mm = root.winfo_screenmmwidth()
        tk_dpi = (px * 25.4 / mm) if mm and mm > 0 else 96.0
    except TclError:
        tk_dpi = 96.0
    try:
        tk_dpi = max(tk_dpi, root.winfo_fpixels('1i'))
    except TclError:
        pass
    return tk_dpi


def _panel_dpi():
    """True panel DPI from xrandr's mode + physical size, or None.

    Used only as a *floor* for very dense panels.  Tk's own numbers are
    unusable under XWayland (it reports 96 DPI and a fictitious screen size
    such as 1016x635 mm), whereas xrandr reports the real panel, e.g.
    ``3840x2400 ... 340mm x 220mm`` -> 287 DPI.
    """
    try:
        out = subprocess.run(['xrandr'], capture_output=True,
                             text=True, timeout=2).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    # e.g. "eDP-1 connected primary 3840x2400+0+0 (normal ...) 340mm x 220mm"
    m = re.search(r'^\S+ connected.*?(\d+)x(\d+)\+\d+\+\d+.*?(\d+)mm x (\d+)mm',
                  out, re.M)
    if not m:
        return None
    px, mm = float(m.group(1)), float(m.group(3))
    return px * 25.4 / mm if mm > 0 else None


def _override_scale():
    """Explicit user override from ANUGA_GUI_SCALE, or None."""
    raw = os.environ.get('ANUGA_GUI_SCALE', '').strip()
    if not raw:
        return None
    try:
        v = float(raw)
    except ValueError:
        return None
    return v if v > 0 else None


def detect_scale(root, lo=1.0, hi=3.0):
    """Font/UI scale factor (1.0 == 96 DPI), clamped to [lo, hi].

    ANUGA_GUI_SCALE overrides everything (clamped only to something sane), for
    displays where the automatic answer is not to taste.

    Otherwise the configured Xft.dpi leads, since it reflects the scaling the
    user asked their desktop for.  On very dense panels that setting can still
    leave the UI physically small -- a 287 DPI laptop panel set to Xft.dpi 192
    gets 2.0x, where the panel itself implies 3.0x -- so the true panel density
    is used to lift the scale halfway towards it.
    """
    override = _override_scale()
    if override is not None:
        return round(min(4.0, max(0.5, override)), 2)

    scale = detect_dpi(root) / 96.0

    panel = _panel_dpi()
    if panel:
        panel_scale = panel / 96.0
        if panel_scale > scale:
            # Halfway: honour the user's setting but do not ignore the panel.
            scale = (scale + panel_scale) / 2.0

    return round(min(hi, max(lo, scale)), 2)


def apply_hidpi_scaling(root, lo=1.0, hi=3.0):
    """Scale Tk's fonts and widget geometry to suit the display; return scale.

    A no-op (returns 1.0) on normal-DPI displays, so it is safe to call
    unconditionally just after creating the root window.
    """
    scale = detect_scale(root, lo, hi)
    if scale <= 1.0:
        return scale

    for name in _NAMED_FONTS:
        try:
            f = tkfont.nametofont(name)
        except TclError:
            continue
        size = f.cget('size')
        if not size:        # 0 == "unspecified"; leave it to Tk
            continue
        # Tk sizes: positive = points, negative = pixels. Preserve the sign.
        sign = 1 if size > 0 else -1
        f.configure(size=sign * max(8, int(round(abs(size) * scale))))

    # Fonts alone leave every *non-text* element the same physical size:
    # padding, borders, ttk element geometry and anything sized in points all
    # derive from Tk's points-per-pixel, which stays at the 96 DPI value unless
    # told otherwise.  That is why scaled fonts alone still look cramped.
    try:
        root.tk.call('tk', 'scaling', (96.0 * scale) / 72.0)
    except TclError:
        pass

    return scale
