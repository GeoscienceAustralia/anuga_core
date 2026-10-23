"""`import anuga` must not require matplotlib.

anuga/__init__.py used to import the plotting helpers eagerly, and those
modules raise at import when matplotlib is missing. That made a PLOTTING
library a prerequisite of running the solver: a nightly build failed on
Windows because matplotlib's _c_internal_utils extension could not load
inside an mpiexec subprocess ("DLL load failed ... The handle is invalid"),
which took `import anuga` down with it -- in a load-balance test that never
plots anything.

Run in a subprocess with matplotlib blocked, because anuga is already
imported by the time this test runs.
"""

import subprocess
import sys
import textwrap

import pytest


def _run(body):
    """Execute body in a fresh interpreter where matplotlib cannot import."""
    script = textwrap.dedent('''
        import builtins
        _real = builtins.__import__
        def _blocked(name, *a, **k):
            if name.startswith('matplotlib'):
                raise ImportError('matplotlib is blocked for this test')
            return _real(name, *a, **k)
        builtins.__import__ = _blocked
    ''') + textwrap.dedent(body)
    return subprocess.run([sys.executable, '-c', script],
                          capture_output=True, text=True)


def test_anuga_imports_without_matplotlib():
    r = _run("import anuga; print('ok')")
    assert r.returncode == 0, \
        'import anuga failed without matplotlib:\n%s' % (r.stderr[-2000:],)
    assert 'ok' in r.stdout


def test_the_solver_runs_without_matplotlib():
    r = _run('''
        import anuga
        d = anuga.rectangular_cross_domain(4, 4)
        b = anuga.Reflective_boundary(d)
        d.set_boundary({t: b for t in d.get_boundary_tags()})
        d.store = False
        for _ in d.evolve(yieldstep=0.5, finaltime=1.0):
            pass
        print('ok')
    ''')
    assert r.returncode == 0, r.stderr[-2000:]
    assert 'ok' in r.stdout


def test_reading_an_sww_does_not_need_matplotlib():
    """plot_utils is the sww reader; only its plotting helpers need matplotlib."""
    r = _run("import anuga.utilities.plot_utils as util; print('ok')")
    assert r.returncode == 0, r.stderr[-2000:]
    assert 'ok' in r.stdout


def test_a_plotting_class_still_says_what_is_wrong():
    """Degrading quietly would be worse than the original hard failure."""
    r = _run('''
        import anuga
        try:
            anuga.SWW_plotter('nonexistent.sww')
        except ImportError as e:
            print('MATPLOTLIB' if 'matplotlib' in str(e) else 'OTHER')
    ''')
    assert r.returncode == 0, r.stderr[-2000:]
    assert 'MATPLOTLIB' in r.stdout, \
        'the error no longer mentions matplotlib: %r' % r.stdout
