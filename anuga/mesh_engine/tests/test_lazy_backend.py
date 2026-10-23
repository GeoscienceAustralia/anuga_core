"""The triangulation backend (meshpy / triangle) is imported on first use,
not by `import anuga` (issue #238).

Both checks run in a subprocess with the backends blocked, because the
in-process test session has long since imported meshpy through other tests.
"""
import subprocess
import sys

import pytest

_BLOCK = ("import sys; "
          "sys.modules['meshpy'] = None; sys.modules['meshpy.triangle'] = None; "
          "sys.modules['triangle'] = None; ")


def _run(code):
    return subprocess.run([sys.executable, '-c', _BLOCK + code],
                          capture_output=True, text=True, timeout=300)


def test_import_anuga_needs_no_triangulator():
    r = _run("import anuga; from anuga.mesh_engine import mesh_engine; "
             "print(mesh_engine.TRILIB)")
    assert r.returncode == 0, r.stderr
    # Only the last line is ours: `import anuga` can print first (mpi4py's
    # UCX layer reports unusable transports on some CI runners, to stdout).
    assert r.stdout.strip().splitlines()[-1] == 'None'   # nothing loaded until first use


def test_generate_mesh_without_a_backend_says_what_to_install():
    r = _run("from anuga.mesh_engine.mesh_engine import generate_mesh; "
             "generate_mesh(points=[[0, 0], [1, 0], [0, 1]], segments=[[0, 1], [1, 2], [2, 0]])")
    assert r.returncode != 0
    assert 'pip install meshpy' in r.stderr
    assert 'only required to create meshes' in r.stderr


def test_backend_is_loaded_on_first_use():
    pytest.importorskip('meshpy')
    from anuga.mesh_engine import mesh_engine
    mesh_engine.generate_mesh(points=[[0, 0], [1, 0], [0, 1]],
                              segments=[[0, 1], [1, 2], [2, 0]], mode='pzq')
    assert mesh_engine.TRILIB in ('meshpy', 'triangle')
