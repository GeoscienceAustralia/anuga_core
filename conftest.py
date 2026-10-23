"""Root pytest configuration.

Fast-test filtering
-------------------
By default **all** tests are run.

Pass ``--run-fast`` to skip slow tests and get a quick feedback loop::

    pytest --run-fast

Tests are considered slow if they carry ``@pytest.mark.slow`` **or** live
under ``anuga/parallel/tests/`` (MPI-based tests that spawn subprocesses).

To run *only* slow tests::

    pytest -m slow
"""

import sys
import pathlib
import pytest

# ---------------------------------------------------------------------------
# Ensure the installed package takes precedence over the source tree.
#
# When pytest loads this conftest.py it inserts the repo root into sys.path
# so that conftest.py itself is importable.  That causes `import anuga` to
# resolve to the source tree, which on a fresh clone lacks the meson-
# generated _version.py and raises:
#     ModuleNotFoundError: No module named 'anuga._version'
#
# Moving the repo root to the END of sys.path lets site-packages win while
# still keeping the repo root on the path (needed for conftest resolution).
# Editable installs are unaffected: their import hooks still point at src.
# ---------------------------------------------------------------------------
_repo_root = str(pathlib.Path(__file__).parent.resolve())
if _repo_root in sys.path:
    sys.path.remove(_repo_root)
    sys.path.append(_repo_root)


# ---------------------------------------------------------------------------
# CLI option
# ---------------------------------------------------------------------------

def pytest_addoption(parser):
    parser.addoption(
        '--run-fast',
        action='store_true',
        default=False,
        help='Skip slow tests (parallel + individually marked) for a quick run.',
    )


# ---------------------------------------------------------------------------
# Auto-mark and conditional skip
# ---------------------------------------------------------------------------

def pytest_collection_modifyitems(config, items):
    # Auto-mark every test under anuga/parallel/tests/ as slow.
    for item in items:
        if 'parallel/tests/' in str(item.fspath).replace('\\', '/'):
            item.add_marker(pytest.mark.slow)

    # Only skip slow tests when --run-fast is explicitly requested.
    if not config.getoption('--run-fast'):
        return

    skip_slow = pytest.mark.skip(reason='slow test — omit --run-fast to include')
    for item in items:
        if 'slow' in item.keywords:
            item.add_marker(skip_slow)


# ---------------------------------------------------------------------------
# SWW writer locks
#
# Every SWW file a test creates takes a sidecar '<file>.sww.lock' that is
# normally removed at interpreter exit. A pytest process that dies without
# running atexit handlers (SIGKILL, a timeout kill, an abort in a native
# runtime) leaves every lock taken so far in the session behind, in whatever
# directory the tests wrote to -- usually the cwd. Releasing them after each
# test bounds that to the test that was running.
# ---------------------------------------------------------------------------

@pytest.fixture(autouse=True)
def _release_sww_locks():
    yield
    from anuga.file.sww import _release_all_sww_locks
    _release_all_sww_locks()
