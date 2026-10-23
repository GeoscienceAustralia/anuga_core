"""_git_version.py, the build-time version helper (issue #247).

An unpacked sdist has no .git, so `meson dist` writes the version into
_version_static.txt next to the script and the script reads that before
trying git. The script lives at the repository root, outside the package, so
these tests copy it into a temporary directory; they skip when run from an
installed package where the source tree is not available.
"""
import os
import shutil
import subprocess
import sys

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPT = os.path.normpath(os.path.join(_HERE, '..', '..', '..', '_git_version.py'))


@pytest.fixture
def script_dir(tmp_path):
    if not os.path.isfile(_SCRIPT):
        pytest.skip('_git_version.py is only present in a source checkout')
    shutil.copy(_SCRIPT, tmp_path / '_git_version.py')
    return tmp_path


def _run(cwd, env=None, *args):
    e = {k: v for k, v in os.environ.items() if k != 'ANUGA_VERSION'}
    if env:
        e.update(env)
    # cwd is a fresh temp dir: no .git, so git describe fails there
    return subprocess.run([sys.executable, '_git_version.py', *args],
                          cwd=cwd, env=e, capture_output=True, text=True)


def test_static_file_is_used_when_there_is_no_git(script_dir):
    (script_dir / '_version_static.txt').write_text('4.0.1\n')
    r = _run(script_dir)
    assert r.stdout.strip() == '4.0.1'


def test_environment_override_wins_over_the_static_file(script_dir):
    (script_dir / '_version_static.txt').write_text('4.0.1\n')
    r = _run(script_dir, {'ANUGA_VERSION': '9.9.9'})
    assert r.stdout.strip() == '9.9.9'


def test_fallback_without_git_or_static_file(script_dir):
    r = _run(script_dir, {'GIT_DIR': str(script_dir / 'nowhere')})
    assert r.stdout.strip() == '0.0.0+unknown'


def test_write_dist_bakes_the_version_into_the_dist_root(script_dir, tmp_path):
    dist_root = tmp_path / 'dist'
    dist_root.mkdir()
    r = _run(script_dir, {'ANUGA_VERSION': '4.0.1', 'MESON_DIST_ROOT': str(dist_root)},
             '--write-dist')
    assert r.returncode == 0, r.stderr
    assert (dist_root / '_version_static.txt').read_text().strip() == '4.0.1'


def test_write_dist_refuses_outside_meson_dist(script_dir):
    r = _run(script_dir, {}, '--write-dist')
    assert r.returncode != 0
    assert 'MESON_DIST_ROOT' in r.stderr
