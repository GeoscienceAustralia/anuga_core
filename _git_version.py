"""
Helper script called by meson.build to produce a PEP 440-compliant version
string from `git describe` output.

Sources, in order:
  1. ANUGA_VERSION in the environment (explicit override);
  2. _version_static.txt next to this file, written into the source
     distribution by `meson dist` (see the dist script below), so a build
     from an unpacked sdist, which has no .git, still knows its version;
  3. `git describe`;
  4. the 0.0.0+unknown fallback.

Run with --write-dist from a dist script: writes the version into
$MESON_DIST_ROOT/_version_static.txt.

git describe format:  TAG[-N-gSHA][-dirty]
PEP 440 mapping:
  3.2.0                          -> 3.2.0
  3.2.0-4-gabcdef                -> 3.2.0.dev4+gabcdef
  3.2.0-4-gabcdef-dirty          -> 3.2.0.dev4+gabcdef.dirty
  3.2.0-dirty                    -> 3.2.0+dirty
  abcdef  (no tag reachable)     -> 0.0.0+gabcdef
"""

import os
import re
import subprocess
import sys


def git_version():
    # Allow an explicit override (PEP 440 string) for builds where git metadata
    # is unavailable — e.g. an sdist, or a Docker build whose context excludes
    # .git. Set ANUGA_VERSION to what `git describe` would have produced.
    override = os.environ.get('ANUGA_VERSION', '').strip()
    if override:
        return override

    static = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_version_static.txt')
    if os.path.isfile(static):
        with open(static) as f:
            value = f.read().strip()
        if value:
            return value

    # Ask git about the tree this script lives in, whatever the current
    # directory is (meson runs dist scripts from the build directory).
    try:
        result = subprocess.run(
            ['git', 'describe', '--tags', '--dirty', '--always'],
            capture_output=True, text=True,
            cwd=os.path.dirname(os.path.abspath(__file__))
        )
    except OSError:          # git not installed at all
        return '0.0.0+unknown'
    if result.returncode != 0:
        return '0.0.0+unknown'

    raw = result.stdout.strip()

    # Full describe: TAG-N-gSHA[-dirty]
    m = re.match(
        r'^(\d+\.\d+(?:\.\d+)?(?:[a-z]\w*)?)'  # base tag
        r'(?:-(\d+)-g([0-9a-f]+))?'             # optional distance + sha
        r'(-dirty)?$',                           # optional dirty flag
        raw
    )
    if not m:
        # No reachable tag — just a bare SHA (possibly dirty)
        sha = re.sub(r'-dirty$', '', raw)
        dirty = raw.endswith('-dirty')
        return f'0.0.0+g{sha}' + ('.dirty' if dirty else '')

    base, distance, sha, dirty = m.groups()
    if distance:
        version = f'{base}.dev{distance}+g{sha}'
        if dirty:
            version += '.dirty'
    elif dirty:
        version = f'{base}+dirty'
    else:
        version = base

    return version


def write_dist_version():
    """Dist-script entry point: bake the version into the staged sdist tree."""
    dist_root = os.environ.get('MESON_DIST_ROOT')
    if not dist_root:
        sys.exit('--write-dist must be run by meson dist (MESON_DIST_ROOT is unset)')
    version = git_version()
    with open(os.path.join(dist_root, '_version_static.txt'), 'w') as f:
        f.write(version + '\n')
    print(f'_git_version.py: wrote {version} into the sdist')


if __name__ == '__main__':
    if '--write-dist' in sys.argv[1:]:
        write_dist_version()
    else:
        print(git_version())
