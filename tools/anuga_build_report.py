#!/usr/bin/env python
"""Report what this ANUGA install actually is, and whether it fits this machine.

Run it after an install, or paste its output into a bug report::

    python tools/anuga_build_report.py

A GPU build names the architectures it was compiled for, and contains *only*
those.  Asking for the wrong one is not a build error -- it compiles cleanly and
then fails at every kernel launch, which surfaces as the GPU tests reporting
CRASH.  ``--check`` exits non-zero on that mismatch (and on a broken install) so
an install script can stop and say why:

    python tools/anuga_build_report.py --check
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys


def _run(cmd, timeout=10):
    """Return stdout of *cmd*, or '' if it cannot be run."""
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return out.stdout.strip() if out.returncode == 0 else ''
    except (OSError, subprocess.SubprocessError):
        return ''


def gpu_compute_caps():
    """Compute capabilities of the GPUs present, e.g. ['8.6']."""
    out = _run(['nvidia-smi', '--query-gpu=compute_cap', '--format=csv,noheader'])
    return [line.strip() for line in out.splitlines() if line.strip()]


def gpu_names():
    out = _run(['nvidia-smi', '--query-gpu=name', '--format=csv,noheader'])
    return [line.strip() for line in out.splitlines() if line.strip()]


def _find_cuobjdump():
    exe = shutil.which('cuobjdump')
    if exe:
        return exe
    pattern = '/opt/nvidia/hpc_sdk/Linux_x86_64/*/cuda/bin/cuobjdump'
    found = sorted(glob.glob(pattern))
    return found[-1] if found else None


def built_sm_architectures(so_path):
    """SM architectures embedded in *so_path*, e.g. ['sm_80', 'sm_86'].

    Returns None when it cannot be determined (no cuobjdump available), which
    is different from an empty list (a CPU build with no device code).
    """
    cuobjdump = _find_cuobjdump()
    if not cuobjdump or not so_path or not os.path.isfile(so_path):
        return None
    out = _run([cuobjdump, '--list-elf', so_path], timeout=60)
    if not out:
        return []
    archs = {tok for line in out.splitlines() for tok in line.split()
             if tok.startswith('sm_')}
    # also catch "sm_80" appearing inside a filename like foo.sm_80.cubin
    for line in out.splitlines():
        for part in line.replace('.', ' ').split():
            if part.startswith('sm_') and part[3:].isdigit():
                archs.add(part)
    return sorted(archs, key=lambda a: int(a[3:]))


def collect():
    """Gather everything worth knowing; never raises."""
    info = {}
    info['python'] = '%d.%d.%d' % sys.version_info[:3]
    info['executable'] = sys.executable
    info['platform'] = _run(['uname', '-sr']) or sys.platform
    info['conda_env'] = os.environ.get('CONDA_DEFAULT_ENV', '(none)')
    info['compute_mode_env'] = os.environ.get('ANUGA_DEFAULT_COMPUTE_MODE',
                                              '(unset -> legacy)')

    try:
        import anuga
    except Exception as exc:                      # noqa: BLE001 - report anything
        info['anuga'] = 'IMPORT FAILED: %s: %s' % (type(exc).__name__, exc)
        return info

    info['anuga'] = getattr(anuga, '__version__', '(unknown)')
    info['anuga_path'] = getattr(anuga, '__file__', '(unknown)')
    # An editable install imports from the source tree and keeps its compiled
    # extensions in build/cp<ver>; deleting that directory breaks the install.
    info['editable'] = 'site-packages' not in (info['anuga_path'] or '')

    try:
        info['gpu_offload_build'] = bool(anuga.gpu_offload_enabled())
    except Exception:
        info['gpu_offload_build'] = None

    ext_path = None
    try:
        import anuga.shallow_water.sw_domain_gpu_ext as ext
        ext_path = getattr(ext, '__file__', None)
        info['ext_path'] = ext_path
        try:
            info['mpi'] = bool(ext.gpu_has_mpi())
        except Exception:
            info['mpi'] = None
    except Exception as exc:                      # noqa: BLE001
        info['ext_path'] = 'IMPORT FAILED: %s' % exc
        info['mpi'] = None

    info['built_for'] = built_sm_architectures(ext_path)
    info['gpu_caps'] = gpu_compute_caps()
    info['gpu_names'] = gpu_names()
    return info


def mismatch(info):
    """(problem, advice) if the build cannot run on this machine, else None."""
    if str(info.get('anuga', '')).startswith('IMPORT FAILED'):
        advice = 'Reinstall: the package does not import.'
        path = info.get('anuga_path') or ''
        if 'build/cp' in str(info.get('anuga', '')) or 'No such file' in str(info.get('anuga', '')):
            advice = ('An editable install lost its build directory. '
                      'Reinstall with pip install --no-build-isolation -e . '
                      '(plus your -Dgpu_* options).')
        return ('ANUGA does not import', advice)

    built = info.get('built_for')
    caps = info.get('gpu_caps') or []
    # Only meaningful for a GPU build on a machine that has a GPU and where the
    # architectures could actually be read.
    if not built or not caps or not info.get('gpu_offload_build'):
        return None

    # Compatibility is one-directional.  nvc embeds PTX alongside the SASS, and
    # the driver JIT-compiles PTX *forward*, so a build for an older
    # architecture runs on a newer GPU (verified: a cc86 build runs on sm_120).
    # The reverse is impossible -- PTX for compute_120 cannot target sm_86 --
    # so the build is unusable only when EVERY architecture in it is newer than
    # the GPU present.
    built_nums = sorted(int(''.join(c for c in a[3:] if c.isdigit()))
                        for a in built)
    cap_nums = sorted(int(c.split('.')[0]) * 10 + int(c.split('.')[1])
                      for c in caps if '.' in c)
    if not built_nums or not cap_nums:
        return None

    newest_gpu = max(cap_nums)
    if min(built_nums) > newest_gpu:
        return ('this build targets only GPUs newer than the one present '
                '(built for %s, GPU is sm_%d)'
                % (' '.join(built), newest_gpu),
                'Rebuild for this GPU: GPU_ARCH=cc%d bash '
                'tools/install_anuga_nvc.sh' % newest_gpu)
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--check', action='store_true',
                        help='exit non-zero if this build cannot run here')
    args = parser.parse_args()

    info = collect()

    def row(label, value):
        print('#   %-22s %s' % (label + ':', value))

    print('#' + '=' * 62)
    print('# ANUGA build report')
    print('#' + '=' * 62)
    row('anuga', info.get('anuga'))
    row('imported from', info.get('anuga_path'))
    if info.get('editable'):
        row('install type', 'editable (needs its build/cp* directory kept)')
    else:
        row('install type', 'copied into site-packages')
    row('python', '%s  (%s)' % (info.get('python'), info.get('conda_env')))
    row('platform', info.get('platform'))
    row('compute mode default', info.get('compute_mode_env'))

    gpu_build = info.get('gpu_offload_build')
    row('GPU offload build', {True: 'yes', False: 'no (CPU multicore)',
                              None: 'unknown'}[gpu_build])
    if info.get('mpi') is not None:
        row('MPI in GPU ext', 'yes' if info['mpi'] else 'no (sequential stubs)')

    built = info.get('built_for')
    if built is None:
        row('built for', '(cuobjdump not available - cannot tell)')
    elif built:
        row('built for', ' '.join(built))

    names, caps = info.get('gpu_names') or [], info.get('gpu_caps') or []
    if caps:
        for name, cap in zip(names or ['GPU'] * len(caps), caps):
            row('GPU present', '%s (compute capability %s -> sm_%s)'
                % (name, cap, cap.replace('.', '')))
    else:
        row('GPU present', 'none detected (nvidia-smi unavailable)')
    print('#' + '=' * 62)

    problem = mismatch(info)
    if problem:
        what, advice = problem
        print('#')
        print('# PROBLEM: %s' % what)
        print('#')
        print('# %s' % advice)
        print('#' + '=' * 62)
        if args.check:
            return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
