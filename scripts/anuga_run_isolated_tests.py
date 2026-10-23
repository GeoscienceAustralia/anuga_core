#!/usr/bin/env python3
"""Run each collected test in its own fresh pytest process.

Why
---
On a GPU-offload build (nvc, ``-Dgpu_offload=true``) the NVHPC OpenMP-target
runtime aborts — and occasionally hangs — once many mode-2 GPU domains have been
created in a single process (see ``claude/KNOWN_ISSUES.md``). Per-test process
isolation sidesteps the accumulation entirely: every test runs in a fresh
interpreter with a clean CUDA context, so one test can never poison the next.
A per-test timeout also converts a genuine hang into a reported failure instead
of a stuck run.

``pytest --forked`` does **not** work here: CUDA contexts are fork-unsafe, so
forking from a GPU-initialised parent poisons every child. This harness spawns
fresh ``python -m pytest <nodeid>`` subprocesses instead.

On a CPU build the harness still works (just slower than a single in-process
run); it is only *required* on a GPU build.

Usage
-----
    anuga_run_isolated_tests [TARGET ...] [options]      # installed command
    python scripts/anuga_run_isolated_tests.py [TARGET ...] [options]  # source tree

TARGET defaults to the installed ``anuga.shallow_water.tests.test_DE_gpu_omp``
file. It may be any pytest target: a file, a directory, a ``path::Class::test``
node id, or (with ``--pyargs``) a dotted module like ``anuga.shallow_water``.

Options
    --timeout SECONDS   per-test wall-clock limit (default 180; a test exceeding
                        it is recorded as TIMEOUT and its process is killed)
    -k EXPR             pytest -k expression to select tests
    --pyargs            interpret TARGETs as importable modules
    -cm, --compute-mode {legacy,unified}
                        default per-domain compute mode for every child process
                        (sets ANUGA_DEFAULT_COMPUTE_MODE); omit to inherit the
                        current environment
    -x / --exitfirst    stop at the first non-pass
    -v / --verbose      echo each failing test's captured output
    -j N / --jobs N     run N tests concurrently (default 1; keep 1 on a GPU —
                        concurrent tests contend for the same device)

Exit status is 0 only if every test passed or was skipped.

Examples
    anuga_run_isolated_tests
    anuga_run_isolated_tests --pyargs anuga.shallow_water -k riverwall
    anuga_run_isolated_tests --pyargs anuga.shallow_water -cm unified
    anuga_run_isolated_tests --pyargs anuga.shallow_water --compute-mode legacy
"""

import argparse
import concurrent.futures
import os
import re
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _default_target():
    """The GPU test file this harness exists for.

    Resolve it from the *installed* anuga package so the command works after a
    meson/pip install (when this script lives in bindir, far from the tests),
    falling back to the in-repo path for an uninstalled source checkout.
    """
    try:
        import importlib.util
        spec = importlib.util.find_spec(
            "anuga.shallow_water.tests.test_DE_gpu_omp")
        if spec and spec.origin:
            return spec.origin
    except Exception:
        pass
    return str(HERE.parent / "anuga" / "shallow_water" / "tests"
               / "test_DE_gpu_omp.py")


DEFAULT_TARGET = _default_target()

# Disable pytest-isolate in the spawned subprocesses if it is installed: it forks
# (fork-unsafe with CUDA — that is exactly what this harness exists to avoid) and
# its 'Failed to get GPU count using pynvml' warning is just noise here. Harmless
# no-op when the plugin is absent.
try:
    import pytest_isolate  # noqa: F401
    _EXTRA_ARGS = ["-p", "no:isolate"]
except Exception:
    _EXTRA_ARGS = []

# Match pytest's own summary counts ("1 passed", "2 failed", ...), anchored on a
# leading integer so incidental words (e.g. the 'Failed to get GPU count'
# warning) are NOT misread as a result.
_COUNT_RE = re.compile(
    r'(\d+)\s+(passed|failed|errors?|skipped|xfailed|xpassed|deselected)\b')


def _find_rootdir(start):
    """Nearest ancestor of *start* with a project/config marker — the directory
    pytest treats as rootdir and emits node ids relative to."""
    p = Path(start).resolve()
    if not p.is_dir():
        p = p.parent
    for d in (p, *p.parents):
        if any((d / m).exists() for m in
               ("pyproject.toml", "setup.py", "setup.cfg", "tox.ini", ".git")):
            return d
    return p


# Run every child pytest process from a stable rootdir and use absolute paths, so
# the harness works no matter which directory it was launched from (e.g.
# sandpit/). pytest emits node ids relative to the rootdir; running the children
# anywhere else makes them "file or directory not found". Seed the search from the
# current working directory rather than this script's own location: now that the
# script is installed to bindir (far from both the repo and the tests), the cwd is
# the reliable anchor — inside a checkout it walks up to the repo root, and from an
# installed command run anywhere it falls back to the cwd.
ROOTDIR = _find_rootdir(Path.cwd())


def _abs_target(tok):
    """Make a file/dir target absolute (relative to the *current* cwd, before we
    switch the children to ROOTDIR). Node-id suffixes and bare names pass through."""
    path, sep, rest = tok.partition("::")
    p = Path(path)
    if not p.is_absolute() and p.exists():
        return str(p.resolve()) + sep + rest
    return tok


# Directories to resolve relative node ids against, most specific first.
# Populated in main() from the resolved targets; ROOTDIR is the last resort.
#
# WHY THIS EXISTS: pytest emits node ids relative to ITS OWN rootdir, which is
# the common ancestor of the targets -- the test file's directory for a single
# file, the package directory for a --pyargs package. That is NOT the same as
# this script's ROOTDIR (derived from the cwd). Inside a source checkout the two
# coincide, so resolving against ROOTDIR alone worked; for an installed package
# outside a checkout (a container, a wheel install) it does not, and every child
# was handed an unresolvable id like 'test_DE_gpu_omp.py::Test::test_x',
# reported "file or directory not found", and was classified CRASH.
_BASE_DIRS = []


def _register_base_dirs(targets, pyargs):
    """Record where pytest will have rooted its node ids for these targets."""
    dirs = []
    for tok in targets:
        path = tok.partition("::")[0]
        if pyargs:
            # A module/package name: resolve it to its location on disk.
            try:
                import importlib.util
                spec = importlib.util.find_spec(path)
            except Exception:
                spec = None
            if spec is None:
                continue
            if spec.submodule_search_locations:      # package
                dirs += [Path(d) for d in spec.submodule_search_locations]
            elif spec.origin:                        # module
                dirs.append(Path(spec.origin).parent)
        else:
            p = Path(path)
            dirs.append(p if p.is_dir() else p.parent)
    seen, out = set(), []
    for d in dirs + [ROOTDIR]:
        d = d.resolve()
        if d not in seen and d.is_dir():
            seen.add(d)
            out.append(d)
    _BASE_DIRS[:] = out


def _abs_nodeid(nodeid):
    """Resolve a collected node id to an absolute one.

    Node ids already absolute pass through. A relative id is resolved against
    ROOTDIR first (preserving in-checkout behaviour exactly), then against the
    target-derived base directories. A candidate is only accepted when it names
    a real FILE -- a bare directory match would fabricate a nonsense id."""
    path, sep, rest = nodeid.partition("::")
    if not path:
        return nodeid
    p = Path(path)
    if p.is_absolute():
        return str(p) + sep + rest
    for base in (ROOTDIR, *_BASE_DIRS):
        cand = base / path
        if cand.is_file():
            return str(cand) + sep + rest
    return nodeid

# Status -> short label used in the live log and summary
PASS, FAIL, SKIP, ERROR, CRASH, TIMEOUT, NOTESTS = (
    "PASS", "FAIL", "SKIP", "ERROR", "CRASH", "TIMEOUT", "NOTESTS")
BAD = {FAIL, ERROR, CRASH, TIMEOUT}


# Default per-domain compute mode for the child processes ('legacy' | 'unified'),
# or None to inherit the caller's ANUGA_DEFAULT_COMPUTE_MODE. Set from --compute-mode
# in main().
COMPUTE_MODE = None


def _base_env():
    env = dict(os.environ)
    # Opt back in to the GPU tests, which skip themselves in a normal in-process
    # run on a GPU build (see the module-level skip in test_DE_gpu_omp.py).
    env["ANUGA_GPU_TESTS_ISOLATED"] = "1"
    env.setdefault("OMP_NUM_THREADS", "1")
    # --compute-mode overrides the default per-domain compute path for every child;
    # left unset, the caller's existing ANUGA_DEFAULT_COMPUTE_MODE is inherited.
    if COMPUTE_MODE is not None:
        env["ANUGA_DEFAULT_COMPUTE_MODE"] = COMPUTE_MODE
    return env


def collect(targets, pyargs, k_expr):
    """Return the list of test node ids pytest would run for the targets."""
    cmd = [sys.executable, "-m", "pytest", "--collect-only", "-q",
           "-p", "no:cacheprovider", *_EXTRA_ARGS]
    if pyargs:
        cmd.append("--pyargs")
    if k_expr:
        cmd += ["-k", k_expr]
    cmd += targets
    proc = subprocess.run(cmd, capture_output=True, text=True,
                          env=_base_env(), cwd=str(ROOTDIR))
    ids = []
    for line in proc.stdout.splitlines():
        line = line.strip()
        # node-id lines contain '::'; the trailing summary lines ("123 tests
        # collected", "no tests ran", warnings, etc.) do not.
        if "::" in line and not line.startswith(("<", "=", "-")):
            ids.append(line.split(" ")[0])
    if not ids and proc.returncode not in (0, 5):
        sys.stderr.write(proc.stdout + proc.stderr)
    return ids


def classify(returncode, output):
    """Map a single-test pytest run to a status label from its exit code and
    pytest's own summary counts (parsed via _COUNT_RE, so incidental 'failed'
    text in warnings is ignored)."""
    if returncode is None:
        return TIMEOUT
    if returncode < 0:                       # killed by a signal (segfault/abort)
        return CRASH
    counts = {}
    for n, word in _COUNT_RE.findall(output):
        key = "error" if word.startswith("error") else word
        counts[key] = counts.get(key, 0) + int(n)
    if counts.get("failed"):
        return FAIL
    if counts.get("error"):
        return ERROR
    if returncode == 5:
        return NOTESTS
    if counts.get("passed"):
        return PASS
    if counts.get("skipped"):
        return SKIP
    if returncode == 0:
        return PASS
    # Non-zero exit but pytest printed no result counts — the process died
    # before summarising (the NVHPC OpenMP-target abort is exit 1, no summary).
    return CRASH


def run_one(nodeid, timeout):
    cmd = [sys.executable, "-m", "pytest", nodeid, "-p", "no:cacheprovider",
           "-q", "--no-header", "-o", "addopts=", *_EXTRA_ARGS]
    start = time.time()
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              env=_base_env(), timeout=timeout, cwd=str(ROOTDIR))
        rc, out = proc.returncode, proc.stdout + proc.stderr
    except subprocess.TimeoutExpired as e:
        rc, out = None, (e.stdout or "") + (e.stderr or "")
        if isinstance(out, bytes):
            out = out.decode("utf-8", "replace")
    return classify(rc, out), out, time.time() - start


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("targets", nargs="*", default=[DEFAULT_TARGET])
    ap.add_argument("--timeout", type=float, default=180.0)
    ap.add_argument("-k", dest="k_expr", default=None)
    ap.add_argument("--pyargs", action="store_true")
    ap.add_argument("-cm", "--compute-mode", choices=("legacy", "unified"),
                    default=None,
                    help="default per-domain compute mode for every child "
                         "(sets ANUGA_DEFAULT_COMPUTE_MODE); omit to inherit the "
                         "current environment")
    ap.add_argument("-x", "--exitfirst", action="store_true")
    ap.add_argument("-v", "--verbose", action="store_true")
    ap.add_argument("-j", "--jobs", type=int, default=1)
    args = ap.parse_args(argv)
    targets = args.targets or [DEFAULT_TARGET]
    if not args.pyargs:
        targets = [_abs_target(t) for t in targets]

    global COMPUTE_MODE
    COMPUTE_MODE = args.compute_mode
    mode = COMPUTE_MODE or os.environ.get("ANUGA_DEFAULT_COMPUTE_MODE", "legacy")
    print(f"Compute mode: {mode}"
          f"{'' if COMPUTE_MODE else ' (inherited)'}", flush=True)
    print(f"Collecting tests from: {' '.join(targets)}", flush=True)
    _register_base_dirs(targets, args.pyargs)
    ids = [_abs_nodeid(n) for n in collect(targets, args.pyargs, args.k_expr)]
    if not ids:
        print("No tests collected.", file=sys.stderr)
        return 2

    # An id that still is not an existing file cannot be run by any child. That
    # is a harness bug, not a test failure: say so once, instead of spawning N
    # children that each report "file or directory not found" as a CRASH.
    unresolved = [n for n in ids if not Path(n.partition("::")[0]).is_file()]
    if unresolved:
        print(f"ERROR: {len(unresolved)} of {len(ids)} collected node ids could "
              f"not be resolved to a file on disk.\n"
              f"  first: {unresolved[0]}\n"
              f"  searched: {', '.join(str(d) for d in (ROOTDIR, *_BASE_DIRS))}\n"
              f"This is a bug in the runner, not a test failure.", file=sys.stderr)
        return 2
    print(f"Running {len(ids)} tests, one per process "
          f"(timeout {args.timeout:g}s, jobs {args.jobs}).\n", flush=True)

    counts = {s: 0 for s in (PASS, FAIL, SKIP, ERROR, CRASH, TIMEOUT, NOTESTS)}
    bad = []
    width = len(str(len(ids)))

    def _short(nodeid):
        path, sep, rest = nodeid.partition("::")
        try:
            path = os.path.relpath(path, ROOTDIR)
        except ValueError:
            pass
        return path + sep + rest

    def record(i, nodeid, status, out, dt):
        counts[status] += 1
        mark = status if status in BAD else status.lower()
        print(f"[{i:>{width}}/{len(ids)}] {mark:<7} {dt:5.1f}s  {_short(nodeid)}",
              flush=True)
        if status in BAD:
            bad.append((nodeid, status))
            if args.verbose and out.strip():
                print("    " + out.strip().replace("\n", "\n    "), flush=True)

    if args.jobs <= 1:
        for i, nodeid in enumerate(ids, 1):
            status, out, dt = run_one(nodeid, args.timeout)
            record(i, nodeid, status, out, dt)
            if args.exitfirst and status in BAD:
                print("\nStopping at first failure (-x).", flush=True)
                break
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
            futs = {ex.submit(run_one, nid, args.timeout): (i, nid)
                    for i, nid in enumerate(ids, 1)}
            for fut in concurrent.futures.as_completed(futs):
                i, nodeid = futs[fut]
                status, out, dt = fut.result()
                record(i, nodeid, status, out, dt)

    print("\n" + "=" * 60)
    print("  ".join(f"{s.lower()}={counts[s]}" for s in
                    (PASS, SKIP, FAIL, ERROR, CRASH, TIMEOUT, NOTESTS)
                    if counts[s]))
    if bad:
        print(f"\n{len(bad)} test(s) did not pass:")
        for nodeid, status in bad:
            print(f"  {status:<7} {_short(nodeid)}")
        return 1
    print("All tests passed (or were skipped).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
