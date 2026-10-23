"""
run_auto_validation_tests.py

Runs two groups of validation tests:

1. validate_*.py  — simulation + accuracy checks (analytical / experimental)
   Run directly as  python validate_*.py -alg DE1 -np 1

2. test_regression_*.py  — pytest-regressions baseline checks (behaviour_only)
   Run via  pytest test_regression_*.py
   On first run use --force-regen to generate the baseline .csv files.
"""

import collections
import os
import sys
import time
import subprocess
import argparse

import anuga
anuga_args = anuga.get_args()


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

dirs_to_skip = ['.']       # avoid infinite recursion
dirs_to_skip += ['patong'] # requires downloaded data, takes many hours

# Directories whose SUBTREES are skipped entirely. A finished run leaves a
# copy of the case's own scripts beside its results -- towradgi snapshots
# code/ into OUTPUT/RUN_<timestamp>/code/ -- so walking into them collects
# the same validate_*.py once per past run: 86 copies against 47 real cases
# on a machine that had run it a few times. OUTPUT is gitignored, so these
# are purely local artefacts and never appear in CI, which is why this went
# unnoticed.
prune_dirs = {'OUTPUT', '__pycache__', '.git', '.svn', '.pytest_cache',
              '.ipynb_checkpoints', 'build', 'dist'}

# Long-running cases, skipped by default to keep the routine run fast; pass
# -l/--long to include them. The HEC-RAS bridge/weir behaviour cases take
# >100 s each; towradgi_simulation builds the catchment mesh from the DEM,
# where the plain towradgi case reuses the cached towradgi.tsh and comes in
# under a minute; van_rijn_trench runs three 15 m flume experiments to
# steady state.
long_dirs = ['bridge_hecras', 'bridge_hecras2', 'lateral_weir_hecras',
             'towradgi_simulation', 'van_rijn_trench']
if not anuga_args.long:
    dirs_to_skip += long_dirs

# (dirpath, filename, runner)  where runner is 'python' or 'pytest'
all_tests = []

for dirpath, dirnames, filenames in os.walk('.'):

    # Prune in place, so os.walk does not descend. The old code skipped a
    # directory's own files with `continue` but still walked its children,
    # so a skipped case's subdirectories were searched anyway.
    dirnames[:] = [d for d in dirnames
                   if d not in prune_dirs and d not in dirs_to_skip]

    dirname = os.path.split(dirpath)[-1]
    if dirname in dirs_to_skip:
        continue

    for filename in sorted(filenames):
        if filename.startswith('validate_') and filename.endswith('.py'):
            all_tests.append((dirpath, filename, 'python'))
        elif filename.startswith('test_regression_') and filename.endswith('.py'):
            all_tests.append((dirpath, filename, 'pytest'))

# De-duplicate the behaviour_only cases: a test_regression_*.py already runs the
# simulation (via ensure_sww) and checks it, so the sibling validate_*.py in the
# same directory would run the same sim a second time. Drop those validate_
# scripts when a regression test covers the directory.
regression_dirs = {p for p, f, r in all_tests if r == 'pytest'}
all_tests = [t for t in all_tests
             if not (t[2] == 'python' and t[0] in regression_dirs)]

all_tests.sort()

# Separate the two groups for display
validate_tests   = [(p, f) for p, f, r in all_tests if r == 'python']
regression_tests = [(p, f) for p, f, r in all_tests if r == 'pytest']


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------

print()
print(80*'=')
print('Running all validation tests - some may take many minutes')
print('and some may require memory in the order of 8-16GB       ')
if not anuga_args.long:
    print('(skipping long cases: %s; pass -l/--long to include)'
          % ', '.join(long_dirs))
print(80*'=')

print('Simulation + accuracy tests (validate_*.py):')
for path, filename in validate_tests:
    print('    ', os.path.join(path, filename))

if regression_tests:
    print()
    print('Regression baseline tests (test_regression_*.py):')
    for path, filename in regression_tests:
        print('    ', os.path.join(path, filename))

print()

t0 = time.time()
parentdir = os.getcwd()

# Results: list of (label, returncode, elapsed, runner)
results = []


# ---------------------------------------------------------------------------
# Run validate_*.py tests
# ---------------------------------------------------------------------------

for path, filename in validate_tests:
    label = os.path.join(path, filename)
    os.chdir(path)

    cmd = [sys.executable, filename, '-alg', str(anuga_args.alg), '-np', str(anuga_args.np)]
    if anuga_args.verbose:
        cmd.append('-v')

    print()
    print(80*'=')
    print(' '.join(cmd))
    print(80*'=')
    sys.stdout.flush()

    t_start = time.time()
    proc = subprocess.run(cmd)
    elapsed = time.time() - t_start

    results.append((label, proc.returncode, elapsed, 'python'))
    status = 'PASSED' if proc.returncode == 0 else f'FAILED (exit code {proc.returncode})'
    print(f'  --> {status}  ({elapsed:.1f} s)')

    os.chdir(parentdir)


# ---------------------------------------------------------------------------
# Run test_regression_*.py tests via pytest
# ---------------------------------------------------------------------------

for path, filename in regression_tests:
    label = os.path.join(path, filename)
    os.chdir(path)

    cmd = ['pytest', filename, '-v']
    if anuga_args.verbose:
        cmd.append('-s')

    print()
    print(80*'=')
    print(' '.join(cmd))
    print(80*'=')
    sys.stdout.flush()

    t_start = time.time()
    proc = subprocess.run(cmd)
    elapsed = time.time() - t_start

    results.append((label, proc.returncode, elapsed, 'pytest'))
    status = 'PASSED' if proc.returncode == 0 else f'FAILED (exit code {proc.returncode})'
    print(f'  --> {status}  ({elapsed:.1f} s)')

    os.chdir(parentdir)


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def format_elapsed(seconds):
    """Seconds, plus a minutes/hours reading once the number stops being
    readable on its own."""
    if seconds < 60.0:
        return '%.1f seconds' % seconds
    m, s = divmod(int(round(seconds)), 60)
    if m < 60:
        return '%.1f seconds (%d m %02d s)' % (seconds, m, s)
    h, m = divmod(m, 60)
    return '%.1f seconds (%d h %02d m %02d s)' % (seconds, h, m, s)


total_elapsed = time.time() - t0
n_passed = sum(1 for _, rc, _, _ in results if rc == 0)
n_failed = len(results) - n_passed

print()
print('VALIDATION SUMMARY')
print()

# The case column names the directory, not the whole path: a case's script
# is nearly always validate_<directory>.py, so repeating it costs 44
# characters a row and says nothing. A directory that holds more than one
# test keeps its filenames, since there the name is what tells them apart.
n_per_dir = collections.Counter(os.path.dirname(label) for label, _, _, _ in results)

rows = []
for label, rc, elapsed, runner in results:
    case = os.path.relpath(os.path.dirname(label), '.')
    if n_per_dir[os.path.dirname(label)] > 1:
        case = os.path.join(case, os.path.basename(label))
    rows.append((case, runner,
                 'PASSED' if rc == 0 else f'FAILED ({rc})',
                 f'{elapsed:.1f} s'))

header = ('Case', 'Runner', 'Result', 'Time')
width = [max(len(r[i]) for r in rows + [header]) for i in range(4)]
rule = '+' + '+'.join('-' * (w + 2) for w in width) + '+'


def print_row(cells, pad=' '):
    print('|' + '|'.join(
        pad + (c.rjust(w) if i == 3 else c.ljust(w)) + pad
        for i, (c, w) in enumerate(zip(cells, width))) + '|')


print(rule)
print_row(header)
print(rule)
group = None
for row in rows:
    # A rule wherever the top-level group changes, so 40-odd rows do not read
    # as one block. The rows are sorted, so each group is contiguous.
    if group is not None and row[0].split(os.sep)[0] != group:
        print(rule)
    group = row[0].split(os.sep)[0]
    print_row(row)
print(rule)

print()
if n_failed == 0:
    print(f'All {n_passed} validation tests PASSED in {format_elapsed(total_elapsed)}.')
else:
    print(f'{n_passed} PASSED,  {n_failed} FAILED, in {format_elapsed(total_elapsed)}.')
    print()
    print('Failed tests:')
    for label, rc, _, runner in results:
        if rc != 0:
            print(f'  * [{runner}] {label}  (exit code {rc})')

print()
if regression_tests:
    print('To regenerate regression baselines after an intentional change:')
    print('  cd validation_tests/behaviour_only')
    print('  pytest test_regression_*.py --force-regen')

print()

sys.exit(0 if n_failed == 0 else 1)
