#!/usr/bin/env bash
# Run test_DE_gpu_omp.py with one fresh process per test class.
#
# Why: on a GPU-offload build (nvc, -Dgpu_offload=true) the NVIDIA OpenMP-target
# runtime cannot churn many mode-2 GPU domains in a single process — after a
# handful of domains the device present table / CUDA context is exhausted and
# the process aborts (silent exit, see claude/KNOWN_ISSUES.md). Each test CLASS
# passes in its own fresh process, so this script runs them separately and
# aggregates the result. pytest's --forked does NOT work here: CUDA contexts are
# fork-unsafe, so forking from a GPU-touched parent poisons the children.
#
# On a CPU build (-Dgpu_offload=false) the whole file runs fine in one process;
# this script is only needed for GPU-offload builds, but is harmless on CPU too.
#
# Usage:
#   bash anuga/shallow_water/tests/run_gpu_tests_isolated.sh
set -u

TESTFILE="$(dirname "$0")/test_DE_gpu_omp.py"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
# Opt back in to the GPU tests: the file skips itself in a normal in-process run
# on a GPU-offload build (see the module-level skip in test_DE_gpu_omp.py).
export ANUGA_GPU_TESTS_ISOLATED=1

mapfile -t CLASSES < <(python -m pytest "$TESTFILE" --collect-only -q -p no:cacheprovider 2>/dev/null \
  | sed -n 's/\(.*::[A-Za-z_0-9]*\)::test.*/\1/p' | sort -u)

if [ "${#CLASSES[@]}" -eq 0 ]; then
  echo "No test classes collected (is the GPU extension importable?)." >&2
  exit 2
fi

fail=0
for cls in "${CLASSES[@]}"; do
  # NB: `rc=$?` after a pipeline reads the LAST command's status (tail, always
  # 0), which silently reported every failing class as ok. Capture the full
  # output first, then take pytest's own status.
  out=$(python -m pytest "$cls" -p no:cacheprovider -q 2>&1)
  rc=$?
  line=$(printf '%s\n' "$out" | tail -1)
  if [ "$rc" -ne 0 ]; then
    echo "FAIL  $cls  -> $line"
    printf '%s\n' "$out" | tail -40
    fail=1
  else
    echo "ok    $cls  -> $line"
  fi
done

if [ "$fail" -ne 0 ]; then
  echo "=== Some GPU test classes failed ==="
  exit 1
fi
echo "=== All GPU test classes passed ==="
