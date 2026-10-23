#!/usr/bin/env bash
# Static analysis of ANUGA's hand-written C sources with gcc -fanalyzer
# (audit issue #342). Run from the repository root:
#
#     bash scripts/c_static_analysis.sh            # analyse, fail on any finding
#     bash scripts/c_static_analysis.sh --list     # just list the files it would analyse
#     CC=gcc-13 bash scripts/c_static_analysis.sh  # a specific gcc (CI pins one; see lint.yml)
#
# Different gcc versions report different things, so a finding that only one
# version raises is normal: fix it if it is real, or add it to EXCLUDE below
# with a note if it is not.
#
# Every *.c under anuga/ that is not a Cython-generated *_ext.c is compiled to
# /dev/null with -fanalyzer and the include paths meson uses. Findings inside
# third-party headers (uthash.h) are filtered out. Anything else fails the
# script, so a new unchecked malloc, leak on an error path, or NULL
# dereference is caught before review. Needs gcc >= 12 and numpy importable
# from the python on PATH.
set -uo pipefail
cd "$(dirname "$0")/.."

NUMPY_INC=$(python -c "import numpy; print(numpy.get_include())") || exit 2
FILES=$(find anuga -name '*.c' -not -name '*_ext.c' | sort)

if [ "${1:-}" = "--list" ]; then echo "$FILES"; exit 0; fi

# -Winfinite-recursion: a function that unconditionally calls itself. Not an
# analyser check, but it is the one that would have caught gpu_set_error
# echoing through itself, so it rides along.
CFLAGS=(-fanalyzer -Winfinite-recursion -fopenmp -DCPU_ONLY_MODE
        -I anuga/utilities -I anuga/shallow_water -I anuga/shallow_water/gpu -I "$NUMPY_INC")
# Third-party headers whose findings are not ours to fix, plus one known
# false positive: the analyser reports the rate-operator cache malloc in
# gpu_rate_operator.c as leaked, but it is stored in the operator slot that
# outlives the call and freed in gpu_rate_operator_finalize.
EXCLUDE='uthash\.h|gpu_rate_operator\.c:[0-9]+:[0-9]+: warning: leak of .malloc'

status=0
for f in $FILES; do
    out=$("${CC:-gcc}" "${CFLAGS[@]}" -c "$f" -o /dev/null 2>&1)
    rc=$?
    findings=$(printf '%s\n' "$out" | grep -E 'warning:|error:' | grep -Ev "$EXCLUDE" || true)
    if [ $rc -ne 0 ] || [ -n "$findings" ]; then
        echo "== $f"
        if [ $rc -ne 0 ]; then printf '%s\n' "$out" | grep -E 'error' | head -5; fi
        [ -n "$findings" ] && printf '%s\n' "$findings"
        status=1
    fi
done

if [ $status -eq 0 ]; then
    echo "$("${CC:-gcc}" --version | head -1): -fanalyzer found nothing in $(echo "$FILES" | wc -l) files"
fi
exit $status
