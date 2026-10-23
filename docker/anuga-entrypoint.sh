#!/usr/bin/env bash
#
# ANUGA container entrypoint — headless/batch friendly with optional S3 I/O.
#
# Behaviour:
#   * If INPUT_S3 is set, its contents are synced INTO the workdir ($ANUGA_WORKDIR,
#     default /work) in place — so your run script and its data land right where
#     the command runs (e.g. `python run.py` finds run.py at /work/run.py).
#   * The command passed to the container ("$@") is run from $ANUGA_WORKDIR.
#   * If OUTPUT_S3 is set, results are synced up to it after the run:
#       - by default the whole workdir (captures output wherever the script
#         writes it, e.g. MODEL_OUTPUTS/), which also re-sends the inputs;
#       - set ANUGA_OUTPUT_DIR to a subdir (e.g. /work/MODEL_OUTPUTS) to upload
#         only that and avoid re-sending large inputs.
#     Uses `aws s3 sync`, so re-runs only upload changed/new files.
#   * For AWS Batch array jobs, AWS_BATCH_JOB_ARRAY_INDEX is appended to
#     OUTPUT_S3 so each index lands in its own prefix (.../0/, .../1/, ...).
#   * With no command, drops to an interactive bash shell (local debugging).
#
# Credentials for aws come from the environment / instance role (e.g. the AWS
# Batch job IAM role, or an EC2 instance profile). None are baked into the image.
set -euo pipefail

WORKDIR="${ANUGA_WORKDIR:-/work}"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

log() { echo "[anuga-entrypoint] $*" >&2; }

if [ -n "${INPUT_S3:-}" ]; then
    log "Staging input: ${INPUT_S3} -> ${WORKDIR}"
    aws s3 sync "${INPUT_S3}" "${WORKDIR}"
fi

OUT="${OUTPUT_S3:-}"
if [ -n "$OUT" ] && [ -n "${AWS_BATCH_JOB_ARRAY_INDEX:-}" ]; then
    OUT="${OUT%/}/${AWS_BATCH_JOB_ARRAY_INDEX}"
fi

stage_out() {
    local rc=$1
    [ -z "$OUT" ] && return 0
    local src="${ANUGA_OUTPUT_DIR:-$WORKDIR}"
    if [ "$rc" -eq 0 ] || [ "${STAGE_OUT_ON_FAILURE:-0}" = "1" ]; then
        log "Uploading results: ${src} -> ${OUT}"
        aws s3 sync "${src}" "${OUT}" || log "WARNING: upload to ${OUT} failed"
    else
        log "Command failed (rc=$rc); skipping upload (set STAGE_OUT_ON_FAILURE=1 to override)"
    fi
}

# No command -> interactive shell for debugging on a local/EC2 box.
if [ "$#" -eq 0 ]; then
    log "No command given; starting interactive shell in ${WORKDIR}"
    exec /bin/bash
fi

log "Running: $*"
set +e
"$@"
rc=$?
set -e
stage_out "$rc"
exit "$rc"
