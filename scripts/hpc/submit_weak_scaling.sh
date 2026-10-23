#!/bin/bash
# =============================================================================
# Submit ANUGA GPU weak scaling sweep to SLURM
# =============================================================================
#
# Submits one SLURM job per rank count in RANK_COUNTS.
# Each job uses N GPUs and a mesh scaled so elements/GPU ≈ constant.
#
# Usage:
#   bash scripts/hpc/submit_weak_scaling.sh [OPTIONS]
#
# Options:
#   --base-nx N      Base grid size per rank (default 500, ~500K tris/rank)
#   --ranks LIST     Comma-separated rank counts (default: 1,2,4,8,16,32,64)
#   --mode M         GPU mode: 2=GPU offload (default), 1=CPU OpenMP
#   --finaltime T    Simulation end time in seconds (default 30)
#   --yieldstep T    Output interval in seconds (default 10)
#   --outdir DIR     Directory for JSON results (default: benchmarks/results/weak_scaling)
#   --dry-run        Print sbatch commands without submitting
#   --after-all      Submit analysis job after all scaling jobs complete
#
# Example:
#   # Standard SC26 weak scaling sweep (500K tris/GPU, 1→64 GPUs)
#   bash scripts/hpc/submit_weak_scaling.sh --base-nx 500
#
#   # Quick test: small mesh, fewer ranks
#   bash scripts/hpc/submit_weak_scaling.sh --base-nx 100 --ranks 1,2,4
#
#   # Dry run to check sbatch commands
#   bash scripts/hpc/submit_weak_scaling.sh --dry-run
#
# =============================================================================

set -euo pipefail

# Locate repo root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SLURM_SCRIPT="${SCRIPT_DIR}/weak_scaling.slurm"

# Defaults
BASE_NX=500
BASE_LEN=1000
RANK_COUNTS="1,2,4,8,16,32,64"
GPU_MODE=2
FINALTIME=30
YIELDSTEP=10
OUTDIR="${REPO_ROOT}/benchmarks/results/weak_scaling"
DRY_RUN=0
AFTER_ALL=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --base-nx)   BASE_NX="$2";      shift 2 ;;
        --base-len)  BASE_LEN="$2";     shift 2 ;;
        --ranks)     RANK_COUNTS="$2";  shift 2 ;;
        --mode)      GPU_MODE="$2";     shift 2 ;;
        --finaltime) FINALTIME="$2";    shift 2 ;;
        --yieldstep) YIELDSTEP="$2";    shift 2 ;;
        --outdir)    OUTDIR="$2";       shift 2 ;;
        --dry-run)   DRY_RUN=1;         shift ;;
        --after-all) AFTER_ALL=1;       shift ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1 ;;
    esac
done

mkdir -p "${OUTDIR}"
mkdir -p "${REPO_ROOT}/logs"

# Convert comma-separated to array
IFS=',' read -ra RANKS_ARRAY <<< "${RANK_COUNTS}"

echo "================================================"
echo "ANUGA GPU weak scaling sweep"
echo "  base_nx=${BASE_NX}  base_len=${BASE_LEN}m  mode=${GPU_MODE}"
echo "  finaltime=${FINALTIME}s  yieldstep=${YIELDSTEP}s"
echo "  ranks: ${RANK_COUNTS}"
echo "  outdir: ${OUTDIR}"
echo "  dry_run: ${DRY_RUN}"
echo "================================================"
echo ""

JOB_IDS=()

for N in "${RANKS_ARRAY[@]}"; do
    N=$(echo "$N" | tr -d '[:space:]')  # trim whitespace

    JOBNAME="anuga_ws_NP$(printf '%04d' ${N})"

    CMD=(sbatch
        --job-name="${JOBNAME}"
        --ntasks="${N}"
        --export="ALL,NRANKS=${N},BASE_NX=${BASE_NX},BASE_LEN=${BASE_LEN},GPU_MODE=${GPU_MODE},FINALTIME=${FINALTIME},YIELDSTEP=${YIELDSTEP},OUTDIR=${OUTDIR}"
        "${SLURM_SCRIPT}"
    )

    if [[ ${DRY_RUN} -eq 1 ]]; then
        echo "[DRY RUN] ${CMD[*]}"
    else
        JOB_OUT=$("${CMD[@]}" 2>&1)
        JOB_ID=$(echo "${JOB_OUT}" | grep -oP '(?<=Submitted batch job )\d+' || echo 'unknown')
        JOB_IDS+=("${JOB_ID}")
        echo "  Submitted NP=${N}  job_id=${JOB_ID}  name=${JOBNAME}"
    fi
done

echo ""

# Optionally submit analysis job that runs after all scaling jobs complete
if [[ ${AFTER_ALL} -eq 1 && ${DRY_RUN} -eq 0 && ${#JOB_IDS[@]} -gt 0 ]]; then
    DEPEND_STR=$(printf ":%s" "${JOB_IDS[@]}")
    DEPEND_STR="afterok${DEPEND_STR}"
    sbatch \
        --job-name="anuga_ws_analyse" \
        --ntasks=1 \
        --cpus-per-task=1 \
        --mem=4G \
        --time=00:05:00 \
        --dependency="${DEPEND_STR}" \
        --wrap="python ${REPO_ROOT}/benchmarks/run_weak_scaling.py --analyse ${OUTDIR}"
    echo "  Submitted analysis job (runs after all scaling jobs)"
fi

if [[ ${DRY_RUN} -eq 0 ]]; then
    echo ""
    echo "Monitor progress:"
    echo "  squeue -u \$USER --name=anuga_ws"
    echo ""
    echo "Analyse results when done:"
    echo "  python benchmarks/run_weak_scaling.py --analyse ${OUTDIR}"
fi
