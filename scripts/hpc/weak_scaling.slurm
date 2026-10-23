#!/bin/bash
# =============================================================================
# ANUGA GPU Weak Scaling SLURM Job Template
# =============================================================================
#
# Submit via submit_weak_scaling.sh, or directly:
#
#   sbatch --ntasks=4 --export=ALL,NRANKS=4 scripts/hpc/weak_scaling.slurm
#
# Required environment variables (set by submit_weak_scaling.sh or --export):
#
#   NRANKS      Number of MPI ranks (= number of GPUs)
#   BASE_NX     Base grid size per rank (default 500, ~500K tris/rank)
#   BASE_LEN    Physical domain length per rank in metres (default 1000).
#               Scaled by sqrt(NRANKS) alongside BASE_NX so cell size =
#               BASE_LEN/BASE_NX is constant and the CFL timestep is
#               invariant across all rank counts.
#   GPU_MODE    multiprocessor_mode: 2=GPU (default), 1=CPU OpenMP
#   FINALTIME   Simulation end time in seconds (default 30)
#   YIELDSTEP   Output interval in seconds (default 10)
#   OUTDIR      Directory for JSON result files
#               (default: $REPO/benchmarks/results/weak_scaling)
#
# =============================================================================
# SLURM directives — adjust for your cluster
# =============================================================================

# Job name includes rank count (set by submit_weak_scaling.sh via --job-name)
#SBATCH --job-name=anuga_ws_NP${NRANKS:-1}
#SBATCH --output=logs/weak_scaling_NP%04a_%j.out
#SBATCH --error=logs/weak_scaling_NP%04a_%j.err

# Resources: one GPU per MPI rank
#SBATCH --ntasks=${NRANKS:-1}
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=4          # cores for host-side Python + OpenMP
#SBATCH --mem-per-gpu=32G          # GPU device VRAM is separate; this is host RAM

# Time limit — scale with NRANKS (weak scaling: constant per-rank work)
# For base_nx=500 (~500K tris/rank, finaltime=30s): ~5 min wall per config
#SBATCH --time=00:30:00

# Partition / QOS — edit for your cluster
##SBATCH --partition=gpu
##SBATCH --qos=gpu

# =============================================================================
# Environment setup
# =============================================================================

# Fail on first error
set -euo pipefail

# Default parameter values (override via --export)
NRANKS=${NRANKS:-1}
BASE_NX=${BASE_NX:-500}
BASE_LEN=${BASE_LEN:-1000}
GPU_MODE=${GPU_MODE:-2}
FINALTIME=${FINALTIME:-30}
YIELDSTEP=${YIELDSTEP:-10}

# Locate repository root (script lives in scripts/hpc/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

OUTDIR=${OUTDIR:-"${REPO_ROOT}/benchmarks/results/weak_scaling"}
mkdir -p "${OUTDIR}"
mkdir -p logs

# Load modules — edit for your HPC environment
# Example for systems using Lmod:
##module purge
##module load gcc/12.3   cuda/12.3   openmpi/4.1-cuda   python/3.11
##module load anaconda3
##conda activate anuga_env

# Or, if conda is already active:
echo "Python: $(which python)"
echo "ANUGA:  $(python -c 'import anuga; print(anuga.__version__)' 2>/dev/null || echo 'not found')"
echo "mpi4py: $(python -c 'import mpi4py; print(mpi4py.__version__)' 2>/dev/null || echo 'not found')"

# =============================================================================
# OpenMP / offload settings
# =============================================================================

# One host thread per rank (GPU does the parallelism)
export OMP_NUM_THREADS=1

# Abort if GPU offload target not found (catches CPU fallback silently)
export OMP_TARGET_OFFLOAD=mandatory

# Prevent memory fragmentation on NVIDIA GPUs
export CUDA_DEVICE_MAX_CONNECTIONS=1

echo "========================================================"
echo "Job:    ${SLURM_JOB_ID:-local}"
echo "Ranks:  ${NRANKS}  (one GPU each)"
echo "BASE_NX=${BASE_NX}  BASE_LEN=${BASE_LEN}m  cell_size=$(python3 -c "print(f'{${BASE_LEN}/${BASE_NX}:.2f}')" 2>/dev/null || echo '?')m  GPU_MODE=${GPU_MODE}"
echo "FINALTIME=${FINALTIME}s  YIELDSTEP=${YIELDSTEP}s"
echo "Host:   $(hostname)"
echo "Date:   $(date)"
echo "========================================================"

# =============================================================================
# Run benchmark
# =============================================================================

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTFILE="${OUTDIR}/weak_scaling_NP$(printf '%04d' ${NRANKS})_${TIMESTAMP}.json"

srun --ntasks="${NRANKS}" \
     python "${REPO_ROOT}/benchmarks/run_weak_scaling.py" \
        --base-nx "${BASE_NX}" \
        --base-len "${BASE_LEN}" \
        --mode "${GPU_MODE}" \
        --finaltime "${FINALTIME}" \
        --yieldstep "${YIELDSTEP}" \
        --output "${OUTFILE}" \
        --verbose

echo "========================================================"
echo "Result: ${OUTFILE}"
echo "Done:   $(date)"
echo "========================================================"
