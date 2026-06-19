# CuPy/CUDA archive

These files implement an earlier GPU acceleration approach using CuPy as a
NumPy drop-in replacement, alongside a raw CUDA C kernel (`cuda_anuga.cu`).

They were moved out of `anuga/shallow_water/` to reduce clutter.  The active
GPU/OpenMP offloading work (for SC26) lives in the main `develop` branch under
`anuga/shallow_water/sw_domain_openmp_ext.pyx` and the associated C headers.

## Contents

| File | Description |
|------|-------------|
| `sw_domain_cupy.py` | Domain subclass using CuPy arrays on the GPU |
| `sw_domain_cuda.py` | Domain subclass driving the raw CUDA kernel |
| `cuda_anuga.cu` | Raw CUDA C kernel for flux/update computations |
| `tests/` | CuPy-specific test scripts and pytest file |

## Status

Archived 2026-04-05.  Not wired into meson, not imported by `anuga/__init__.py`.
Preserved in case the CuPy approach is revisited.
