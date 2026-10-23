// NVTX profiling macros for ANUGA GPU kernels
//
// When compiled with -DNVTX_ENABLED (set automatically by meson when
// nvToolsExt is found and -Duse_nvtx=true is passed), these macros wrap
// NVIDIA Tools Extension (NVTX) range markers so that Nsight Systems
// (nsys) and Nsight Compute (ncu) can attribute time to named regions.
//
// Without -DNVTX_ENABLED every macro expands to ((void)0), adding zero
// overhead to release builds.
//
// Usage in C source files:
//
//   #include "gpu_nvtx.h"
//
//   void gpu_compute_fluxes(struct gpu_domain *GD) {
//       NVTX_PUSH("gpu_compute_fluxes");
//       // ... kernel work ...
//       NVTX_POP();
//   }
//
// To profile with Nsight Systems:
//
//   OMP_NUM_THREADS=1 OMP_TARGET_OFFLOAD=mandatory
//       nsys profile --trace=openmp,nvtx,cuda
//       python my_script.py
//
// Both NVTX v1 (nvToolsExt.h from CUDA toolkit) and NVTX v3 (nvtx3 headers,
// standalone NVIDIA package) are supported:
//   - v1: include nvToolsExt.h  (CUDA_HOME/include)
//   - v3: include nvtx3/nvToolsExt.h  (standalone install)
//
// The meson build probes for v3 first, then v1.

#ifndef GPU_NVTX_H
#define GPU_NVTX_H

#ifdef NVTX_ENABLED

#  ifdef HAVE_NVTX3
#    include <nvtx3/nvToolsExt.h>
#  else
#    include <nvToolsExt.h>
#  endif

// Push a named range — visible in Nsight timeline as a coloured bar
#  define NVTX_PUSH(name)  nvtxRangePushA(name)
// Pop the innermost range
#  define NVTX_POP()       nvtxRangePop()

#else  /* NVTX_ENABLED not defined */

#  define NVTX_PUSH(name)  ((void)0)
#  define NVTX_POP()       ((void)0)

#endif /* NVTX_ENABLED */

#endif /* GPU_NVTX_H */
