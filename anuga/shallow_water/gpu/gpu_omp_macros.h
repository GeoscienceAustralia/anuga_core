// OpenMP pragma macros for GPU vs CPU execution
//
// When CPU_ONLY_MODE is defined (via -Dgpu_offload=false), use regular OpenMP
// Otherwise, use OpenMP target for GPU offloading
//
// This allows the same code to run on both GPU and CPU multicore

#ifndef GPU_OMP_MACROS_H
#define GPU_OMP_MACROS_H

// Helper macro to stringify pragma arguments (allows macro expansion before stringification)
#define DO_PRAGMA(x) _Pragma(#x)

#ifdef CPU_ONLY_MODE

// ============================================================================
// CPU MULTICORE MODE - Regular OpenMP, no device offloading
// ============================================================================

// Parallel loops with SIMD vectorization
#define OMP_PARALLEL_LOOP _Pragma("omp parallel for simd")
#define OMP_PARALLEL_LOOP_SIMD _Pragma("omp parallel for simd")

// Reductions - use DO_PRAGMA to allow variable name expansion
// Note: simd is omitted here; reduction loops cannot be SIMD-vectorized and
// clang warns (-Wpass-failed=transform-warning) when simd is requested.
#define OMP_PARALLEL_LOOP_REDUCTION_PLUS(var) DO_PRAGMA(omp parallel for reduction(+:var))
#define OMP_PARALLEL_LOOP_REDUCTION_MIN(var) DO_PRAGMA(omp parallel for reduction(min:var))
#define OMP_PARALLEL_LOOP_REDUCTION_MAX(var) DO_PRAGMA(omp parallel for reduction(max:var))

// Data mapping - no-op on CPU (data already in host memory)
#define OMP_TARGET_ENTER_DATA_MAP_TO(...)
#define OMP_TARGET_ENTER_DATA_MAP_ALLOC(...)
#define OMP_TARGET_EXIT_DATA_MAP_DELETE(...)

// Data transfer - no-op on CPU
#define OMP_TARGET_UPDATE_TO(...)
#define OMP_TARGET_UPDATE_FROM(...)

// Device pointer clause - no-op on CPU
#define OMP_IS_DEVICE_PTR(ptr)

// ============================================================================
// CPU MULTICORE MODE stubs for OpenMP target-offloading API
//
// The standard libomp on macOS (and minimal Linux builds) does NOT ship the
// OpenMP 4.5 target-allocation API (omp_target_alloc, omp_target_free,
// omp_target_memcpy, omp_target_is_present).  In CPU_ONLY_MODE there is no
// GPU device, so we provide static inline stubs and redirect every call site
// via macros.  The redirect must happen here (after <omp.h> is already
// included by the .c files) so the macros suppress external symbol references
// without conflicting with omp.h's extern declarations.
// ============================================================================

#include <stdlib.h>
#include <string.h>

static inline void *_anuga_omp_target_alloc(size_t size, int dev) {
    (void)dev; return malloc(size);
}
static inline void _anuga_omp_target_free(void *ptr, int dev) {
    (void)dev; free(ptr);
}
static inline int _anuga_omp_target_memcpy(
        void *dst, const void *src, size_t length,
        size_t dst_offset, size_t src_offset,
        int dst_dev, int src_dev) {
    (void)dst_dev; (void)src_dev;
    memcpy((char *)dst + dst_offset, (const char *)src + src_offset, length);
    return 0;
}
static inline int _anuga_omp_target_is_present(const void *ptr, int dev) {
    (void)ptr; (void)dev; return 1;
}
static inline int _anuga_omp_get_initial_device(void) { return 0; }

#define omp_target_alloc        _anuga_omp_target_alloc
#define omp_target_free         _anuga_omp_target_free
#define omp_target_memcpy       _anuga_omp_target_memcpy
#define omp_target_is_present   _anuga_omp_target_is_present
#define omp_get_initial_device  _anuga_omp_get_initial_device

#else

// ============================================================================
// GPU MODE - OpenMP target offloading
// ============================================================================

// Parallel loops on device
#define OMP_PARALLEL_LOOP _Pragma("omp target teams loop")
#define OMP_PARALLEL_LOOP_SIMD _Pragma("omp target teams loop")

// Reductions on device - use DO_PRAGMA to allow variable name expansion
// Note: Using distribute parallel for for better reduction support
#define OMP_PARALLEL_LOOP_REDUCTION_PLUS(var) DO_PRAGMA(omp target teams distribute parallel for reduction(+:var))
#define OMP_PARALLEL_LOOP_REDUCTION_MIN(var) DO_PRAGMA(omp target teams loop reduction(min:var))
#define OMP_PARALLEL_LOOP_REDUCTION_MAX(var) DO_PRAGMA(omp target teams distribute parallel for reduction(max:var))

// Data mapping to device
// Note: These need to be used with care - the actual map clause arguments
// must be specified in the calling code since macros can't handle variable arguments well
// For now, these are placeholders - actual mapping is done inline
#define OMP_TARGET_ENTER_DATA_MAP_TO(...) _Pragma("omp target enter data map(to:" #__VA_ARGS__ ")")
#define OMP_TARGET_ENTER_DATA_MAP_ALLOC(...) _Pragma("omp target enter data map(alloc:" #__VA_ARGS__ ")")
#define OMP_TARGET_EXIT_DATA_MAP_DELETE(...) _Pragma("omp target exit data map(delete:" #__VA_ARGS__ ")")

// Data transfer
#define OMP_TARGET_UPDATE_TO(...) _Pragma("omp target update to(" #__VA_ARGS__ ")")
#define OMP_TARGET_UPDATE_FROM(...) _Pragma("omp target update from(" #__VA_ARGS__ ")")

// Device pointer clause
#define OMP_IS_DEVICE_PTR(ptr) is_device_ptr(ptr)

#endif // CPU_ONLY_MODE

#endif // GPU_OMP_MACROS_H
