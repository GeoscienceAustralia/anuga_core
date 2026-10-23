// GPU-accelerated shallow water solver
// Split from sw_domain_gpu.c for maintainability

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
// MPI (or single-process stubs) come in via gpu_domain.h
#include "gpu_domain.h"

// FLOP counters for Gordon Bell performance profiling
#include "gpu_device_helpers.h"

// ============================================================================
// FLOP Counter Functions (Gordon Bell Performance Profiling)
// ============================================================================

void gpu_flop_counters_init(struct gpu_domain *GD) {
    // Initialize all FLOP counters to zero
    memset(&GD->flops, 0, sizeof(struct flop_counters));
    GD->flops.enabled = 1;  // Enable by default
}

void gpu_flop_counters_reset(struct gpu_domain *GD) {
    // Reset counters but keep enabled state
    int enabled = GD->flops.enabled;
    memset(&GD->flops, 0, sizeof(struct flop_counters));
    GD->flops.enabled = enabled;
}

void gpu_flop_counters_enable(struct gpu_domain *GD, int enable) {
    GD->flops.enabled = enable;
}

void gpu_flop_counters_start_timer(struct gpu_domain *GD) {
    GD->flops.start_time = MPI_Wtime();
}

void gpu_flop_counters_stop_timer(struct gpu_domain *GD) {
    GD->flops.elapsed_time = MPI_Wtime() - GD->flops.start_time;
}

uint64_t gpu_flop_counters_get_total(struct gpu_domain *GD) {
    // Recalculate total from individual counters
    GD->flops.total_flops =
        GD->flops.extrapolate_flops +
        GD->flops.compute_fluxes_flops +
        GD->flops.update_flops +
        GD->flops.protect_flops +
        GD->flops.manning_flops +
        GD->flops.backup_flops +
        GD->flops.saxpy_flops +
        GD->flops.rate_operator_flops +
        GD->flops.ghost_exchange_flops;
    return GD->flops.total_flops;
}

double gpu_flop_counters_get_flops(struct gpu_domain *GD) {
    // Return FLOP/s (floating point operations per second)
    if (GD->flops.elapsed_time <= 0.0) {
        return 0.0;
    }
    return (double)gpu_flop_counters_get_total(GD) / GD->flops.elapsed_time;
}

void gpu_flop_counters_print(struct gpu_domain *GD) {
    uint64_t total = gpu_flop_counters_get_total(GD);
    double gflops = (double)total / 1.0e9;
    double elapsed = GD->flops.elapsed_time;
    double gflops_per_sec = (elapsed > 0.0) ? gflops / elapsed : 0.0;

    printf("\n");
    printf("============================================================\n");
    printf("FLOP Counter Summary (Gordon Bell Profiling)\n");
    printf("============================================================\n");
    printf("Kernel                       |      FLOPs |     Calls |  FLOPs/call\n");
    printf("-----------------------------|------------|-----------|------------\n");

    if (GD->flops.extrapolate_calls > 0) {
        printf("extrapolate_second_order     | %10lu | %9lu | %10lu\n",
               (unsigned long)GD->flops.extrapolate_flops,
               (unsigned long)GD->flops.extrapolate_calls,
               (unsigned long)(GD->flops.extrapolate_flops / GD->flops.extrapolate_calls));
    }
    if (GD->flops.compute_fluxes_calls > 0) {
        printf("compute_fluxes               | %10lu | %9lu | %10lu\n",
               (unsigned long)GD->flops.compute_fluxes_flops,
               (unsigned long)GD->flops.compute_fluxes_calls,
               (unsigned long)(GD->flops.compute_fluxes_flops / GD->flops.compute_fluxes_calls));
    }
    if (GD->flops.update_calls > 0) {
        printf("update_conserved_quantities  | %10lu | %9lu | %10lu\n",
               (unsigned long)GD->flops.update_flops,
               (unsigned long)GD->flops.update_calls,
               (unsigned long)(GD->flops.update_flops / GD->flops.update_calls));
    }
    if (GD->flops.protect_calls > 0) {
        printf("protect                      | %10lu | %9lu | %10lu\n",
               (unsigned long)GD->flops.protect_flops,
               (unsigned long)GD->flops.protect_calls,
               (unsigned long)(GD->flops.protect_flops / GD->flops.protect_calls));
    }
    if (GD->flops.manning_calls > 0) {
        printf("manning_friction             | %10lu | %9lu | %10lu\n",
               (unsigned long)GD->flops.manning_flops,
               (unsigned long)GD->flops.manning_calls,
               (unsigned long)(GD->flops.manning_flops / GD->flops.manning_calls));
    }
    if (GD->flops.saxpy_calls > 0) {
        printf("saxpy_conserved_quantities   | %10lu | %9lu | %10lu\n",
               (unsigned long)GD->flops.saxpy_flops,
               (unsigned long)GD->flops.saxpy_calls,
               (unsigned long)(GD->flops.saxpy_flops / GD->flops.saxpy_calls));
    }
    if (GD->flops.rate_operator_calls > 0) {
        printf("rate_operator_apply          | %10lu | %9lu | %10lu\n",
               (unsigned long)GD->flops.rate_operator_flops,
               (unsigned long)GD->flops.rate_operator_calls,
               (unsigned long)(GD->flops.rate_operator_flops / GD->flops.rate_operator_calls));
    }

    printf("-----------------------------|------------|-----------|------------\n");
    printf("TOTAL                        | %10lu |\n", (unsigned long)total);
    printf("============================================================\n");
    printf("Total GFLOPs:     %.3f\n", gflops);
    printf("Elapsed time:     %.3f s\n", elapsed);
    printf("Performance:      %.3f GFLOP/s\n", gflops_per_sec);
    printf("============================================================\n\n");
}

// Per-kernel FLOP getters
uint64_t gpu_flop_counters_get_extrapolate(struct gpu_domain *GD) {
    return GD->flops.extrapolate_flops;
}

uint64_t gpu_flop_counters_get_compute_fluxes(struct gpu_domain *GD) {
    return GD->flops.compute_fluxes_flops;
}

uint64_t gpu_flop_counters_get_update(struct gpu_domain *GD) {
    return GD->flops.update_flops;
}

uint64_t gpu_flop_counters_get_protect(struct gpu_domain *GD) {
    return GD->flops.protect_flops;
}

uint64_t gpu_flop_counters_get_manning(struct gpu_domain *GD) {
    return GD->flops.manning_flops;
}

uint64_t gpu_flop_counters_get_backup(struct gpu_domain *GD) {
    return GD->flops.backup_flops;
}

uint64_t gpu_flop_counters_get_saxpy(struct gpu_domain *GD) {
    return GD->flops.saxpy_flops;
}

uint64_t gpu_flop_counters_get_rate_operator(struct gpu_domain *GD) {
    return GD->flops.rate_operator_flops;
}

uint64_t gpu_flop_counters_get_ghost_exchange(struct gpu_domain *GD) {
    return GD->flops.ghost_exchange_flops;
}

// ============================================================================
// MPI Reduction for Multi-GPU FLOP Counters (Gordon Bell)
// ============================================================================

uint64_t gpu_flop_counters_get_global_total(struct gpu_domain *GD) {
    // Get local total first
    uint64_t local_total = gpu_flop_counters_get_total(GD);

    // MPI_Allreduce to sum across all ranks
    // MPI_UNSIGNED_LONG_LONG maps to uint64_t
    uint64_t global_total = 0;
    MPI_Allreduce(&local_total, &global_total, 1, MPI_UNSIGNED_LONG_LONG,
                  MPI_SUM, GD->comm);

    return global_total;
}

double gpu_flop_counters_get_global_flops(struct gpu_domain *GD) {
    // Get global total FLOPs
    uint64_t global_total = gpu_flop_counters_get_global_total(GD);

    // Use local elapsed time (should be same across ranks if synchronized)
    if (GD->flops.elapsed_time <= 0.0) {
        return 0.0;
    }
    return (double)global_total / GD->flops.elapsed_time;
}

void gpu_flop_counters_print_global(struct gpu_domain *GD) {
    // Gather per-kernel FLOPs from all ranks
    uint64_t local_flops[9];
    local_flops[0] = GD->flops.extrapolate_flops;
    local_flops[1] = GD->flops.compute_fluxes_flops;
    local_flops[2] = GD->flops.update_flops;
    local_flops[3] = GD->flops.protect_flops;
    local_flops[4] = GD->flops.manning_flops;
    local_flops[5] = GD->flops.backup_flops;
    local_flops[6] = GD->flops.saxpy_flops;
    local_flops[7] = GD->flops.rate_operator_flops;
    local_flops[8] = GD->flops.ghost_exchange_flops;

    uint64_t global_flops[9];
    MPI_Reduce(local_flops, global_flops, 9, MPI_UNSIGNED_LONG_LONG,
               MPI_SUM, 0, GD->comm);

    // Get global total and max elapsed time
    uint64_t local_total = gpu_flop_counters_get_total(GD);
    uint64_t global_total = 0;
    MPI_Reduce(&local_total, &global_total, 1, MPI_UNSIGNED_LONG_LONG,
               MPI_SUM, 0, GD->comm);

    double local_elapsed = GD->flops.elapsed_time;
    double max_elapsed = 0.0;
    MPI_Reduce(&local_elapsed, &max_elapsed, 1, MPI_DOUBLE,
               MPI_MAX, 0, GD->comm);

    // Only rank 0 prints
    if (GD->rank != 0) return;

    double gflops = (double)global_total / 1.0e9;
    double gflops_per_sec = (max_elapsed > 0.0) ? gflops / max_elapsed : 0.0;

    printf("\n");
    printf("============================================================\n");
    printf("GLOBAL FLOP Counter Summary (Gordon Bell - %d GPUs)\n", GD->nprocs);
    printf("============================================================\n");
    printf("Kernel                       |      GFLOPs | %% of total\n");
    printf("-----------------------------|-------------|------------\n");

    const char* names[] = {
        "extrapolate_second_order",
        "compute_fluxes",
        "update_conserved_quantities",
        "protect",
        "manning_friction",
        "backup_conserved_quantities",
        "saxpy_conserved_quantities",
        "rate_operator_apply",
        "ghost_exchange"
    };

    for (int i = 0; i < 9; i++) {
        if (global_flops[i] > 0) {
            double kernel_gflops = (double)global_flops[i] / 1.0e9;
            double pct = 100.0 * (double)global_flops[i] / (double)global_total;
            printf("%-28s | %11.3f | %9.1f%%\n", names[i], kernel_gflops, pct);
        }
    }

    printf("-----------------------------|-------------|------------\n");
    printf("TOTAL                        | %11.3f | 100.0%%\n", gflops);
    printf("============================================================\n");
    printf("Number of GPUs:   %d\n", GD->nprocs);
    printf("Total GFLOPs:     %.3f\n", gflops);
    printf("Elapsed time:     %.3f s\n", max_elapsed);
    printf("Performance:      %.3f GFLOP/s\n", gflops_per_sec);
    printf("Per-GPU average:  %.3f GFLOP/s\n", gflops_per_sec / GD->nprocs);
    printf("============================================================\n\n");
}
