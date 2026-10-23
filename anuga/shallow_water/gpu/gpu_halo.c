// GPU-accelerated shallow water solver
// Split from sw_domain_gpu.c for maintainability

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
// MPI (or single-process stubs) come in via gpu_domain.h
#include "gpu_domain.h"
#include "gpu_omp_macros.h"
#include "gpu_nvtx.h"

// Halo exchange setup and MPI ghost exchange

// ============================================================================
// Halo Exchange Setup
// ============================================================================

// Number of doubles exchanged per halo element.
//
// The halo carries the CONSERVED centroid quantities: stage, xmom, ymom, and
// then one slot per tracer for m = h*c. Derived quantities are not exchanged --
// height is recomputed from the exchanged stage in extrapolation Step 1, and c
// is recomputed from the exchanged m in exactly the same place, for ghost cells
// as well as owned ones.
//
// The stride is DERIVED from number_of_tracers rather than stored in
// struct halo_exchange, because that struct is embedded in struct gpu_domain,
// which embeds struct domain D -- adding a member there risks the silent offset
// aliasing of HANDOVER.md 2.1. Deriving it is safe because the count cannot
// change between init and exchange: gpu_halo_init runs after
// get_domain_pointers has set number_of_tracers, and add_tracer() tears the
// whole GPU interface down, forcing a re-init.
//
// At Ns=0 this is exactly 3, so the no-tracer path is byte-for-byte unchanged.
static inline int gpu_halo_stride(const struct gpu_domain *GD) {
    return 3 + (int)GD->D.number_of_tracers;
}

int gpu_halo_init(struct gpu_domain *GD,
                  int num_neighbors,
                  int *neighbor_ranks,
                  int *send_counts,
                  int *recv_counts,
                  int *flat_send_indices,
                  int *flat_recv_indices) {
    struct halo_exchange *H = &GD->halo;
    int stride;
#ifdef GPU_AWARE_MPI
    int dev;
#endif

    H->num_neighbors = num_neighbors;

    // Start from a clean slate so gpu_halo_finalize() can be used as the
    // single cleanup path if any allocation below fails.
    H->neighbor_ranks = NULL;
    H->send_counts = NULL;
    H->recv_counts = NULL;
    H->send_offsets = NULL;
    H->recv_offsets = NULL;
    H->flat_send_indices = NULL;
    H->flat_recv_indices = NULL;
    H->send_buffer = NULL;
    H->recv_buffer = NULL;
    H->host_send_buffer = NULL;
    H->host_recv_buffer = NULL;
    H->requests = NULL;

    if (num_neighbors == 0) {
        // No communication needed
        return 0;
    }

    // Allocate and copy neighbor info
    H->neighbor_ranks = (int *)malloc(num_neighbors * sizeof(int));
    H->send_counts = (int *)malloc(num_neighbors * sizeof(int));
    H->recv_counts = (int *)malloc(num_neighbors * sizeof(int));
    H->send_offsets = (int *)malloc((num_neighbors + 1) * sizeof(int));
    H->recv_offsets = (int *)malloc((num_neighbors + 1) * sizeof(int));
    if (!H->neighbor_ranks || !H->send_counts || !H->recv_counts ||
        !H->send_offsets || !H->recv_offsets) {
        gpu_set_error(GD, "ERROR: malloc failed for halo neighbour tables");
        goto fail;
    }

    memcpy(H->neighbor_ranks, neighbor_ranks, num_neighbors * sizeof(int));
    memcpy(H->send_counts, send_counts, num_neighbors * sizeof(int));
    memcpy(H->recv_counts, recv_counts, num_neighbors * sizeof(int));

    // Compute total sizes and offsets
    H->total_send_size = 0;
    H->total_recv_size = 0;
    H->send_offsets[0] = 0;
    H->recv_offsets[0] = 0;

    for (int ni = 0; ni < num_neighbors; ni++) {
        H->total_send_size += send_counts[ni];
        H->total_recv_size += recv_counts[ni];
        H->send_offsets[ni + 1] = H->total_send_size;
        H->recv_offsets[ni + 1] = H->total_recv_size;
    }

    // Allocate and copy flattened index arrays
    H->flat_send_indices = (int *)malloc(H->total_send_size * sizeof(int));
    H->flat_recv_indices = (int *)malloc(H->total_recv_size * sizeof(int));
    if (!H->flat_send_indices || !H->flat_recv_indices) {
        gpu_set_error(GD, "ERROR: malloc failed for halo index arrays");
        goto fail;
    }

    memcpy(H->flat_send_indices, flat_send_indices, H->total_send_size * sizeof(int));
    memcpy(H->flat_recv_indices, flat_recv_indices, H->total_recv_size * sizeof(int));

    // Allocate communication buffers
    // stride quantities per element: stage, xmom, ymom, then one m per tracer
    stride = gpu_halo_stride(GD);
#ifdef GPU_AWARE_MPI
    // Device buffers for GPU pack/unpack kernels
    dev = omp_get_default_device();
    H->send_buffer = (double *)omp_target_alloc(stride * H->total_send_size * sizeof(double), dev);
    H->recv_buffer = (double *)omp_target_alloc(stride * H->total_recv_size * sizeof(double), dev);
    if (!H->send_buffer || !H->recv_buffer) {
        gpu_set_error(GD, "ERROR: omp_target_alloc failed for halo buffers");
        goto fail;
    }
    // Host staging buffers for MPI calls.
    // Some UCX transports (e.g. uct_mm shared-memory used intra-node) cannot
    // access omp_target_alloc device pointers, causing a SIGSEGV in MPI_Isend.
    // We always stage through host memory; the overhead is small because halos
    // are tiny compared to the full domain.
    H->host_send_buffer = (double *)malloc(stride * H->total_send_size * sizeof(double));
    H->host_recv_buffer = (double *)malloc(stride * H->total_recv_size * sizeof(double));
    if (!H->host_send_buffer || !H->host_recv_buffer) {
        gpu_set_error(GD, "ERROR: malloc failed for halo host staging buffers");
        goto fail;
    }
#else
    H->send_buffer = (double *)malloc(stride * H->total_send_size * sizeof(double));
    H->recv_buffer = (double *)malloc(stride * H->total_recv_size * sizeof(double));
    H->host_send_buffer = NULL;
    H->host_recv_buffer = NULL;
    if (!H->send_buffer || !H->recv_buffer) {
        gpu_set_error(GD, "ERROR: malloc failed for halo buffers");
        goto fail;
    }
#endif

    // Allocate MPI request array
    H->requests = (MPI_Request *)malloc(2 * num_neighbors * sizeof(MPI_Request));
    if (!H->requests) {
        gpu_set_error(GD, "ERROR: malloc failed for halo MPI request array");
        goto fail;
    }

    if (GD->rank == 0) {
        printf("GPU halo exchange initialized:\n");
        printf("  Neighbors: %d\n", num_neighbors);
        printf("  Total send: %d elements\n", H->total_send_size);
        printf("  Total recv: %d elements\n", H->total_recv_size);
    }

    return 0;

fail:
    // Release whatever was acquired before the failure. Every pointer is
    // either NULL or a live allocation at this point, which is what
    // gpu_halo_finalize expects.
    gpu_halo_finalize(GD);
    return -1;
}

void gpu_halo_finalize(struct gpu_domain *GD) {
    struct halo_exchange *H = &GD->halo;

    if (H->neighbor_ranks) free(H->neighbor_ranks);
    if (H->send_counts) free(H->send_counts);
    if (H->recv_counts) free(H->recv_counts);
    if (H->send_offsets) free(H->send_offsets);
    if (H->recv_offsets) free(H->recv_offsets);
    if (H->flat_send_indices) free(H->flat_send_indices);
    if (H->flat_recv_indices) free(H->flat_recv_indices);
#ifdef GPU_AWARE_MPI
    int dev = omp_get_default_device();
    if (H->send_buffer) omp_target_free(H->send_buffer, dev);
    if (H->recv_buffer) omp_target_free(H->recv_buffer, dev);
    if (H->host_send_buffer) free(H->host_send_buffer);
    if (H->host_recv_buffer) free(H->host_recv_buffer);
#else
    if (H->send_buffer) free(H->send_buffer);
    if (H->recv_buffer) free(H->recv_buffer);
#endif
    if (H->requests) free(H->requests);

    H->num_neighbors = 0;
    H->neighbor_ranks = NULL;
    H->send_counts = NULL;
    H->recv_counts = NULL;
    H->send_offsets = NULL;
    H->recv_offsets = NULL;
    H->flat_send_indices = NULL;
    H->flat_recv_indices = NULL;
    H->send_buffer = NULL;
    H->recv_buffer = NULL;
    H->host_send_buffer = NULL;
    H->host_recv_buffer = NULL;
    H->requests = NULL;
}


// ============================================================================
// Ghost Exchange - Key MPI Function
// ============================================================================

// Exchange ghost cell data between MPI ranks
// Adapted from miniapp_mpi.c exchange_halo()
void gpu_exchange_ghosts(struct gpu_domain *GD) {
    NVTX_PUSH("gpu_exchange_ghosts");
    struct halo_exchange *H = &GD->halo;

    if (H->num_neighbors == 0) {
        NVTX_POP();
        return;
    }

    int send_size = H->total_send_size;
    int recv_size = H->total_recv_size;

    const int stride = gpu_halo_stride(GD);
    const int ns = (int)GD->D.number_of_tracers;
    const anuga_int n = GD->D.number_of_elements;

    double *stage = GD->D.stage_centroid_values;
    double *xmom = GD->D.xmom_centroid_values;
    double *ymom = GD->D.ymom_centroid_values;
    // Loaded at function scope, never inside the loop: these loops are 'omp
    // target' regions on a GPU build and D is not mapped to the device.
    double *t_cons = GD->D.tracer_conserved_values;
    double *send_buf = H->send_buffer;
    double *recv_buf = H->recv_buffer;
    int *flat_send = H->flat_send_indices;
    int *flat_recv = H->flat_recv_indices;

    // Pack send buffer on GPU
#ifdef GPU_AWARE_MPI
    // send_buf/recv_buf are omp_target_alloc'd device pointers — use is_device_ptr
    #pragma omp target teams distribute parallel for is_device_ptr(send_buf)
#else
    OMP_PARALLEL_LOOP
#endif
    for (int idx = 0; idx < send_size; idx++) {
        int k = flat_send[idx];  // Local element index
        send_buf[stride*idx + 0] = stage[k];
        send_buf[stride*idx + 1] = xmom[k];
        send_buf[stride*idx + 2] = ymom[k];
        for (int s_i = 0; s_i < ns; s_i++) {
            send_buf[stride*idx + 3 + s_i] = t_cons[(anuga_int)s_i * n + k];
        }
    }

#ifdef GPU_AWARE_MPI
    // GPU_AWARE_MPI path: device buffers are used for GPU pack/unpack, but MPI
    // communication uses host staging buffers.  Some UCX transports (e.g.
    // uct_mm, the intra-node shared-memory transport) cannot access
    // omp_target_alloc device pointers and segfault inside MPI_Isend if we
    // pass device pointers directly.  Staging through host is safe for ALL
    // transports, and halos are small enough that the D2H/H2D cost is minimal.
    {
        double *host_send = H->host_send_buffer;
        double *host_recv = H->host_recv_buffer;

        // Copy packed send buffer from device to host staging buffer
        int host = omp_get_initial_device();
        int dev  = omp_get_default_device();
        omp_target_memcpy(host_send, send_buf,
                          stride * send_size * sizeof(double),
                          0, 0, host, dev);

        int req_count = 0;
        int send_offset = 0, recv_offset = 0;

        // Post all receives first (into host staging buffer)
        for (int ni = 0; ni < H->num_neighbors; ni++) {
            int partner = H->neighbor_ranks[ni];
            int count = H->recv_counts[ni];
            MPI_Irecv(&host_recv[stride*recv_offset], stride*count, MPI_DOUBLE,
                      partner, 0, GD->comm, &H->requests[req_count++]);
            recv_offset += count;
        }

        // Post all sends (from host staging buffer)
        for (int ni = 0; ni < H->num_neighbors; ni++) {
            int partner = H->neighbor_ranks[ni];
            int count = H->send_counts[ni];
            MPI_Isend(&host_send[stride*send_offset], stride*count, MPI_DOUBLE,
                      partner, 0, GD->comm, &H->requests[req_count++]);
            send_offset += count;
        }

        // Wait for all communication to complete
        MPI_Waitall(req_count, H->requests, MPI_STATUSES_IGNORE);

        // Copy received halo data from host staging buffer to device
        omp_target_memcpy(recv_buf, host_recv,
                          stride * recv_size * sizeof(double),
                          0, 0, dev, host);
    }

#else
    // Non-GPU-aware MPI path: transfer halo buffers through host
    // This is still efficient because halo is much smaller than full domain

    // Copy packed send buffer from GPU to host
    #pragma omp target update from(send_buf[0:stride*send_size])

    // MPI communication on host
    int req_count = 0;
    int send_offset = 0, recv_offset = 0;

    // Post all receives first
    for (int ni = 0; ni < H->num_neighbors; ni++) {
        int partner = H->neighbor_ranks[ni];
        int count = H->recv_counts[ni];
        MPI_Irecv(&recv_buf[stride*recv_offset], stride*count, MPI_DOUBLE,
                  partner, 0, GD->comm, &H->requests[req_count++]);
        recv_offset += count;
    }

    // Post all sends
    for (int ni = 0; ni < H->num_neighbors; ni++) {
        int partner = H->neighbor_ranks[ni];
        int count = H->send_counts[ni];
        MPI_Isend(&send_buf[stride*send_offset], stride*count, MPI_DOUBLE,
                  partner, 0, GD->comm, &H->requests[req_count++]);
        send_offset += count;
    }

    // Wait for all communication to complete
    MPI_Waitall(req_count, H->requests, MPI_STATUSES_IGNORE);

    // Copy received halo data from host to GPU
    #pragma omp target update to(recv_buf[0:stride*recv_size])
#endif

    // Unpack receive buffer on GPU
#ifdef GPU_AWARE_MPI
    #pragma omp target teams distribute parallel for is_device_ptr(recv_buf)
#else
    OMP_PARALLEL_LOOP
#endif
    for (int idx = 0; idx < recv_size; idx++) {
        int k = flat_recv[idx];  // Local ghost element index
        stage[k] = recv_buf[stride*idx + 0];
        xmom[k] = recv_buf[stride*idx + 1];
        ymom[k] = recv_buf[stride*idx + 2];
        for (int s_i = 0; s_i < ns; s_i++) {
            t_cons[(anuga_int)s_i * n + k] = recv_buf[stride*idx + 3 + s_i];
        }
    }
    NVTX_POP();
}


// Exchange the BED of ghost cells (#424): each ghost takes its owner's bed
// centroid, and its three edge values shift by the same amount (the DE bed is
// discontinuous; the sediment kernels move a cell's edges with its centroid).
// The per-step exchange above carries the conserved quantities and tracers
// but not elevation, so under bedload -- whose update of a cell reads its
// neighbours, and an outer ghost lacks some -- ghost beds drifted from their
// owners'. Reuses the halo buffers, one value per cell (their stride is >= 3),
// with its own tag. Called after the sediment step when the bed evolves.
void gpu_exchange_bed_ghosts(struct gpu_domain *GD) {
    struct halo_exchange *H = &GD->halo;
    if (H->num_neighbors == 0) return;

    int send_size = H->total_send_size;
    int recv_size = H->total_recv_size;
    double *bed_cv = GD->D.bed_centroid_values;
    double *bed_ev = GD->D.bed_edge_values;
    double *send_buf = H->send_buffer;
    double *recv_buf = H->recv_buffer;
    int *flat_send = H->flat_send_indices;
    int *flat_recv = H->flat_recv_indices;
    const int tag = 2;   // 0 is the conserved-quantity exchange

#ifdef GPU_AWARE_MPI
    #pragma omp target teams distribute parallel for is_device_ptr(send_buf)
#else
    OMP_PARALLEL_LOOP
#endif
    for (int idx = 0; idx < send_size; idx++) {
        send_buf[idx] = bed_cv[flat_send[idx]];
    }

#ifdef GPU_AWARE_MPI
    {
        double *host_send = H->host_send_buffer;
        double *host_recv = H->host_recv_buffer;
        int host = omp_get_initial_device();
        int dev  = omp_get_default_device();
        omp_target_memcpy(host_send, send_buf, send_size * sizeof(double),
                          0, 0, host, dev);
        int req_count = 0, send_offset = 0, recv_offset = 0;
        for (int ni = 0; ni < H->num_neighbors; ni++) {
            int count = H->recv_counts[ni];
            MPI_Irecv(&host_recv[recv_offset], count, MPI_DOUBLE,
                      H->neighbor_ranks[ni], tag, GD->comm, &H->requests[req_count++]);
            recv_offset += count;
        }
        for (int ni = 0; ni < H->num_neighbors; ni++) {
            int count = H->send_counts[ni];
            MPI_Isend(&host_send[send_offset], count, MPI_DOUBLE,
                      H->neighbor_ranks[ni], tag, GD->comm, &H->requests[req_count++]);
            send_offset += count;
        }
        MPI_Waitall(req_count, H->requests, MPI_STATUSES_IGNORE);
        omp_target_memcpy(recv_buf, host_recv, recv_size * sizeof(double),
                          0, 0, dev, host);
    }
#else
    #pragma omp target update from(send_buf[0:send_size])
    {
        int req_count = 0, send_offset = 0, recv_offset = 0;
        for (int ni = 0; ni < H->num_neighbors; ni++) {
            int count = H->recv_counts[ni];
            MPI_Irecv(&recv_buf[recv_offset], count, MPI_DOUBLE,
                      H->neighbor_ranks[ni], tag, GD->comm, &H->requests[req_count++]);
            recv_offset += count;
        }
        for (int ni = 0; ni < H->num_neighbors; ni++) {
            int count = H->send_counts[ni];
            MPI_Isend(&send_buf[send_offset], count, MPI_DOUBLE,
                      H->neighbor_ranks[ni], tag, GD->comm, &H->requests[req_count++]);
            send_offset += count;
        }
        MPI_Waitall(req_count, H->requests, MPI_STATUSES_IGNORE);
    }
    #pragma omp target update to(recv_buf[0:recv_size])
#endif

#ifdef GPU_AWARE_MPI
    #pragma omp target teams distribute parallel for is_device_ptr(recv_buf)
#else
    OMP_PARALLEL_LOOP
#endif
    for (int idx = 0; idx < recv_size; idx++) {
        int k = flat_recv[idx];
        double shift = recv_buf[idx] - bed_cv[k];
        bed_cv[k] = recv_buf[idx];
        bed_ev[3 * k + 0] += shift;
        bed_ev[3 * k + 1] += shift;
        bed_ev[3 * k + 2] += shift;
    }
}
