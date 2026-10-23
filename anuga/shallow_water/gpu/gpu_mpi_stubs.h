// Single-process MPI stubs for building the GPU kernels without MPI/mpi4py.
//
// Included by gpu_domain.h only when HAVE_MPI is NOT defined. The GPU kernels
// already guard every real MPI call behind `nprocs > 1`, `num_neighbors == 0`,
// or rank-role checks, so in single-process mode these symbols only need to
// *compile* - they are never executed. The exceptions are the (unguarded)
// flop-counter reductions and MPI_Wtime, which are given correct single-process
// behaviour below.
#ifndef GPU_MPI_STUBS_H
#define GPU_MPI_STUBS_H

#include <stddef.h>
#include <string.h>
#include <omp.h>

// Opaque handles - only ever stored in structs, never dereferenced when
// running single-process.
typedef int MPI_Comm;
typedef int MPI_Request;

#define MPI_COMM_WORLD      0
#define MPI_STATUSES_IGNORE NULL

// Reduction ops are unused single-process; values are irrelevant.
#define MPI_SUM 0
#define MPI_MAX 0
#define MPI_MIN 0

// Datatype "handles" encode the element byte size, so the reduction stubs
// below can copy the right number of bytes for any type.
#define MPI_DOUBLE             sizeof(double)
#define MPI_UNSIGNED_LONG_LONG sizeof(unsigned long long)

// Single-process reductions: the global result is just the local input.
// (For MPI_Reduce, root == self, so the copy is always correct here.)
#define MPI_Allreduce(sendbuf, recvbuf, count, dtype, op, comm) \
    memcpy((recvbuf), (sendbuf), (size_t)(count) * (size_t)(dtype))
#define MPI_Reduce(sendbuf, recvbuf, count, dtype, op, root, comm) \
    memcpy((recvbuf), (sendbuf), (size_t)(count) * (size_t)(dtype))

// Point-to-point calls are only reached when num_neighbors > 0, i.e. never in
// single-process mode. Provide no-op stubs so the code compiles.
#define MPI_Isend(buf, count, dtype, dest, tag, comm, request)   (0)
#define MPI_Irecv(buf, count, dtype, source, tag, comm, request) (0)
#define MPI_Waitall(count, requests, statuses)                   (0)

// Wall-clock timer for the flop counters (gpu_flop.c already includes <omp.h>).
#define MPI_Wtime() omp_get_wtime()

#endif // GPU_MPI_STUBS_H
