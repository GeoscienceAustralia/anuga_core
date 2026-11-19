# The `anuga-gpu` wheel

A binary wheel of ANUGA for NVIDIA GPUs, so that a GPU run needs neither the
NVIDIA HPC SDK nor a source build:

```bash
pip install anuga_gpu-<version>-cp313-cp313-manylinux_2_28_x86_64.whl
```

```python
domain.set_compute_mode('unified')    # the kernels run on the GPU
```

It is the whole `anuga` package, built with `nvc` (`-Dgpu_offload=true`), under
the distribution name `anuga-gpu`. It installs the same files as `anuga`, so
install one or the other, not both.

| | |
|---|---|
| Platform | Linux x86-64, glibc 2.28 or later |
| GPUs | T4, A100, A10G, L4, H100, RTX 50-series (`cc75` to `cc120`) |
| Needs on the machine | an NVIDIA driver; no CUDA toolkit, no SDK |
| MPI | none: one process, one GPU |
| Size | 17 MB to download, 65 MB installed |

The wheel has not been published anywhere. It has been run on an RTX 5070 and
on a Google Colab T4.

## Building

From the repository root, once per Python version:

```bash
docker build -f tools/gpu_wheel/Dockerfile \
  --build-arg ANUGA_VERSION="$(python _git_version.py)" \
  --build-arg PYTAG=cp313-cp313 \
  --target out --output type=local,dest=dist-gpu .
```

The compiler is copied out of NVIDIA's SDK container
(`nvcr.io/nvidia/nvhpc:25.9-devel-cuda_multi-ubuntu24.04`, about 49 GB) into a
`manylinux_2_28` container, where the wheel is built. CUDA 12.9 is used rather
than 13 so that older drivers still work.

GitHub's runners have no GPU, so a build there can only be tested on the host
fallback. Test a wheel on a real GPU before releasing it:

```bash
OMP_TARGET_OFFLOAD=MANDATORY python my_run.py      # fails if no GPU is used
anuga_run_isolated_tests -cm unified               # the GPU test file
```

## The NVIDIA run-time libraries

An `nvc` build needs the SDK's OpenMP run time (`libnvomp`, `libacchost` and
seven others). We read the SDK licence (section 2A iii) as allowing an
application built with the SDK to carry the files in the SDK's `REDIST`
directories, on the conditions of section 2B. The build is designed around
those conditions:

* `bundle_runtime.py` copies only files found under `REDIST`, byte for byte,
  into `anuga/_nvhpc_runtime/`, and stops if a library the wheel needs is not
  there. Only ANUGA's own extension modules are changed, to find that
  directory. (`auditwheel repair` is not used: it renames and rewrites the
  libraries it bundles, and we do not read the licence as allowing them to
  be modified.)
* `NOTICE_NVIDIA_RUNTIME.txt` is meant to pass the conditions on to the user: no further
  distribution, no modification or reverse engineering, and the warranty and
  liability disclaimers on NVIDIA's behalf. It and the SDK licence are in the
  wheel's metadata and beside the libraries.
* The wheel's licence is `Apache-2.0 AND LicenseRef-NVIDIA-HPC-SDK`.
* The name has no NVIDIA trademark in it.

The GPU container images (`docker/Dockerfile.gpu.slim`, `Dockerfile.gpu.mpi`)
place the same licence and notice in `/opt/nvhpc-redist/`.

This is the project's reading of the licence. It is not legal advice and has
not been confirmed by NVIDIA. Whoever builds and distributes the wheel is the
licensee: read sections 2B and 7 (indemnity) of the licence before publishing
one.
