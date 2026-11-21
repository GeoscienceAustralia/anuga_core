# The anuga-gpu wheel: ANUGA built with nvc for NVIDIA GPUs (OpenMP target
# offload), single GPU, no MPI, with the NVIDIA HPC SDK run-time libraries it
# needs. See README.md in this directory.
#
# Build, with the repository root as the context:
#   docker build -f tools/gpu_wheel/Dockerfile \
#     --build-arg ANUGA_VERSION="$(python _git_version.py)" \
#     --build-arg PYTAG=cp313-cp313 \
#     --target out --output type=local,dest=dist-gpu .
#
# The compiler comes from the SDK's container; the wheel is built on
# manylinux_2_28 so that it installs on any Linux with glibc 2.28 or later.
# CUDA 12.9 is used rather than 13 so that older NVIDIA drivers still work.
ARG NVHPC_IMAGE=nvcr.io/nvidia/nvhpc:25.9-devel-cuda_multi-ubuntu24.04
FROM ${NVHPC_IMAGE} AS nvhpc

FROM quay.io/pypa/manylinux_2_28_x86_64 AS build
ARG NVHPC_VER=25.9
ARG CUDA_VER=12.9
ENV NVROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/${NVHPC_VER}
COPY --from=nvhpc ${NVROOT}/compilers ${NVROOT}/compilers
COPY --from=nvhpc ${NVROOT}/cuda/${CUDA_VER} ${NVROOT}/cuda/${CUDA_VER}
# REDIST holds links into compilers/lib: the redistributable run-time files.
COPY --from=nvhpc ${NVROOT}/REDIST/compilers ${NVROOT}/REDIST/compilers
ENV PATH=${NVROOT}/compilers/bin:$PATH
# Point nvc at this container's gcc, headers and linker.
RUN makelocalrc -x ${NVROOT}/compilers/bin && nvc --version

ARG PYTAG=cp312-cp312
ENV PY=/opt/python/${PYTAG}/bin
RUN $PY/pip install -q "numpy>=2.0" cython "meson>=1.4" meson-python ninja pybind11 build wheel

WORKDIR /src
COPY . .

# The distribution is anuga-gpu; the package it installs is still anuga.
# Its licence is ANUGA's plus NVIDIA's for the run-time files.
RUN cp ${NVROOT}/compilers/license/LICENSE.txt LICENSE_NVIDIA_HPC_SDK.txt && \
    cp tools/gpu_wheel/NOTICE_NVIDIA_RUNTIME.txt . && \
    sed -i \
      -e "s/^name = 'anuga'$/name = 'anuga-gpu'/" \
      -e 's/^license-files = \["LICENSE.txt"\]$/license-files = ["LICENSE.txt", "NOTICE_NVIDIA_RUNTIME.txt", "LICENSE_NVIDIA_HPC_SDK.txt"]/' \
      -e 's/^license = "Apache-2.0"$/license = "Apache-2.0 AND LicenseRef-NVIDIA-HPC-SDK"/' \
      -e "s/^description = .*/description = 'ANUGA built for NVIDIA GPUs: tsunami and flood modelling'/" \
      pyproject.toml && \
    grep -q "^name = 'anuga-gpu'$" pyproject.toml && \
    grep -q "LicenseRef-NVIDIA-HPC-SDK" pyproject.toml && \
    grep -q "NOTICE_NVIDIA_RUNTIME.txt" pyproject.toml

# -tp=haswell: nvc otherwise targets the build machine's CPU.
ARG GPU_ARCH=cuda12.9,cc75,cc80,cc86,cc89,cc90,cc120
ARG ANUGA_VERSION=
RUN set -o pipefail && \
    ANUGA_VERSION="${ANUGA_VERSION}" CC=nvc CFLAGS="-tp=haswell" PATH=$PY:$PATH \
    $PY/python -m build --wheel --no-isolation -o /wheel-raw \
      -Csetup-args=-Dgpu_offload=true -Csetup-args=-Dgpu_arch="${GPU_ARCH}" \
      2>&1 | tee /build.log | grep -E "GPU extension will|GPU arch|Successfully built" && \
    grep -q "GPU offloading enabled" /build.log

RUN PATH=$PY:$PATH $PY/python tools/gpu_wheel/bundle_runtime.py /wheel-raw/*.whl \
      ${NVROOT}/REDIST/compilers/lib LICENSE_NVIDIA_HPC_SDK.txt \
      tools/gpu_wheel/NOTICE_NVIDIA_RUNTIME.txt /wheel && ls -la /wheel

FROM scratch AS out
COPY --from=build /wheel /
COPY --from=build /build.log /
