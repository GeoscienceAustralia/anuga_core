"""Tiny ANUGA smoke-test run for the Docker images.

Writes an SWW to /work/output so the S3 stage-out path has something to upload.
Set ANUGA_DEFAULT_COMPUTE_MODE=unified (and run with --gpus all on the GPU
image) to exercise the GPU offload path — you'll see the "GPU Domain" banner.

    docker run --rm -v "$PWD/docker/data:/work" anuga:cpu python /work/example_run.py
    docker run --rm --gpus all -e ANUGA_DEFAULT_COMPUTE_MODE=unified \
        -v "$PWD/docker/data:/work" anuga:gpu python /work/example_run.py
"""
import os
import anuga

print("ANUGA", anuga.__version__)
# GPU-offload API only exists on the unified/develop line (the source-built GPU
# image), not in the released PyPI wheel (the CPU image) — query defensively.
_supported = getattr(anuga, "gpu_offload_supported", None)
_enabled = getattr(anuga, "gpu_offload_enabled", None)
print("gpu_offload_supported:", _supported() if _supported else "n/a (release build)",
      "| gpu_offload_enabled:", _enabled() if _enabled else "n/a")

outdir = os.path.join(os.environ.get("ANUGA_WORKDIR", "/work"), "output")
os.makedirs(outdir, exist_ok=True)

points, vertices, boundary = anuga.rectangular_cross(40, 40, len1=100.0, len2=100.0)
domain = anuga.Domain(points, vertices, boundary)
domain.set_name("example")
domain.set_datadir(outdir)
domain.set_quantity("elevation", lambda x, y: -0.05 * x)
domain.set_quantity("stage", 0.5)
domain.set_quantity("friction", 0.03)

Br = anuga.Reflective_boundary(domain)
domain.set_boundary({"left": Br, "right": Br, "top": Br, "bottom": Br})

_mode = getattr(domain, "get_compute_mode", lambda: "n/a")
print("compute_mode:", _mode(), "| omp_num_threads:", getattr(domain, "omp_num_threads", "n/a"))

for t in domain.evolve(yieldstep=1.0, finaltime=5.0):
    print(domain.timestepping_statistics())

print("Wrote", os.path.join(outdir, "example.sww"))
