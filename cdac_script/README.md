# JORGE's secret sauce 

Ok, buckle up. 

## convert.py 

Easy script to transform a tch to msh file thereby reducing storage. You can get a "help" by running the script without args:

```
(anuga_env_3.13) [jlv900@gadi-gpu-v100-0035 cdac_script]$ python convert.py

Convert .tsh (ASCII) mesh files to .msh (NetCDF binary) format

Binary .msh files load 10-20x faster than ASCII .tsh files.

Usage:
    python convert_tsh_to_msh.py input.tsh output.msh
    python convert_tsh_to_msh.py input.tsh  # Creates input.msh
```

For further use, make sure the msh file gets redirected as `${resolution}sqm.msh`

## partition_mahanadi.py

This basically invokes metis to partition the domain and INITIALIZE it. It needs to be run from the `300_data` directory or the paths inside the script
need to be amended. 


```
usage: partition_mahanadi.py [-h] --mesh-size {100,300,900} --nprocs NPROCS [--output OUTPUT] [--chilka-stage CHILKA_STAGE] [--verbose]

Partition Mahanadi Delta mesh

options:
  -h, --help            show this help message and exit
  --mesh-size {100,300,900}
                        Mesh resolution: 100, 300, or 900 sqm
  --nprocs NPROCS       Number of partitions to create
  --output OUTPUT       Output directory for partition files
  --chilka-stage CHILKA_STAGE
                        Initial stage inside Chilka polygon (default: 0.15)
  --verbose             Verbose output
```

Things will get stored in a partitions/ directory (around 5.7 GB total for 900sqm).


## run_model_3_partitioned.py  

Loads the prepartitioned data from `partitions/` and gets going asap. IT has no nice CLI arguments, 
you'll have to go in an ponit to the partitions/ and also the grid size: 

```
MESH_SIZE = 900  # <-- CHANGE THIS for different mesh sizes
```

This is why it is nice to call things `${resolution}sqm.msh` so that I can just change that name and boom. 


## Building ANUGA GPU

Requirements:

- Preferably the conda_env_3.13 environment, 3.12 was evil (memory leak)
- A supported compiler with OpenMP support (see table below)
- CUDA for GPU code generation (NVIDIA GPUs)
- Depending on the behaviour of your system `module load your-mpi` mpi4py seems to come with one, but be careful!

### Supported Compilers

| Compiler | CPU Multicore | GPU Offload | Notes |
|----------|---------------|-------------|-------|
| GCC | ✅ | ✅ (nvptx) | `-fopenmp`, GPU via libgomp offloading |
| Clang/LLVM | ✅ | ✅ (NVIDIA/AMD) | `-fopenmp`, GPU via LLVM offloading |
| NVIDIA HPC SDK | ✅ | ✅ (NVIDIA) | `-mp=multicore` / `-mp=gpu`, best NVIDIA support |
| Intel oneAPI | ✅ | ❌ | `-qopenmp` or `-fiopenmp` |
| Apple Clang | ✅ | ❌ | Requires libomp from Homebrew |

For GPU offloading:
- **NVIDIA GPUs**: Use NVIDIA HPC SDK (recommended) or Clang with LLVM offloading
- **AMD GPUs**: Use Clang/AOMP with ROCm

### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `gpu_offload` | `false` | Target GPU (`true`) or CPU multicore (`false`) |
| `gpu_aware_mpi` | `false` | Enable GPU-aware MPI for direct device communication |
| `gpu_arch` | `cc70` | GPU architecture: `cc70` (V100), `cc80` (A100), `cc90` (H100), or AMD `gfx*` |

### CPU-only build (no CUDA required)

If using NVHPC compiler but don't want GPU support (e.g., no CUDA installed):

```bash
pip install -e . --no-build-isolation
```

This uses multicore OpenMP only, no GPU flags.

### GPU build (V100 - default architecture)

```bash
pip install -e . --no-build-isolation -Csetup-args=-Dgpu_offload=true
```

### Specifying GPU architecture

Use the `-Dgpu_arch` flag to target different GPUs:

```bash
# V100 (default - cc70)
pip install -v -e . --no-build-isolation \
    -Csetup-args=-Dgpu_offload=true \
    -Csetup-args=-Dgpu_arch=cc70

# A100 (cc80)
pip install -v -e . --no-build-isolation \
    -Csetup-args=-Dgpu_offload=true \
    -Csetup-args=-Dgpu_arch=cc80

# H100 (cc90)
pip install -v -e . --no-build-isolation \
    -Csetup-args=-Dgpu_offload=true \
    -Csetup-args=-Dgpu_arch=cc90
```

### GPU build with GPU-aware MPI

```bash
# V100 (default - cc70)
pip install -v -e . --no-build-isolation \
    -Csetup-args=-Dgpu_offload=true \
    -Csetup-args=-Dgpu_arch=cc70 \
    -Csetup-args=-Dgpu_aware_mpi=true

# A100 (cc80)
pip install -v -e . --no-build-isolation \
    -Csetup-args=-Dgpu_offload=true \
    -Csetup-args=-Dgpu_arch=cc80 \
    -Csetup-args=-Dgpu_aware_mpi=true

# H100 (cc90)
pip install -v -e . --no-build-isolation \
    -Csetup-args=-Dgpu_offload=true \
    -Csetup-args=-Dgpu_arch=cc90 \
    -Csetup-args=-Dgpu_aware_mpi=true
```

### Supported GPU architectures

AMD GPUs have not been tested as of now (2026.02.05)

| GPU | Architecture | Flag |
|-----|--------------|------|
| NVIDIA V100 | cc70 | `-Dgpu_arch=cc70` (default) |
| NVIDIA A100 | cc80 | `-Dgpu_arch=cc80` |
| NVIDIA H100 | cc90 | `-Dgpu_arch=cc90` |
| AMD MI100 | gfx908 | `-Dgpu_arch=gfx908` (requires clang/AOMP) |
| AMD MI210/MI250 | gfx90a | `-Dgpu_arch=gfx90a` (requires clang/AOMP) |
| AMD MI300X | gfx942 | `-Dgpu_arch=gfx942` (requires clang/AOMP) |

### Runtime settings

To ensure GPU execution: `export OMP_TARGET_OFFLOAD=mandatory`

All unit tests pass, using `pytest --pyargs anuga`. It will be slow-ish because there's many of them.

### Multiprocessor Modes

| Mode | RK2 Loop | Use Case |
|------|----------|----------|
| `set_multiprocessor_mode(1)` | Python (`use_c_rk_loop=False`) | Debugging, flexibility, Python callbacks |
| `set_multiprocessor_mode(2)` | C (`use_c_rk_loop=True`) | Performance, GPU (data stays on device) |

Both modes use the same unified `core_kernels.c`. The only difference is:
- **Mode 1**: RK2 time-stepping orchestrated by Python (more flexible)
- **Mode 2**: RK2 time-stepping in C (faster, data stays on GPU)

*To enable the GPU code*: `domain.set_multiprocessor_mode(2)`

The GPU and CPU code share the same base kernels, you can find these in `anuga/shallow_water/gpu/core_kernels.*`, they
contain the main math and logic needed for the fluxes, extrapolation, etc. Data transfer kernels are defined in the GPU
specific files. 

The main difference between using the RK2 loop in Python or in C is that in C the GPU aware MPI implementation can do more
of its own work. I could not get mpi4py to work well with GPU aware MPI due to some pointer addressing issues.

### Runtime GPU Detection

The domain tracks whether GPU offload is actually active via `domain.gpu_offload_active`:
- Checks `OMP_TARGET_OFFLOAD` environment variable at runtime
- If `OMP_TARGET_OFFLOAD=disabled`, prints warning and suppresses GPU-specific messages

| Build Flag | Runtime Env | Result |
|------------|-------------|--------|
| `gpu_offload=true` | (default) | GPU offload active |
| `gpu_offload=true` | `OMP_TARGET_OFFLOAD=disabled` | CPU execution (warning shown) |
| `gpu_offload=false` | (any) | CPU execution (no warning) |

Example output when GPU offload is disabled:
```
+==============================================================================+
| WARNING: GPU mode enabled but OMP_TARGET_OFFLOAD=disabled                   |
| Running on CPUs with OMP_NUM_THREADS=4
+==============================================================================+
```

### Script Order (IMPORTANT!)

**Boundaries MUST be set BEFORE enabling GPU mode.** If you call `set_multiprocessor_mode(2)` before `set_boundary()`, you will get a RuntimeError:

```
RuntimeError: GPU mode requires boundaries to be set before calling set_multiprocessor_mode(2).
Please call domain.set_boundary({...}) BEFORE domain.set_multiprocessor_mode(2).
```

This is something I aim to fix by doing some initialization checks. TODO JORGE. 


Correct pattern:

```python
# 1. Create and distribute domain
domain = anuga.distribute(domain)

# 2. Set runtime parameters
domain.set_flow_algorithm('DE1')

# 3. Set boundaries FIRST!
Br = anuga.Reflective_boundary(domain)
Bt = anuga.Transmissive_n_momentum_zero_t_momentum_set_stage_boundary(domain, function=tide_function)
domain.set_boundary({'exterior': Br, 'open': Bt})

# 4. THEN enable GPU mode
domain.set_multiprocessor_mode(2)

# 5. Evolve
for t in domain.evolve(yieldstep=60, finaltime=3600):
    print(domain.timestepping_statistics())
```

You should see something like:
```
GPU halo exchange initialized:
  Neighbors: 2
  Total send: 7662 elements
  Total recv: 7634 elements
GPU domain initialized: 4869676 elements
GPU Domain Info (rank 0/4):
  Elements: 4869676
  GPU initialized: 0
  GPU-aware MPI: 0
  Halo neighbors: 2
  Total send: 7662, Total recv: 7634
+==============================================================================+
| GPU interface initialized using OpenMP target offloading                    |
+==============================================================================+
GPU arrays mapped to device 0
```

## Multi GPU execution

The code is GPU aware by design, `mpirun -np X python run_my_sim.py` will use X GPUs. If you ask for more processes than there are GPUs you will 
oversubscribe and that could lead to performance degradation. Be wary. I have not tested this. 

To run on multi node: `mpirun -np 64 --bind-to core --map-by ppr:4:node python run_model_3_partitioned.py` 

For example that line would run a 16 node, 64 GPU job for the partitioned scheme. You would need a prepartitioned mesh on 64 processes.


## GPU Gotchas & Troubleshooting

### Wrong Results Checklist

If GPU mode produces different results than CPU mode:

1. **Check boundary initialization**: Look for messages like `Reflective boundary arrays mapped to GPU: XX edges`. If you see `num_edges=0`, boundaries weren't initialized properly.

2. **Did you set boundaries first?**: GPU mode now requires `set_boundary()` before `set_multiprocessor_mode(2)`. If the order is wrong, you'll get a RuntimeError.

3. **Verify with OMP_TARGET_OFFLOAD=disabled**: This runs GPU code on CPU. If results are still wrong, it's a logic bug not a GPU memory issue.

### Performance Tips

1. **Don't mix modes**: Once in GPU mode, stay in GPU mode. Mode switching is expensive.

2. **Minimize yieldsteps**: GPU→CPU sync only happens at yieldsteps. Longer yieldsteps = fewer syncs.

3. **Use prepartitioned meshes**: Loading partitioned data is much faster than partitioning at runtime.

4. **Consolidate Rate_operators into a single Quantity** (see below).

### Rate_operator Optimization (CRITICAL for GPU)

**The Problem**: Many scripts create multiple Rate_operators for rainfall, one per polygon:

```python
# BAD for GPU - creates 26+ operators, each with a file_function
for filename in os.listdir(Rainfall_Gauge_directory):
    polygon = anuga.read_polygon(Gaugefile)
    rainfall = anuga.file_function(Rainfile, quantities='rate')
    op = Rate_operator(domain, rate=rainfall, polygon=polygon, ...)
```

Each operator has `rate_type='t'` (time-dependent scalar), causing:
- 26+ Python→GPU round trips per RK2 step
- 26+ GPU kernel launches
- 26+ `file_function(t)` evaluations

**The Solution**: Consolidate into a single Quantity-based Rate_operator:

```python
# GOOD for GPU - single operator with per-cell array
from anuga.geometry.polygon import inside_polygon

# Step 1: Pre-load all rainfall data and polygon masks
rainfall_data = []
for filename in os.listdir(Rainfall_Gauge_directory):
    polygon = anuga.read_polygon(Gaugefile)
    rainfall_func = anuga.file_function(Rainfile, quantities='rate')

    # Create mask for cells inside this polygon
    centroids = domain.get_centroid_coordinates(absolute=True)
    mask = inside_polygon(centroids, polygon)
    if len(mask) > 0:
        rainfall_data.append((mask, rainfall_func))

# Step 2: Create a Quantity to hold per-cell rainfall rates
rainfall_quantity = anuga.Quantity(domain, name='rainfall_rate')
rainfall_quantity.set_values(0.0)

# Step 3: Single Rate_operator using the Quantity
rainfall_operator = Rate_operator(domain, rate=rainfall_quantity, factor=1.0e-3,
                                  default_rate=0.0, label='Combined_Rainfall')

# Step 4: Function to update rainfall (call at yieldsteps)
_last_update_time = [-1.0]

def update_rainfall_quantity(t):
    if abs(t - _last_update_time[0]) < 60.0:  # Only update every 60s sim time
        return
    _last_update_time[0] = t

    rainfall_quantity.centroid_values[:] = 0.0
    for mask, rainfall_func in rainfall_data:
        try:
            rate = rainfall_func(t)
            if hasattr(rate, '__len__'):
                rate = rate[0]
            rainfall_quantity.centroid_values[mask] = rate
        except:
            pass

    # Signal GPU cache refresh
    rainfall_operator._gpu_rate_array_cache = None
    rainfall_operator._gpu_rate_changed = True

# Step 5: Call in evolve loop
for t in domain.evolve(yieldstep=60, finaltime=3600):
    update_rainfall_quantity(t)
    print(domain.timestepping_statistics())
```

**Result**:
- `rate_type='quantity'` uses the array GPU kernel
- Array cached on GPU, only transfers when `_gpu_rate_changed=True`
- 1 kernel launch instead of 26+ per RK2 step
- `file_function` evaluated at yieldsteps only, not every RK2 step

This optimization eliminated Rate_operator from the CPU profile entirely.

### Debugging

Enable verbose GPU output:
```python
domain.set_multiprocessor_mode(2)
# Check domain.gpu_interface.gpu_dom for internal state
```

Run with CPU fallback to compare:
```bash
OMP_TARGET_OFFLOAD=disabled python your_script.py
```

## GPU Feature Status

### Fully Supported on GPU

- **Unified kernel architecture**: CPU (mode 1) and GPU (mode 2) share ALL kernel code in `core_kernels.c`
- **Riverwall/weir support**: Weir discharge adjustments work on GPU (Villemonte submergence, blending)
- **Boundary flux tracking**: For `boundary_flux_integral_operator`
- **MPI ghost cell handling**: `tri_full_flag` checks for MPI domains
- **All friction models**: Manning friction (flat and sloped, semi-implicit)
- **Second-order extrapolation**: Edge-based gradient limiting

### Supported Boundary Types (GPU-native)

These boundary types run entirely on GPU (no CPU fallback):
- `Reflective_boundary`
- `Dirichlet_boundary`
- `Transmissive_boundary`
- `Transmissive_n_momentum_zero_t_momentum_set_stage_boundary`
- `Time_boundary`

If you use an unsupported boundary type, ALL boundaries fall back to CPU evaluation (with a warning).

### Known Issues

mpi4py failure due to lack of GPU awareness, fix by: `CFLAGS="-noswitcherror" CC=mpicc pip install mpi4py --no-cache-dir --no-binary mpi4py ` rebuilding mpi4py locally with a GPU aware MPI install

On Gadi sometimes it does not find the right python, gotta remove some pkgconfig:

```
export PKG_CONFIG_PATH=$(echo $PKG_CONFIG_PATH | tr ':' '\n' | grep -v '/half-root' | paste -sd:)
```
