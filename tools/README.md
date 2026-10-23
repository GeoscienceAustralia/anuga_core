# tools/

Helper scripts for installing ANUGA and a few developer utilities.

> **Authoritative install source:** the conda environment files in
> [`../environments/`](../environments/) and the build/test instructions in
> [`../CLAUDE.md`](../CLAUDE.md). The scripts here are conveniences and can drift;
> when in doubt, follow those.

## Recommended install

For most users/developers:

```bash
# 1. Miniforge + a conda env from environments/environment_<pyver>.yml
bash tools/install_miniforge.sh          # Linux/macOS
tools\install_miniforge_windows.bat      # Windows

# 2. Build ANUGA (editable)
conda activate anuga_env_3.14
pip install --no-build-isolation -e .
```

`install_miniforge*.{sh,bat}` are the only scripts here referenced by the docs
(`docs/source/installation/install_anuga_developers.rst`).

## GPU build

```bash
bash tools/install_anuga_nvc.sh          # NVIDIA HPC SDK (nvc) GPU build
```

Builds the GPU-offload extension (`-Dgpu_offload=true`) with `nvc`, runs the
isolated GPU test suite, and clears a stale build dir so the compiler switch
takes effect. Env vars: `PY` (default 3.12, same as `install_miniforge.sh`; ignored when a
conda env is already active), `GPU_ARCH` (default: **detected from this
machine's GPU** via `nvidia-smi`; a wrong value builds a package that
compiles but crashes on every kernel launch), `NVHPC_ROOT`. See `../CLAUDE.md` → "Testing a GPU-offload (nvc) build".

## Contents

| Script | Purpose |
|--------|---------|
| `install_miniforge.sh` / `install_miniforge_windows.bat` | **Recommended** — install Miniforge and create the conda env. Referenced by the docs. |
| `install_anuga_nvc.sh` | GPU-offload build via the NVIDIA HPC SDK (nvc). |
| `install_ubuntu.sh` | Install ANUGA into a Python venv on Ubuntu (22.04 / 24.04 / 26.04) via apt + pip. Detects the release and picks numpy/gdal pins. The conda path (`install_miniforge.sh`) is more robust. |
| `count_lines.py` | Report total lines of code in ANUGA (dev utility). |
| `clear_for_master.sh` | Remove generated `*.c`/`*.so` build artifacts (dev utility). |

### Legacy / unreferenced (candidates for removal)

Kept for now but referenced by nothing in CI or docs — older conda/apt install
variants superseded by `install_miniforge*` + `environments/*.yml`:

- `install_conda_ubuntu.sh`, `install_conda_ubuntu_22_04.sh`

## Notes

- Minimum Python is **3.10** (`pyproject.toml`: `requires-python = ">=3.10, <3.15"`).
- CI uses GitHub Actions (`../.github/workflows/`), not the scripts here.
