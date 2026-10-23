# Upgrading to ANUGA 4.0.0

Most scripts need **no changes**. The default compute path is unchanged and GPU
offload is opt-in. This guide covers the three things that can affect you.

---

## 1. `anuga.culvert_flows` has been removed

The package was superseded by the structure operators years ago.

```python
# before (4.0.0: ModuleNotFoundError)
from anuga.culvert_flows.culvert_class import Culvert_flow

# after
from anuga import Boyd_box_operator

Boyd_box_operator(domain,
                  end_points=[[x0, y0], [x1, y1]],
                  width=2.0, height=2.0,
                  losses=1.5, manning=0.013,
                  apron=0.5)
```

Pipes use `Boyd_pipe_operator(..., diameter=...)`; weirs use
`Weir_orifice_trapezoid_operator`. A worked example is in
`examples/structures/run_open_slot_wide_bridge.py`.

## 2. Forcing classes: deprecated in 4.0, removed in 4.1

In 4.0.x the classes `Rainfall`, `Inflow`, `Wind_stress`, `Barometric_pressure`
(and the `_fast` variants) in `anuga.shallow_water.forcing` still work but emit
`DeprecationWarning`, and they are **silently skipped in the `'unified'`
compute mode**. **In 4.1 they are gone**: importing them raises `ImportError`.
Migrate to the operators, which work in both compute modes:

```python
# before
from anuga.shallow_water.forcing import Rainfall, Inflow, Wind_stress
domain.forcing_terms.append(Rainfall(domain, rate=10.0))     # mm/s
domain.forcing_terms.append(Inflow(domain, rate=5.0))        # m^3/s
domain.forcing_terms.append(Wind_stress(s=10.0, phi=45.0))

# after
from anuga import Rate_operator, Wind_stress_operator
Rate_operator.rainfall(domain, rate=36000.0)   # mm/hr  (10 mm/s == 36000 mm/hr)
Rate_operator.inflow(domain, rate=5.0)         # m^3/s
Wind_stress_operator(domain, speed=10.0, phi=45.0)
```

Operators attach themselves to the domain when constructed — there is no
`forcing_terms.append(...)`.

**Watch the units.** `Rainfall` took **mm/s**; `Rate_operator.rainfall()` takes
**mm/hr** — multiply by 3600 when porting. `Inflow` and `Rate_operator.inflow()`
both take m³/s, so those carry over unchanged.
`Barometric_pressure` → `Barometric_pressure_operator`.

**Wind and pressure fields read from a file** (`Wind_stress(F,
use_coordinates=False)`, `Wind_stress_fast`, `Barometric_pressure_fast`) are
now handled by the operators with the same `use_coordinates=False` switch; the
`file_function` is interpolated at the centroids for wind and at the nodes for
pressure:

```python
F = anuga.file_function('wind.sww', domain,
                        quantities=['wind_speed', 'wind_angle'],
                        interpolation_points=domain.get_centroid_coordinates())
Wind_stress_operator(domain, F, use_coordinates=False)

P = anuga.file_function('pressure.sww', domain,
                        quantities=['barometric_pressure'],
                        interpolation_points=domain.get_nodes())
Barometric_pressure_operator(domain, P, use_coordinates=False)
```

## 3. Structures on sloping ground give different results

Structure operators used to write a uniform *depth* across each inlet, which
tilted a level water surface onto a sloping bed. They now level the *stage*.

* **Flat bed under the inlet** — results are bit-identical. Nothing to do.
* **Sloping bed** — results change, in proportion to the bed elevation range
  across the inlet. On Towradgi (median inlet bed spread 0.77 m) the peak stage
  difference was 0.32 m.

If you have a calibrated model with structures on sloping ground, re-run and
re-check it against your observations. The new behaviour is the correct one: a
lake at rest is now undisturbed by an idle structure, which was not true before.

To see how exposed a model is, measure the bed spread across its inlets:

```python
z = domain.quantities['elevation'].centroid_values
for op in domain.fractional_step_operators:
    for i, inlet in enumerate(getattr(op, 'inlets', []) or []):
        if inlet is None or not len(inlet.triangle_indices):
            continue
        zz = z[inlet.triangle_indices]
        print(f'{getattr(op, "label", type(op).__name__)} inlet {i}: '
              f'bed range {zz.max() - zz.min():.3f} m')
```

A range near zero means this change does not affect that structure.

---

## Optional: trying the unified compute mode

Nothing below is required to run 4.0.0.

```python
domain.set_compute_mode('unified')     # per domain
```

or process-wide:

```bash
export ANUGA_DEFAULT_COMPUTE_MODE=unified
```

**This is not a GPU switch.** `'unified'` selects the unified C kernels, which
run CPU-multicore on an ordinary build. Offloading to a GPU is a separate,
process-wide decision:

```python
anuga.set_gpu_offload(True)            # needs a build made with nvc
print(anuga.gpu_offload_supported())   # whether this build can offload at all
```

Caveats:

* The deprecated forcing classes are not applied in 'unified' (a warning is
  issued) — use the operators.
* `protect_against_isolated_degenerate_timesteps` is not implemented in
  'unified' (default-off; a warning is issued if enabled).
* Under MPI with GPU offload, the number of ranks must match the number of
  visible GPUs.

`domain.set_multiprocessor_mode(2)` remains as a thin wrapper over
`set_compute_mode('unified')`, but new code should prefer the named form.

Container images with a working GPU toolchain are published to GHCR — see
`docker/README.md`.
