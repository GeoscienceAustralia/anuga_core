# Riverwall Throughflow — Design and Implementation Plan

Written 2026-04-01.

---

## Motivation

The current riverwall models flow **over** the wall crest using the Villemonte
(1947) weir formula. There is no mechanism for flow **through** the body of the
wall. This is needed to represent:

- Levees or embankments with a culvert or pipe at the base
- Porous embankments (seepage)
- Walls with a partial gap or gate below the crest

The new feature allows a separate discharge component to pass through the wall,
driven by the stage difference between the two sides and limited by the depth
of water available on each side below the wall crest.

---

## Physics

### Overtopping flow (existing)

The existing Villemonte weir formula computes discharge using heads measured
**above the wall crest** (`h_left = stage_left - z_wall`, `h_right = stage_right - z_wall`).
It only activates when water overtops the wall.

### Throughflow (new)

Throughflow passes through the wall body at or below the crest. The appropriate
model is **submerged orifice flow**:

```
Q_through / unit_width = Cd_through * h_eff * sqrt(2 * g * |Δstage|) * sign(Δstage)
```

where:

```
Δstage      = stage_left - stage_right          [driving head]

h_left_sub  = max(min(stage_left,  z_wall) - z_bed_left,  0)   [depth below wall on left]
h_right_sub = max(min(stage_right, z_wall) - z_bed_right, 0)   [depth below wall on right]
h_eff       = h_left_sub  if stage_left >= stage_right          [upstream submerged depth]
            = h_right_sub otherwise
```

`h_eff` is the submerged depth on the **upstream** (driving) side of the wall —
the depth of water pressing through the wall from the high-stage side. This naturally:

- Gives positive throughflow when the downstream side is dry (the common case
  of a wet river against a dry floodplain): h_eff = h_upstream > 0 even when
  the downstream side has no water
- Scales with the upstream water depth (deeper upstream → larger effective
  seepage/orifice cross-section)
- Works equally whether water is above or below the crest on either side
- Reduces to zero when Cd_through = 0 (backward compatible default)

Note: using `min(h_left_sub, h_right_sub)` was considered but rejected — it
gives zero flow when the downstream side is dry, which is physically wrong for
both culvert and seepage models. Using `max` is equivalent to using the upstream
depth when bed levels are similar, but `h_upstream_sub` is more physically
precise for unequal bed elevations on the two sides.

### Interaction with overtopping flow

The two components are **additive**:

```
total_flux = flux_over_wall  (Villemonte, existing)
           + flux_through_wall  (orifice, new)
```

The existing code scales `edgeflux` to match the Villemonte discharge. The
throughflow term is then **added** to `edgeflux[0]` (mass flux) with a
corresponding update to the momentum fluxes. This preserves full backward
compatibility when `Cd_through = 0`.

### Momentum for the throughflow component

The throughflow discharge `Q_through` is a mass flux per unit edge length.
Momentum is carried with it at the velocity implied by orifice flow:

```
u_through = sqrt(2 * g * |Δstage|)   [orifice exit velocity]
```

The momentum contribution is small relative to the mass flux term and can be
approximated as zero for the first implementation (set `edgeflux[1]` and
`edgeflux[2]` contributions from throughflow to zero). A refined version could
add `Q_through * u_through` in the normal direction. Start simple.

---

## New Parameters

Following the existing constraint — **append to the end** of
`hydraulic_variable_names`, never reorder — add:

```python
# riverwall.py
hydraulic_variable_names = ('Qfactor', 's1', 's2', 'h1', 'h2', 'Cd_through')
```

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `Cd_through` | `0.0` | Discharge coefficient for flow through the wall. `0` = impermeable (current behaviour). Typical values for a culvert opening: 0.5–0.8. |

A single parameter keeps the interface simple. Future extensions could add
`through_depth` (cap the effective opening height) or `through_width_fraction`
(fraction of wall length that has an opening), but one parameter is enough for
a first implementation.

---

## Files to Change

### 1. `anuga/structures/riverwall.py`

**a) Add `Cd_through` to default parameters and variable names:**

```python
self.default_riverwallPar = {
    'Qfactor':    1.0,
    's1':         0.9,
    's2':         0.95,
    'h1':         1.0,
    'h2':         1.5,
    'Cd_through': 0.0,    # NEW — impermeable by default
}

self.hydraulic_variable_names = ('Qfactor', 's1', 's2', 'h1', 'h2', 'Cd_through')
```

**b) `ncol_hydraulic_properties` is computed as `len(hydraulic_variable_names)`** — no change needed.

**c) Update docstring** to document the new parameter and the physics.

### 2. `anuga/shallow_water/gpu/gpu_device_helpers.h`

Add a new device function (keeping `gpu_adjust_edgeflux_with_weir` unchanged
for clarity):

```c
// Apply throughflow (orifice) discharge through the wall body.
// Called AFTER gpu_adjust_edgeflux_with_weir.
// Adds Q_through to edgeflux[0]; momentum contribution set to zero (first implementation).
//
// Parameters:
//   edgeflux: [3] mass+momentum fluxes (modified in place)
//   stage_left, stage_right: cell-centroid stage values
//   z_bed_left, z_bed_right: bed elevations at left/right centroids
//   z_wall: riverwall crest elevation
//   g: gravitational acceleration
//   Cd_through: discharge coefficient (0 = no throughflow)
//   max_speed_local: updated if throughflow is significant
static inline void gpu_adjust_edgeflux_with_throughflow(
    double * restrict edgeflux,
    double stage_left, double stage_right,
    double z_bed_left, double z_bed_right,
    double z_wall, double g,
    double Cd_through,
    double * restrict max_speed_local) {

    if (Cd_through <= 0.0) return;

    // Depth of water on each side below the wall crest
    double h_left_sub  = fmax(fmin(stage_left,  z_wall) - z_bed_left,  0.0);
    double h_right_sub = fmax(fmin(stage_right, z_wall) - z_bed_right, 0.0);
    // h_eff = upstream submerged depth (driving side).
    // Using min() was rejected: gives zero flow when downstream is dry, which
    // is wrong for both culvert and seepage models.
    double h_eff = (stage_left >= stage_right) ? h_left_sub : h_right_sub;

    if (h_eff <= 0.0) return;

    double delta_stage = stage_left - stage_right;
    double abs_delta   = fabs(delta_stage);
    if (abs_delta <= 0.0) return;

    double Q_through = Cd_through * h_eff * sqrt(2.0 * g * abs_delta);
    if (stage_right > stage_left) Q_through = -Q_through;

    edgeflux[0] += Q_through;
    // Momentum: zero for first implementation

    // Update max_speed estimate for timestep
    double speed_est = Cd_through * sqrt(2.0 * g * abs_delta);
    if (speed_est > *max_speed_local) {
        *max_speed_local = speed_est;
    }
}
```

### 3. `anuga/shallow_water/gpu/core_kernels.c`

In the riverwall block (around line 871), after the existing
`gpu_adjust_edgeflux_with_weir` call, add:

```c
// Read the new Cd_through (column index 5)
double Cd_through = (ncol_riverwall_hp > 5)
    ? riverwall_hydraulic_properties[ii + 5]
    : 0.0;

// Apply throughflow (additive to the overtopping flux)
gpu_adjust_edgeflux_with_throughflow(
    edgeflux,
    stage_cv[k], stage_cv[neighbour],
    zl, zr,
    zwall, g, Cd_through, &max_speed_local);
```

The `ncol_riverwall_hp > 5` guard ensures old SWW files that don't carry the
6th column still work correctly.

### 4. CPU Cython path

The CPU flux computation uses the equivalent C code (not the GPU path). Locate
the `adjust_edgeflux_with_weir` call in:

- `anuga/shallow_water/sw_domain.pyx` or
- `anuga/shallow_water/sw_domain_ext.c`

Apply the same additive throughflow formula immediately after the existing weir
call. The logic is identical — no GPU-specific code involved.

### 5. `anuga/shallow_water/sw_domain_gpu_ext.pyx` (GPU Cython bridge)

The `riverwall_hydraulic_properties` array is already passed through to C via
the domain struct (`D->riverwall_hydraulic_properties`). The new column is
automatically included because the array is sized by `ncol_hydraulic_properties`.
No change needed in the Cython bridge — the C code reads column 5 directly.

---

## Test Plan

### Unit test — throughflow formula

In `anuga/structures/tests/test_riverwall*.py` (or a new file):

```python
def test_throughflow_zero_by_default(self):
    """Cd_through=0 (default) produces no change vs current behaviour."""
    # Create domain with riverwall, no Cd_through specified
    # Compare evolve output to reference (should match existing behaviour exactly)

def test_throughflow_direction(self):
    """Throughflow flows from high stage to low stage."""
    # Set Cd_through > 0, stage_left > stage_right < z_wall
    # Verify edgeflux[0] > 0 (left to right)

def test_throughflow_dry_side(self):
    """Zero throughflow when either side is dry below the wall."""
    # h_left_sub = 0 → Q_through = 0

def test_throughflow_plus_overtopping(self):
    """Both components active when water overtops on one side."""
    # stage_left > z_wall > stage_right > z_bed
    # Verify total flux > Villemonte-only flux

def test_throughflow_additive(self):
    """Throughflow adds to (not replaces) the overtopping flux."""
    # Run with and without Cd_through; check difference equals expected Q_through

def test_throughflow_evolve(self):
    """End-to-end: water equilibrates faster with throughflow enabled."""
    # Two basins separated by a wall below both water levels
    # Cd_through > 0 → stages converge; Cd_through = 0 → stages stay separated
    # (mark @pytest.mark.slow if runtime > 5 s)
```

### Regression test — backward compatibility

Existing riverwall tests must pass without modification (Cd_through defaults to
0, so no change in behaviour).

---

## Implementation Order

1. **riverwall.py** — add parameter (30 min, no C changes, immediately testable)
2. **gpu_device_helpers.h** — add `gpu_adjust_edgeflux_with_throughflow` (1 h)
3. **core_kernels.c** — call the new function (30 min)
4. **CPU Cython path** — mirror the same logic (1 h)
5. **Tests** — unit + regression (half day)
6. **Docstring / user docs update** (1 h)

Total estimated effort: **1–2 days**.

---

## Open Questions

1. **Momentum flux for throughflow**: Setting momentum to zero is conservative
   (safe, stable). A more accurate treatment would add `Q_through * u_through`
   in the edge-normal direction. Defer to after validation.

2. **`through_depth` cap**: If the user wants to model a specific culvert
   diameter `d`, they currently cannot — `h_eff` uses the full submerged depth.
   Adding an optional `through_depth` parameter (cap on `h_eff`) is a natural
   follow-on, appended as column 6 in `hydraulic_variable_names`.

3. **Interaction near the crest**: When `stage ≈ z_wall` on both sides, both
   the overtopping and throughflow terms are small and approaching zero
   simultaneously — this is the correct limiting behaviour and requires no
   special treatment.

4. **MPI / parallel**: No parallel-specific changes needed. The throughflow
   calculation is purely local (per-edge, uses centroid values that are already
   available after halo exchange).

5. **GPU path**: The GPU path uses `gpu_device_helpers.h` which is already
   compiled with `#pragma omp declare target` — the new function goes in the
   same block and is automatically available on device. No extra mapping needed.
