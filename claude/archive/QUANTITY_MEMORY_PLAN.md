# Quantity Memory Reduction — Design Plan (v4.0.0)

Written 2026-04-02.

---

## Current Situation

Every `Quantity` unconditionally allocates 9 arrays in `__init__`:

| Array | Shape | Bytes |
|-------|-------|-------|
| `vertex_values` | (N, 3) | 24N |
| `centroid_values` | (N,) | 8N |
| `edge_values` | (N, 3) | 24N |
| `x_gradient` | (N,) | 8N |
| `y_gradient` | (N,) | 8N |
| `phi` | (N,) | 8N |
| `explicit_update` | (N,) | 8N |
| `semi_implicit_update` | (N,) | 8N |
| `centroid_backup_values` | (N,) | 8N |
| **Total** | | **104N bytes** |

For a typical 8-quantity domain with 1M triangles: **832 MB**, all unconditionally allocated.

---

## v4.0.0 Change: Elevation is Centroid-Primary

In v4.0.0 elevation is defined by centroid values, not vertex values.
- Vertex values of elevation can now be lazy (computed on demand for output).
- Edge values of elevation are computed once at initialisation from centroids
  via the gradient extrapolation — after that elevation is static and the
  gradient arrays are no longer needed.

This changes the elevation category from "static vertex-primary" to
"static centroid-primary", enabling significantly more savings than pre-4.0.0.

---

## What Each Quantity Type Actually Needs

| Array | Evolved (stage/xmom/ymom) | Height | Elevation | Friction | Velocity (xvel/yvel) |
|-------|:---:|:---:|:---:|:---:|:---:|
| `centroid_values` | ✅ | ✅ | ✅ primary | ✅ | ✅ |
| `edge_values` | ✅ | ✅ | ✅ (set as stage_e − height_e) | ❌ | ❌ |
| `vertex_values` | 🔶 lazy | 🔶 lazy | 🔶 lazy | 🔶 lazy | 🔶 lazy |
| `x_gradient` | ✅ | ✅ | ❌ | ❌ | ❌ |
| `y_gradient` | ✅ | ✅ | ❌ | ❌ | ❌ |
| `phi` | ✅ | ✅ | ❌ | ❌ | ❌ |
| `explicit_update` | ✅ | ❌ | ❌ | ❌ | ❌ |
| `semi_implicit_update` | ✅ xmom/ymom; ❌ stage | ❌ | ❌ | ❌ | ❌ |
| `centroid_backup_values` | ✅ | ❌ | ❌ | ❌ | ❌ |
| **Minimum bytes** | **80N** | **56N** | **32N** | **8N** | **8N** |

Key reasoning per type:

### Evolved: stage, xmomentum, ymomentum (80N each)
- Primary state is centroid values; edge values computed by extrapolation each step.
- `vertex_values` are only needed to write `.sww` output — lazy allocation.
- `stage.semi_implicit_update` is always zero (Manning acts on momentum only) —
  can omit for stage, saving 8N bytes (optional micro-optimisation).

### Height (56N)
- Centroid values computed as `stage_c − elevation_c` each timestep.
- Edge values are extrapolated from height centroids via gradient + limiter —
  height IS actively extrapolated, so `x_gradient`, `y_gradient`, and `phi` are needed.
- No update arrays or backup (height is not a conserved quantity; it is
  recomputed from stage and elevation after each RK2 step, not backed up).
- `vertex_values` — lazy.
- Minimum: `centroid` (8N) + `edge` (24N) + `x_grad` (8N) + `y_grad` (8N) + `phi` (8N) = **56N**.

### Elevation (32N — v4.0.0)
- Centroid values are the primary input.
- Edge values are set during the extrapolate phase as `stage_edge − height_edge`
  — elevation is never independently extrapolated, so it needs **no gradient arrays**.
- `vertex_values` — lazy.
- No update arrays, no phi, no backup.
- Minimum: `centroid_values` (8N) + `edge_values` (24N) = **32N**.

### Friction (8N)
- Manning formula reads `friction.centroid_values` to compute the coefficient,
  then writes the result into `xmomentum.semi_implicit_update` and
  `ymomentum.semi_implicit_update` — never into `friction.semi_implicit_update`.
- No edge values, no gradients, no phi, no update arrays of any kind, no backup.
- Minimum: `centroid_values` (8N) = **8N**.

### xvelocity, yvelocity (8N each)
- Computed diagnostically from momentum / height at centroid level.
- Edge values not needed for flux computation (flux kernel uses momentum edges directly).
- `vertex_values` — lazy.
- Minimum: `centroid_values` (8N) = **8N**.

---

## Memory Savings (1M triangles)

### Current
| Quantities | Count | Per-quantity | Total |
|------------|-------|-------------|-------|
| stage, xmom, ymom | 3 | 104 MB | 312 MB |
| elevation | 1 | 104 MB | 104 MB |
| friction | 1 | 104 MB | 104 MB |
| height, xvel, yvel | 3 | 104 MB | 312 MB |
| **Total** | **8** | | **832 MB** |

### Optimised
| Quantities | Count | Per-quantity | Total |
|------------|-------|-------------|-------|
| stage, xmom, ymom | 3 | 80 MB | 240 MB |
| height | 1 | 56 MB | 56 MB |
| elevation | 1 | 32 MB | 32 MB |
| friction | 1 | 8 MB | 8 MB |
| xvel, yvel | 2 | 8 MB | 16 MB |
| **Total** | **8** | | **352 MB** |

**Saving: 480 MB — 58% reduction** for a typical 1M-triangle domain.

For a 10M triangle domain (SC26-scale): saving ~4.8 GB.

### Additional: Shared Gradient Workspace (optional)
`x_gradient` and `y_gradient` are recomputed from scratch every timestep
(C kernel `_compute_gradients()` overwrites them). They serve purely as
workspace between gradient computation and edge extrapolation within one call.
Replacing per-quantity storage with two domain-level scratch arrays would save
a further 3 × 16N = **48 MB** for evolved quantities — but requires changing
the C extension call signatures.

---

## Implementation Strategy

### Phase 1 — No C changes required (easy, ~2 days)

**1a. Strip unused arrays from non-evolved quantities**

Introduce an `_allocate_arrays(quantity_type)` helper called from `__init__`:

```python
# quantity.py
_EVOLVED    = 'evolved'    # stage, xmom, ymom
_STATIC     = 'static'     # elevation (centroid-primary in v4.0.0)
_FORCING    = 'forcing'    # friction
_DIAGNOSTIC = 'diagnostic' # height, xvel, yvel
```

Each type allocates only what it needs (see table above). The domain
constructor passes `quantity_type` when creating each quantity.

**1b. Lazy `vertex_values` via property**

Replace the unconditional `vertex_values` allocation with a property:

```python
@property
def vertex_values(self):
    if self._vertex_values is None:
        self._vertex_values = num.zeros((self.N, 3), float)
    return self._vertex_values

@vertex_values.setter
def vertex_values(self, value):
    self._vertex_values = num.array(value, float)
```

First access allocates; no code changes needed at call sites. Existing
`q.vertex_values = ...` assignment still works via the setter.

**1c. Elevation gradient freed after init**

After `domain.set_quantity('elevation', ...)` finalises and edge values are
computed, set:
```python
domain.quantities['elevation'].x_gradient = None
domain.quantities['elevation'].y_gradient = None
```

The gradient arrays will be garbage-collected. Any subsequent access raises
`AttributeError` which would catch misuse.

### Phase 2 — Shared gradient workspace (medium, ~1 week)

Replace per-quantity `x_gradient`/`y_gradient` with domain-level scratch:

```python
# generic_domain.py
self._gradient_workspace_x = num.zeros(N, float)
self._gradient_workspace_y = num.zeros(N, float)
self._phi_workspace = num.zeros(N, float)
```

The C extension `extrapolate_second_order_and_limit_by_edge()` receives these
as explicit arguments instead of reading from `quantity.x_gradient`. Evolved
quantities no longer store gradient arrays.

Additional saving: 3 × 16N = 48 MB (1M triangles) + 3 × 8N = 24 MB (phi) = **72 MB**.

---

## Risk Register

| Risk | Severity | Mitigation |
|------|----------|------------|
| User code accessing `q.vertex_values` directly on a non-output step causes unexpected allocation | Low | Lazy property is transparent; allocation just happens on first access |
| `elevation.x_gradient = None` breaks code that reads gradients after init | Medium | Grep all access sites; none expected in core — gradients only used inside `extrapolate_*` |
| Diagnostic quantities missing `edge_values` breaks `get_values(location='edges')` | Medium | Compute on-demand in `get_values()` from centroid interpolation |
| `stage.semi_implicit_update` omission breaks custom forcing terms that write to it | Low | Keep for stage by default; document that stage semi-implicit is unused in core |
| Phase 2 C signature change breaks GPU/OpenMP path | Medium | GPU path reads quantity arrays via `sw_domain.h` struct pointers — update struct to use workspace pointers |
| Parallel domain: ghost quantities also allocate full arrays | Low | Parallel quantities inherit from Quantity — same change applies automatically |

---

## Testing

- Existing test suite must pass unchanged (backward compatibility via lazy property).
- Add memory-reporting test: for N=10000, assert each quantity type uses ≤ expected bytes.
- Add regression test: 10 s evolution with optimised quantities matches reference centroid values.
- Test `get_values(location='vertices')` and `get_values(location='edges')` on diagnostic quantities post-optimisation.

---

## Implementation Order

| Step | Change | Saving | Effort |
|------|--------|--------|--------|
| 1a | Strip update arrays from elevation, friction, diagnostic | 376 MB | 1 day |
| 1b | Lazy `vertex_values` property (all types) | 72 MB | 1 day |
| 1c | Free elevation gradients after init | 16 MB | half day |
| 2 | Shared gradient workspace (C changes) | 72 MB | 1 week |
| **Total** | | **~536 MB** | **~2 weeks** |
