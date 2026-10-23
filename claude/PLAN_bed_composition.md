# Plan: per-fraction bed composition (active layer)

Created 2026-10-06. Branch `feat/bed-composition` (worktree `~/anuga_core_bedcomp`,
conda env `anuga_bedcomp`, so the main checkout can keep running jobs).

## Why

The sediment bed is one elevation field over an optional floor `sediment_z_base`.
The erodible thickness `z - z_base` is shared by every fraction, and when it runs
short the [L-5] block scales all fractions' erosion by one factor
(`gpu/core_kernels.c` ~2139-2152; `docs/source/appendices/physics_spec.md`
1201-1208 already says "If bed stratigraphy is ever added, this rule is what
must change").

So any fraction can be entrained wherever the shear allows. In the Delta-X Wax
Lake runs (8 days, 497k cells, de Leeuw / Nghiem 2022) mud is pulled out of
deep, high-shear channel beds that are really sand: in-domain mud reaches
150-400 g/L (inflow 0.12 g/L), 2e6 m^3 is in suspension after a day and the
scoured channels drop the upstream stage 0.1-0.15 m. The coarse mesh did the
same in September, so it is the physics, not the mesh.

The Caltech Wax Lake sediment model (Wang, Salter & Lamb, ORNL DAAC 2309), which
our physics otherwise follows, has a 1 cm active layer with per-fraction
composition: E_k = F_k * E_k,pot. As mud leaves a channel bed the surface goes
sandy and mud entrainment shuts itself off (armouring). That is the missing piece.

## Reference: the Caltech model (`SSC_2D_Model_WLD_MainCode.m`)

- Active layer `L_active = 0.01` m over one finite substrate ("strat"),
  total erodible 1 m; `bedcomp(y,x,k)` (active fractions), `stratcomp`.
- Initial composition from a mud-fraction map (`3_ModelInput/bedcomp_GSDInterp.txt`,
  100 m grid, Fmud 0.17-0.96, Wax Lake delta only); `stratcomp = bedcomp`.
- E_k = min(L/0.01, 1) * F_k * E_k,pot (de Leeuw, Nghiem 2022 constants).
- Hirano exchange each step: aggradation overflows active -> strat at the active
  composition; degradation refills active from strat at the strat composition.
- Per-class erosion capped at that class's content of the column; deposition
  capped at h*c/dt.
- Mixture skin roughness k_s = 3 * 2^(mean psi + var psi).
- They store fractions (class mass not exactly conserved on aggradation) and
  have no porosity. We store amounts and keep porosity.

## Design (phase 1)

State per cell k and fraction s (solid volume per bed area, m):
`bed_active[s,k]`, `bed_sub[s,k]`; scalar `active_layer` H_A. Invariant:
`sum_s (bed_active + bed_sub) / (1 - porosity) == z - z_base`.

In `core_apply_sediment_source` (shared CPU / GPU kernel), only when
composition is on:

1. After E_s is computed (~line 2082): `E_s *= F_s`, F_s = bed_active[s] / sum bed_active
   (0 if the active layer is empty).
2. Replace the shared [L-5] scale by a per-fraction cap: erosion of s in a step
   <= (bed_active[s] + bed_sub[s]) / (dt * morfac) (solid volume rate).
3. In the apply loop (~2155-2180): update bed_active[s] by the net solid
   exchange, then restore the active layer to H_A: overflow to the substrate at
   the active composition when aggrading, refill from the substrate at the
   substrate composition when degrading (all amounts, so exactly conservative).
4. Everything is cell-local: ghosts recompute it identically, no new halo.

Python: `Domain.set_bed_composition(fractions, active_layer=0.01, substrate=None)`
(fractions: {name: value | array | f(x, y)} as volume fractions summing to 1;
substrate defaults to the same). Needs `set_erodible_base` for the thickness.
Not called -> behaviour exactly as now.

Plumbing (the `sediment_source_limited` pattern):
- `sw_domain.h` sediment block (fields; mind struct-offset aliasing, HANDOVER 2.1)
- `sw_domain_openmp_ext.pyx` and `sw_domain_gpu_ext.pyx` externs + pointer binding
- `gpu/gpu_domain_core.c`: map / unmap (mirror exactly) / sync_from_device
- `sediment_summary`, `erodible_thickness`; optional sww output of the active mud
  fraction; checkpoint pickles the host arrays once synced
- invalidation convention: reset `_Domain_C_struct`, `gpu_interface`

Phase 1 limits (warn): bedload (`core_apply_bedload`, summed over classes) and
repose (`core_apply_repose`) still move bulk bed without composition.

Tests:
- one fraction with composition == today's results
- per-class mass conservation in a closed box (suspended + bed)
- no erosion of a fraction absent from the bed
- armouring: a mixed bed under steady shear loses its mud and mud entrainment
  decays to ~0
- CPU (legacy) vs GPU (unified) parity; serial vs MPI
- docs: physics_spec §5 (an unimplemented RDycore multi-layer design is there)

## Phase 2

Bedload and repose carry composition; mixture skin roughness; possibly several
substrate layers (physics_spec §5 design).

## Applying it to WLAD

A whole-domain initial mud-fraction map: the Caltech map where it covers the
Wax Lake delta, a rule elsewhere (e.g. sandy below a channel depth, muddy marsh),
then rerun the 8-day sediment case against the Caltech results.
