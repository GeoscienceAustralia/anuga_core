# Sediment transport examples

Ported from the examples shipped with **anugaSed** (Mariela Perignon, 2016),
which is MIT licensed. `topo.asc`, `topo.prj` and `outline.csv` are their data
files, copied unchanged.

The originals are Python 2 and drive `Sed_transport_operator` through an
`evolved_quantities` list containing `'concentration'`. These ports keep their
geometry, boundary conditions, durations and parameters, and drive our API
instead: `domain.add_sediment_fraction()`, which registers the fractional-step
`Sediment_transport_operator` for you.

## Bed material: both routes are available

anugaSed erodes with the **cohesive** Hanson & Simon route `[E-3]`
(`tau_crit = 0.088` Pa, dimensional). The default here is the **non-cohesive**
Shields route `[E-1]`, the bed material this work targets. Spec 4.1.1 is
explicit that these describe *different sediment*, not competing formulations
of the same physics.

`example_1_channel.py --bed cohesive` runs anugaSed's own configuration.
The difference is not a detail:

| `--bed` | erosion law | bed change over 30 s |
|---------|-------------|----------------------|
| `noncohesive` (default) | `[E-1]` Shields, `tau_c*` = 0.04 | **erodes** 3.1 to 5.8 cm, all 16 cells |
| `cohesive` | `[E-3]` Hanson & Simon, `tau_c` = 0.088 Pa | **accretes** 0.6 to 1.3 mm, all 16 cells |

The sign reverses. A cohesive bed resists erosion strongly enough
(`K_e = 6.74e-7 m3/N/s`) that the sediment-laden inflow deposits instead of
scouring. This is what spec 4.1.1 means by "a physics error, not a tuning
error".

## How far the anugaSed comparison goes

`test_cohesive.py` checks our `[E-3]` against their `erosion()` transcribed
from `sed_transport_operator.py`, and reproduces their `edot` **exactly**
(difference 0.000e+00) across `tau_b` in [0, 20] Pa, with the same
`K_e = 6.742e-7`.

It stops there, deliberately. anugaSed cannot be *run* here -- it is Python 2
and drives a `Sed_transport_operator` this ANUGA does not have -- and a
whole-model comparison would differ anyway for a reason unrelated to erosion:
they compute bed shear from the depth-slope closure `[T-7]`, we use quadratic
drag `[T-1]` (divergence **D1**). Equal erosion laws fed unequal `tau_b` give
unequal answers.

Their other defaults are carried over: porosity 0.3, rho_s 2650, rho_w 1000,
D50 = 65 um (Griffin et al. 2014).

| file | ported from |
|------|-------------|
| `example_1_channel.py` | `run_simple_sed_transport.py` |
| `example_2_raster.py`  | `run_raster_sed_transport.py` |

`example_3_erodible_dambreak.py` and `example_4_dune_collapse.py` are not ports;
see below.

`run_simple_veg.py` is not ported: vegetation drag (spec 8) is Phase 5 and is
not implemented.

## Example 4 -- what slope collapse changes

`example_4_dune_collapse.py` runs the same dune-overtopping event twice, once
with angle-of-repose relaxation off and once with it on, because either run
alone says very little. The comparison is the example.

```
                                       repose off      repose on
steepest bed slope (deg)                    38.38          33.03
surviving crest (m)                         1.780          1.780
budget: suspended + bed (m3)            -2.84e-14      -4.26e-14
```

The dune starts at 27.8 degrees, below the 33 degree limit, so everything
relaxation does here was caused by the scour rather than by an unreasonable
initial condition. Without it, scour leaves a 38.4 degree face -- steeper than
sand stands. With it, the bed is held at the limit.

The crest is *not* where the difference shows: scour takes the top of the dune
either way. The difference is on the face, where relaxation removes a further
0.052 m and lays it down 0.051 m thick just downslope. And it is moved, not
lost -- the budget closes to machine precision with relaxation on exactly as it
does with it off, which is the property that separates this from
`sanddune_erosion_operator`, where the collapsed material simply disappears.

## Example 3 is not a port

Because examples 1 and 2 are faithful ports, they only exercise what anugaSed
does: a single grain size, suspension only, compute mode 1. Four of the
things this module adds were therefore not demonstrated anywhere runnable.
`example_3_erodible_dambreak.py` covers them -- two grain sizes at once,
bedload beside suspension, a bed that moves, and the unified GPU path:

```
$ python example_3_erodible_dambreak.py
```

It prints `sediment_summary()` first, so the run begins with a full statement
of what was configured, then tracks the sediment budget as it evolves. The
budget is the point of the example. Boundaries are reflective, so nothing
leaves the domain, and over 60 s:

```
Sediment budget:  suspended +5.752703e+03  +  bed -5.752703e+03  =  -9.09e-13 m3
```

Every cubic metre in suspension came out of the bed, and `(1-lambda) dz`
accounts for it to machine precision -- across both grain sizes, with erosion,
deposition and bedload all running and the bed moving under them. Bed change
ranges from 0.43 m of scour to 0.11 m of deposition, so this is not a
small-perturbation result.

The two grain sizes separate as they should: fine sand (150 um, `v_s` = 1.7 cm/s)
climbs monotonically to 5.7e3 m3 in suspension, while coarse sand (800 um,
`v_s` = 15 cm/s) peaks near 1.1e2 m3 at t = 30 s and then falls to 6.0e1 as
deposition overtakes entrainment.

Mode 2 selects the unified code path; whether it offloads is a property of
the build, so the example reports `anuga.gpu_offload_enabled()` rather than
inferring a device from the mode. A build without offload runs the same
kernels on the host and the example still works.

## What these examples do and do not cover

Runnable demonstration is not the same as coverage. This table is what the
three examples actually exercise; everything else is covered only by the test
suites in `anuga/shallow_water/tests/test_sediment_*.py`.

| capability | spec | ex 1 | ex 2 | ex 3 | ex 4 | tests |
|---|---|:--:|:--:|:--:|:--:|---|
| suspended transport | `[G-3]` | yes | yes | yes | yes | `test_tracer_ns1`, `test_ns2` |
| multiple grain sizes | 2.2 | -- | -- | **yes** | -- | `test_add_tracer` |
| Shields erosion, non-cohesive | `[E-1]` | yes | yes | yes | yes | `test_sediment_source` |
| Hanson & Simon, cohesive | `[E-3]` | `--bed cohesive` | -- | -- | -- | `test_sediment_cohesive` |
| Partheniades | `[E-4]` | -- | -- | -- | -- | `test_sediment_cohesive` |
| deposition | `[D-1]` | yes | yes | yes | yes | `test_sediment_source` |
| Rouse near-bed ratio | `[S-4]` | -- | -- | **yes** | -- | `test_sediment_rouse` |
| Exner, suspended exchange | `[G-4]` | yes | yes | yes | yes | `test_sediment_exner` |
| bedload and its bed term | `[K-1]`, `[G-5]` | -- | -- | **yes** | -- | `test_sediment_bedload` |
| quadratic drag | `[T-1]` | yes | yes | yes | yes | `test_sediment_shear_closure` |
| depth-slope closure | `[T-7]` | -- | -- | -- | -- | `test_sediment_shear_closure` |
| Wilson / Larsen-Lamb friction | `[T-8..15]` | -- | -- | -- | -- | `test_sediment_friction` |
| limiters | `[L-1..4]` | yes | yes | yes | yes | `test_sediment_source` |
| non-erodible base | `[L-5]` | -- | -- | -- | -- | `test_sediment_erodible_base` |
| erodible region | `[L-5]` | -- | -- | -- | -- | `test_sediment_erodible_base` |
| angle-of-repose relaxation | 7 | -- | -- | -- | **yes** | `test_sediment_repose` |
| external tracer source | 2.6 | -- | -- | -- | -- | `test_sediment_mms_full` |
| compute mode 1 | -- | yes | yes | -- | yes | all |
| compute mode 2 (GPU) | -- | -- | -- | **yes** | -- | `test_mode1_vs_mode2` |
| parallel halo exchange | -- | -- | -- | -- | -- | `run_parallel_tracer.py` |

Gaps are visible in that table and are worth naming rather than leaving to be
noticed: the depth-slope closure `[T-7]`, the Wilson/Larsen-Lamb friction
closures `[T-8..15]`, the non-erodible base and region restriction `[L-5]`, and
the external tracer source all have tests but appear in no example, and nothing
here runs in parallel. None is a correctness concern; they are simply
undemonstrated.

## A note on the data files

`topo.asc` (594 KB) and `outline.csv` are force-added: the repository's root
`.gitignore` excludes `*.asc` and `*.csv` as build artefacts, and without `-f`
example 2 would have been committed without the data it needs.

## Status

These demonstrate the module running end to end. They are **not validation**:
nothing here has been checked against measured data or against a published
result. Verification -- that the equations are solved correctly -- lives in
`anuga/shallow_water/tests/test_sediment_*.py` (MMS convergence, the RDycore
benchmarks, mode 1 against
mode 2), and is a separate claim from validation. The validation rungs of spec 10 -- Rio Puerco (rung 8) and the crater
breach (rung 7) -- have not been attempted.
