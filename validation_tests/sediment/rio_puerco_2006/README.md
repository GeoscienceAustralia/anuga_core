# Lower Rio Puerco, 2006 flood (planned case)

A field-scale sediment validation case that has not been built yet. There is no
`validate_*.py` here, so the validation runner skips this directory. These notes
were section 10.1 of the sediment physics specification
(`docs/source/appendices/physics_spec.md`, "validation rung 8"). They hold what is
needed to build the case.

`[P13]` is Perignon, M.C., Tucker, G.E., Griffin, E.R. & Friedman, J.M. (2013),
JGR Earth Surface 118(3), 1193–1209, doi:10.1002/jgrf.20073. `P14` and `aS16` are
the Perignon thesis and the anugaSed code, as listed in the specification's sources.


`[P13]` specifies this case completely enough to build it. Values are theirs.

| Quantity | Value | Source |
|----------|-------|--------|
| Peak discharge `Q` (Highway 6 bridge) | 625 m³ s⁻¹ | P13 §[55] |
| Mean active flow width `b` | 235 m | P13 §[55] |
| Unit discharge `q` | 2.66 m² s⁻¹ | P13 §[55] |
| Peak duration `t` | 50 400 s (14 h, Bernardo NM gauge) | P13 §[54] |
| Representative grain size | **0.045 mm** (coarse silt) | P13 §[54], after Nordin (1963) |
| Settling velocity `v_s` | 0.00175 m s⁻¹ | P13 §[54], Ferguson & Church (2004) |
| Bed porosity `φ` | 0.30 | P13 §[54], Beard & Weyl (1973) |
| Profile factor `d*` (`p`) | 1 (uniform) | P13 §[54], Nordin (1963) |
| Study reach | 12 km | P13 abstract |
| Deposit thickness (field) | 1–5 cm unconsolidated sand | P13 §[73] |
| Field control | 26 pits, measured 2006 deposit thickness | P13 Figure 4 |

**Volumetric budget** (P13 §[68]–[69], with Vincent et al. 2009):

| Component | Volume |
|-----------|--------|
| Eroded, arroyo bottom (sprayed reach) | 680 000 m³ |
| Eroded, arroyo walls (sprayed reach) | 56 000 m³ |
| **Total eroded, sprayed + study area** | **≈ 1 150 000 m³** |
| Aggraded, floodplains of sprayed reach | 500 000 m³ |
| **Total aggradation, both reaches** | **≈ 1 500 000 m³** |

**What the case must reproduce**, in order of stringency:

1. **Qualitative** — deposition depth correlates with vegetation density, and its
   *variability* correlates with vegetation type (P13 Figure 7). This is the result
   the vegetation operator exists to capture.
2. **Profile** — median deposit thickness decaying exponentially downstream from the
   sediment source, per P13's 1-D model:
   ```
   [P-1]   C(x) = (C₀ − S_P/(p v_s)) · exp(−p v_s x / q) + S_P/(p v_s)     [P13 6]
   ```
   with `S_P` a local source term. Deviations from P-1 are what P13 attribute to
   arroyo morphology and vegetation — so a 2-D model should reproduce the deviations,
   not just the profile.
3. **Budget** — total aggradation within the study area, against the table above.

> **Caution on the budget.** P13 are explicit that the sprayed-reach volumes of
> Vincent et al. (2009) are "not as well constrained" as their own lidar
> differencing. Treat the 1.5×10⁶ m³ total as an order-of-magnitude check; the lidar
> differencing over the study area is the defensible target.

> **Third grain size.** P13 use 0.045 mm, P14 uses 0.1 mm, `aS16` uses 0.065 mm — all
> for the same river, all citing the same body of work. See divergence D8 in the physics specification; none is
> wrong, they are choices for different purposes, but the spec must not silently
> inherit one.
