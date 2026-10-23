# Sediment transport physics specification

This is the specification the sediment implementation was written from. Every
equation is transcribed from a published, citable source, so each line of the
implementation can point at an equation label here and at the paper behind it.
The bracketed labels used throughout the ANUGA documentation -- {speclit}`E-1`,
{speclit}`T-7`, {speclit}`G-4` and the rest -- are defined in this document.

```{note}
This is the reference half of a longer internal document. Two sections are
omitted: an audit of a third-party implementation against this specification,
and the working list of open items. Neither defines any term the
documentation refers to. What remains is the notation, the governing
equations, every closure, the numerical scheme and the parameter defaults.
```

## Sources

| Tag | Source |
|-----|--------|
| **P14** | Perignon, M.C. (2014), *A Rolling Stone Gathers No Moss*, PhD thesis, University of Colorado Boulder, §3.3.2. Equations 3.1–3.11. |
| **FG21** | Fassett, C.I. & Goudge, T.A. (2021), *Modeling the Hydrodynamics, Sediment Transport, and Valley Incision of Outlet-Forming Floods From Martian Crater Lakes*, JGR Planets 126, e2021JE006979. Equations 1–9 + Supporting Information S1–S5. |
| **RDy26** | Feng, D., Tan, Z., Xu, D., Johnson, J. & Bisht, G. (2026), *RDycore-sediment v1.0*, EGUsphere preprint 2026-4859 (CC BY 4.0). Equations 1–13 + Appendix A1–A7. |
| **aSM16** | Perignon, M.C. (2016), *Using the Sediment Transport and Vegetation Operators in ANUGA*, [anugaSed/docs/anugaSed_manual.pdf](https://github.com/mperignon/anugaSed/blob/master/docs/anugaSed_manual.pdf), 4 May 2016. Equations 1–14. **This is the authoritative specification for the shipped code** and supersedes P14 where they differ. |
| **DL09** | Davy, P. & Lague, D. (2009), *Fluvial erosion/transport equation of landscape evolution models revisited*, JGR Earth Surface 114, F03007, doi:10.1029/2008JF001146. Equations 2–8, 19. Source of the E–D framework and of {speclit}`S-4`. |
| **P13** | Perignon, M.C., Tucker, G.E., Griffin, E.R. & Friedman, J.M. (2013), JGR Earth Surface 118(3), 1193–1209, doi:10.1002/jgrf.20073. Rio Puerco lidar differencing — the field data for a planned validation case (§10.1). |
| **W04** | Wilson, L., Ghatan, G.J., Head, J.W. & Mitchell, K.L. (2004), JGR Planets 109, E09003, doi:10.1029/2004JE002281. Eqs 4, 13–17: Darcy–Weisbach `f_c` by bed type. |
| **LL16** | Larsen, I.J. & Lamb, M.P. (2016), Nature 538, 229–232, doi:10.1038/nature19817. Methods: Manning–Strickler roughness closure. **Uses ANUGA.** |
| **EH67** | Engelund, F. & Hansen, E. (1967), *A Monograph on Sediment Transport in Alluvial Streams*, Teknisk Forlag. Eqs 3.1.3, 4.3.5. |
| **aS16** | [`anugaSed`](https://github.com/mperignon/anugaSed) source (Perignon 2016, MIT). |

## Clean-room statement

No code from `ANUGA-erosion` was read in preparing this document. FG21 states that
their implementation of the bedload equations is "nearly literal" from the paper,
so the paper is a sufficient specification. Commit messages implementing sections of
this document should cite the equation label and source tag.

---

## 1. Notation

| Symbol | Meaning | Units |
|--------|---------|-------|
| `h` | water depth | m |
| `u, v` | depth-averaged velocity components | m s⁻¹ |
| `\|v\|` | velocity magnitude √(u²+v²) | m s⁻¹ |
| `z` | bed elevation | m |
| `g` | gravitational acceleration (**domain parameter, not a constant**) | m s⁻² |
| `ρ` | water density (1000) | kg m⁻³ |
| `ρ_s` | sediment particle density | kg m⁻³ |
| `R` | submerged specific gravity, (ρ_s − ρ)/ρ | – |
| `D`, `d_g` | grain diameter | m |
| `ν` | kinematic viscosity of water (1×10⁻⁶) | m² s⁻¹ |
| `κ` | von Kármán constant (0.41) | – |
| `λ`, `φ` | bed porosity | – |
| `n` | Manning coefficient | s m⁻¹ᐟ³ |
| `f` | Darcy–Weisbach friction factor | – |
| `f_c`, `C_D` | friction coefficient, ≡ f/8 | – |
| `τ_b`, `τ` | bed shear stress | Pa |
| `τ*` | Shields (dimensionless) shear stress | – |
| `τ_c*` | critical Shields stress for motion | – |
| `u*` | shear velocity | m s⁻¹ |
| `c_s` | depth-averaged volumetric concentration, class s | – |
| `m_s` | conserved sediment variable, ≡ h·c_s | m |
| `c_b` | near-bed concentration | – |
| `w_s`, `v_s` | settling velocity | m s⁻¹ |
| `Z` | Rouse number, w_s/(κu*) | – |
| `d*`, `D*` | near-bed concentration profile factor, c_b/c | – |
| `E_s`, `Ė` | entrainment (erosion) flux | m s⁻¹ |
| `D_s`, `Ḋ` | deposition flux | m s⁻¹ |
| `q_b` | volumetric bedload flux per unit width | m² s⁻¹ |
| `q_b*` | dimensionless bedload flux (Einstein number) | – |
| `M_ℓ,s` | bed mass per unit area, layer ℓ class s | kg m⁻² |
| `β` | ratio of sediment transport speed to water speed | – |
| `N_s` | number of sediment classes | – |

**Sign convention (fixed for the whole document):** `E` is positive from bed into
the water column; `D` is positive from water column onto the bed. Bed elevation
rises under net deposition.

---

## 2. Governing equations

### 2.1 Hydrodynamics — unchanged

ANUGA already solves the conservative shallow water equations; this is stated only
to fix notation. `[P14 3.1–3.3]`, `[RDy26 3–5]`:

```{math}

\frac{\partial \mathbf{U}}{\partial t} + \frac{\partial \mathbf{E}}{\partial x} + \frac{\partial \mathbf{G}}{\partial y} = \mathbf{S}, \qquad \mathbf{U} = \begin{bmatrix} h & hu & hv \end{bmatrix}^{T}
```

with bed friction slope and bed slope in `S`. Sediment does **not** modify these
equations in any phase of the planned work (see §2.4).

### 2.2 Suspended sediment transport

The prognostic sediment variable is **mass per unit area**, not concentration
`[RDy26 1]`:

```{index} single: physics label; [G-1]
```
(spec-g-1)=
```{math}
:nowrap:

\begin{align}
m_s = h\,c_s, \qquad s = 1 \dots N_s \tag{G-1}
\end{align}
```

State vector `[RDy26 2]`:

```{index} single: physics label; [G-2]
```
(spec-g-2)=
```{math}
:nowrap:

\begin{align}
\mathbf{U} = \begin{bmatrix} h & hu & hv & m_1 & m_2 & \dots & m_{N_s} \end{bmatrix}^{T} \tag{G-2}
\end{align}
```

Conservative transport, one equation per class `[RDy26 6]`, which is `[DL09 4]`
written per class — DL09 give it in Lagrangian form,
`D(c_S h)/Dt = ∂(c_S h)/∂t + div(q_S) = ė − ḋ`:

```{index} single: physics label; [G-3]
```
(spec-g-3)=
```{math}
:nowrap:

\begin{align}
\frac{\partial m_s}{\partial t} + \frac{\partial (u\,m_s)}{\partial x} + \frac{\partial (v\,m_s)}{\partial y} = E_s - D_s + S_{m_s} \tag{G-3}
\end{align}
```

`S_ms` is an optional external source (hillslope yield, tributary load, rainfall
washoff).

> **Reconciliation note.** P14 Eq 3.4 writes the same balance as
> `∂(Ch)/∂t = Ė − Ḋ − (∂q_sx/∂x + ∂q_sy/∂y)` with `q_s(x,y) = β C q(x,y)`, where
> `q` is specific discharge and **β is the ratio of sediment transport speed to
> water speed**. {speclit}`G-3` is the β = 1 case. β is retained as an exposed parameter —
> see the note on `beta` in the ANUGA sediment documentation.

### 2.3 Bed evolution (Exner)

`[P14 3.5]`, `[FG21 6]`:

```{index} single: physics label; [G-4], single: physics label; [G-5]
```
(spec-g-4)=
(spec-g-5)=
```{math}
:nowrap:

\begin{align}
\frac{\partial z}{\partial t} &= \frac{D - E}{1 - \lambda} && \text{suspended contribution} \tag{G-4} \\
\frac{\partial z}{\partial t} &= -\frac{1}{1 - \lambda}\,\nabla\cdot\mathbf{q}_b && \text{bedload contribution} \tag{G-5}
\end{align}
```

Both act on the same `z`; when both operators are active the contributions sum.

A morphological acceleration factor `M` (`morphological_factor`, default 1) multiplies
both contributions per step, with the erodible-base cap [L-5] divided by `M`; the water
column is not scaled. Bed and water-column sediment budgets then differ by exactly `M`.

### 2.4 Coupling stages

| Stage | Bed elevation | Bed → flow | Sediment → momentum | Reference model |
|-------|---------------|-----------|---------------------|-----------------|
| Phase 3 | fixed | no | no | RDy26 v1.0 |
| Phase 4 | evolving via {speclit}`G-4`,{speclit}`G-5` | yes | no | FG21 |
| not scoped | evolving | yes | yes | — |

Neither reference model implements the momentum feedback. FG21 flags this as the
assumption most likely to fail for very energetic flows that "could easily have
debris flow-like behaviors". Fluid density is held at ρ regardless of `c_s` in all
scoped phases.

---

## 3. Bed shear stress — the reconciled closure

**This is the most important section in the document.** Every erosion and deposition
rate depends on `τ_b`. The shipped `anugaSed` uses a *different but documented*
closure — see §3.4 — which this specification recommends superseding.

### 3.1 The quadratic drag law

FG21, RDy26 and P14 all specify the same quadratic drag law, in different notation:

```{index} single: physics label; [T-1]
```
(spec-t-1)=
```{math}
:nowrap:

\begin{align}
\tau_b = \rho\, f_c\, |\mathbf{v}|^2, \qquad f_c \equiv f/8 \tag{T-1}
\end{align}
```

| Source | As written | Equivalent to |
|--------|-----------|---------------|
| FG21 Eq 1 | `τ = ρ f_c v²`, "f_c is the Darcy–Weisbach friction factor divided by 8" | {speclit}`T-1` |
| P14 Eq 3.7 | `τ* = u*²/(R g d_g)` with `u* = \|v\|√(f/8)` | {speclit}`T-1` via {speclit}`T-3` |
| RDy26 A5–A7 | `C_D = g n² h⁻¹ᐟ³` ; `τ_b = ρ_w C_D (u² + v²)` | {speclit}`T-1` with `f_c = C_D` |
| aSM16 Eq 6 | `τ_b = ρ_w u*²` | {speclit}`T-1` via {speclit}`T-2`, but `u*` from {speclit}`T-7` not {speclit}`T-2` |

Shear velocity:

```{index} single: physics label; [T-2]
```
(spec-t-2)=
```{math}
:nowrap:

\begin{align}
u_* = \sqrt{\tau_b/\rho} = |\mathbf{v}|\sqrt{f_c} = |\mathbf{v}|\sqrt{f/8} \tag{T-2}
\end{align}
```

Shields stress `[FG21 2]`, `[P14 3.7]`:

```{index} single: physics label; [T-3]
```
(spec-t-3)=
```{math}
:nowrap:

\begin{align}
\tau^{*} = \frac{\tau_b}{(\rho_s - \rho)\, g\, D} = \frac{u_*^2}{R\, g\, D} \tag{T-3}
\end{align}
```

Excess Shields stress `[FG21 3]`:

```{index} single: physics label; [T-4]
```
(spec-t-4)=
```{math}
:nowrap:

\begin{align}
\tau_x = \tau^{*} - \tau_c^{*} \tag{T-4}
\end{align}
```

> **Typographic note on P14 Eq 3.7.** The thesis prints `u* = √(f/8)`, omitting the
> velocity. `u*` must be a velocity and `√(f/8)` is dimensionless, so the intended
> form is {speclit}`T-2`. FG21 Eq 1 and RDy26 A7 both confirm the intent.

### 3.2 Velocity regularisation

Velocity must be recovered from momentum with the standard ANUGA depth-limiting
form to avoid division by near-zero depth. RDy26 A22–A23 explicitly adopts
"the ANUGA-style form (Mungkasi & Roberts, 2012)":

```{index} single: physics label; [T-5]
```
(spec-t-5)=
```{math}
:nowrap:

\begin{align}
u = \frac{(uh)\,h}{h^2 + h_\epsilon^2}, \qquad v = \frac{(vh)\,h}{h^2 + h_\epsilon^2} \tag{T-5}
\end{align}
```

with `h_ε` a regularisation depth (`aS16` uses 1×10⁻⁶). Set `u = v = c_s = 0`
below the wet/dry threshold.

### 3.3 Friction closures

`f_c` varies **per cell, per timestep** as a function of the evolving flow depth.
It must be recomputed inside the sediment kernel, not read once from Manning's `n`.
The Manning ↔ Darcy–Weisbach bridge `[FG21 9]` is:

```{index} single: physics label; [T-6]
```
(spec-t-6)=
```{math}
:nowrap:

\begin{align}
f = \frac{8\, g\, n^2}{h^{1/3}} \quad\Longleftrightarrow\quad f_c = C_D = g\, n^2\, h^{-1/3} \tag{T-6}
\end{align}
```

which is exactly RDy26 A5. Three selectable closures `[FG21 §2.2.3]`:

| Mode | Definition | `n` varies? |
|------|-----------|-------------|
| `constant` | user supplies `n` | no |
| `larsen_lamb` | constant `n` from bed roughness parameters `k_s`, `r_d`, `r_br`, `σ_br` — Larsen & Lamb (2016) Eqs 3–4 | no |
| `wilson` | `f` as an empirical function of grain size and depth — Wilson et al. (2004) Eqs 13–15; `n` then follows from {speclit}`T-6` | yes, spatially and temporally |

**In all three modes `f_c` still varies per cell per timestep**, because {speclit}`T-6`
depends on `h`. This is the coupling most easily missed.

### 3.3.1 ⚠ Friction-factor conventions — a factor-of-2 and factor-of-8 trap

Three conventions appear across these sources, all describing the same physics.
**Getting this wrong silently scales every transport rate.**

| Convention | Definition | Relation |
|------------|-----------|----------|
| Darcy–Weisbach `f` | `τ₀ = ρ f V²/8` | — |
| **This spec, `f_c`** (= FG21, RDy26 `C_D`) | `τ₀ = ρ f_c V²` | `f_c = f/8` |
| Engelund–Hansen `f_EH` `[EH67 3.1.3]` | `I = f_EH V²/(2gD)` | `f_EH = f/4 = 2 f_c` |

W04 explicitly warns of this: "some theoretical treatments adopt factors other than 8
in the above expression and thus quote systematically different, but still
dimensionless, values for the friction factor." Convert before substituting.

### 3.3.2 `wilson` mode — `f_c` by bed type  `[W04 13–17]`

W04 give `f_c` for five bed types. `R` is hydraulic radius (`R ≈ h` for wide
channels); `D50`/`D84`/`D90` are bed clast percentiles.

```{index} single: physics label; [T-8], single: physics label; [T-9], single: physics label; [T-10], single: physics label; [T-11], single: physics label; [T-12]
```
(spec-t-8)=
(spec-t-9)=
(spec-t-10)=
(spec-t-11)=
(spec-t-12)=
```{math}
:nowrap:

\begin{align}
\text{Sand bed} & \quad (8/f_c)^{1/2} = 8.46\,(R/D_{50})^{0.1005} \tag{T-8} \\
\text{Gravel bed} & \quad (8/f_c)^{1/2} = 5.75\,\log_{10}(R/D_{84}) + 3.514 \tag{T-9} \\
\text{Boulder bed} & \quad (8/f_c)^{1/2} = 5.62\,\log_{10}(R/D_{84}) + 4.0 \tag{T-10} \\
\text{Steep pool-fall} & \quad (8/f_c)^{1/2} = 4.60\,\log_{10}(d_s/D_{90}) + 4.203 \tag{T-11} \\
\text{Fixed roughness} & \quad (8/f_c)^{1/2} = 5.657\,\log_{10}(R/r) + 6.6303 \tag{T-12}
\end{align}
```

FG21 cite Eqs 13–15, i.e. {speclit}`T-8`–{speclit}`T-10`. {speclit}`T-11` uses `d_s` = depth of water
*plus* sediment; {speclit}`T-12` uses roughness element size `r`. Velocity follows from
`[W04 4]`: `U = ((8gRS)/f_c)^½`, which is {speclit}`T-1` rearranged.

> ### ⚠⚠ `f_c` ABOVE IS **NOT** THIS SPEC'S `f_c`.
>
> {speclit}`T-8`–{speclit}`T-12` reproduce W04's symbol verbatim, and **W04's `f_c` is the
> Darcy–Weisbach `f`, not `f/8`.** Their Eq 4, quoted immediately above, is
> `U = ((8gRS)/f_c)^½` — the standard Darcy–Weisbach velocity, which carries `f`.
> This spec's `f_c` is defined in §3.3.1 as `f/8`. The two collide on one symbol,
> which is precisely the trap §3.3.1 exists to warn about — and §3.3.2 walked
> into it by adopting W04's notation without converting.
>
> **The conversion, and a convenient cancellation.** Writing `X` for the
> right-hand side of {speclit}`T-8`–{speclit}`T-12`:
>
> ```
>     X = (8/f_W04)^½   with  f_W04 = f      ⟹   f = 8/X²
>     f_c(this spec) = f/8                   ⟹   f_c = 1/X²
> ```
>
> The eights cancel: **`f_c = 1/X²`**. Taking W04's `f_c` as ours would make
> `f_c` — and therefore `τ_b`, and therefore every transport rate — **8× too
> large**.
>
> **Confirmed against W04's own headline result.** They derive
> `n = 0.0545 s m⁻¹ᐟ³` for Martian channels. Under {speclit}`T-6`,
> `n = √(f_c h^⅓ / g)` with *this spec's* `f_c`. Taking a sand bed at
> `R/D50 = 1000` and Martian `g = 3.71 m s⁻²`, `X = 8.46 × 1000^0.1005 = 16.94`:
>
> | reading | `f_c` | depth giving `n = 0.0545` |
> |---------|------:|--------------------------:|
> | **correct**, `f_c = 1/X²` | 0.00349 | **31.6 m** — an outflow channel |
> | literal, `f_c = 8/X²` | 0.02788 | 0.062 m — not an outflow channel |
>
> Martian outflow channels are tens of metres deep, so the correct reading is
> the one that reproduces W04's own number.

These are simplifications of W04's fuller forms, obtained by fixing weakly-varying
terms at Martian averages (slope `S`, clast spread `σ_g`, and the `a` term in the
gravel formula). W04 quantify the resulting error as ~2%, <5% and ~13% respectively.
**On Earth those fixed values may not apply** — for terrestrial work prefer W04's
unsimplified Eqs 7a/7b, 8–9, 10, 11 (also in §14's source).

### 3.3.3 `larsen_lamb` mode — Manning–Strickler  `[LL16 Methods]`

```{index} single: physics label; [T-13], single: physics label; [T-14], single: physics label; [T-15]
```
(spec-t-13)=
(spec-t-14)=
(spec-t-15)=
```{math}
:nowrap:

\begin{align}
\frac{u}{u_*} &= 8.1\,\left(\frac{h}{k_s}\right)^{1/6} && \text{Manning–Strickler} \tag{T-13} \\
n &= \frac{k_s^{1/6}}{8.1\,\sqrt{g}} && \text{rearranged with the Manning equation} \tag{T-14} \\
k_s &= r_d\, r_{br}\, \sigma_{br} && \text{bedrock roughness length scale} \tag{T-15}
\end{align}
```

`r_d` and `r_br` are hydraulic roughness scaling parameters and `σ_br` is one
standard deviation of bedrock bed elevation. LL16 set `r_d = 2`, `r_br = 2`, and
measured `σ_br ≈ 5 m` across five reaches, giving `k_s = 20 m` and:

```{math}

n = \frac{20^{1/6}}{8.1\,\sqrt{9.81}} = 0.0649 \approx 0.065 \qquad \text{matches LL16's stated value}
```

`n` is then **spatially and temporally uniform**; only `f_c` varies, through {speclit}`T-6`.

> **LL16 used ANUGA.** Their Methods state "ANUGA implements bed friction with
> Manning's roughness coefficient (n)" and that ANUGA generated their triangular
> mesh from a USGS 10 m DEM. This is a *Nature* paper built on this codebase — worth
> knowing both as provenance for {speclit}`T-13`–{speclit}`T-15` and as a citable precedent for
> ANUGA in megaflood work. Note `σ_br` is a **site-measured** quantity: 5 m is Moses
> Coulee, not a universal default.

> **Both non-constant closures are megaflood parameterisations, not general-purpose
> terrestrial flood closures.** Wilson et al. (2004) is a *Mars
> outflow channel* study — its `n = 0.0545 s m⁻¹ᐟ³` is derived for Martian channels
> under Martian gravity, and is FG21's constant-`n` default for exactly that reason.
> Larsen & Lamb (2016) is the Channeled Scablands, a terrestrial *megaflood*. Neither
> is calibrated for ordinary river or urban flood modelling, which is ANUGA's main
> use. **For standard flood work, `constant` with an `n` from conventional tables
> remains the right default**; treat `wilson` and `larsen_lamb` as specialist modes
> for planetary and megaflood applications and label them as such in the API.

> **Open item F1.** Wilson et al. (2004) Eqs 13–15 and Larsen & Lamb (2016) Eqs 3–4
> are not reproduced in FG21's main text (they are in Supporting Information S3).
> Transcribe from the primary papers — both are now fully identified in §14 — before
> implementing those modes. `constant` is sufficient for Phases 1–3.

### 3.4 The depth–slope closure used by anugaSed

`aSM16` Eq 7 specifies a different shear velocity — the **depth–slope product**:

```{index} single: physics label; [T-7], single: physics label; [T-7e]
```
(spec-t-7)=
(spec-t-7e)=
```{math}
:nowrap:

\begin{align}
\tau_b &= \rho\, g\, h\, S \quad\Longrightarrow\quad u_* = \sqrt{g\, S\, h} && \text{[aSM16 6–7]} \tag{T-7} \\
\tau_b &= \rho\, g\, h\, S_f && S_f \text{ the free-surface slope} \tag{T-7e}
\end{align}
```

{speclit}`T-7e` is ANUGA's variant, not in the sources: the same depth–slope
product with the energy (free-surface) slope in place of the bed slope, which
removes the normal-flow assumption from the slope while keeping the rest of
the closure. It is selected with `set_shear_closure('energy_slope')`.

This is a documented closure, not undocumented drift. It is the
**steady uniform (normal) flow** approximation:
it assumes the energy slope equals the bed slope and that the flow is locally in
equilibrium.

**It should nonetheless be superseded by {speclit}`T-1`, for three reasons:**

1. **Validity.** Normal-flow equilibrium is exactly what does not hold in a
   dam-breach, outburst or flash flood — the cases this add-on exists to model.
   {speclit}`T-1` makes no equilibrium assumption.
2. **`S` is taken from the bed, not the energy grade line.** `aSM16` says "`S` is the
   local slope of the bed". The energy slope `S_f = f_c|v|²/(gh)` is the correct
   quantity, and substituting it into {speclit}`T-7` recovers {speclit}`T-1` identically. The
   commented-out energy-slope block in `aS16` shows this was recognised.
3. **The clamp is not in the manual.** `aS16` applies `S ← min(S, mean(S)/2)`, a
   domain-global rescaling with no counterpart in `aSM16` and no physical basis.
   Under {speclit}`T-1` it is unnecessary.

Retain {speclit}`T-7` behind a `legacy_depth_slope` flag for reproducing published
anugaSed results. Default to {speclit}`T-1`.

---

## 4. Suspended sediment

### 4.1 Entrainment — two published formulations

The sources do **not** agree here, and the disagreement is physical rather than
typographic. Both formulations should be implemented behind a selector.

**(a) Saturating reference-concentration form** — Smith & McLean (1977), Parker
(1998), as used by `[FG21 7]`:

```{index} single: physics label; [E-1], single: physics label; [E-2]
```
(spec-e-1)=
(spec-e-2)=
```{math}
:nowrap:

\begin{align}
E^{*} &= \frac{0.65\,\gamma_0\, S}{1 + \gamma_0\, S}, \qquad S = \frac{\tau^{*}}{\tau_c^{*}} - 1 \tag{E-1} \\
E &= v_s\, E^{*} \tag{E-2}
\end{align}
```

with `γ₀ = 0.0024` (empirical) and `τ_c* = 0.04` (FG21's choice for suspension).
`E*` is a dimensionless near-bed reference concentration; 0.65 is the maximum
packing fraction, so `E*` saturates rather than growing without bound. `E` is
recovered as a flux by multiplying by the settling velocity {speclit}`E-2`.

**(b) Cohesive / fine-grained form** — Hanson (1990), Hanson & Simon (2001), as
specified by `[aSM16 3–5]` and implemented in `aS16`:

```{index} single: physics label; [E-3], single: physics label; [E-5]
```
(spec-e-3)=
(spec-e-5)=
```{math}
:nowrap:

\begin{align}
\dot{E} &= K_e\,(\tau_b - \tau_c) && \text{dimensional stress, Pa} \tag{E-3} \\
K_e &= \frac{0.2\times 10^{-6}}{\tau_c^{0.5}} && [\mathrm{m^3\,N^{-1}\,s^{-1}}], \nonumber \\
\tau_c &= \tau_c^{*}\,(\rho_s - \rho)\, g\, D_{50} && \tag{E-5}
\end{align}
```

{speclit}`E-5` is the jet-test erodibility relation of Hanson & Simon (2001) for cohesive
streambeds, `k_d = 0.2 τ_c^(−0.5)` in cm³ N⁻¹ s⁻¹; the `10⁻⁶` converts cm³ to m³.

**(c) Partheniades form on dimensional stress** — as used by `[RDy26 A8]`, per bed
layer ℓ:

```{index} single: physics label; [E-4]
```
(spec-e-4)=
```{math}
:nowrap:

\begin{align}
E_\ell^{\mathrm{pot}} = \begin{cases} K_{e,\ell}\, \dfrac{\tau_b - \tau_{c,\ell}}{\tau_{c,\ell}} & \tau_b > \tau_{c,\ell} \\[4pt] 0 & \text{otherwise} \end{cases} \tag{E-4}
\end{align}
```

Note {speclit}`E-4` and {speclit}`E-3` both use *dimensional* stress; {speclit}`E-1` uses Shields stress.

**de Leeuw et al. (2020)** `[dL20 26a]` -- a non-cohesive alternative to {speclit}`E-1`,
driven by the skin-friction shear velocity and the Froude number, with no
Shields threshold:

```{index} single: physics label; [E-6]
```
(spec-e-6)=
```{math}
:nowrap:

\begin{align}
E^{*} = \frac{A\,X^{\beta}}{1 + 3A\,X^{\beta}}, \qquad
X = \left(\frac{u_{*,\mathrm{sk}}}{v_s}\right)^{\alpha}\mathrm{Fr} - 0.015 \tag{E-6}
\end{align}
```

`E = v_s E*` as in {speclit}`E-2`; `E*` is the near-bed concentration at `0.1 h`, capped at
1/3. `u*_sk` from Manning-Strickler on a skin roughness `k_s` (`[dL20 7]`:
`U/u*_sk = 8.1 (H_sk/k_s)^{1/6}`, `u*_sk = sqrt(g H_sk S)`, `k_s = 3 D84`),
with `S = tau_b/(rho g h)` from the active shear closure, `H_sk <= h`, and
`k_s` floored at the viscous sublayer `nu/(8 u*)`. Constants: `A = 4.74e-4,
alpha = 1.5, beta = 1.18` (`[dL20 26a]`, sand and gravel), or `A = 7.04e-4,
alpha = 0.945, beta = 1.81` as the Delta-X Wax Lake model uses them for mud as
flocculated bed material load and for sand (`[Wan23]`, after Nghiem et al.
2022).

> **On the dimensions of `K_e`.**
> Reading P14 alone suggests a dimensional inconsistency in `K_e`. It is an
> artefact of that reading. P14 Eq 3.6 applies
> `K_e` to *Shields* stress (dimensionless), which cannot balance; `aSM16` Eq 3
> applies it to *dimensional* stress, which balances exactly:
>
> ```
>   [m³ N⁻¹ s⁻¹] × [N m⁻²]  =  m s⁻¹   ✓
> ```
>
> The shipped code implements {speclit}`E-3`/{speclit}`E-5` correctly and with proper provenance.
> **No action needed.** P14 Eq 3.8 (`K_e = K_e* d_g √(Rgd_g)`, `K_e* = 12` after
> Wilson 1966) is a different, non-cohesive parameterisation that the manual
> abandoned; do not implement it.

### 4.1.1 Which erosion law? — a regime choice, not a preference

{speclit}`E-1` and {speclit}`E-3` are **not competing formulations of the same physics**. They
describe different sediment:

| | {speclit}`E-3` cohesive | {speclit}`E-1`/{speclit}`K-1` non-cohesive |
|---|---|---|
| Sediment | silt, clay, cohesive bank material | sand, gravel, boulders |
| Threshold physics | inter-particle cohesion, jet-test calibrated | grain weight, Shields entrainment |
| Calibration source | Hanson & Simon (2001) | Smith & McLean (1977) / Wong & Parker (2006) |
| Planned field validation (§10.1) | Rio Puerco | Mars crater breach |
| Typical `D` | < 0.1 mm | mm to cm |

Selecting between them is a **statement about the bed material**, and the API should
present it that way (`bed_material='cohesive'|'noncohesive'`) rather than as an
opaque numerical switch. Getting this wrong is a physics error, not a tuning error.

### 4.2 Settling velocity

Two published options; they differ in provenance, not in intent.

**Ferguson & Church (2004)** — used by `[P14 3.10]` and implemented in `aS16`:

```{index} single: physics label; [S-1]
```
(spec-s-1)=
```{math}
:nowrap:

\begin{align}
v_s = \frac{R\, g\, d_g^2}{C_1\,\nu + \sqrt{0.75\, C_2\, R\, g\, d_g^3}} \tag{S-1}
\end{align}
```

with `C₁ = 18`, `C₂ = 0.4` for smooth spheres (1.0 and 1.1 respectively for natural
grains — expose both). This form is smooth across the Stokes → turbulent transition
and is cheap, branch-free and GPU-friendly.

**Dietrich (1982) Eqs 5–9** — used by `[FG21 §2.2.2]` via Supporting Information S1.
A polynomial fit in the dimensionless particle size; more accurate for natural
irregular grains, but requires transcription from Dietrich or FG21 S1.

**Recommendation:** implement {speclit}`S-1` as the default. It is fully specified here,
dimensionally clean, and already the basis of the existing code. Add Dietrich as an
option only if a crater-breach validation case needs it to reproduce FG21. Note DL09 themselves
take settling velocities from Dietrich (1982).

> **{speclit}`S-1` numerically verified.** P13 §[54] reports `v_s = 0.00175 m s⁻¹` for
> 0.045 mm quartz using Ferguson & Church (2004). Evaluating {speclit}`S-1` at
> `D = 4.5×10⁻⁵ m`, `R = 1.65`, `g = 9.81`, `ν = 10⁻⁶`:
>
> | `C₁` | `C₂` | `v_s` | error |
> |------|------|-------|-------|
> | **18** | **0.4** | **0.00176** | **+0.3%** |
> | 18 | 1.0 | 0.00172 | −1.7% |
> | 20 | 1.1 | 0.00155 | −11.2% |
>
> This confirms both the transcription of {speclit}`S-1` and that the smooth-sphere
> constants `C₁ = 18`, `C₂ = 0.4` are the ones in use. Use this as a unit test.

### 4.3 Rouse number and the near-bed concentration profile

```{index} single: physics label; [S-2], single: physics label; [S-3]
```
(spec-s-2)=
(spec-s-3)=
```{math}
:nowrap:

\begin{align}
Z &= \frac{v_s}{\kappa\, u_*}, \qquad \kappa = 0.41 \tag{S-2} \\
c_b &= c\; d^{*}(Z) \tag{S-3}
\end{align}
```

`d*` is the ratio of near-bed to depth-averaged concentration. Physically: at low
`Z` sediment is well mixed through the column (`d* → 1`, washload); at high `Z` it
concentrates near the bed (`d* ≫ 1`).

The sources diverge:

| Source | `d*` treatment |
|--------|---------------|
| P14 Eq 3.9 | `D* = 1` — uniform suspension assumed, citing Nordin (1963), Griffin et al. (2014) |
| **aSM16 Eq 10** | **the defining Rouse–Vanoni integral, from Davy & Lague (2009)** |
| FG21 Fig 2 / S2 | `d*` a function of `Z`, plotted; numeric form only in Supporting Information S2 |
| aS16 | 8-term polynomial in `Z`, plus a linear branch for `Z > 4` — a *fit* to {speclit}`S-4` |

**The defining expression** `[aSM16 10]`, after Davy & Lague (2009), obtained by
requiring that sediment discharge be the depth integral of concentration × velocity
over a Rouse–Vanoni concentration profile and a logarithmic velocity profile:

```{index} single: physics label; [S-4]
```
(spec-s-4)=
```{math}
:nowrap:

\begin{align}
d^{*} = \frac{\displaystyle\int_a^h \ln(z/z_0)\,dz}{\displaystyle\int_a^h \left(\frac{h - z}{h - a}\cdot\frac{a}{z}\right)^{Z} \ln(z/z_0)\,dz} \tag{S-4}
\end{align}
```

with `Z` the Rouse number {speclit}`S-2`, `z₀` the roughness length and `a` the reference
height near the bed.

> **The profile factor is misprinted in DL09 as published.**
> `[DL09 19]` prints it as `((z−a)/(h−a) · a/z)^Z`. That is wrong, and
> DL09's own paper contains the proof.
>
> *Internal contradiction in DL09.* Immediately above their Eq 19 they write the
> sediment flux as
>
> ```
>     q_S = ∫ₐʰ c_s(z) u(z) dz = c_S(a) ∫ₐʰ [ … ]^Z u(z) dz
> ```
>
> Factoring `c_S(a)` out requires the bracket to equal **1** at `z = a`. The
> printed factor `(z−a)/(h−a) · a/z` is **0** there, which would make
> `c_s(a) = 0 ≠ c_S(a)`. The reference height is by definition where the near-bed
> concentration is evaluated, so the profile cannot vanish there.
>
> *The intended expression.* The Rouse–Vanoni profile
> `c(z)/c_a = [(h−z)/z · a/(h−a)]^Z` rearranges **exactly** to
> `((h−z)/(h−a) · a/z)^Z`, which is 1 at `z = a` and 0 at `z = h`. So `(z−a)` is a
> slip for `(h−z)`. The model was never in doubt — only one glyph.
>
> *It is not a small numerical difference.* At `z₀/h = 10⁻⁴`:
>
> | `Z` | as printed `(z−a)`, `a/h`=0.005 | corrected `(h−z)`, `a/h`=0.005 |
> |-----|-------------------------------:|-------------------------------:|
> | 0.10 | 1.28 | 1.29 |
> | 0.50 | 14.3 | 11.1 |
> | 2.00 | **41 227** | **356** |
>
> Both forms satisfy DL09's Figure 4 caption for `Z < 0.1` (they agree closely
> there), so the caption alone does not discriminate. The `Z = 2` column does:
> `d* ≈ 4×10⁴` is not a value that appears on their plotted figure, whereas the
> corrected form stays in a plottable range. **DL09's Figure 4 was evidently
> computed with the correct profile; their Eq 19 as typeset is a typo.**
>
> Everything downstream in this document, and the implementation, uses the
> corrected form.

**{speclit}`S-4` is confirmed as `[DL09 19]`** — the manual's Eq 10 reproduces it exactly.
DL09 derive it by assuming a Rouse–Vanoni concentration profile and a logarithmic
velocity profile, and integrating `c_s(z)u(z)` from a reference height `a` to `h`.

**DL09's own assessment of the expression** (their §4.1, Figure 4):

> "The above expression depends mostly on `Z`, and only slightly on the ratios
> `a/h` and `z₀/h`."

and from the Figure 4 caption:

> "`d*` is always larger than 1 and smaller than 3 for Rouse numbers `Z` smaller
> than 0.1. For `Z` between 0.1 and 2, `d*` increases rapidly mainly as a function
> of `a/h`."

**Documented implementation choices in `aSM16`** — these are approximations, and
each is a candidate for improvement:

- `z₀ = D50/30`  (standard hydraulically-rough wall)
- **flow depth assumed to be 1 m when evaluating the integral**, irrespective of the
  actual local depth
- `d*` recomputed only **once every 10 timesteps**, for speed

Physical expectation stated in `aSM16`: `d*` lies between 1 and 3 for `Z < 0.1`, is
close to 1 for large rivers or fine particles, and becomes ≫ 1 for small rivers with
coarse particles where entrainment is largely bedload.

> **On the shipped polynomial.**
> It is a fit to {speclit}`S-4`. Its fixed 1 m depth is a limitation rather than a
> defect:
>
> - `d*` depends *mostly on `Z`*, and `Z` is computed per cell from the actual `u*`.
>   The fixed depth enters only through the weakly-dependent ratios `a/h` and `z₀/h`.
>   Fitting `d*` as a function of `Z` alone is therefore what DL09's own analysis
>   supports, and the shipped polynomial is a **defensible** engineering choice.
> - The caveat is the mid-range: for `0.1 < Z < 2`, DL09 state `d*` "increases
>   rapidly mainly as a function of `a/h`". A single-depth fit is weakest precisely
>   there. **Quantify `d*` sensitivity to `a/h` across the `Z` range of interest
>   before deciding whether a depth-dependent fit is warranted.**
> - `a`, the near-bed reference height, is a required input to {speclit}`S-4` that `aSM16`
>   does not state. Recover it from the code, document it, and expose it.
>
> Remaining engineering (tracked as **S1a**): {speclit}`S-4` is a 1-D quadrature per cell,
> far too expensive inside a GPU kernel every step, so a fitted form is the right
> answer. Regenerate the fit from {speclit}`S-4`, document its `Z` range, and make
> out-of-range behaviour explicit rather than extrapolating a polynomial. Keep
> `d* = 1` available as the P14/P13 limiting case and the 10-timestep refresh as a
> tunable rather than a hard-coded stride.

> **Implemented and tested**; see
> `test_rouse.py` (12/12) and `core_rouse_d_star()`.
>
> **S1b first, because it changes S1a.** The `a/h` sensitivity, from the
> corrected {speclit}`S-4` at `z₀/h = 10⁻⁴`:
>
> | `Z` | `a/h`=0.01 | 0.02 | 0.05 | 0.10 | spread |
> |-----|-----------:|-----:|-----:|-----:|-------:|
> | 0.05 | 1.268 | 1.226 | 1.174 | 1.139 | 1.11× |
> | 0.25 | 3.064 | 2.598 | 2.109 | 1.821 | 1.68× |
> | 1.00 | 34.237 | 19.654 | 9.967 | 6.270 | **5.46×** |
> | 4.00 | 507.489 | 225.822 | 80.775 | 38.604 | **13.15×** |
>
> `z₀/h` over two decades moves `d*` by 1.01×–1.14× — genuinely weak, exactly as
> DL09 say. **`a/h` is not.** So the Draft-3 conclusion above, that a `Z`-only fit
> is what DL09's analysis supports, holds only for **small `Z`**. Above `Z ≈ 0.25`
> it fails, and DL09's other remark — that for `0.1 < Z < 2`, `d*` "increases
> rapidly mainly as a function of `a/h`" — is the one that governs the range of
> interest. **A `Z`-only fit is not adequate; the fit must take `a/h` too.**
>
> **S1a.** Fitted in both variables. Using the structure of the integral rather
> than fitting blind is what makes it cheap enough for a kernel: near the bed the
> Rouse profile behaves like `z^(−Z)`, so `d*` diverges roughly as `(a/h)^(−Z)`.
> Factoring that out,
>
> ```
>       ln d*  =  −Z ln(a/h)  +  P(Z, ln(a/h))
> ```
>
> leaves a mild remainder that a low-order polynomial captures well:
>
> | form | terms | max err |
> |------|------:|--------:|
> | direct polynomial in `(ln Z, ln(a/h))`, `Z ≤ 4`, `a/h ≥ 0.002` | 36 | 10.0% |
> | direct polynomial, suspension range | 18 | 3.1% |
> | **structural form above** | **15** | **0.84%** |
>
> Fitted range `Z ∈ [0.01, 2.5]`, `a/h ∈ [0.01, 0.15]`. Beyond `Z ≈ 2.5` transport
> is essentially bedload and this ratio is the wrong model. Out-of-range inputs are
> **clamped, not extrapolated**. `d* = 1` remains available as the P14/P13 limiting
> case (`sediment_d_star_mode = 0`, still the default).
>
> **`a` is exposed and remains a judgement call.** `aSM16` require it and never
> state it, and it was not recovered from their code — the default used here
> (`a = 2 d_g`, with the standard van Rijn floor `a ≥ 0.01 h`) is **ours**. Given
> the 13× spread above this is a first-order modelling choice, not a detail, and
> recovering `aSM16`'s actual value remains worth doing before any validation run
> that claims to reproduce them.
>
> The fixed-1 m-depth and 10-timestep-refresh approximations of `aSM16` are both
> now moot: `d*` is evaluated per cell per step from the local `h` and `u*`, and
> costs nothing measurable (Ns=0 gate unchanged, and the evaluation is a Horner
> in `Z` over coefficients quadratic in `ln(a/h)`).

### 4.3.1 Transport length — a scale that constrains model setup

DL09's central quantity is the **disequilibrium (transport) length** `ξ`: the mean
distance a particle travels in the flow before being trapped on the bed.

```{index} single: physics label; [S-5]
```
(spec-s-5)=
```{math}
:nowrap:

\begin{align}
\xi = \frac{q}{d^{*}\, v_s} \qquad \text{[DL09 8]} \tag{S-5}
\end{align}
```

It is the length over which the sediment load relaxes toward capacity. `ξ → 0`
recovers transport-limited (local-capacity) behaviour; `ξ → ∞` recovers
detachment-limited. **The explicit `E − D` balance of {speclit}`G-3` exists precisely
because `ξ` is finite** — the flow is *not* instantaneously at capacity.

Evaluated for the Rio Puerco 2006 flood from P13's own parameters
(`q = 2.66 m² s⁻¹`, `d* = 1`, `v_s = 0.00175 m s⁻¹`):

```{math}

\xi \approx 1520\ \mathrm{m}, \qquad \text{about } 6.5\times \text{ the 235 m mean arroyo width}
```

DL09 note the same generally: "except for coarse sand, the disequilibrium length is
about larger than the river width".

**Two practical consequences:**

1. **Mesh and domain sizing.** A domain shorter than a few `ξ` cannot develop the
   deposition pattern; a cell much larger than `ξ` makes the local-equilibrium error
   irrelevant, while `Δx ≪ ξ` is where the explicit balance earns its keep. Report
   `ξ` at setup as a diagnostic — it is cheap and immediately tells a user whether
   their domain resolves the process they are modelling.
2. **A sanity check on results.** Deposition concentrated within one or two cells of
   a source, when `ξ` spans hundreds of cells, indicates a bug rather than physics.

### 4.4 Deposition

```{index} single: physics label; [D-1]
```
(spec-d-1)=
```{math}
:nowrap:

\begin{align}
D = c_b\, v_s = d^{*}(Z)\, c\, v_s \qquad \text{[P14 3.9], [FG21 §2.2.2]} \tag{D-1}
\end{align}
```

RDy26 A12 uses a different, threshold-based form with a critical *deposition* stress:

```{index} single: physics label; [D-2]
```
(spec-d-2)=
```{math}
:nowrap:

\begin{align}
D_s = \begin{cases} w_s\, c_s\,\left(1 - \dfrac{\tau_b}{\tau_{d,s}}\right) & \tau_{d,s} > 0 \text{ and } \tau_b < \tau_{d,s} \\[4pt] 0 & \text{otherwise} \end{cases} \tag{D-2}
\end{align}
```

Setting `τ_d,s = 0` disables deposition entirely — RDy26 uses this for their passive
transport benchmarks, which is a useful test hook worth preserving.

Both forms take the near-bed concentration to be the *equilibrium* one for the local
flow at every instant. The vertical profile actually adjusts over a time of order
`h/(α w_s)` (Galappatti & Vreugdenhil 1985): grains entrained at the bed must diffuse
up the column before they are carried, and grains high in the column must settle
through it before deposition is felt. Their depth-integrated model relaxes the load
toward the same equilibrium at that rate:

```{index} single: physics label; [D-3]
```
(spec-d-3)=
```{math}
:nowrap:

\begin{align}
E - D &= \alpha\, v_s\, (c_{eq} - c), \tag{D-3} \\
\frac{1}{\alpha} &= \frac{a}{h} + \left(1 - \frac{a}{h}\right)
\exp\!\left[-1.5\left(\frac{a}{h}\right)^{-1/6} \frac{w_s}{u_*}\right] \nonumber
\end{align}
```

(Galappatti & Vreugdenhil 1985; Armanini & Di Silvio 1988), with `c_eq = E*/d*` the concentration at which [E-1]/[E-6] balance [D-1]. In the
implementation both `E` and `D` are scaled by `α/d*`, so every equilibrium of the
instantaneous forms is unchanged and only the transient slows; `α → 1` in the
well-mixed limit `w_s/u* → 0` and `α → h/a` when fully stratified. Off by default;
the evidence for it is van Rijn's pick-up flume (equilibrium reached in ~15 depths
without it, >40 measured) and his migrating trench (fill ~25 % too fast), issue #389.

A second form keeps the *stratification* of the suspension as a state instead of
scaling the rate. Each fraction carries the ratio `r_b = c_b/c` of near-bed to
depth-averaged concentration, advected with the flow and relaxed toward `d*` only when
`d*` has risen above it, which is when the flow has slowed and the profile is collapsing:

```{index} single: physics label; [D-4]
```
(spec-d-4)=
```{math}
:nowrap:

\begin{align}
\frac{\partial (h r_b)}{\partial t} + \nabla\cdot(h r_b \mathbf{u})
  &= h\,\frac{d^{*} - r_b}{T_D} \qquad \text{for } d^{*} > r_b, \tag{D-4} \\
r_b &= d^{*} \qquad \text{otherwise}, \nonumber \\
T_D &= \frac{z_c(Z, a/h) - a}{w_s}, \qquad D = v_s\, r_b\, c \nonumber
\end{align}
```

with `z_c` the centroid height of the equilibrium Rouse profile (a fitted companion of the
`d*` fit). The grains high in the column settle to the bed from about the centroid, so a
parcel entering slower water deposits first at the stratification it brought and works up to
the local one over `T_D`. The lag is one-sided by construction: it never acts in uniform or
quickening flow, where the near-bed concentration is at the bed and responds at once, and it
does not represent the slow filling of the upper column when a bed loads clear water, which
needs a layered suspension. A ratio of zero means "not set" and takes the local `d*`; the
kernel writes the local `d*` on the ratio's boundary edges every step, so inflows carry no
lag in. `T_D → 0` recovers [D-1]. Off by default.

A third form gives the suspension a vertical structure: a near-bed layer of thickness
`h_1 = a` (the reference height, floor included) and the rest of the column, `h_2 = h - a`.
The fraction's own tracer still carries the total mass `m = h c` with its advection, limiters
and bed exchange unchanged; a second tracer carries the upper layer's mass `m_2 = h_2 c_2`,
and the lower layer is the remainder, `m_1 = m - m_2`:

```{index} single: physics label; [T-16]
```
(spec-t-16)=
**[T-16] Bed-shear amplification field.** `set_shear_amplification(factor)` multiplies `f_c`
by a per-centroid factor in both the suspended exchange and the bedload kernel, so `u_*`, the
Shields stress, the Rouse number and every closure built on them see it. Off by default. The
depth-averaged shear is the uniform-flow one and misses the local amplification at obstacles:
at a bridge pier the horseshoe vortex raises the bed shear to two to four times the approach
value within about a diameter (Melville and Raudkivi), and that is what digs the scour hole;
the field is where a case or a structure operator puts it.

```{index} single: physics label; [K-7]
```
(spec-k-7)=
**[K-7] Bedload depth ramp.** `set_bedload(..., min_depth=h_min)` scales the bedload transport
vector by `clamp((h − h_min)/h_min, 0, 1)`: none below `h_min`, full above `2 h_min`. Off by
default. Bedload relations assume a flow many grains deep; at a wetting front over an erodible
bed the film cells carry the front's momentum and a large nominal shear, and unramped they dig
steps that collapse the time step (the Zaragoza pier dam-break case). A few grain diameters is
the physical choice.

```{index} single: physics label; [S-2b]
```
(spec-s-2b)=
**[S-2b] Rouse-number correction.** `rouse_beta='van_rijn'` replaces `Z` by
`Z / β`, `β = min(2, 1 + 2 (w_s/u_*)^2)` (van Rijn 1984b), and `rouse_scale` multiplies the
result. Both enter only the `'rouse'` near-bed fit (`d*`, and the centroid of [D-4]); off by
default.

```{index} single: physics label; [D-5]
```
(spec-d-5)=
```{math}
:nowrap:

\begin{align}
\frac{\partial m_2}{\partial t} + \nabla\cdot(m_2 \mathbf{u}) &= K\,(c_1 - c_2) - v_s\, c_2, \nonumber \\
D &= v_s\, c_1, \qquad E \text{ into layer 1}, \nonumber \\
\rho &= \frac{c_2}{c_1}\Big|_{eq} = \frac{1/d^{*} - a/h}{1 - a/h}, \qquad
K = v_s\,\frac{\rho}{1 - \rho} \tag{D-5}
\end{align}
```

Settling moves sediment down and the exchange `K` up; `K` is set so that the two-layer
equilibrium reproduces the fitted Rouse ratio `c_1/c = d*`, so every equilibrium of [D-1] is
unchanged and only the transient differs. The partition relaxes toward its equilibrium at the
rate `K/h_1 + (K + v_s)/h_2`, integrated exactly per step. Both directions lag: a parcel
entering slower water keeps its upper-layer load and settles it out over about `h_2/v_s`, and
a bed loading clear water fills the lower layer first, so the near-bed concentration and the
deposition lead the depth-averaged load, which grows only as sediment is exchanged up. By
default both layers are advected with the depth-averaged velocity. With `velocity_profile=True`
each layer moves at its log-law mean, `u_1/ū = 1 + ln f_1 / L` and
`u_2/ū = 1 − f_1 ln f_1 / ((1 − f_1) L)` with `L = ln(h/z_0) − 1 = κ/√f_c` and `f_1 = h_1/h`:
the fraction's own tracer (the total) takes the mass-weighted mean of the two and the
upper-layer tracer `u_2`, through a per-cell, per-tracer factor on the advective flux that the
flux kernel applies with the donor cell's value (conservative, since both cells of an edge see
the same product). The sediment transport is then `∫ u c dz < ū h c`, and an inflow that must
carry a given load needs the correspondingly higher concentration. A negative
upper-layer mass means "not set" and takes the equilibrium partition; the kernel writes the
upper-layer tracer's boundary value as the equilibrium partition of the fraction's own
boundary concentration every step. Where `d* > h/a` the lower layer cannot hold the
stratification and the ratio is capped at `h/a`.

**This is the DEFAULT closure**, with `layer_fraction = 0.2` and the velocity split on;
`adaptation='none'` restores the instantaneous exchange {speclit}`D-1`, which was the default
up to ANUGA 4.0. Three flume datasets put the measured adaptation rate at about half the
equilibrium stratification `d*` that {speclit}`D-1` deposits at: Wang & Ribberink's (1986)
settling flume gives 1.44 and 1.55 against `d*` = 2.99, and van Rijn's (1986b) trench and
pick-up flume improve from 3.9, 3.2 and 2.85 cm to 2.45, 2.19 and 1.74 cm, and from a
discrepancy ratio of 0.60 to 0.86.

**Recommendation:** {speclit}`D-1` as default (consistent with the `d*` machinery of §4.3),
{speclit}`D-2` available for the layered bed model and required to reproduce RDy26's
passive-transport validation cases.

⚠ When `d* ≠ 1`, {speclit}`D-1` **must** be evaluated with the near-bed ceiling {speclit}`L-4`
of §4.5(d). Without it the equilibrium `d*` makes the deposition rate diverge as
shear vanishes.

### 4.5 Limiters

Four distinct limiters, easily confused. All four are needed.

**(a) Positivity — physical, mandatory** `[RDy26 9]`:

```{index} single: physics label; [L-1]
```
(spec-l-1)=
```{math}
:nowrap:

\begin{align}
F_s^{\mathrm{net}} = E_s - D_s \;\ge\; -\frac{m_s}{\Delta t} \tag{L-1}
\end{align}
```

Deposition can never remove more suspended sediment than is present. This is a
hard physical constraint, not a stability patch, and it replaces `aS16`'s ad-hoc
concentration clamp.

**(b) Maximum concentration — physical ceiling** `[FG21 §2.2.2]`:

```{index} single: physics label; [L-2]
```
(spec-l-2)=
```{math}
:nowrap:

\begin{align}
c_s \le c_{\max} \qquad \text{FG21: } c_{\max} = 0.30 \text{ by volume} \tag{L-2}
\end{align}
```

FG21 notes this matters mainly for the finest grain sizes (1–2 mm in their runs).
`aS16` uses 0.2. Expose it; do not hard-code.

**(c) Maximum bed change rate — numerical, diagnostic** `[FG21 §2.2.1, Table 2]`:

```{index} single: physics label; [L-3]
```
(spec-l-3)=
```{math}
:nowrap:

\begin{align}
\left|\frac{\partial z}{\partial t}\right| \le \mathrm{max\_dz} \qquad \text{FG21: 5 m/s (deliberately non-physical)} \tag{L-3}
\end{align}
```

FG21 is explicit that this is a stability device set "nonphysically large" so it
"seldom comes into effect", and that when it does fire it acts as an artificial
reduction in transport efficiency. Requirements for our implementation:

- user-settable, default deliberately large;
- **count and report firings** — a run where {speclit}`L-3` fires often is not trustworthy;
- rescale consistently (FG21 rescale erosion in *all* changing cells, preserving
  relative pattern) rather than clipping cell-by-cell.

{speclit}`L-3` is **not implemented in ANUGA**; the label is kept so the gap in
the numbering is not mistaken for an omission. It is unrelated to the `beta`
edge-reconstruction limiter, which the tracers share with the flow variables.

**(d) Near-bed concentration ceiling — physical, mandatory when `d* ≠ 1`**
**[not in the source literature]**:

```{index} single: physics label; [L-4]
```
(spec-l-4)=
```{math}
:nowrap:

\begin{align}
c_b = d^{*}(Z)\, c_s \;\le\; c_{\mathrm{pack}}, \qquad c_{\mathrm{pack}} = 0.65 \tag{L-4}
\end{align}
```

{speclit}`D-1` is `D = c_b v_s`, and **no source bounds `c_b`**. It needs bounding, for a
reason that only appears once {speclit}`S-4` is used inside a transient solver.

`d*` is derived from the **equilibrium** Rouse–Vanoni profile, which presupposes a
suspension maintained by shear. As shear vanishes that premise fails: `u* → 0`, so
`Z = v_s/(κu*) → ∞` {speclit}`S-2`, and `d*` grows without bound — in the fitted
implementation it saturates at its range clamp, `d* ≈ 250` at `a/h = 0.01`. The
deposition rate `d* c v_s` then becomes enormous *precisely where the physical
settling rate should be at its most ordinary*.

Observed before the bound was added: a tilted lake starting from rest, `d_g = 0.1 mm`
(`v_s = 8.0×10⁻³ m s⁻¹`, `h = 1 m`), deposited its **entire** suspended load in under
one second. The physical timescale is `h/v_s ≈ 125 s`. At `t = 1 s` the median Rouse
number was 2.34 while `max|u|` was still 0.11 m s⁻¹ — the flow had barely started.

`c_b` is a concentration, so it cannot exceed maximum packing. `0.65` is the same
constant that bounds `E*` in {speclit}`E-1`, where it appears for exactly this reason, so
{speclit}`L-4` introduces no new parameter. With the bound the same case removes 26.0% in
the first second, matching the packing-limited rate `c_pack v_s / m` to two figures.

Scope and interactions:

- **Inactive in the well-mixed limit.** With `d* = 1`, `c_b = c_s ≤ c_max = 0.30`
  by {speclit}`L-2`, comfortably below `c_pack`. {speclit}`L-4` therefore changes nothing in the
  P14/P13 `d* = 1` configuration, and only engages at moderate-to-high `Z`.
- **Does not conflict with {speclit}`L-1`.** {speclit}`L-4` bounds the deposition *rate*;
  {speclit}`L-1` bounds the *net source* by what is present. {speclit}`L-4` applies first and
  can only reduce a removal, so {speclit}`L-1` remains the binding positivity guarantee.
- **`c_pack` should be exposed**, defaulting to 0.65, alongside `c_max`.

> **Provenance.** {speclit}`L-4` is **not** in P14, FG21, RDy26, DL09 or `aSM16`. It is
> required by the combination of an equilibrium profile ({speclit}`S-4`) with a transient
> solver, which none of those sources does in this form: DL09 derive `d*` for a
> steady sediment-discharge relation, and `aS16` evaluate it at a fixed 1 m depth
> every 10 steps. Flagged here rather than buried in code, because it is a
> modelling addition and future validation against FG21 or `anugaSed` must know it
> is present.

**(e) Non-erodible base — physical, optional**
**[not in the source literature]**:

```{index} single: physics label; [L-5]
```
(spec-l-5)=
```{math}
:nowrap:

\begin{align}
z(x, y, t) \ge z_{\mathrm{base}}(x, y) \qquad z_{\mathrm{base}} \text{ a per-cell field, default } -\infty \tag{L-5}
\end{align}
```

Erosion in {speclit}`G-4`/{speclit}`G-5` is otherwise unbounded below: the bed lowers for as long
as the flow can lift material. That is right for a deep alluvial bed and wrong
wherever the erodible layer is finite — a reach floored by an outcrop, a lined
culvert, a dam apron, a soil layer of known depth over rock. `z_base` is a **field,
not a constant**, because bedrock is a surface.

Applied to the **source**, never by clamping `z`. The erodible thickness
`T = z − z_base` bounds the net removal over a step, so for the suspended exchange

```{math}

\sum_s F_s^{\mathrm{net}} \;\le\; \frac{T\,(1 - \lambda)}{\Delta t}
```

and the erosive part of the source is scaled to satisfy it. Sediment that is not
eroded never enters the water column, so {speclit}`G-3` and {speclit}`G-4` still use the same
limited source and the budget closes exactly, as it does for {speclit}`L-1` and {speclit}`L-2`.
Clamping `z` afterwards would instead leave suspended sediment that came from
nowhere — the same failure mode as the {speclit}`L-1` sign bug (spec 12, D-obs).

Scope and interactions:

- **Shared across classes, proportionally.** `T` is a property of the cell, not of a
  class, so the classes are limited *together*: their erosive sources are scaled by
  one common factor. Serving them in registration order would make the answer depend
  on the order `add_sediment_class` was called, which is not physics. By default
  the bed carries no per-class composition, so no class has a stronger claim on the
  last millimetre; proportional is the only choice that invents nothing. **With bed
  composition (§5.0) this rule changes:** each class is capped at what the cell holds
  of it, in both layers, and the shared scale no longer binds.
- **Deposition is never scaled.** A shortage of bed material does not restrain
  deposition — deposition is what relieves it.
- **Ordered after {speclit}`L-1` and {speclit}`L-2`.** Those bound the source by what the water
  column holds and can hold; {speclit}`L-5` bounds it by what the bed holds. It can only
  reduce an erosive source, so it cannot defeat {speclit}`L-1`'s positivity.
- **Bedload {speclit}`K-3`/{speclit}`G-5` is limited differently, and less strictly.** Bedload is a
  divergence: clipping it cell-by-cell would create bed material, because one side of
  an edge would refuse to give up what the other has already received. The limit is
  therefore applied to the transport vector `q_b` and to whole edges — both of which
  the two cells sharing an edge evaluate identically, so the flux stays antisymmetric
  and {speclit}`K-3` stays exactly conservative. The cost is that the floor is **not exact
  under bedload**: closing an edge for a cell that cannot pay also cancels its
  neighbour's inflow, so the deficit migrates one cell per sweep. Measured overshoot
  5.1×10⁻⁶ m on a 1.0×10⁻² m layer (14 cells of 960). Driving it to zero requires the
  exhaustion flag iterated to a fixed point with double buffering, so the result stays
  independent of thread order; **not implemented**.

**Region restriction.** *Where* the bed may erode is the same constraint as *how
deep*, so it is expressed through the same field rather than as a second
mechanism: a locked cell is one whose erodible thickness is zero, `z_base = z` at
the moment it is locked. It follows that locked means **unscourable, not inert** —
deposition still lands on it, and that new material is erodible again, being above
the base. Under erosive flow such a cell settles at exactly net zero, the limiter
scaling erosion back until it cancels deposition. This matches ANUGA's existing
region-based erosion operators in interface (the same `Region` polygon/circle/
indices arguments) while differing from them in kind: those clamp `z` directly,
which is sound only because they carry no sediment budget to violate.

> **Provenance.** {speclit}`L-5` is **not** in P14, FG21, RDy26, DL09 or `aSM16`, none of
> which limits the erodible depth. It is an engineering requirement rather than a
> closure, and it is off by default, so every configuration in this spec is
> unaffected unless a base is set.

---

## 5. Bed layer model

### 5.0 What is implemented: one active layer over one substrate

`Domain.set_bed_composition` (4.1.0) implements the active-layer part of this model,
with **one** substrate layer rather than a stack. The state is the solid volume per
unit bed area of each class in each layer, `A_s` (active) and `B_s` (substrate)
[m], filling the erodible thickness:

    (1/(1-λ)) Σ_s (A_s + B_s) = z - z_base

Each step, in `core_apply_sediment_source` (both compute modes):

1. Entrainment by share: `E_s = F_s E*_s`, `F_s = A_s / Σ_j A_j` at the start of the
   step; an empty active layer supplies nothing.
2. Per-class cap, after L-1 and L-2: `S_s ≤ (A_s + B_s) / (Δt M)`. This replaces the
   shared proportional scale of L-5 (§4.5).
3. Exchange through the active layer, `A_s ← A_s − M Δt S_s` (a class driven below
   zero takes the rest from its substrate), then restore its solid content
   `L_A = H_A (1−λ)`. Aggradation sends the overflow down at the **active**
   composition (§5.3, with one substrate); degradation refills from the substrate at
   the **substrate** composition (§5.2, with one layer). Every transfer moves an
   amount of one class, so each class is conserved exactly.

Not implemented: the substrate stack (§5.1-5.2 with more than one layer), mass per
area in kg m⁻² (volume is used, with porosity λ), composition carried by bedload
(K-3) and by angle-of-repose relaxation (§7), and a mixture skin roughness. The rest
of this section is the original RDycore-derived specification. User-facing
description: *Sediment transport* setup page, *Bed composition*; equations:
*Sediment physics* appendix.

From `[RDy26 §2.2, A1–A3]`. One **active layer** exchanging with the water column,
over up to 32 **substrate layers**. Per class per layer the state is mass per unit
area `M_ℓ,s` [kg m⁻²].

### 5.1 Initialisation

```{index} single: physics label; [B-1], single: physics label; [B-2], single: physics label; [B-3], single: physics label; [B-4]
```
(spec-b-1)=
(spec-b-2)=
(spec-b-3)=
(spec-b-4)=
```{math}
:nowrap:

\begin{align}
\sum_s f_s &= 1 && \text{class fractions} \tag{B-1} \\
M_{\ell,s}(0) &= f_s\, C_\ell\, H_\ell && C_\ell \text{ layer concentration, } H_\ell \text{ thickness} \tag{B-2} \\
\rho_{\mathrm{grain}}(0) &= \sum_s f_s\, \rho_s && \text{bulk particle density} \tag{B-3} \\
\phi_\ell &= 1 - \frac{C_\ell}{\rho_{\mathrm{grain}}} && \text{layer porosity} \tag{B-4}
\end{align}
```

### 5.2 Erosion through the stack

Layers erode **sequentially**: the active layer first; if it is exhausted within the
timestep, the remaining time is applied to the next substrate layer, and so on.

```{index} single: physics label; [B-5], single: physics label; [B-6], single: physics label; [B-7]
```
(spec-b-5)=
(spec-b-6)=
(spec-b-7)=
```{math}
:nowrap:

\begin{align}
\Delta M_\ell^{\mathrm{ero}} &= \min\!\left(M_\ell,\; E_\ell^{\mathrm{pot}}\, \Delta t_\ell\right) && \Delta t_\ell = \text{time layer } \ell \text{ is exposed} \tag{B-5} \\
\Delta M_{\ell,s}^{\mathrm{ero}} &= \Delta M_\ell^{\mathrm{ero}}\, \frac{M_{\ell,s}}{M_\ell} && \text{partition by composition} \tag{B-6} \\
E_s &= \frac{1}{\Delta t} \sum_\ell \Delta M_{\ell,s}^{\mathrm{ero}} \tag{B-7}
\end{align}
```

with `E_ℓ^pot` from {speclit}`E-4` and `M_ℓ = Σ_s M_ℓ,s` `[RDy26 7]`.

Erosion is disabled in shallow water `[RDy26 A8]`:

```{index} single: physics label; [B-8]
```
(spec-b-8)=
```{math}
:nowrap:

\begin{align}
h < h_{\min,\mathrm{ero}} = \max(0.1\ \mathrm{m},\; 10\, h_{\mathrm{wet}}) \quad\Longrightarrow\quad \text{no erosion} \tag{B-8}
\end{align}
```

where `h_wet` is the solver's wet/dry threshold. **This threshold matters** — it is
far more conservative than `aS16`'s `h > 0.1` gate and prevents spurious erosion in
thin films where velocity is poorly resolved.

### 5.3 Deposition and active-layer update

```{index} single: physics label; [B-9], single: physics label; [B-10], single: physics label; [B-11]
```
(spec-b-9)=
(spec-b-10)=
(spec-b-11)=
```{math}
:nowrap:

\begin{align}
\Delta M_{A,s}^{\mathrm{dep}} &= \max\!\left(0,\; (E_s - F_s^{\mathrm{net}})\,\Delta t\right) \tag{B-9} \\
H_\ell &= \frac{M_\ell}{\bar\rho_\ell\,(1 - \phi_\ell)} \tag{B-10} \\
\bar\rho_\ell &= \sum_s \frac{M_{\ell,s}}{M_\ell}\, \rho_s \tag{B-11}
\end{align}
```

After exchange, if the active layer exceeds its prescribed thickness `H_A`, the
excess is transferred down to the first substrate layer; if thinner, mass is drawn
up until refilled or the substrate is exhausted. This conserves class-specific bed
mass locally and **does not itself change bed elevation**. Active-layer concentration
is updated by mixing with the donor layer `[RDy26 A16]`.

### 5.4 Memory

```{math}

\text{bytes} = N_{\mathrm{classes}} \times N_{\mathrm{layers}} \times N_{\mathrm{cells}} \times 8
```

5 classes × 32 layers × 5×10⁶ cells = **6.4 GB** for the bed alone, before
hydrodynamics. The API must surface this estimate at setup and fail with a clear
message rather than a mid-run allocation failure. Cap `N_layers`.

---

## 6. Bedload

From `[FG21 §2.2.1]`, Eqs 1–6. Uses the shear stress chain of §3 directly.

```{index} single: physics label; [K-1], single: physics label; [K-2], single: physics label; [K-3]
```
(spec-k-1)=
(spec-k-2)=
(spec-k-3)=
```{math}
:nowrap:

\begin{align}
q_b^{*} &= K\, \tau_x^{\,m} && \text{[FG21 4]} \tag{K-1} \\
q_b &= q_b^{*}\, \sqrt{(\rho_s/\rho - 1)\, g}\; D^{1.5} && \text{[FG21 5]} \tag{K-2} \\
\frac{\partial z}{\partial t} &= -\frac{1}{1 - \lambda}\,\nabla\cdot\mathbf{q}_b && \text{[FG21 6]} \tag{K-3}
\end{align}
```

Parameter sets for {speclit}`K-1`:

| Set | `K` | `m` | `τ_c*` | Notes |
|-----|-----|-----|--------|-------|
| Wong & Parker (2006) Eq 24 | 3.97 | 1.5 | 0.0495 | MPM-tradition, reanalysed. Bedload only. Their Eq 23 is the alternative free-exponent fit (`K` = 4.93, `m` = 1.60, `τ_c*` = 0.0470) — confirm which FG21 used. |
| Engelund & Hansen (1967) `[EH67 4.3.5]` | `0.05/f_c` | 2.5 | 0 (no threshold) | **Total load** — includes suspension implicitly. See {speclit}`K-5`. |

> **Critical usage rule** `[FG21 §2.2.2]`: the suspended operator is **not** applied
> when using Engelund & Hansen, because E&H already incorporates suspended
> transport. Running both double-counts. The API must enforce this.

**Engelund & Hansen total load** `[EH67 4.3.5]` —

EH67 write their transport relation as

```{math}

f_{EH}\, \Phi = 0.1\, \theta^{5/2} \qquad \text{[EH67 4.3.5]}
```

where `Φ = q_T/√((s−1)g d³)` is dimensionless **total** sediment discharge, `θ` is
the Shields stress, and `f_EH` is EH67's own friction factor from their Eq 3.1.3.
Converting to this spec's convention using §3.3.1 (`f_EH = 2 f_c`):

```{index} single: physics label; [K-5]
```
(spec-k-5)=
```{math}
:nowrap:

\begin{align}
q_b^{*} = \frac{0.05\, (\tau^{*})^{5/2}}{f_c} \qquad \text{total load, no threshold} \tag{K-5}
\end{align}
```

so in the {speclit}`K-1` power-law form, `K = 0.05/f_c` and `m = 2.5`, with `τ_c* = 0`.

> **Open item K1 — RESOLVED.** `K` is indeed friction-dependent, not constant, and
> the familiar literature coefficient `0.05` is EH67's `0.1` after the factor-of-2
> convention conversion. Two cautions carried forward:
>
> 1. **Do not apply a threshold.** {speclit}`K-5` has none; subtracting `τ_c*` would be a
>    different model. This is also why EH67 must replace, not supplement, the
>    suspended operator.
> 2. **`θ` vs `θ'`.** EH67 derive {speclit}`K-5` via the *effective* (skin-friction) bed
>    shear `θ'` — their Eq 4.2.4, `θ' − 0.06 = 0.4 θ²` — but the final empirical fit
>    4.3.5 is stated in terms of **total** `θ`, which is what {speclit}`K-5` uses. Do not
>    mix the two.

**Grass** —

The law the classic Exner test cases are written in (Hudson & Sweby 2003), and
the one with a closed-form reference for a migrating bed hump:

```{index} single: physics label; [K-6]
```
(spec-k-6)=
```{math}
:nowrap:

\begin{align}
\mathbf{q}_b = A_g\, |\mathbf{u}|^{m-1}\, \mathbf{u} \qquad \text{total load, no threshold} \tag{K-6}
\end{align}
```

with `m = 3` and `A_g` (s²/m, 0 to 1) a calibration with no default. It is
dimensional as written, so {speclit}`K-2` does not apply and no grain size
enters; like {speclit}`K-5` it is total load and replaces the suspended
operator.

**Boundary edges.** A boundary edge carries no bedload unless its tag has been
declared open (`set_bedload(..., open_boundaries=...)`), in which case the
flux across it is the cell's own `q_b · n` (zero gradient): an outflow carries
bedload away at the rate it arrives, an inflow supplies it at the rate the
first cell carries it off. Walls stay closed, which is what keeps a closed
domain exactly conservative; a reach's inflow must be opened or its first cell
exports and never imports.

**Vectorisation** `[FG21 §2.2.1]`: `q_b` is a scalar magnitude and must be
partitioned into components. FG21 follow Parker (1998) Eqs 2.11–2.12: **the sediment
transport vector is parallel to the boundary shear stress vector**. With {speclit}`T-1`
that means parallel to `(u, v)`:

```{index} single: physics label; [K-4]
```
(spec-k-4)=
```{math}
:nowrap:

\begin{align}
q_{bx} = q_b\, \frac{u}{|\mathbf{v}|}, \qquad q_{by} = q_b\, \frac{v}{|\mathbf{v}|} \tag{K-4}
\end{align}
```

> **Known defect to avoid.** FG21 report ~1% sediment mass conservation error in
> their bedload operator, "equivalent to a sediment volume lost by advection out of
> the domain without intermediate deposition", which they were unable to correct.
> Computing {speclit}`K-3` as a **flux difference across shared edges** — the same
> `q_b·n·L` debited from one cell and credited to its neighbour — makes conservation
> structural. Do not compute `∇·q_b` from a reconstructed cell-centred gradient.

---

## 7. Angle-of-repose relaxation

From `[FG21 §2.2.4]`. Where the cell-to-cell bed slope exceeds a threshold angle
(FG21 use 35°), topography is "near-instantaneously diffused until it falls below
the threshold".

FG21 are explicit that this is **a numerical heuristic, not physics** — real bed
slope failures are advective. It exists to stop other parts of the model breaking on
over-steep slopes, and it has a side effect worth recording: it limits the steepness
of canyon walls and knickpoints, suppressing knickpoint retreat that may be real.

Implementation requirements:

- iterative neighbour diffusion — **the only non-cell-local sediment kernel**;
- MPI halo exchange **per sweep**;
- hard iteration cap, with a report when the cap is hit;
- conserve mass: material removed from over-steep cells is deposited on neighbours,
  never discarded.

### 7.1 Implementation status

Implemented as `core_apply_repose`, enabled by `domain.set_angle_of_repose(angle)`
and **off by default** — it is a heuristic, and it suppresses knickpoint retreat.

Three of the four requirements are met as written. Mass conservation is exact
(measured drift 2.3×10⁻¹³ m³ on 4.0×10² m³ of bed) and structural rather than
incidental: the transfer is computed per **edge** from data both cells share, so
both compute the identical volume,

```{math}

V = \frac{\mathrm{relax}\,\left(|\Delta z| - \tan\theta_c\, d\right)}{3\,\left(1/A_k + 1/A_{nb}\right)}
```

and the pair balances by construction. This is the same device that keeps {speclit}`K-3`
conservative, and it is where this differs from ANUGA's `sanddune_erosion_operator`,
which lowers an over-steep cell and lets the material vanish.

Sweeps are **Jacobi**, not Gauss–Seidel: each sweep reads elevation and writes
increments to a separate array, so no cell sees a neighbour that has already moved.
Reading live elevation would make the result depend on thread order and the two
compute modes would disagree.

Two properties found by measurement rather than assumed, both recorded because they
are not obvious from FG21's one-sentence description:

- **The `/3` above is a stability limit.** A cell has three edges and can be the
  donor on all of them, so without it a relaxation factor of 1 moves up to three
  times what any single edge intended — an explicit diffusion step past its
  stability limit. It overshoots, creates fresh over-steep edges, and oscillates.
- **Convergence is asymptotic and slow.** An over-steep cone (36.8°) needed **793
  sweeps** to reach a 30° limit from cold. Convergence is therefore declared within
  a 10⁻³ relative tolerance on the threshold slope; a strict test never terminates.
  In a running model the bed is already near-relaxed and a step needs a handful of
  sweeps, so the cap is for the pathological case. Hitting it is not fatal —
  progress carries over between timesteps — but it is reported, as required.

**Deviation from the requirements: the halo exchange is per TIMESTEP, not per
sweep.** The sweep loop lives inside the kernel so that it stays on the device;
lifting it into Python to exchange between sweeps would put a host round trip in
the middle of every sweep and drop the run off the GPU path entirely. The
consequence in parallel is that relaxation crosses a subdomain boundary one sweep
per timestep rather than one per sweep, so a slump spanning a boundary relaxes more
slowly there than it would in serial. Serial and single-subdomain results are
unaffected. Resolving this properly needs a device-side halo exchange callable from
within the sweep loop.

---

## 8. Vegetation drag

From `[aSM16 §4.2–4.3, Eqs 12–14]`, after Kean & Smith (2004) and Nepf (1999).
`aSM16` is fuller than P14 Eq 3.11 and should be preferred.

### 8.1 Drag force

```{index} single: physics label; [V-1], single: physics label; [V-2]
```
(spec-v-1)=
(spec-v-2)=
```{math}
:nowrap:

\begin{align}
F_D &= \tfrac{1}{2}\, \rho\, C_D\, \alpha\, U_{\mathrm{ref}}^2 && \text{[aSM16 12]} \tag{V-1} \\
\alpha &= \frac{d_s}{\lambda^2} && \text{projected plant area per unit volume } [\mathrm{m^{-1}}] \tag{V-2}
\end{align}
```

`U_ref` is the flow velocity in the absence of vegetation, `d_s` the stem diameter,
`λ` the mean stem spacing.

> P14 Eq 3.11 omits `ρ` from {speclit}`V-1`. `aSM16` Eq 12
> includes it. The manual's form is dimensionally correct and is the one to implement.

### 8.2 Velocity reduction

```{index} single: physics label; [V-3]
```
(spec-v-3)=
```{math}
:nowrap:

\begin{align}
U = U_{\mathrm{ref}} - F_D\, \Delta t \qquad \text{[aSM16 13]} \tag{V-3}
\end{align}
```

**Constraint stated in `aSM16`:** the drag force may reduce the flow velocity to
zero but must **not reverse its direction**. Clamp at zero; do not allow sign change.
This is a hard requirement, and an obvious source of instability if missed.

### 8.3 Drag coefficient as a function of stem density

For a sparse regular array of emergent cylinders `C_D = 1.2`. `aSM16` follows
Nepf (1999) to capture the effect of stem population density, via the force balance
`[aSM16 14]`:

```{index} single: physics label; [V-4]
```
(spec-v-4)=
```{math}
:nowrap:

\begin{align}
(1 - \alpha d)\, C_B\, U^2 + \tfrac{1}{2}\, C_D\, \alpha\, d\, \frac{h}{d}\, U^2 = g\, h\, \frac{\partial h}{\partial x} \tag{V-4}
\end{align}
```

where `C_B` is the bed drag coefficient. `aSM16` fits `∂h/∂x` as a cubic in `αd`
matched to Figure 6 of Nepf (1999), and solves {speclit}`V-4` to obtain:

```{index} single: physics label; [V-5]
```
(spec-v-5)=
```{math}
:nowrap:

\begin{align}
C_D = \begin{cases} 1.2, & x \le 0.006, \\[4pt] 56.11\,x^2 - 15.28\,x + 1.3 - 5.465\times 10^{-4}\,x^{-1}, & x > 0.006, \end{cases}
\tag{V-5}
\end{align}
```

with `x = αd`. `C_D` is evaluated per cell whenever the stem-diameter or stem-spacing fields change.

> **Caveat.** {speclit}`V-5` is a fit to a fit — a cubic fit to a digitised figure, fed
> through a force balance. It is reasonable within the range of Nepf's experiments
> and should not be extrapolated. Note the `(αd)⁻¹` term diverges as `αd → 0`, so the
> `αd ≤ 0.006` branch is doing real work, not just simplification. Guard the boundary.

Field values from P14 (Rio Puerco, sandbar willow and young tamarisk, after Griffin
et al. 2005): `d_s = 0.02 m`, `λ = 0.13 m`.

This operator is cell-local and drops directly into the L2 kernel pattern.

---

## 9. Numerical scheme

### 9.1 Tracer flux — upwind on the water flux

This is the L1 design, and `[RDy26 §2.3, 12–13]` specifies it exactly:

```{index} single: physics label; [N-1], single: physics label; [N-2]
```
(spec-n-1)=
(spec-n-2)=
```{math}
:nowrap:

\begin{align}
c_s^{\mathrm{up}} &= \begin{cases} c_{s,L} & \hat F_{\mathrm{water}} \ge 0 \\ c_{s,R} & \hat F_{\mathrm{water}} < 0 \end{cases} \tag{N-1} \\
\hat F_s &= \hat F_{\mathrm{water}}\; c_s^{\mathrm{up}} \tag{N-2}
\end{align}
```

RDy26 states plainly: *"This upwind sediment-transport scheme is conservative for
the transported sediment mass m_s = h c_s."*

Key ordering constraint `[RDy26 A4]`: the hydrodynamic flux is computed **first**,
from `h, u, v` only; sediment concentrations are used only afterwards. This
guarantees sediment cannot alter the hydrodynamics, which is what makes the one-way
coupling of Phase 3 exact rather than approximate.

For ANUGA the same principle applies to whichever Riemann solver is in use — take
the already-computed edge water flux and multiply by the upwind concentration. The
flux is shared between the two cells sharing the edge, so mass balance is exact
regardless of cell size ratio.

### 9.2 Semi-discrete update

`[RDy26 11]`, matching ANUGA's existing finite-volume form:

```{index} single: physics label; [N-3]
```
(spec-n-3)=
```{math}
:nowrap:

\begin{align}
\frac{d\mathbf{U}_i}{dt} = -\frac{1}{A_i} \sum_{e\in\partial i} \sigma_{i,e}\, L_e\, \hat{\mathbf F}_e + \mathbf{S}_i \tag{N-3}
\end{align}
```

### 9.3 Boundary conditions

`[RDy26 §2.3]`: apply hydrodynamic BCs to `h, hu, hv`; apply sediment BCs through
the **conserved** tracer `h·c_s`. A prescribed inflow concentration is converted to
`h·c_s` using the boundary depth. Dirichlet (prescribed concentration) is the
minimum required set.

### 9.4 Wet/dry

Below the wet/dry threshold set `u = v = c_s = 0` `[RDy26 A21–A23]`. Note the
separate, more conservative erosion threshold {speclit}`B-8`.

### 9.5 Well-balancedness

RDy26 added hydrostatic reconstruction specifically because sediment amplifies
hydrodynamic error: *"small hydrodynamic errors can be amplified through
shear-stress-dependent erosion and deposition"* — erosion goes as `τ^m` with
`m` up to 2.5. Since `τ ∝ v²`, transport goes as `v^2m` — at `m = 2.5` that is `v⁵`,
so a 10% velocity error becomes a **61%** transport error (1.1⁵ ≈ 1.61). Verify
ANUGA's existing reconstruction is adequate under an erodible bed as part of Phase 4.

`LM15` is the constructive counterpart to that warning, and on the same mesh
topology (triangular, unstructured). It shows that once `Z` becomes
time-varying the bed-slope term needs its **own** discretisation, matched to the
reconstruction, for the lake-at-rest state to stay exact — a fixed-bed
well-balanced scheme does not remain well-balanced for free when the bed starts
to move. Their Example 2 (sand deposition into quiescent water) is the test that
isolates this: the bed aggrades while the free surface must remain flat and the
velocity exactly zero. **Use it as the Phase 4 acceptance test** — it is cheap,
needs no field data, and fails loudly if the coupling breaks well-balancedness.

Note the scope limit: `LM15` derives this for a *fully coupled* system, whereas
§9.1 here advects `c` with the already-computed water flux. The bed-slope
requirement is a property of a moving `Z`, not of the coupling, so it applies
either way — but their specific discretisation is not transferable as-is.

---

## 10. Parameter defaults

Defaults are FG21 Table 2 / P14 §3.3.2 where stated. **All must be settable.**

| Parameter | Symbol | Default | Source |
|-----------|--------|---------|--------|
| Water density | ρ | 1000 kg m⁻³ | FG21 T2 |
| Sediment density | ρ_s | 2650 (P14) / 3000 (FG21, Mars) kg m⁻³ | P14, FG21 T2 |
| Porosity | λ | 0.2 (FG21) / 0.3 (aS16) | FG21 T2 |
| Kinematic viscosity | ν | 1×10⁻⁶ m² s⁻¹ | P14 |
| von Kármán | κ | 0.41 | FG21 |
| Grain size | D | 0.1 mm (P14 suspension) … 8 mm–1 cm (FG21 bedload) | P14, FG21 T2 |
| Settling C₁, C₂ | | 18, 0.4 (smooth spheres) | P14 3.10 |
| Suspension critical stress | τ_c* | 0.04 | FG21 |
| Smith–McLean constant | γ₀ | 0.0024 | FG21 |
| Bedload critical stress | τ_c* | 0.0495 | Wong & Parker 2006 |
| Bedload coefficient | K, m | 3.97, 1.5 | Wong & Parker 2006 |
| Max concentration | c_max | 0.30 | FG21 |
| Max bed change rate | max_dz | large; report firings | FG21 T2 |
| Angle of repose | | 35° | FG21 T2 |
| Vegetation drag | C_D | 1.2 | P14 3.11 |
| Sediment speed ratio | β | 1.0 | P14 3.4 |
| Active layer thickness | H_A | user | RDy26 |
| Substrate layers | | ≤ 32 | RDy26 |
| Gravity | g | **from domain** (9.81 Earth, 3.711 Mars) | FG21 T2 |

---

## 10.1 Validation

This specification says what the code must compute. Whether those are the right
equations for a real flow is a separate question. The cases in
`validation_tests/sediment/` answer it, and they run nightly with the rest of the
validation suite (see the *Validation Test Suite* page for each case and how to run it):

- **Against closed-form solutions:**
  - `sediment_settling`: deposition {speclit}`D-1`, the mass balance and the bed update.
  - `sediment_erosion`: Smith–McLean entrainment {speclit}`E-1`, its threshold, and
    the depth-slope closure {speclit}`T-7`.
  - `sediment_settling_basin`: the steady concentration profile, the bed-rise rate and
    the sediment budget.
  - `sediment_equilibrium_flow`: clear water approaching the equilibrium
    concentration in normal flow.
  - `sediment_bed_hump`: Grass bedload and the Exner equation, against the
    characteristic solution.
- **Against laboratory measurements:**
  - `van_rijn_pickup_flume`: entrainment {speclit}`E-6` and the near-bed ratio.
  - `van_rijn_trench`: the settling lag, deposition and the Exner coupling over a
    dredged trench.
  - `wang_ribberink_settling`: deposition alone, with nothing else acting. It scores
    the two-layer closure {speclit}`D-5` and runs {speclit}`D-3` and {speclit}`D-4`
    for comparison.
  - `zaragoza_pier_dambreak`: transient bedload over a finite sand layer, the angle of
    repose, and scour at a pier.

No field-scale case has been built. The original specification planned two: the
Rio Puerco 2006 flood `[P13]` and a Mars crater breach `[FG21]`. The Rio Puerco
values `P13` gives, which are enough to build that case, are in
`validation_tests/sediment/rio_puerco_2006/README.md`.

---

## 11. Equation citation map

Implementation should cite these labels. Suggested mapping to the layer plan:

| Layer | Sections |
|-------|----------|
| **L1** core tracer kernel | §2.2 `[G-1..G-3]`, §9 `[N-1..N-3]` |
| **L2** device bed operator | §3 `[T-1..T-6]`, §4 `[E-*, S-*, D-*, L-*]`, §5 `[B-*]`, §6 `[K-*]`, §8 `[V-*]` |
| **L2** separate kernel | §7 angle of repose (non-cell-local) |
| **L3** Python API | §10 parameters, closure selection, §4.3 pluggable `d*` |

---

## 14. References

The full reference list, with one entry per short label, is kept in one place:
{ref}`References <sediment_references>` on the companion page.
The sources that only this specification draws on, and that the companion
page does not cite, are:

- Ganti, V. et al. (2014) — cited by FG21 for the `d*` concept.
- Griffin, E.R. et al. (2005, 2014) — Rio Puerco vegetation and grain size.
- Hanson, G.J. (1990) — excess-shear entrainment for cohesive soils.
- Kean, J.W. & Smith, J.D. (2004) — vegetation drag.
- Nepf, H.M. (1999), *Drag, turbulence, and diffusion in flow through emergent vegetation*, Water Resources Research 35(2), 479–489. Figure 6; source of {speclit}`V-4`,{speclit}`V-5`.
- Nordin, C.F. (1963) — uniform suspension observations.
- Parker, G. (2005); Kleinhans, M.G. (2005) — transport law reviews.
- Simpson, G. & Castelltort, S. (2006) — coupled flow/sediment/topography framework cited by aSM16.
