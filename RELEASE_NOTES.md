# ANUGA 4.1.0

ANUGA 4.1.0 adds two new capabilities to the solver: **passive tracers**,
concentrations carried with the water, and **sediment transport** built on
them, with suspended load, bedload and an evolving bed. Both run in both
compute modes, offload to a GPU, and work in parallel. Vegetation drag joins
the friction, and every release now ships its validation reports.

Most existing scripts run unchanged. The exceptions are scripts that still use
the forcing classes deprecated in 4.0.0, which are now removed, and runs with
the `DE_ader2` flow algorithm, whose results change (see *Breaking changes*).

## Highlights

### Passive tracers

A tracer is a depth-averaged concentration `c` carried as the conserved
quantity `h c` through the same flux kernel as the water: dye, salt, a
pollutant, or a sediment fraction.

```python
domain.add_tracer('salinity')
domain.set_tracer('salinity', 35.0)                    # initial field
domain.set_tracer_boundary('salinity', 'ocean', 35.0)  # inflow values per tag
domain.set_tracer_boundary('salinity', 'river', 0.0)
```

Tracers have their own boundary values, are written to the SWW as
`<name>_c`, are distributed and reordered with the mesh in parallel, and come
with mass diagnostics (`get_tracer_mass`,
`get_tracer_boundary_flux_integral`, `check_tracer_conservation`). The
operators that add or remove water carry them too: inlets, structures
(culverts, weirs, bridges) and rate operators (see *Selected fixes*). See the
*Passive Tracers* page of the setup guide.

### Sediment transport

One or more sediment fractions, each a tracer, with suspended load
(entrainment, settling and the adaptation of the suspension), bedload, and
bed evolution through the Exner equation, limited by an erodible base and an
angle of repose. One line turns it on; the rest is optional:

```python
domain.add_sediment_fraction('sand', diameter=2.0e-4)
```

It offers a choice of entrainment laws (Smith-McLean by default, de Leeuw et
al. 2020), deposition and adaptation closures, bed-shear closures, and
bedload laws including Grass. A morphological factor reaches long
timescales, and the bed shear sees vegetation and can be amplified locally
at structures. It can all be set up from a TOML scenario's `[sediment]`
table. The physics is written out with references in the *Sediment physics*
appendix, and nine validation cases, five analytical and four laboratory
flumes, have their own chapter in the validation report. See the *Sediment
transport* page of the setup guide.

### Bed composition

`domain.set_bed_composition({'mud': f, 'sand': 1 - f}, active_layer=0.01)`
tracks how much of each sediment fraction the bed holds, in a thin active
layer over a substrate (Hirano). A fraction is entrained only in proportion to
its share of the active layer and never beyond what the cell holds, and a bed
that loses its fines armours. Without it every fraction comes out of one shared
bed, so mud is picked up from sand-bedded channels for as long as the bed
lasts: in the Delta-X Wax Lake case suspended mud ran away within a day, and
with composition the eroded area fell 15-fold and the median suspended mud
stayed at the inflow level. Each fraction is conserved exactly, in both
compute modes and in parallel; it is off by default. See *Bed composition* on
the *Sediment transport* page.

### Vegetation drag

`domain.set_vegetation_drag(density, diameter, height)` adds the drag of a
stem field, emergent or submerged, with the vegetated Chezy coefficient of
Baptist et al. (2007), in the friction kernel of both compute modes (see
*Smaller improvements*).

### Validation reports with every release

Publishing a release now builds the validation reports from the tagged
commit: one combined report that shows every case with the `DE0`, `DE1` and
`DE_ader2` flow algorithms side by side, and the full report for each. They
are attached to the GitHub release and archived on Zenodo with a DOI. The
validation suite also runs nightly in CI, and `rundown_mild_slope` now shows
the roll-wave instability of fast shallow flow on a slope, compared across
the flow algorithms.

## Breaking changes

* **The legacy forcing-function classes have been removed**, as announced in
  4.0.0: `Rainfall`, `Inflow`, `Wind_stress`, `Wind_stress_fast`,
  `Barometric_pressure`, `Barometric_pressure_fast` and `General_forcing` are
  gone from `anuga.shallow_water.forcing`, and `anuga.Inflow`, `anuga.Rainfall`
  and `anuga.Wind_stress` are no longer exported. Use `Rate_operator.rainfall()`
  / `Rate_operator.inflow()`, `Wind_stress_operator` and
  `Barometric_pressure_operator`, which work in both compute modes. Note the
  rainfall units change from mm/s to mm/hr. See `UPGRADING_TO_4.0.md` §2.
  `Cross_section` stays in `anuga.shallow_water.forcing`.

* **The suspended-sediment deposition closure has a new default.**
  `set_deposition()` now selects the two-layer suspension with the log-law
  velocity split (`adaptation='two_layer'`, spec [D-5]/[D-5v]), where before
  it deposited instantaneously at the equilibrium near-bed ratio `d*`
  ([D-1]). `set_deposition(adaptation='none')` restores the old behaviour
  exactly.

  The old default used `d*`, the *equilibrium* stratification of the
  suspension, as the *rate* at which it adapts. Those are different numbers,
  and the ratio between them grows with the settling velocity relative to the
  shear: 1.2 at `w_s/u* = 0.05`, 2.0 at 0.2, 4.1 at 1.0. Three flume datasets
  agree. Wang and Ribberink's (1986) settling flume, which has deposition and
  nothing else, measures a decay coefficient of 1.44 and 1.55 where their own
  theory gives 1.47 and ANUGA gave 2.99.

  | validation case | old default | new default | measured |
  |---|---|---|---|
  | Wang & Ribberink, decay coefficient | 2.99 | 1.50 | 1.44, 1.55 |
  | van Rijn trench, bed RMS at 15 h (cm) | 3.9, 3.2, 2.85 | 2.45, 2.19, 1.74 | |
  | van Rijn pick-up flume, discrepancy ratio | 0.60 | 0.86 | 1.0 |

  **Results change** for any model whose near-bed ratio differs from 1, which
  means any run using `near_bed='rouse'`. They are unchanged where `d* = 1`,
  the well-mixed limit that `near_bed='constant'` gives by default, because
  the two-layer partition reduces to the old closure exactly there; the five
  analytical sediment validation cases all sit in that limit and are
  untouched. Each fraction gains one tracer, `<name>_upper`, registered at
  the first `evolve` and appearing in the SWW output.

* **`DE_ader2` results change** (#404). The ADER-2 predictor had no friction,
  so on friction-dominated flow it over-drove the state the fluxes were taken
  from: on a 1:10 slope with 3.9 mm of water it carried 0.45 of the
  prescribed discharge where DE0 and DE1 are exact. Manning drag is now
  applied in the predictor, semi-implicitly. The default limiter beta for
  `DE_ader2` rises from 0.5 to 1.0, because beta 0.5 flattens the
  reconstructed bed on a linear slope and damps propagating waves (32% of
  the forced amplitude kept across the `deep_wave` domain, against 100% at
  1.0). The slope case now carries 0.9997 of the discharge. `DE0`, `DE1` and
  `DE2` are bit-identical to 4.0. `domain.set_beta(0.5)` restores the old
  limiter but not the old predictor.

## Selected fixes

* Culvert hydraulics. Four transcription errors, each shared by the Python
  operator and its C port (which also runs serial mode-1 culverts, so CPU and
  GPU agreed while both were wrong), are fixed. Results change for the
  culverts each one touches.
  - `Boyd_pipe_operator`: the critical depth used `Q/√g·D^2.5` where
    `Q/(√g·D^2.5)` was meant, throttling pipes under about 0.8 m (a 0.375 m
    pipe passed a fifth of its flow) and forcing pipes over about 1.2 m to run
    full. It is now an exact solve, per barrel (#419, thanks to Dave
    Kennewell).
  - `Boyd_box_operator` with several barrels took the wetted perimeter of one
    wide box, overstating outlet-controlled flow by up to 20% (#420), and
    skipped the barrel-friction limit once the head difference reached the
    driving energy, so flow jumped up to 2.9× on long culverts with a free
    outlet (#421).
  - `Weir_orifice_trapezoid_operator`: barrels and blockage were applied
    twice in the C orifice area; several barrels acted as one wide trapezoid;
    the full-flow perimeter left out the roof; the friction limit had the
    same gap as the box. Its tests now reproduce the independent spreadsheet
    they were written from (#422).
* Sediment mass across MPI ranks. In mode 1 the sediment step ran before the
  ghost exchange and elevation was not exchanged, so ghost beds drifted from
  their owners' and about 0.5-1% of the eroded volume was lost at partition
  boundaries (#423, #425); bedload leaked in both modes (#424, #426). Ghost
  cells now take their owners' bed every step, and a partitioned run
  conserves sediment to round-off and matches the serial one.

* Inlet operators now carry tracers. `Inlet_operator` moved only water: an
  outlet (negative `Q`) took the water out and left its tracer behind, so
  the tracer concentrated without bound in the outlet cells as they drained,
  and an inflow always came in clean. Water removed now leaves at the inlet
  pool's mean concentration (the pool spans ranks in parallel), and water
  added carries `tracer_concentrations={'name': c or f(t)}`, zero by default
  as before. Both are exact, in both compute modes and in parallel; each
  operator keeps its running totals in `op.tracers.total_in` / `total_out`.
  Results change wherever tracers or sediment meet an extracting inlet.
* Structures carry tracers. Culverts, weirs, bridges and internal boundaries
  (`Boyd_box_operator`, `Boyd_pipe_operator`,
  `Weir_orifice_trapezoid_operator`, `Internal_boundary_operator`) moved
  water only, so the tracer stayed at the upstream end and clean water
  arrived downstream. Each cell that loses water now keeps its
  concentration, and the tracer it gives up arrives with the water, spread
  over the outflow cells in proportion to the water each gains: exact, and a
  structure that moves no water changes nothing. In both compute modes (the
  batched GPU culvert kernel included) and in parallel, where a culvert
  whose ends are on different ranks takes the host MPI path when tracers are
  present. Results change wherever a structure passes tracer or sediment.
* `Rate_operator` (and its `rainfall` / `inflow` factories) takes
  `tracer_concentrations={'name': c or f(t)}` for the water a positive rate
  adds, and `tracer_extraction='retain' | 'carry'` (or a dict per tracer)
  for what a negative rate does: `'retain'` leaves the tracer behind, the old
  behaviour and right for salt under evaporation; `'carry'` removes it with
  the water at the cell's concentration, as for a drain or pump. Exact per
  cell, in both compute modes. An operator with neither set behaves exactly
  as before and costs nothing extra.

* Structure operators (`Boyd_box_operator`, `Boyd_pipe_operator`,
  `Weir_orifice_trapezoid_operator`, `Internal_boundary_operator`) built their
  inlets from different regions depending on whether mpi4py was installed:
  the sequential class used the exchange line extended back by the apron,
  the parallel class (which the factory returns whenever mpi4py imports, on a
  serial domain too) computed the same polygon and then used the bare line
  (#367). Both now use the line plus apron, and the processors allocated to
  an inlet under MPI are chosen from that same region. Discharge and drained
  volume change for structure users whose install had mpi4py, by the apron's
  share of the exchange area.
* Line regions (`Region(line=...)`, `Set_quantity(line=...)`, the exchange
  lines of the structure operators, and the edges of an expanded polygon
  region) selected a compiler-dependent set of triangles when the line lay on
  mesh lines: the same script gave a 68% larger inlet area under one compiler
  than another (#231). The segment-triangle test is now tolerance-aware
  clipping with one rule: a triangle is selected when the segment passes
  through or along it with positive length, and not when it only touches a
  vertex or ends on an edge from outside. Selections change wherever a line
  touched triangles at a point: a 10 m exchange line on a 5 m mesh line now
  selects the 8 triangles along it rather than 32 (or 54), and an
  `Inlet_operator` line that passes through mesh vertices now feeds only the
  triangles it crosses, not every triangle meeting those vertices, so the
  same inflow enters over a smaller area and the local stage at the inlet
  is higher (the HEC-RAS bridge regression references were regenerated for
  this; the gauge agreement was unchanged). A zero-length segment, such as
  the side edges of a structure's apron polygon when `apron = 0`, selects
  nothing.
* Two runs writing the same SWW file used to interleave their frames on the
  shared time dimension, giving a file that opened fine and held plausible
  data with a non-monotonic time variable (#232). The writer now holds a
  sidecar `<name>.sww.lock` for the life of the run and a second run that
  would create the same file stops with `SWWFileInUseError`, naming the
  process and host that hold the lock (#399). A stale lock left by a process
  that is no longer running is taken over silently (logged at DEBUG, #395). It also refuses to append a
  frame earlier than the last one on file that rewrites no existing frame
  (`SWWTimeOrderError`); checkpoint resumes, which rewrite frames, still work.
* Sediment transport: the bed slope for the `'depth_slope'` shear closure
  (and the surface slope for `'energy_slope'`) is now the least-squares
  gradient of the centroid values over the cell and its neighbours. It was
  read from the edge values, which the DE algorithms rebuild every step
  through the hydrodynamic limiter: cells along a reflective wall saw no
  slope and never eroded, and the slope was reduced wherever the limiter
  engaged. Results change for those closures near boundaries, kinks and
  steps; `'quadratic_drag'` (the default) is unaffected.
* The C extension audit's correctness items (#345, #349): a gradient on a
  degenerate triangle falls back to flat instead of producing NaN; both
  conjugate-gradient solvers report breakdown instead of dividing by zero;
  the arrays handed to the C domain are checked for size; allocation
  failures in the GPU halo and inlet setup, the quad tree and the sparse
  matrices are reported instead of exiting or crashing; and GPU errors that
  were only printed to stderr now raise in Python with the message.
* `anuga_toml_run` read `[[culverts]]` and `[[weirs]]` and counted them, but
  never built them, so a TOML scenario ran without its structures; the
  scripts written by `--emit-script` had the same gap. It also skipped
  `[sediment]`. All three are now built, with the parallel operators under
  MPI (#431).
* `anuga_sww_gui` and `SWW_plotter`: generating a maximum-value frame (max
  depth, speed, ...) failed in the GUI with a Tk `TclError` about the
  toolbar, and in a notebook a plot made after a save came back closed
  (#444); animating a tracer over several frames failed in each worker process
  with a `KeyError` (#441); and an SWW file with no `timezone` attribute, as
  older versions wrote, would not open. It is now read as UTC, as the SWW
  writer assumes (#443).

## Smaller improvements

* Strict compute mode (#346): `set_compute_mode(..., strict=True)`,
  `set_gpu_offload(..., strict=True)` or `ANUGA_STRICT_COMPUTE_MODE=1` turn
  every silent fall-back to the host, for example a boundary type the GPU
  path does not support, into an error at setup, so a run that should be on
  the GPU cannot quietly run on the CPU.
* `Domain.set_plotter()` returns the plotter it attaches (#398).
* Mesh resolution from a size field (#436): `create_domain_from_regions`,
  `create_pmesh_from_regions` and `Pmesh.generate_mesh` take `size_field`,
  the maximum triangle area as a function of position, either a function
  `f(x, y)` or an `anuga.Size_field` raster. Unlike `interior_regions` it
  adds no segments to the mesh, so the outline of a resolution zone can no
  longer force the small triangles that set the timestep.
  `anuga.size_field_from_regions` builds a graded field from the same
  `(polygon, maximum area)` list. The new `relax_passes` argument, and
  `Pmesh.relax_worst_cells`, move the vertices of the smallest cells to
  enlarge them. On a 1,407 km² catchment mesh the two together lengthened
  the timestep 2.4 to 9.9 times at a similar triangle count. In tests on
  three river-flood models (#436), they cut the number of steps by 26% and
  34% on two of them, with the hydraulics essentially unchanged (on one,
  against HEC-RAS on the same terrain), and by 3% on the third, whose
  smallest cells sit against buildings meshed as holes, whose vertices stay
  fixed.
* `SWW_plotter` reads a quantity a timestep at a time, or one triangle's
  values over time for a timeseries, instead of holding the whole run in
  memory several times over, so `anuga_sww_gui` opens long runs that used to
  need many gigabytes (#442).
* The validation suite: the reports can be rebuilt from the SWW files already
  on disk (`-nr`, #403), the runner prints a summary table (#400, #401), and
  a failing case now makes the report build fail instead of being skipped
  silently. The combined three-algorithm report and the release workflow are
  described under *Highlights*.
* Documentation: the sediment physics specification is published as an
  appendix with its equations typeset (#318, #322), the `Domain` reference
  lists every public method with a guard against drift (#320), the flow
  algorithms page warns that `DE0` damps propagating waves (#405), and the
  visualisation chapter has a *Choosing a viewer* table and a page on viewing
  SWW files in ParaView and in a web browser with the separate
  `anuga-paraview` package (#445).

* The sediment kernel sees the vegetation: in a vegetated cell its bed shear
  is the bed's share of the Baptist resistance by default,
  `tau_b/rho = g u_v^2 / Cb^2` with the canopy velocity `u_v = U Cv_r / Cv`,
  or the whole vegetated resistance with
  `set_vegetation_drag(..., sediment_shear='total')`, as the Delta-X Wax Lake
  sediment model has it (`'ignore'` keeps the cell's own closure). Before,
  a vegetated class with Manning n = 0 gave the sediment no shear at all.
  Entrainment, the Rouse profile and bedload all see it.

* Closures for the adaptation of a suspension, all off by default except the
  two-layer one above. `set_deposition(adaptation=...)` takes:
  `'armanini'` ([D-3]), Galappatti and Vreugdenhil's depth-integrated lag in
  Armanini and Di Silvio's closed form, which scales the exchange rate
  without moving its equilibrium; `'carried'` ([D-4]), which carries the
  near-bed ratio as a tracer and relaxes it one-sidedly over a settling time
  when the flow decelerates; and `'two_layer'` ([D-5]), a near-bed layer and
  the rest of the column exchanging at a rate chosen so the equilibrium
  reproduces `d*` exactly, tuned by `layer_fraction` and `exchange_factor`.
  `velocity_profile=True` ([D-5v]) advects the two layers at their log-law
  mean velocities, so the transport is the profile integral of `u c` rather
  than `u h c`; a case that sets an inflow concentration to carry a known
  load must divide it by that factor.

* `set_deposition(rouse_beta='van_rijn', rouse_scale=...)` corrects the Rouse
  number of the `'rouse'` fit ([S-2b]), dividing it by van Rijn's
  `1 + 2 (w_s/u*)^2` and then by a scale. Van Rijn's measured concentration
  profiles sit at an effective Rouse number about two thirds of the plain
  one.

* `set_bedload(min_depth=h)` ramps bedload off below a depth ([K-7]), none
  below `h` and full above `2h`. At a wetting front over an erodible bed the
  film cells carry the front's momentum and a large nominal shear, and
  without the ramp they dig steps that collapse the time step. A few grain
  diameters is the physical choice.

* `set_shear_amplification(factor)` multiplies the bed shear the sediment
  sees by a per-cell factor ([T-16]), taken as a scalar, an array or a
  function of `x, y`. It exists for the structures a depth-averaged model
  cannot resolve: at a bridge pier the horseshoe vortex raises the bed shear
  to several times the approach value within about a diameter, and that is
  what digs the scour hole.

* `set_sediment_parameters(morphological_factor=M)` multiplies every step's
  bed change, from the suspended exchange and from bedload alike, leaving the
  water column untouched: the usual way to reach morphological timescales.
  `set_bedload(supply={tag: q_b})` prescribes a bedload inflow across a
  boundary in m²/s per unit width, for a flume fed at a set rate, in place of
  the zero-gradient import, which equals the inflow cell's own export and can
  run away.

* Two new experimental sediment validation cases, both from open data. Wang
  and Ribberink's (1986) settling flume, a suspension crossing onto
  perforated plate through which everything that reaches the bed is lost, so
  it tests deposition alone; and the Zaragoza pier dam-break
  (Segovia-Burillo et al. 2026, data doi:10.5281/zenodo.17777387), a
  dam-break over a sand bed around a bridge pier, bedload-dominated and
  transient, with the bed scanned after each of three events. The nine
  sediment cases, analytical and experimental, are now gathered in
  `validation_tests/sediment/` with their own chapter in the validation
  report, and each experimental case's report opens with a scale diagram of
  the flume it reproduces.

* A second non-cohesive entrainment law: de Leeuw et al. (2020),
  `set_bed_material('noncohesive', entrainment='de_leeuw')`, spec [E-6].
  `E* = A X^beta / (1 + 3 A X^beta)` with
  `X = (u*_skin / v_s)^alpha Fr - 0.015`, the skin-friction shear velocity
  from Manning-Strickler on a roughness `k_s`. Two constant sets:
  `de_leeuw_fit='de_leeuw_2020'` (their Eq 26a, sand and gravel) and
  `'nghiem_2022'`, as the Delta-X Wax Lake Delta sediment model uses it for
  mud as flocculated bed material load and for sand. Smith-McLean stays the
  default. Also in the TOML `[sediment]` table (`entrainment`,
  `de_leeuw_fit`, `skin_roughness`).

* Vegetation drag: `domain.set_vegetation_drag(density, diameter, height)`
  adds the drag of a stem field to the friction, with the vegetated Chezy
  coefficient of Baptist et al. (2007), emergent or submerged, applied
  semi-implicitly alongside Manning friction. The three fields are the
  quantities `veg_density`, `veg_diameter` and `veg_height`; cells with no
  stems are untouched. It runs in the friction kernel of both compute modes,
  so it offloads on a GPU build, which is what the Delta-X Wax Lake Delta
  model (1.5 million cells, a week of tides) needs to run in mode 2. See the
  vegetation page of the setup guide.
* `set_shear_closure()` takes two bounds for the slope closures:
  `max_slope` caps the per-cell slope (a stated replacement for anugaSed's
  undocumented global clamp), and `freeze_slope=True` (depth-slope only)
  takes the slope from the bed at setup and keeps it, so the anugaSed
  closure can run on an evolving bed without feeding on the roughness it
  creates. `Domain.bed_slope_magnitude()` returns the slope the kernel uses.
  Both are available from the TOML `[sediment]` table.
* Bedload has the Grass law, `set_bedload('grass', K=A_g)`: total load,
  `q_b = A_g |u|^(m-1) u` with `m` = 3, no threshold and no grain size, the
  law the classic Exner test cases are written in. And bedload can now pass
  through a boundary: `set_bedload(..., open_boundaries=('inflow',
  'outflow'))` makes the flux across those edges the cell's own `q_b . n`
  (zero gradient), so an outflow carries bedload away at the rate it arrives
  and an inflow supplies it at the rate the first cell removes it. Every
  boundary was closed before, so a reach's inflow cell exported bedload it
  never received and dug a hole that travelled downstream. Walls stay
  closed, and a closed domain is still exactly conservative. Both are in the
  TOML `[sediment]` table as `bedload = "grass"`, `bedload_K` and
  `bedload_open_boundaries`.
* The bedload edge flux now carries Rusanov dissipation scaled by a bound
  on the bed-wave speed, derived per formula. The bare centred flux was
  unstable on a migrating bed form: in the new bed-hump validation case it
  grew a scour hole at the upstream toe and an overshoot at the crest that
  fed back into the flow. The flux is still antisymmetric, so bedload is
  still exactly conservative, and the dissipation vanishes on a flat bed;
  results change wherever bedload moves a bed form. The bed update is
  applied in a pass of its own so that no cell reads a bed its neighbour
  has already changed.
* Four sediment validation cases join the automated suite under
  `validation_tests/analytical_exact/`: `sediment_erosion` (Smith-McLean
  entrainment on an evolving bed with the slope frozen, checked against the
  per-cell ODE), `sediment_settling_basin` (a settling basin with flow,
  whose steady concentration decays as `c0 exp(-x/L_s)` with
  `L_s = q/(d* v_s)`, checked cell by cell together with the bed-rise rate
  and the sediment budget between boundaries, water column and bed) and
  `sediment_equilibrium_flow` (clear water in normal flow on a Manning
  slope, where quadratic drag gives `tau_b = rho g h S` exactly, picking up
  its load as `c_eq (1 - exp(-x/L_s))`) and `sediment_bed_hump` (the
  Hudson-Sweby Exner test, a hump migrating under Grass bedload against its
  characteristic solution).
* A ready-made `anuga.Region` is accepted as `region=` by
  `Quantity.set_values()` (and so `Domain.set_quantity()` and friends), the
  erosion operators and `Set_w_uh_vh_operator`, alongside the existing
  polygon, circle and index forms, so one region definition serves every
  place it is needed (#15).
* The TOML scenario interface configures sediment transport: a `[sediment]`
  table for the domain-wide choices (porosity, bed evolution, shear closure,
  bed material, deposition law, friction mode, bedload, angle of repose,
  erodible base and regions) and `[[sediment.fractions]]` entries for the
  grain fractions, with inflow concentrations per boundary tag. Validated
  like the other sections, applied by `setup_sediment` in the generated run
  script, counted in the scenario summary.
* TOML scenarios take a `[project] starttime`, e.g. `"-3*3600"` for a 3 h
  burn-in before the forcing starts, on the same clock as `finaltime` (#433).
  Inlet, stage-boundary and CSV-rainfall series hold their first value before
  they begin (rain is 0), instead of raising.
* An inlet whose enquiry point falls in one of its own triangles now warns
  once per run with `anuga.EnquiryPointWarning`, instead of once per inlet;
  `verbose=True` on the structure still names each point (#417).
* GeoPackage in and out (#39, #40): `gpkg2polygons()` / `polygons2gpkg()`
  read and write polygon and polyline layers with attributes and a CRS as
  the `[x, y]` lists ANUGA's regions, buildings and breaklines take, and
  `polygon_csv_files2gpkg()` / `gpkg2polygon_csv_files()` and
  `building_csv2gpkg()` / `gpkg2building_csv()` convert the two CSV
  conventions ANUGA already reads. Needs `fiona` and `shapely`
  (`anuga[data]`), imported only on call.
* `Domain.get_water_volume()` and `compute_total_volume()` take `region=` (an
  `anuga.Region`) or `indices=` to sum the water in part of the domain, as
  `Quantity.get_integral()` already did (#22). Regional queries do not enter
  `volume_history`.
* `Wind_stress_operator` and `Barometric_pressure_operator` accept
  `use_coordinates=False` with a `file_function`, replacing the file-driven
  wind and pressure fields of the removed `_fast` classes. The file's
  precomputed time series is interpolated for every point at once
  (`anuga.utilities.function_utils.evaluate_file_function_all_points`).

## Requirements

* Python 3.10 – 3.15 (3.15 added, #397), numpy ≥ 2.0.
* GeoPackage support needs `fiona` and `shapely` (`pip install anuga[data]`),
  imported only when used.
* GPU offload requires the NVIDIA HPC SDK (`nvc`), as in 4.0; the standard
  build needs no GPU.

## Known issues

* The `'depth_slope'` shear closure takes its slope from the bed, so on an
  evolving bed it can amplify itself: a local scour steepens the slope,
  which raises the shear and scours faster. Bound it with `max_slope` or
  `freeze_slope=True`, or use the default `'quadratic_drag'` closure.
* Bedload is first order in the bed.
* `DE_ader2`'s steady state depends slightly on the time step, so in a
  physically unstable flow (roll waves on a steep slope) it amplifies its own
  time-step changes where `DE1` and `DE2` stay exactly steady. The effect is
  about 10⁻³ of the depth or less. See `claude/KNOWN_ISSUES.md`.

## Thanks

Thanks to Dave Kennewell, for the fix to the Boyd pipe culvert's critical
depth (#419) and for taking up the call for testing (#411); and to
@PMacAddress, who tested the size field and worst-cell relaxation on three
river-flood models, one against HEC-RAS (#436), and reported on viewing
large results remotely (anuga-community/anuga-viewer#2), which led to the
`SWW_plotter` memory fix (#442).

---

# ANUGA 4.0.0

First release on the 4.x line, and the first release of the work that has been
accumulating on `develop` since 3.3.10 — 993 commits across 935 files.

The headline is a **GPU/multicore solver**: the shallow-water kernels now run
through a single unified C implementation that executes on the CPU or offloads
to a GPU, with mode-1 and mode-2 agreeing bit-for-bit on the paths that matter.
Alongside that are a decade-old correctness fix in the structure operators, a
substantial set of parallel and wet/dry fixes, and container images.

**Existing scripts keep working.** The default compute mode is unchanged
(`'legacy'`), and both the unified mode and GPU offload are opt-in. What makes this a major version is the
removal of long-deprecated code, not a change in day-to-day behaviour.

---

## Highlights

### The unified compute mode

The solver and its operators now have a single C implementation that runs
CPU-multicore, and can offload to a GPU on a suitably built install.

```python
domain.set_compute_mode('unified')     # per domain; CPU-multicore by default
```

or process-wide with `ANUGA_DEFAULT_COMPUTE_MODE=unified`.

* **No GPU required.** 'unified' is a compute mode, not a GPU switch: on an
  ordinary build it runs the same unified C kernels on the CPU. GPU offload is a
  separate, process-wide opt-in — `anuga.set_gpu_offload(True)` on a build made
  with the NVIDIA HPC SDK. `gpu_offload_supported()` reports whether a build can.
* **One implementation of the physics.** `culvert_compute_one()` is shared by
  both paths, so 'legacy' and 'unified' give bit-identical culvert results.
* On GPU: fluxes, extrapolation, boundaries, riverwalls, culverts and operators
  all execute on device, with multi-GPU via MPI and device-side halo exchange.
* New controls: `set_gpu_offload()`, `gpu_offload_enabled()`,
  `gpu_offload_supported()`, `set_omp_num_threads()`, `get_omp_num_threads()`.

`set_multiprocessor_mode(1|2)` still works as a thin wrapper, but
`set_compute_mode('legacy'|'unified')` is the preferred API for new code.

### New timestepping method

* **`DE_ader2`** — a fused ADER-2 predictor/extrapolation step, measured **1.75×
  faster than DE1** at equivalent accuracy. Select with
  `domain.set_flow_algorithm('DE_ader2')`.

### Rainfall and forcing

* `Raster_rate_operator` — apply gridded rainfall rasters directly.
* Australian Rainfall & Runoff support: `ARR_rate_operator`, `Arr_hub_rain`,
  `Arr_ifd_rain`, `ARR_point_rainfall_patterns`, `Arr_grd`.
* `Wind_stress_operator` and `Barometric_pressure_operator` — operator
  equivalents of the legacy forcing classes (see Deprecations).

### Boundaries

* `Absorbing_wave_boundary` and `Characteristic_wave_boundary`.

### Meshes and memory

* `sequential_mesh_dump` / `sequential_mesh_load` — partition once, run many
  scenarios; NetCDF, self-describing, safe to share.
* `uniform_refine_domain`, `sequential_mesh_refine`, `create_parallel_mesh`.
* ~58% reduction in quantity memory; `quantity_memory_stats()`,
  `domain_memory_stats()` to inspect it.

### Tooling

* Container images: CPU, GPU (NVHPC/nvc), and GPU+CUDA-aware-MPI, published to
  GHCR on release.
* GUIs: `anuga_sww_gui`, `anuga_animate_gui`.
* `anuga_toml_run` for the TOML scenario interface (which itself shipped in
  3.3.x; new here are `[[mesh.interior_holes]]` and `[[erosion]]` sections).
* `anuga_run_isolated_tests` — per-test process isolation, needed on GPU builds.

---

## Breaking changes

* **`anuga.culvert_flows` is removed** (7 790 lines). It was superseded years ago.
  Migrate to the structure operators — `Boyd_box_operator`, `Boyd_pipe_operator`,
  `Weir_orifice_trapezoid_operator`. `examples/structures/run_open_slot_wide_bridge.py`
  shows the equivalent setup.
* **Local-timestepping attributes removed**: `max_flux_update_frequency`,
  `flux_update_frequency`, `update_next_flux`, `update_extrapolation`,
  `edge_timestep`, `allow_timestep_increase`. The 3.1.9 implementation was not
  GPU-compatible and had been dead code; see `claude/FUTURE_WORK.md` P3.1 for the
  redesign.
* **Structure operators now write back a level water surface** (see below). This
  changes results for any structure whose inlet sits on a sloping bed.

### The structure write-back change (#229)

Structure operators wrote their transfer back as a **uniform depth** across each
inlet. On a sloping bed that tilts a level water surface onto the bed, so a lake
at rest was disturbed by roughly half the bed elevation range across the inlet —
every timestep, at zero discharge.

The write-back now applies the volume change by **levelling the stage** (water
finds its level), clamping cells at their bed.

What to expect:

| case | change |
|---|---|
| flat bed under the inlet | **none** — bit-identical by construction |
| sloping bed | results move; magnitude scales with the bed range across the inlet |

Measured on the Towradgi catchment (22 culverts on real terrain, inlet bed
spread median 0.77 m): peak stage difference **0.32 m**, 99th percentile 1.5 mm.
Three of four flat-bed validation cases were bit-identical.

If you calibrated a model against pre-4.0.0 results with structures on sloping
ground, re-check it.

---

## Deprecations

The legacy forcing classes in `anuga.shallow_water.forcing` now all emit
`DeprecationWarning` and **will be removed in 4.1**:

| deprecated | use instead |
|---|---|
| `Rainfall` | `Rate_operator.rainfall()` |
| `Inflow` | `Rate_operator.inflow()` |
| `Wind_stress`, `Wind_stress_fast` | `Wind_stress_operator` |
| `Barometric_pressure`, `Barometric_pressure_fast` | `Barometric_pressure_operator` |

These are **silently skipped** in the `'unified'` compute mode, which applies
forcing in C and handles only Manning friction; a warning is issued when that
happens. Migrating to the operators is required to use them with 'unified'.

`manning_friction_semi_implicit` is *not* deprecated — it is an in-step
semi-implicit term, not an operator.

---

## Selected fixes

* **Parallel inlet mass balance** (#193) and **startup mass loss** (#200) — both
  affected conserved volume in ordinary runs.
* **Culvert stack overflow** with more than 64 culverts (#217) — buffers were
  sized by a constant that was only the initial capacity.
* **GPU culvert inlet cap** (#225) — an inlet was limited to 64 triangles.
* **Riverwall crest changes reach the device** (#224) — runtime
  `set_elevation()` was silently ignored under GPU offload.
* **DE0 non-GPU boundary fallback** — Euler steps silently ignored
  Python-evaluated boundary types in the 'unified' mode.
* Degenerate-timestep protection now warns rather than silently not running
  in the 'unified' mode (#189).

---

## Smaller improvements

* Importing anuga from an unbuilt source tree now explains itself. It used to
  fail with a bare `ModuleNotFoundError: No module named 'anuga._version'`,
  which pointed at neither cause: `_version.py` is generated at build time, and
  because Python searches the working directory first, an unbuilt tree also
  shadows a correctly installed anuga whenever you run from the repository
  root. The error now names the directory it imported from and says how to fix
  it (#237).

## Requirements

* Python 3.10 – 3.14, numpy ≥ 2.0
* GPU offload requires the NVIDIA HPC SDK (`nvc`); the standard build is
  unaffected and needs no GPU.

## Known issues

* Windows CI pins the conda-forge mingw sysroot chain to build 10. Build 11
  (2026-08-20) produces a toolchain whose binaries abort in the stack protector;
  tracked at conda-forge/m2w64-sysroot-feedstock#23. This pins only the CI
  environment — released wheels are unaffected.
* A structure sitting in near-still water amplifies roundoff; see
  `claude/KNOWN_ISSUES.md`.

## Thanks

Stephen Roberts, Jorge Luis Gálvez Vallejo, wangshuo, Samir Shaikh, and everyone
who reported issues against 3.3.x.

The GPU and multicore solver was developed in collaboration with the Centre for
Development of Advanced Computing (C-DAC), India, as part of a catchment flood
prediction project aimed at two-day forecasts, and the National Computational
Infrastructure (NCI), Australia, whose Gadi supercomputer provided the GPU
systems it was developed, tested and benchmarked on. Thanks to Samir Shaikh,
Shweta Das and Rutvik Gulhane (C-DAC) and Jorge Luis Gálvez Vallejo (NCI).
