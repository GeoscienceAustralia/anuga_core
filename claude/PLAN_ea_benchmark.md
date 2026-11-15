# Plan: add the UK Environment Agency 2D benchmark suite to validation tests

2026-11-01. Issue: #476. Living version: https://claude.ai/artifact/LXoAMu62G4RBMn71ozerY8

## The benchmark suite

The Environment Agency's SC120002 report (Néelz and Pender, completed August 2013) ran 19 packages through ten tests. ANUGA 1.1beta took part in 2009 and ran Tests 1–6B. Most tests have no analytical answer: a package is judged against the spread of the others at fixed output points. The inputs are a public 51 MB ZIP under the Open Government Licence v3.0.

| Test | What it checks | Domain | Specified grid | Run to | Manning n | Reference | ANUGA in 2009 | EA recommendation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 Flooding a disconnected water body | wetting over a bank into a basin, final level | 700 × 100 m | 10 m | 20 h | 0.03 | other packages | ran | remove (trivial) |
| 2 Filling of floodplain depressions | momentum vs simplified models; levels in 16 depressions | 2 × 2 km | 20 m | 48 h | 0.03 | other packages | ran, 1,130 s | keep |
| 3 Momentum over a small obstruction | momentum conservation down a 1:200 slope | 300 × 100 m | 5 m | 15 min | 0.01 | other packages | ran | keep |
| 4 Flood propagation over an extended floodplain | front speed on a flat plain, 20 m³/s breach inflow | 1 × 2 km | 5 m | 5 h | 0.05 | other packages | ran, 3,650 s | keep; boundary spec clarified |
| 5 Valley flooding | dam-failure flood down a real valley: peak level, speed, arrival | 0.8 × 17 km | 50 m | 30 h | 0.04 | other packages | ran, 4,160 s; arrival ~1.5 h late and level ~0.5 m low at point 5; spurious shallow-flow velocities | keep |
| 6A Dam break, laboratory | flow round a building (Soares-Frazão and Zech 2002 flume) | flume | 0.1 m | 2 min | 0.01 | **measured** | ran | remove |
| 6B Dam break, field scale | 6A × 20: shock capturing, high-velocity flow | ×20 flume | 2 m | 30 min | 0.05 | other packages | ran, 1,390 s; high peak velocities | keep; geometry as GIS polygons |
| 7 River and floodplain linking | 1D river linked to 2D floodplain | river reach | 20 m | 72 h | 0.028 river | other packages | not run (no 1D–2D link) | keep; clearer set-up |
| 8A Rainfall and point source, urban | rain on grid and a point inflow over streets and buildings | 0.388 km² (Glasgow) | 2 m | 5 h | 0.02 roads | other packages | not run | merge with 8B |
| 8B Surcharging sewer, urban | 1D sewer surcharging onto the 2D surface | as 8A | 2 m | 5 h | 0.02 roads | other packages | not run | merge with 8A |

The report flagged ANUGA's peak-velocity maps for investigation (Tests 5 and 6B). It could not say whether the cause was the scheme or the gridding of its output. Today's `DE_ader2` and the gridding tools differ from 2009's code, so rerunning these tests would show what has changed.

## What validation_tests/ already covers

None of the ten tests is in `validation_tests/` today, and no case cites Néelz and Pender or Soares-Frazão and Zech. The nearest existing cases cover the same physics in idealised form:

- **Dam breaks (6A, 6B):** `analytical_exact/dam_break_dry`, `dam_break_wet` and `other_references/radial_dam_break_*` are 1D or radial with no obstacle. None has a building in the flow.
- **Real-terrain floods (5, 8A):** `case_studies/merewether` (urban, houses, observed levels) and `towradgi` are the closest. Neither has a reference from other packages.
- **Momentum over obstacles (3):** `subcritical_over_bump` and `transcritical_*` test a bump in a channel, against an analytical answer.
- **Package comparison:** the `behaviour_only/*_hecras` cases already compare ANUGA with HEC-RAS results at fixed points. That is the pattern the EA tests need.

The EA ZIP holds the inputs only: a DEM (`.asc`), boundary tables (`.csv`), output-point lists and shapefiles for Tests 1–7, plus cross-section workbooks for Tests 5 and 6B. **It holds no Test 8A/8B data and none of the packages' results.** A pass envelope therefore has to come from somewhere else (see Open questions).

## Plan

Add the tests as a new `validation_tests/ea_benchmark/` group, one directory per test, run by `run_auto_validation_tests.py` like the others. Each case is checked two ways:

- **Against the other packages:** at each output point, ANUGA's level (and speed, where the report gives it) must lie inside the envelope of the full shallow-water packages, plus a stated tolerance. Test 6A is checked against the measured gauges instead.
- **Against itself:** a stored ANUGA baseline, as in `behaviour_only/`, so that a later change that moves the answer is caught even while it stays inside the envelope.

Order of work:

1. **Shared pieces.** A fetch step (done: `validation_tests/ea_benchmark/data_download.py`, unpacks to `data/TestN/`) for the gov.uk ZIP with a checksum, a reader for the DEM, boundary and output-point files, and one plot of ANUGA against the envelope at each point.
2. **Tests 3, 4 and 6B.** Small, quick and cleanly specified. 6B also rechecks the 2009 peak-velocity issue.
3. **Tests 2 and 5.** Larger runs (48 h and 30 h). Test 5 rechecks the 2009 late arrival at point 5.
4. **Test 6A.** Needs the flume gauge records from Soares-Frazão and Zech (2002).
5. **Test 1.** Cheap and near-trivial, so it goes last; it completes the set.
6. **Later, if wanted: Tests 7 and 8A.** Test 7 can be run fully 2D, with the river in the mesh and the bridge and culvert as structure operators. Test 8A needs its data from the EA, since the ZIP lacks it.

All runs use `DE_ader2` at the report's specified grid size, then once more on a finer mesh to show convergence.

## Open questions

- [ ] **Where does the envelope come from?** The ZIP has no package results. Options: digitise the report's time-series figures, ask the EA or Heriot-Watt for the submitted results, or use published HEC-RAS 2D and SRH-2D runs of the same tests.
- [ ] **Commit the inputs or download them?** The OGL allows redistribution with attribution. Tests 1, 3, 4 and 6 are under 100 kB each and could be committed. Test 2 (3.3 MB), Test 5 (1.1 MB) and Test 7 (48 MB) are better fetched.
- [ ] **Test 6A gauges:** are the Soares-Frazão and Zech (2002) measurements available as data, or only as figures?
- [ ] **Test 8A/8B data:** is it still available from the EA? It is not in the public ZIP.
- [ ] **Nightly runtime:** in 2009, single-core ANUGA took 1,130 s (Test 2) to 4,160 s (Test 5). Do Tests 2 and 5 fit the nightly budget, or do they belong in the `long` set?
- [ ] **Tolerance:** how far outside the envelope counts as a failure, given that the specified grids are coarse (50 m in Test 5)?

## Sources

- [2D benchmarking: evaluating the latest generation of the hydraulic models for FCRM (SC120002)](https://gov.uk/flood-and-coastal-erosion-risk-management-research-reports/2d-benchmarking-evaluating-the-latest-generation-of-the-hydraulic-models-for-fcrm): the project page, with the OGL v3.0 licence.
- [Néelz and Pender (2013), Benchmarking the latest generation of 2D hydraulic modelling packages, technical report](https://assets.publishing.service.gov.uk/media/6033a943d3bf7f721f4b0d49/_SC120002_Benchmarking_2D_hydraulic_models_Report.pdf): the test specifications (Appendix A), results, ANUGA's 2009 run times (Table 4.13) and recommendations (Section 5).
- [Benchmarking model data (ZIP, 51 MB)](https://assets.publishing.service.gov.uk/media/6033a9bbd3bf7f72182e99ca/Benchmarking_Model_Data.zip): the inputs for Tests 1–7, downloaded and listed on 1 November 2026.
- [Validation of AquaTwin Sewer 2D against the EA benchmark suite](https://aquanuity.com/resources/validation-of-aquatwin-sewer-2d-against-the-uk-environment-agency-2d-hydraulic-model-benchmark-suite/): an example of a package judged against the min–max envelope of 21 reference models.
