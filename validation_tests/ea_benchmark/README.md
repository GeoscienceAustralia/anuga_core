# Environment Agency 2D benchmark (SC120002)

Cases from the UK Environment Agency's benchmark of 2D flood models: Néelz, S.
and Pender, G. (2013), *Benchmarking the latest generation of 2D hydraulic
modelling packages*, Report SC120002, Environment Agency, Bristol.
Plan and status: issue #476 and `claude/PLAN_ea_benchmark.md`.

The input data (DEMs, boundary tables, output points) are published under the
Open Government Licence v3.0. `data_download.py` fetches the 51 MB ZIP from
gov.uk, checks its SHA-256 and unpacks it into `data/TestN/`. Each case calls
it on first use, so there is no separate step.

| Directory | Test | Run time (8 threads) | Checked against |
| --- | --- | --- | --- |
| `test1_disconnected` | 1: flooding a disconnected water body | 10 s | peak and final levels, flat pond, sheet flow (Section 4.2); baseline |
| `test2_depressions` | 2: filling floodplain depressions | 50 s | volume; sill levels (DEM); final levels (TUFLOW 2017); transient peaks, wet/dry depressions (Section 4.3.3); baseline |
| `test3_obstruction` | 3: momentum over a small obstruction | 2 s | volume; report text (Section 4.4.3); baseline |
| `test4_floodplain` | 4: flood propagation over a flat floodplain | 3 min | volume; front order and ~1 h travel time (Section 4.5.4); baseline |
| `test5_valley` | 5: dam-failure flood down a 17 km valley | 40 s | volume; arrival order; arrival at point 5, peak speeds (Section 4.6.4); final pond level (TUFLOW 2017); baseline |
| `test6a_dambreak` | 6A: dam break past a building, laboratory flume | 30 s (report: 70 s, with the limiter study) | volume; bore order; gate-flow period (Section 4.7.1); speed bound; baseline |
| `test6b_dambreak` | 6B: dam break past a building, field scale | 1 min | volume; bore order; speed bound; level end state; baseline |

The report publishes the other packages' results only as raster figures (and
reserves its rights in them), so the cases do not test against their envelope.
Each case checks what the report states in its text, physical constraints, and
a stored baseline per flow algorithm (`baseline_<test>_<alg>.csv`), so that a
change that moves the answer fails. Tests 7 and 8A/8B are not included yet
(issues #477 and #478).

Each case directory holds:

| File | Role |
| --- | --- |
| `run_testN.py` | the simulation; writes `testN_gauges.csv`, `testN_volume.txt` and `parameters.tex` |
| `checks.py` | the checks, as a list of (name, value, criterion, pass) |
| `validate_testN.py` | runs the case and asserts the checks (nightly and `run_auto_validation_tests.py`) |
| `produce_results.py`, `plot_results.py` | the validation report: set-up, time-series and checks figures; `plot_results.py` exits non-zero if a check fails |
| `results.tex`, `report.tex` | the case's section of the full report, and a stand-alone report |

The cases sit in their own chapter of the full validation report
(`reports/validations_report_body.tex`). To build one case's report:

```bash
cd test6b_dambreak
python produce_results.py -alg DE_ader2      # report.pdf
```

To write or rewrite the baseline for one algorithm:

```bash
cd test4_floodplain
EA_REGENERATE_BASELINE=1 python validate_test4.py -alg DE1
```

The 6A report also runs `run_limiter_study.py`: the first 40 s with the slope
limiter's beta set to 0, 0.5 and 1, for its figures of how the limiter moves
the hydraulic jump (#479). The nightly validation does not run it.

Tests 6A and 6B share their set-up (`ea_dambreak.py`). They take two
settings for studies, both off by default: `EA_6A_MAX_AREA`/`EA_6B_MAX_AREA`
(triangle area; finer meshes are best run in `unified` mode on a GPU) and
`EA_6A_BUILDING`/`EA_6B_BUILDING` (`hole`, the default; `block`, the fitted
rectangle as raised floor; or `dem`, the DEM's own building cells). Tests 2 and 5
take `EA_2_MAX_AREA` and `EA_5_MAX_AREA`. The outputs are then named `test<N>_a<area>_*` and
`test<N>_<building>_*`.

The cases run in serial: gauges are sampled from the domain during
`evolve`, which needs the whole mesh in one process.
