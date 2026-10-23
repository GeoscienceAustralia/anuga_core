# ANUGA Code & Documentation Improvement Progress

Last updated: 2026-10-02 (session 58)
Branch: `develop` (all feature branches merged)

Full history of completed work: `claude/PROGRESS_ARCHIVE.md`
Future work recommendations: `claude/FUTURE_WORK.md`

---

## Overview

| Area | Total | Done | Remaining |
|------|-------|------|-----------|
| Code improvements (original list) | 60 | 53 | 7 |
| Documentation improvements | 20 | 20 | 0 |
| Additional enhancements | 57 | 57 | 0 |
| Hydrata Phases 0–4 | 21 | 21 | 0 |
| GPU Phase 1–3 | 15 | 15 | 0 |
| GPU Phase 4 — SC26 paper | 3 | 0 | 3 |
| Riverwall throughflow | 6 | 6 | 0 |
| Quantity memory reduction | 7 | 7 | 0 |
| Domain memory reduction | 3 | 3 | 0 |
| Benchmark suite | 2 | 2 | 0 |
| Bug fixes | 7 | 7 | 0 |
| Kinematic viscosity parallelisation | 4 | 4 | 0 |
| Mode-2 test triage + isolated runner tooling | 4 | 4 | 0 |
| Documentation overhaul (session 47) | 15 | 15 | 0 |
| GPU mode-2 Time_boundary fix | 2 | 2 | 0 |
| Operator-timing fixes (DE1 rk2, DE2 rk3) | 2 | 2 | 0 |
| Validation reports per release (session 58) | 6 | 5 | 1 |
| **Total** | **234** | **225** | **9** |

---

## Remaining / Deferred Items

### Code improvements

- [ ] **3.3** `anuga/scenario/` tests — deferred; depends on compiled `spatialInputUtil` and real test data
- [ ] **4.1** Reduce parameter counts (`gauge.py`, `generic_domain.py`, `boyd_box_operator.py`) — deferred
- [ ] **5.2** `polygon.intersection()` — not a confirmed hotspot; deferred
- [ ] **5.x** `util.py:301` (`csv2timeseries_graphs`) — dominated by matplotlib; deferred
- [ ] **1.3** `anuga/file/urs.py:29` — intentionally skipped: file handle stored as `self.mux_file` for iterator lifecycle

### GPU Phase 4 — SC26 paper preparation (months 4–6)

- [ ] **G4.1** Gordon Bell metrics — per-kernel timing (not just totals), roofline model comparison, peak theoretical FLOP/s
- [ ] **G4.2** Physical benchmark validation — Thacker paraboloid, dam break (Ritter), tide gauge comparison in GPU mode
- [ ] **G4.3** Multi-node strong scaling — 20 M triangles, 1→64 GPUs; demonstrate ~50× runtime reduction

### Validation reports per release (session 58)

- [x] Roll-wave + stable cases in `rundown_mild_slope`, solver comparison (PR #407)
- [x] DE_ader2 dt-dependent steady state documented (KNOWN_ISSUES)
- [x] Release workflow: per-algorithm reports, release assets, Zenodo draft (PR #408)
- [x] Zenodo sandbox dry run (create + new version), tokens, `anuga` community
- [x] Combined three-algorithm report (PR #409); CITATION.cff valid and
  checked in CI (#410, #412); 4.1.0 release notes (#413)
- [ ] After the next develop→main merge, run the workflow by hand with
  zenodo = sandbox to test its wiring; add the Zenodo community co-managers;
  after the first release is archived, add the software record to the
  `anuga` community.

### GPU correctness

- [x] **Time_boundary substep** (issue #170): option B, routing domains with
  Python-evaluated boundaries to the Python-orchestrated GPU loop (PR #171), is
  the kept design — decided 2026-09-19, issue closed. Option A (per-substep
  values passed into the single-call C RK loop) was judged not worth its cost
  for a measured ≤4% routing overhead. Reopen only if a tidal/file-boundary GPU
  workload is dominated by it.

---

## Remaining Work (priority order)

### Short term — SC26 prerequisites (needs GPU hardware)
1. **G2.1** GPU benchmark suite (actual GPU runs — `benchmarks/run_gpu_benchmarks.py` is ready)
2. **G2.4** Weak scaling experiment (1→64 GPUs on HPC cluster — scripts in `scripts/hpc/`)
3. **G4.1** Gordon Bell metrics — per-kernel timing, roofline model comparison
4. **G4.2** Physical benchmark validation — Thacker paraboloid, dam break (Ritter) in GPU mode
5. **G4.3** Multi-node strong scaling — 20 M triangles, 1→64 GPUs

### Medium effort — see `claude/FUTURE_WORK.md`

Top P2 recommendations: P2.6 (fast-suite coverage), P2.7 (gauge module modernisation), P2.4 (culvert compute_rates deduplication).

### Long-term / opportunistic
- **H4.1** Modernise test patterns — convert `unittest.TestCase` to plain pytest functions incrementally when files are touched. Not worth a dedicated pass.
- **Coverage** — Full suite ~58%, fast suite ~55%. `fail_under=57` in `.coveragerc`. Add tests opportunistically when touching files.
- **Local-timestepping** — `compute_flux_update_frequency` is a `pass` stub. When implementing, allocate work arrays on demand in `set_local_time_stepping()`. See P3.1 in `FUTURE_WORK.md` for GPU-compatible redesign notes.
