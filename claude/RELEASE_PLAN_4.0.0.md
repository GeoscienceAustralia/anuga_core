# ANUGA 4.0.0 — Release Plan

Drafted 2026-08-21 (session 57), expanding R1 of the engineering review.
Owner: Stephen Roberts. Target: **tag within 2 weeks of the go decision.**

---

## The decision, and why 4.0.0 directly (no 3.4 bridge)

`main` is at **3.3.10**; `develop` is **976 commits ahead** (933 files, +400k/−173k
lines). The review's R1 offered two routes; the facts point at one:

- **Both branches already share the same foundation** — meson-python, numpy ≥ 2,
  Python 3.10–3.14, the same wheel workflow. There is no packaging cliff to bridge.
- **develop still defaults to `legacy` compute mode** (mode 1;
  `shallow_water_domain.py:643` — `unified` is opt-in via
  `ANUGA_DEFAULT_COMPUTE_MODE`). Existing user scripts run unchanged on develop.
  The scary flag-day (mode-2 by default) is NOT part of this release.
- What makes it a major version is **removal and architecture**, not behaviour:
  `culvert_flows/` deleted, local-timestepping infrastructure deleted, forcing
  classes deprecated (removal in 4.1), and the GPU/mode-2 solver as the headline.
- A 3.4 bridge would mean curating ~976 commits into a partial cherry-pick — weeks
  of work to ship *less*, with a second release to follow anyway.

**Precondition check (from ROADMAP.md):** "SC26 paper submitted / merge sp26" — done;
sp26 has been merged into develop since 2026-04-01. The stated gate is met.

---

## Scope

**In (everything on develop at the freeze commit), headlines for the notes:**

| Theme | Items |
|---|---|
| New solver capability | GPU/mode-2 (`multiprocessor_mode=2`, OpenMP target offload); `DE_ader2` algorithm (1.75× DE1); unified CPU/GPU culvert kernel (mode-1 == mode-2 bit-for-bit) |
| Correctness fixes | parallel-inlet mass balance (#193), startup mass loss (#200), culvert stack overflow (#217), riverwall crest→device (#224), inlet cap (#225), structure well-balancedness (#229), DE0 non-GPU-boundary fallback |
| New interfaces | TOML scenario system (`anuga_run_toml`, incl. interior holes + erosion, #214); `RiverWall.set_elevation()` runtime API; `Rate_operator.rainfall()/inflow()` factories; type hints on the public API |
| Performance/memory | ~58% quantity-memory reduction; OpenMP tuning; benchmark suite |
| Removals (BREAKING) | `anuga/culvert_flows/` (migration: `Boyd_box_operator` — example in `run_open_slot_wide_bridge.py`); local-timestepping dead attrs (`flux_update_frequency` etc.) |
| Deprecations (warn, remove in 4.1) | `Inflow`, `Rainfall`, `Wind_stress`, `Barometric_pressure` forcing classes → operators |
| Tooling | Docker images (CPU/GPU/GPU-MPI) + GHCR publishing; hardened installers; `anuga_run_isolated_tests`; SWW GUI improvements |

**Explicitly out (deferred to 4.1+):**
- Flipping the default compute mode to `unified` (PLAN_default_mode2_cpu.md)
- Removing the deprecated forcing classes (P2.10 — one release of warning first)
- cibuildwheel migration (#141), OpenACC backend (#188)

---

## Phase 0 — Freeze & decision (Day 0–1)

- [ ] Team says **go** (this plan is the proposal).
- [ ] Freeze commit chosen on `develop`; only release-blocking fixes land after.
- [ ] Merge-window discipline: hold open PRs (#148 etc.) until after the tag.
- [ ] Decide the branch-protection question **before** the release PR: either
      enforce review on `main`/`develop` or remove the rule. The release itself
      should not need `--admin`.

## Phase 1 — Pre-flight verification (Day 1–5)

Each gate is a named command with a recorded result; a red gate blocks the tag.

- [ ] **Full suite, CPU build** (as CI runs it):
      `cd sandpit && OMP_NUM_THREADS=1 pytest -rs --pyargs anuga` — expect ~2 937 pass.
- [ ] **Unified-mode suite, CPU build** (one process, documented all-green config):
      `ANUGA_DEFAULT_COMPUTE_MODE=unified pytest --pyargs anuga --run-fast`.
- [ ] **GPU build, isolated runner** (local RTX 5070 + one AWS g6):
      `bash anuga/shallow_water/tests/run_gpu_tests_isolated.sh` and
      `anuga_run_isolated_tests --pyargs anuga.shallow_water -cm unified`.
- [ ] **Validation suite**: `python validation_tests/run_auto_validation_tests.py`.
      **Expected delta:** #229 moves results for structures on sloping beds —
      record before/after for the affected cases and put the numbers in the
      release notes rather than being surprised by a user report.
- [ ] **Towradgi case study**, mode 1 vs mode 2, archived to S3 with the build report.
- [ ] **MPI smoke**: `pytest anuga/parallel/tests/` (auto-slow) on 2 and 4 ranks.
- [ ] **Wheel smoke**: install the CI-built wheel in a fresh venv on
      linux/macos/windows (the 37-check matrix on PR #230 is recent evidence this
      is green); `import anuga; anuga.test()` + one example script.
- [ ] **Fresh-environment install**: `environments/environment_3.10.yml` and 3.14,
      `pip install --no-build-isolation -e .`, run an example.
- [ ] **Docs build** clean; add `docs/source/reference/generated/` to `.gitignore`.

## RESOLVED: Windows CI (2026-08-21/22) — Phase 2 hold LIFTED

Fixed in PR #233, merged as `6d5e7c4b`. All 37 checks green (1 skip: Publish to
PyPI, releases only). The hold below is kept as the record of what happened.

Every Windows job in **both** workflows (`conda-setup.yml` and
`python-publish-pypi.yml` — they share an install line) has failed since
2026-08-20 at meson's compiler sanity check, before any ANUGA code compiles:

```
meson.build:9:0: ERROR: Executables created by c compiler
  .../x86_64-w64-mingw32-cc are not runnable.
```

**Established by measurement:**

* Not our code — commits touching only markdown fail identically.
* Last green 08-20T13:04; first red 08-21T00:50; red on every run since.
* `gcc_win-64` 16.1.0, `binutils_win-64` 2.46.1 and the runner image
  (`windows-2025-vs2026`, 20260729.566) are **identical** either side of the break.
* A full package diff shows **exactly four** differences, same version,
  build 10 -> 11: `m2w64-sysroot_win-64`,
  `mingw-w64-ucrt-x86_64-{crt,headers,winpthreads}-git`.
* Those come from conda-forge `m2w64-sysroot-feedstock` PR #21 ("finish v1
  transition"), merged 2026-08-20 — the same day.

**Falsified:** that the newly split-out `libwinpthread` (the package holding
`libwinpthread-1.dll`) was missing from the env and caused it. Installing it
explicitly changed nothing; the error was byte-identical. Do not re-try this.

**Actual cause**, found by running meson's sanity check by hand:

```
printf 'int main(void){return 0;}' > sanity.c
x86_64-w64-mingw32-cc sanity.c -o sanity.exe   ->  compile: OK
./sanity.exe   ->  *** stack smashing detected ***: terminated  (exit 127)
```

`__stack_chk_fail` fires on a function with no locals and no buffers: the
stack-protector ABI between `gcc_impl_win-64` 16.1.0 and the build-11 CRT does
not line up, so *any* mingw-built binary in that env aborts on startup.

**Fix shipped:** pin those four packages to build 10 in both workflows, marked
TEMPORARY with a link to conda-forge/m2w64-sysroot-feedstock#23.

**Tracking:** PR #233 (draft, ours) · conda-forge/m2w64-sysroot-feedstock#23
(upstream). Fallback if the cause stays unclear: pin the sysroot chain to
build 10 (`m2w64-sysroot_win-64=*=*_10`) as a temporary hold.

**Effect on Phase 1:** gate 7 is satisfied again — all 5 Windows wheel builds
and all 5 Example jobs pass. Phase 1 is complete, nine gates green.

**Carry into the release notes:** while the pin stands, Windows CI holds the
mingw sysroot at build 10. It pins only the CI environment, not anything a user
receives in a wheel. Remove it when upstream fixes build 11.

---

## Phase 2 — Release candidate (Day 5–9)

- [ ] **Release notes** (`docs/` + GitHub Release body). Method: walk
      `git log 3.3.10..HEAD` segmented by the session summaries in
      `claude/SESSION_GUIDE.md` §21–56 — they are already a categorised changelog.
      Structure: Highlights / Breaking changes / Deprecations / Fixes / Thanks.
- [ ] **Upgrade guide** (one page): `culvert_flows` → `Boyd_box_operator` mapping;
      forcing classes → operators table; note that mode 2 is opt-in and how to try it
      (`domain.set_multiprocessor_mode(2)` / `ANUGA_DEFAULT_COMPUTE_MODE=unified`).
- [x] **Version/citation hygiene**: CITATION.cff bumped to 4.0.0; the stale
      2022 `commit:` pin dropped and `date-released` omitted rather than left
      wrong. Docs take the version from `anuga.__version__` (no hardcoding) and
      README badges are version-agnostic, so neither needs a change.
- [ ] **RC on TestPyPI** (optional but cheap): `4.0.0rc1` tag, workflow_dispatch the
      wheel build, `pip install --index-url test.pypi.org` smoke.
- [x] **RC published**: `4.0.0rc1` on TestPyPI (21 artifacts, cp310–cp314,
      Linux/macOS arm64+x86_64/Windows, plus sdist), built via the
      `ANUGA_VERSION` override with no rc tag. Dispatch:
      `gh workflow run python-publish-pypi.yml -f publish_to=testpypi -f version=4.0.0rc1`
- [x] **Send the announcement** — draft at
      `claude/archive/RELEASE_4.0.0rc1_ANNOUNCEMENT.md`, addressed to Ole, Rudy, Petar,
      David and Jorge, with a specific ask for each. Then 3–4 days of soak.
      **Wait for Petar on Towradgi before Phase 3**: it is his case study, it is
      the evidence behind the #229 numbers in the release notes, and if he says
      the changed result looks wrong for that catchment the physics decision
      reopens rather than the wording.
- [x] **PETAR'S VERDICT — 2026-08-30: the culvert stage reconstruction on slopes
      is OK.** This was the Phase-3 gate. The #229 physics decision stands, the
      release-notes delta table stands as written, and the soak is closed.
      During the review he asked for the Collins St culvert logs and the
      upstream approach flow; the peak was extracted from the existing 24 h run
      (RUN_20260810_195336) as 32.4 m³/s at t=38 700 s, harness preserved in
      `sandpit/culvert_issue/`.

### Commits after 4.0.0rc1

The RC is not byte-identical to what 4.0.0 will be. Anything landing after it
goes here so the difference is deliberate and reviewable, rather than
discovered at tag time:

| commit | change | in the RC? |
|---|---|---|
| `32577739` | #237 — explain the missing `anuga/_version.py` instead of a bare ModuleNotFoundError | no |
| `9c5d1215` | release-plan bookkeeping (this table) | no |
| `13195657` | `claude/PPA_FEASIBILITY.md` for issue #25 — notes only | no |
| `64f53a6a` | **Fix the build under Cython 3.3.0** (dict views are not lists) | **no** |
| `9dc87f0f` | Towradgi reproduction recipe from the soak — notes only | no |

Judgement: no solver or API change among them, so the soak stands.

`64f53a6a` is the one that matters and it is worth stating plainly: Cython 3.3.0
landed on conda-forge on 23 Aug, a few hours *after* the RC wheels were built.
**The RC on TestPyPI therefore cannot be built from its sdist today** — the same
defect that left `main` unbuildable 23–30 Aug and that the released 3.3.10 tag
still carries. 4.0.0 fixes it. The difference between the soaked artifact and
the tag is in the safe direction (broken → fixed), and it is a build-time fix
that cannot move a simulation result, so this is a row rather than an RC rebuild.

Before cutting the release PR, confirm `origin/develop` is current: as of
2026-08-30 the local branch was 2 commits ahead (`bdbc9d8f` session notes,
`db7230c4` aws_run_gpu.sh) and unpushed.

## Phase 3 — Ship (Day 10–12)

Procedure per the 3.3.8 runbook (SESSION_GUIDE §"Release procedure"):

```bash
# 1. Release PR: develop → main (review per the Phase-0 policy decision)
gh pr create --repo anuga-community/anuga_core --base main --head develop \
  --title "Release 4.0.0"
# 2. Set the release date in CITATION.cff, then tag.
#    (version is already 4.0.0; date-released is deliberately absent until now)
git checkout main && git pull origin main
sed -i "s/^version: 4.0.0$/version: 4.0.0\ndate-released: '$(date +%F)'/" CITATION.cff
git commit -am "CITATION.cff: record the 4.0.0 release date"
#    annotated tag, BARE version (no v prefix), on the merge commit
git tag -a 4.0.0 -m "ANUGA 4.0.0"
git push origin 4.0.0
# 3. The GitHub Release is what triggers PyPI upload (a bare tag does NOT)
gh release create 4.0.0 --verify-tag --title "ANUGA 4.0.0" --notes-file RELEASE_NOTES.md
```

- [ ] Watch `python-publish-pypi.yml` → wheels + sdist on PyPI; `pip install anuga==4.0.0` smoke.
- [ ] `docker-publish.yml` fires on the Release → CPU image to GHCR automatically;
      trigger the GPU image via `workflow_dispatch` (`build_gpu=true`).
- [ ] conda-forge: feedstock bot opens the PR automatically; review it (new
      build deps are unlikely — meson stack unchanged since 3.3.2) and merge.
- [ ] Close the shipped issues milestone; announce (list, README news, Hydrata).

## Phase 4 — After (Day 12+)

- [ ] `main` becomes the 4.0.x patch base (same cherry-pick model as 3.3.x).
- [ ] Retire the "no develop→main" rule from ROADMAP.md/CLAUDE.md; normal cadence
      resumes — aim for a release every 2–3 months so a backlog like this cannot
      re-form.
- [ ] Open the 4.1 milestone seeded with: forcing-class removal (P2.10),
      default-mode decision (PLAN_default_mode2_cpu.md), cibuildwheel (#141),
      C-audit P1 items (review R3).

---

## Risks

| Risk | Likelihood | Mitigation |
|---|---|---|
| #229 moves sloping-bed structure results; users see changed numbers | certain (it's the fix) | Measure in Phase 1, publish the deltas in the notes with the physics rationale |
| conda-forge feedstock friction | low | Build stack unchanged since 3.3.2; review the bot PR same-day |
| GPU build regressions invisible to CI | medium | Phase-1 manual GPU gates on two architectures; known limitation, R2 of the review is the durable fix |
| A 976-commit release has an unknown regression | medium | RC soak + downstream pings; 4.0.1 within days is cheap once the patch line exists |
| Release stalls again on perfection | medium | The gates above are the complete list; nothing else blocks the tag |
