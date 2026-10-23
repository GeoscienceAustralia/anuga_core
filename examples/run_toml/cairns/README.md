# Cairns scenario — a real DEM tsunami run

The full real-world `anuga_toml_run` example: a tsunami inundation model of
Cairns, Queensland, driven by a real Digital Elevation Model and a shapefile
domain outline. Where [`../simple/`](../simple/) and [`../complex/`](../complex/)
use tiny hand-made CSVs so they run anywhere in seconds, this one uses an actual
dataset and demonstrates every major TOML section.

It exercises:

- **elevation from a DEM raster** (`cairns.asc`) rather than a point CSV,
- a **shapefile bounding polygon** (`cairnsextent.shp`) that carries its boundary
  tags as feature attributes — no `[[mesh.boundary_tags]]` block needed,
- a sinusoidal tsunami wave on the ocean boundary and mean-sea-level elsewhere,
- one of each hydraulic structure: a **bridge**, a **culvert**, a **pipe
  culvert**, a **weir**, and a **pumping station**,
- the projection given as an EPSG code (`"EPSG:32755"`, UTM zone 55 south), which
  propagates into the SWW georeferencing.

## Data layout

The heavy DEM is **shared**, not duplicated: it lives once under
[`../../data/cairns/`](../../data/cairns/) and is referenced by this scenario
(and by the legacy [`../../cairns_toml_excel/`](../../cairns_toml_excel/) Excel
example) via a relative path. The small, scenario-specific config files are
local:

| Path | Role |
|------|------|
| `../../data/cairns/cairns.asc` | the DEM raster (≈9 MB, shared) |
| `cairns_example.toml` | the scenario definition (heavily commented — documents every section) |
| `cairns_mesh/cairnsextent.shp` (+ `.shx/.dbf/.prj`) | domain outline with tagged boundary edges |
| `cairns_mesh/*_exchange_*.csv`, `*_deck.csv`, `*_rating_curve.csv` | structure exchange lines / decks |
| `cairns_mesh/pump1_*.csv` | pumping-station basin and wet-well/discharge lines |
| `cairns_boundarycond/sine_30m.csv` | the incoming tsunami wave (`time,stage`) |
| `cairns_boundarycond/msl_boundary.csv` | mean-sea-level boundary (`time,stage`) |
| `user_functions.py` | optional evolve-loop reporting hooks |

## Run it

```bash
# serial
anuga_toml_run cairns_example.toml

# parallel (e.g. 6 processes)
mpiexec -np 6 anuga_toml_run cairns_example.toml
```

In a source checkout without the installed command, use
`python <repo>/scripts/anuga_toml_run.py cairns_example.toml`.

> **This is a real run.** The default `finaltime` is 21600 s (6 hours of model
> time). For a quick smoke test, lower `finaltime` in `[project]` (e.g. to a few
> hundred seconds) — the mesh is coarse (`default_res = 1e6`) so setup is fast.

## What you get

A new `OUTPUT/RUN_<timestamp>_cairns/` directory with `cairns.sww`, the
`Simulation_logfile.log`, a `code/` copy of the inputs, and peak-quantity
rasters. The SWW carries the EPSG:32755 georeferencing so it overlays correctly
in GIS.

See [`../README.md`](../README.md) for the TOML file-format reference, and
[`cairns_example.toml`](cairns_example.toml) itself for fully documented
examples of every section — including the structure blocks that
[`../complex/`](../complex/) only points at.
