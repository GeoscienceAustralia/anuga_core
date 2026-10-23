# Simple scenario — dam break in a flat box

The smallest useful `anuga_toml_run` scenario. A 100 m × 100 m flat basin with
solid (reflective) walls on all four sides. A 4 m deep column of water fills the
western 40 m — "the reservoir" — and the rest starts dry. When the run begins
the column collapses and sloshes back and forth off the walls.

Everything it needs is in this directory; no external data.

## Files

| File | Role |
|------|------|
| `dam_break.toml` | the scenario definition (read this first — it is heavily commented) |
| `bounding_polygon.csv` | `x,y` polygon for the 100 m × 100 m domain outline (no header) |
| `dam.csv` | `x,y` polygon for the initial reservoir (western 40 m, no header) |

## Run it

```bash
# serial
anuga_toml_run dam_break.toml

# parallel (e.g. 4 processes)
mpiexec -np 4 anuga_toml_run dam_break.toml
```

In a source checkout without the installed command, use
`python <repo>/scripts/anuga_toml_run.py dam_break.toml`.

## What you get

A new `OUTPUT/RUN_<timestamp>_dam_break/` directory containing:

- `dam_break.sww` — the result (open in the ANUGA SWW viewer or
  `anuga.SWW_plotter`),
- `Simulation_logfile.log` — the full run log,
- `code/` — a copy of the inputs used for this run,
- peak-quantity rasters (stage/depth/velocity maxima).

The run is tiny (~500 triangles, 30 s simulated) and finishes in a few seconds.

> A final `GeoTif creation failed: ...` message is harmless — it is a known
> issue in the optional GeoTIFF export step and does not affect the `.sww`
> output.

## Things to try

- Change `flow_algorithm` from `"DE1"` to `"DE0"` (faster, first-order).
- Move the reservoir: edit the `dam.csv` polygon, or its `value` (water depth)
  in `[initial_conditions.stage]`.
- Open one wall to the sea: change the `east` boundary from `Reflective` to a
  `Stage` type (you will need a small `time,stage` CSV — see the `complex/`
  scenario for an example).

See [`../README.md`](../README.md) for the file-format reference and how this
relates to the full Cairns example.
