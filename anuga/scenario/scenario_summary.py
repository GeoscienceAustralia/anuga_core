#!/usr/bin/python
"""
Render an anuga_toml_run scenario as a self-contained HTML summary.

A "dry run" for a TOML scenario: parse the config (no mesh build, no evolve),
summarise the domain / mesh / friction / structures / forcing, and — for
rainfall — draw a catchment-mean hyetograph (15-minute intensity bars plus a
cumulative-depth curve) read straight from the gauge timeseries. The output is
a single self-contained HTML file (inline CSS + inline SVG, no external assets)
that can be opened in a browser.

Public API:
    build_summary_html(config_path, base_dir=None) -> str          # HTML text
    write_scenario_summary(config_path, output_html=None,
                           base_dir=None, open_browser=False) -> str  # -> path

Used by `anuga_toml_run --dry-run`.
"""

import os
import glob
import html
import collections


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def _load_toml(path):
    try:
        import tomllib
    except ImportError:  # Python < 3.11
        import tomli as tomllib
    with open(path, 'rb') as fh:
        return tomllib.load(fh)


def _read_timeseries(path):
    """Return (time_s, rate_mm_per_s) for a rainfall timeseries file.

    Accepts an ANUGA ``.tms`` file (a NetCDF 'rate' quantity, mm/s) or a CSV
    ``[time_s, rate_mm_per_hr]`` with an optional header (converted to mm/s).
    Returns ``None`` if the file cannot be read.
    """
    import numpy as np
    if not os.path.exists(path):
        return None
    try:
        if path.lower().endswith('.tms'):
            from anuga.file.netcdf import NetCDFFile
            fd = NetCDFFile(path, 'r')
            t = np.asarray(fd.variables['time'][:], dtype=float)
            r = np.asarray(fd.variables['rate'][:], dtype=float)  # mm/s
            fd.close()
            return t, r
        # CSV: [time_s, mm/hr] with optional single header row
        with open(path) as fh:
            first = fh.readline()
        skip = 1 if any(c.isalpha() for c in first) else 0
        arr = np.genfromtxt(path, delimiter=',', skip_header=skip)
        if arr.ndim != 2 or arr.shape[1] < 2:
            return None
        return arr[:, 0].astype(float), arr[:, 1].astype(float) / 3600.0
    except Exception:
        return None


def _read_two_col(path):
    """Return (x, y) float arrays from a 2-column CSV with an optional header.

    Used for inlet discharge hydrographs ``[time_s, discharge_m3s]``. Returns
    ``None`` if the file is missing or unreadable.
    """
    import numpy as np
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path) as fh:
            first = fh.readline()
        skip = 1 if any(c.isalpha() for c in first) else 0
        arr = np.genfromtxt(path, delimiter=',', skip_header=skip)
        if arr.ndim != 2 or arr.shape[1] < 2:
            return None
        return arr[:, 0].astype(float), arr[:, 1].astype(float)
    except Exception:
        return None


def _is_constant(y):
    import numpy as np
    y = np.asarray(y, dtype=float)
    return float(y.max() - y.min()) <= 1e-9 * max(1.0, abs(float(y.max())))


def _line_chart_svg(t, y, unit, accent='water'):
    """Compact line chart of *y* (in *unit*) versus time in hours."""
    import numpy as np
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    th = t / 3600.0
    T = th[-1] if th[-1] > 0 else 1.0
    ymax = float(y.max()) if y.max() > 0 else 1.0
    W, H, padL, padR, padT, padB = 720, 150, 40, 12, 14, 24
    pw, ph = W - padL - padR, H - padT - padB
    ybase = padT + ph

    def X(v):
        return padL + v / T * pw

    def Y(v):
        return padT + (1 - v / ymax) * ph

    pts = ' '.join(f'{X(a):.1f},{Y(b):.1f}' for a, b in zip(th, y))
    step = 6 if T > 14 else (2 if T > 4 else 1)
    ticks = []
    h = 0
    while h <= T + 1e-6:
        x = X(h)
        ticks.append(
            f'<line class="ax" x1="{x:.1f}" y1="{ybase:.1f}" x2="{x:.1f}" y2="{ybase + 4:.1f}"/>'
            f'<text class="axl" x="{x:.1f}" y="{H - 6}" text-anchor="middle">{h}h</text>')
        h += step
    axes = (
        f'<text class="hy-axl hy-{accent}" x="{padL - 6}" y="{padT + 4:.0f}" text-anchor="end">{ymax:g}</text>'
        f'<text class="hy-axl hy-{accent}" x="{padL - 6}" y="{ybase:.0f}" text-anchor="end">0</text>'
        f'<text class="axl" x="{padL - 6}" y="{padT - 3:.0f}" text-anchor="end">{esc(unit)}</text>')
    return (
        f'<svg viewBox="0 0 {W} {H}" class="hyeto" role="img" '
        f'aria-label="Inlet discharge over time">'
        f'<line class="grid" x1="{padL}" y1="{ybase:.1f}" x2="{padL + pw}" y2="{ybase:.1f}"/>'
        f'<polyline class="ln-{accent}" points="{pts}"/>{"".join(ticks)}{axes}</svg>')


# ---------------------------------------------------------------------------
# Rainfall aggregation
# ---------------------------------------------------------------------------

def _aggregate_rainfall(rain_entries, base_dir):
    """Catchment-mean hyetograph across all rainfall entries.

    Returns a dict with grid (s), mean_rate (mm/s), and summary scalars, or
    ``None`` when no timeseries could be read.
    """
    import numpy as np
    series = []
    for r in rain_entries:
        ts = r.get('timeseries_file')
        if not ts:
            continue
        got = _read_timeseries(os.path.join(base_dir, ts))
        if got is not None and len(got[0]) > 1:
            mult = float(r.get('multiplier', 1.0))
            # A .tms carries the raw gauge rate; the mm-based ones use a 1e-3
            # unit factor in the run. Only apply a multiplier that is NOT that
            # unit-conversion factor, so intensities read in real mm/hr.
            if abs(mult - 1.0e-3) > 1e-12:
                got = (got[0], got[1] * mult)
            series.append(got)
    if not series:
        return None

    tmax = max(float(t[-1]) for t, _ in series)
    if tmax <= 0:
        return None
    grid = np.arange(0.0, tmax + 1.0, 300.0)          # 5-minute grid
    acc = np.zeros_like(grid)
    for t, r in series:
        acc += np.interp(grid, t, r, left=0.0, right=0.0)
    mean_rate = acc / len(series)                     # mm/s
    inten = mean_rate * 3600.0                        # mm/hr
    try:
        depth = float(np.trapezoid(mean_rate, grid))  # mm
    except AttributeError:                            # numpy < 2.0
        depth = float(np.trapz(mean_rate, grid))
    wet = grid[inten > 0.5]
    return dict(
        grid=grid, mean_rate=mean_rate,
        n_series=len(series),
        peak=float(inten.max()),
        peak_t_h=float(grid[int(inten.argmax())]) / 3600.0,
        depth_mm=depth,
        duration_h=(float(wet.max() - wet.min()) / 3600.0) if len(wet) else 0.0,
        total_h=tmax / 3600.0,
    )


def _hyetograph_svg(rain):
    """SVG: 15-minute intensity bars (left axis, mm/hr) + cumulative-depth line
    (right axis, mm)."""
    import numpy as np
    grid, mean_rate = rain['grid'], rain['mean_rate']
    inten = mean_rate * 3600.0
    ipeak = rain['peak'] or 1.0
    total = rain['depth_mm'] or 1.0
    dt = grid[1] - grid[0]
    cum = np.cumsum(mean_rate * dt)
    T = grid[-1] or 1.0

    W, H, padL, padR, padT, padB = 720, 216, 34, 40, 26, 26
    pw, ph = W - padL - padR, H - padT - padB
    ybase = padT + ph

    def X(t):
        return padL + t / T * pw

    def YL(v):
        return padT + (1 - v / ipeak) * ph

    def YR(v):
        return padT + (1 - v / total) * ph

    binw = 900.0                                    # 15-minute bars
    bars = []
    for b in range(int(np.ceil(T / binw))):
        lo, hi = b * binw, min((b + 1) * binw, T)
        sel = (grid >= lo) & (grid <= hi) if hi >= T else (grid >= lo) & (grid < hi)
        if not sel.any():
            continue
        v = float(inten[sel].mean())
        x0, x1 = X(lo), X(hi)
        bars.append(
            f'<rect class="hy-bar" x="{x0 + (x1 - x0) * 0.09:.1f}" y="{YL(v):.1f}" '
            f'width="{max((x1 - x0) * 0.82, 0.6):.1f}" '
            f'height="{max(ybase - YL(v), 0):.1f}"/>')
    cum_pts = ' '.join(f'{X(t):.1f},{YR(c):.1f}' for t, c in zip(grid, cum))

    # x-axis: hour ticks at a sensible spacing for the record length
    total_h = T / 3600.0
    step = 6 if total_h > 14 else (2 if total_h > 4 else 1)
    ticks = []
    h = 0
    while h * 3600.0 <= T + 1:
        x = X(h * 3600.0)
        ticks.append(
            f'<line class="ax" x1="{x:.1f}" y1="{ybase:.1f}" x2="{x:.1f}" y2="{ybase + 4:.1f}"/>'
            f'<text class="axl" x="{x:.1f}" y="{H - 6}" text-anchor="middle">{h}h</text>')
        h += step
    axes = (
        f'<text class="hy-axl hy-l" x="{padL - 6}" y="{padT + 4:.0f}" text-anchor="end">{ipeak:.0f}</text>'
        f'<text class="hy-axl hy-l" x="{padL - 6}" y="{ybase:.0f}" text-anchor="end">0</text>'
        f'<text class="hy-axl hy-r" x="{W - padR + 6}" y="{padT + 4:.0f}" text-anchor="start">{total:.0f}</text>'
        f'<text class="hy-axl hy-r" x="{W - padR + 6}" y="{ybase:.0f}" text-anchor="start">0</text>')
    legend = (
        f'<g class="hy-legend" transform="translate({padL + 6},6)">'
        f'<rect class="hy-bar" x="0" y="1" width="12" height="9"/>'
        f'<text x="17" y="9">intensity (mm/hr, 15-min)</text>'
        f'<line class="hy-cum" x1="185" y1="6" x2="205" y2="6"/>'
        f'<text x="210" y="9">cumulative depth (mm)</text></g>')
    return (
        f'<svg viewBox="0 0 {W} {H}" class="hyeto" role="img" '
        f'aria-label="Rainfall intensity bars with cumulative depth over time">'
        f'<line class="grid" x1="{padL}" y1="{ybase:.1f}" x2="{padL + pw}" y2="{ybase:.1f}"/>'
        f'{"".join(bars)}<polyline class="hy-cum" points="{cum_pts}"/>'
        f'{"".join(ticks)}{axes}{legend}</svg>')


# ---------------------------------------------------------------------------
# Friction categorisation (generic — by roughness magnitude)
# ---------------------------------------------------------------------------

def _friction_tier(n):
    """(css-category, tier-label) for a Manning's n by magnitude."""
    if n < 0.02:
        return 'water', 'smooth / channel'
    if n < 0.06:
        return 'stone', 'open / paved'
    if n < 1.0:
        return 'reed', 'vegetated / rough'
    return 'silt', 'building / blockage'


# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------

esc = html.escape


def _badge(label, value):
    return (f'<span class="badge"><span class="badge-k">{esc(str(label))}</span>'
            f'<span class="badge-v">{esc(str(value))}</span></span>')


def _stat(v, label, note=''):
    n = f'<div class="stat-n">{esc(note)}</div>' if note else ''
    return (f'<div class="stat"><div class="stat-v">{esc(str(v))}</div>'
            f'<div class="stat-l">{esc(label)}</div>{n}</div>')


def _glob_count(patterns, base_dir):
    n = 0
    for pat in patterns or []:
        n += len(glob.glob(os.path.join(base_dir, pat)))
    return n


def build_summary_html(config_path, base_dir=None):
    """Return the self-contained HTML summary for a TOML scenario as a string."""
    base_dir = base_dir or os.path.dirname(os.path.abspath(config_path)) or '.'
    cfg = _load_toml(config_path)
    p = cfg.get('project', {})
    m = cfg.get('mesh', {})
    ic = cfg.get('initial_conditions', {})
    bcs = cfg.get('boundary_conditions', {}).get('boundaries', [])
    inlets = cfg.get('inlets', [])
    rain_entries = cfg.get('rainfall', [])
    culverts = cfg.get('culverts', [])
    weirs = cfg.get('weirs', [])
    bridges = cfg.get('bridges', [])
    pumps = cfg.get('pumping_stations', [])
    erosion = cfg.get('erosion', [])
    sediment = cfg.get('sediment', {})
    irs = m.get('interior_regions', [])

    scenario = str(p.get('scenario', os.path.splitext(os.path.basename(config_path))[0]))
    compute = p.get('compute_mode')
    if compute is None and 'multiprocessor_mode' in p:
        compute = {1: 'legacy', 2: 'unified'}.get(int(p['multiprocessor_mode']), 'legacy')

    # ---- friction --------------------------------------------------------
    fr = ic.get('friction', [])
    fr_numeric = [f for f in fr if isinstance(f.get('value'), (int, float))]
    fr_files = [f for f in fr if not isinstance(f.get('value'), (int, float))]
    fr_counts = collections.Counter(round(float(f['value']), 4) for f in fr_numeric)

    # ---- rainfall --------------------------------------------------------
    rain = _aggregate_rainfall(rain_entries, base_dir) if rain_entries else None

    # ---- structures ------------------------------------------------------
    def _enabled(lst):
        return [s for s in lst if s.get('enabled', True)]
    culverts, weirs, bridges, pumps, erosion = map(
        _enabled, (culverts, weirs, bridges, pumps, erosion))
    n_struct = len(culverts) + len(weirs) + len(bridges)

    # ---- riverwalls / breaklines ----------------------------------------
    n_riverwalls = _glob_count(m.get('riverwall_csv_files', []), base_dir)

    # ---- stat cards ------------------------------------------------------
    stats = []
    if 'default_res' in m:
        stats.append(_stat(f"{float(m['default_res']):g}", 'm² max triangle', 'background'))
    if irs:
        res = [float(r.get('resolution', 0)) for r in irs if 'resolution' in r]
        stats.append(_stat(len(irs), 'interior regions',
                           f'to {min(res):g} m²' if res else ''))
    if n_riverwalls:
        stats.append(_stat(n_riverwalls, 'riverwalls', 'breaklines'))
    if fr_numeric or fr_files:
        stats.append(_stat(len(fr), 'friction zones', f'{len(fr_counts)} n-values'))
    if n_struct:
        stats.append(_stat(n_struct, 'structures',
                           f'{len(culverts)} culv / {len(weirs)} weir / {len(bridges)} bridge'))
    if rain_entries:
        note = f"{rain['n_series']} read" if rain else 'timeseries'
        stats.append(_stat(len(rain_entries), 'rainfall inputs', note))
    if inlets:
        stats.append(_stat(len(inlets), 'inlets', 'line sources'))
    if pumps:
        stats.append(_stat(len(pumps), 'pumping stations'))
    if erosion:
        stats.append(_stat(len(erosion), 'erosion operators'))
    if sediment:
        fractions = sediment.get('fractions', [])
        note = ', '.join(str(f.get('name', '?')) for f in fractions)
        stats.append(_stat(len(fractions), 'sediment fractions', note))

    # ---- badges ----------------------------------------------------------
    badges = []
    if 'projection_information' in p:
        badges.append(_badge('CRS', p['projection_information']))
    if 'flow_algorithm' in p:
        badges.append(_badge('algorithm', p['flow_algorithm']))
    if compute:
        badges.append(_badge('compute', compute))
    # Times may be quoted arithmetic ("24*3600"), as the parser allows.
    from anuga.scenario.parse_input_data_toml import _num
    if p.get('starttime'):
        badges.append(_badge('start time', f"{_num(p['starttime']):g} s"))
    if 'finaltime' in p:
        badges.append(_badge('final time', f"{_num(p['finaltime']):g} s"))
    if 'yieldstep' in p:
        badges.append(_badge('yield', f"{_num(p['yieldstep']):g} s"))

    sections = []

    # ---- domain & mesh ---------------------------------------------------
    kv = []
    if 'default_res' in m:
        kv.append(('Background resolution', f"{float(m['default_res']):g} m²"))
    elev = ic.get('elevation', [])
    if elev:
        ev = elev[0].get('value')
        kv.append(('Elevation source',
                   os.path.basename(str(ev)) if isinstance(ev, str) else f'constant {ev}'))
    if n_riverwalls:
        kv.append(('Riverwalls', f'{n_riverwalls} breaklines'))
    if bcs:
        kv.append(('Boundary tags', ', '.join(str(b.get('tag', '?')) for b in bcs)))
    kv_html = ''.join(
        f'<div class="kv"><span class="k">{esc(k)}</span><span class="v">{esc(str(v))}</span></div>'
        for k, v in kv)
    ir_html = ''.join(
        f'<li><code>{esc(os.path.basename(str(ir.get("polygon", "?"))))}</code>'
        f'<span class="ir-res">{float(ir.get("resolution", 0)):g} m²</span></li>'
        for ir in irs)
    mesh_right = (f'<div class="card"><p class="eyebrow" style="margin-bottom:.6rem">'
                  f'Refined interior regions</p><ul class="regions">{ir_html}</ul></div>'
                  if irs else '')
    if kv_html or mesh_right:
        grid_cls = 'grid2' if mesh_right else ''
        sections.append(
            f'<section><h2>Domain &amp; mesh</h2><div class="{grid_cls}">'
            f'<div class="card">{kv_html}</div>{mesh_right}</div></section>')

    # ---- boundaries ------------------------------------------------------
    if bcs:
        chips = ''.join(
            f'<div class="bc"><span class="bc-tag">{esc(str(b.get("tag", "?")))}</span>'
            f'<span class="bc-type">{esc(str(b.get("type", "?")))}</span></div>' for b in bcs)
        sections.append(f'<section><h2>Boundary conditions</h2><div class="bcs">{chips}</div></section>')

    # ---- friction --------------------------------------------------------
    if fr_counts:
        fmax = max(fr_counts.values())
        rows = []
        for val, cnt in sorted(fr_counts.items()):
            cat, tier = _friction_tier(val)
            rows.append(
                f'<div class="fbar"><span class="fbar-n">{val:g}</span>'
                f'<span class="fbar-track"><span class="fbar-fill c-{cat}" '
                f'style="width:{cnt / fmax * 100:.1f}%"></span></span>'
                f'<span class="fbar-c">{cnt}</span><span class="fbar-l">{esc(tier)}</span></div>')
        extra = (f'<p class="note">+ {len(fr_files)} spatial friction input(s) from file</p>'
                 if fr_files else '')
        sections.append(
            f'<section><h2>Friction — Manning\'s n · {len(fr)} zones</h2>'
            f'<div class="card"><div class="fbars">{"".join(rows)}</div>{extra}</div></section>')

    # ---- structures ------------------------------------------------------
    def _size(s):
        if s.get('type') == 'boyd_pipe' or 'diameter' in s:
            return f"⌀ {s.get('diameter', '?')} m"
        w, h = s.get('width', '?'), s.get('height', s.get('width', '?'))
        return f"{w} × {h} m"
    srows = []
    for s in culverts + weirs:
        is_bridge = 'bridge' in str(s.get('label', '')).lower()
        kind = 'bridge' if is_bridge else ('pipe' if s.get('type') == 'boyd_pipe'
                                           else ('weir' if s in weirs else 'box'))
        srows.append(
            f'<tr><td class="cv-label">{esc(str(s.get("label", "?")).replace("_", " "))}</td>'
            f'<td><span class="chip chip-{kind}">{kind}</span></td>'
            f'<td class="num">{esc(_size(s))}</td>'
            f'<td class="num">{s.get("manning", "—")}</td>'
            f'<td class="num">{s.get("apron", "—")}</td>'
            f'<td class="num">{s.get("enquiry_gap", "—")}</td></tr>')
    for b in bridges:
        srows.append(
            f'<tr><td class="cv-label">{esc(str(b.get("label", "?")).replace("_", " "))}</td>'
            f'<td><span class="chip chip-bridge">bridge</span></td>'
            f'<td class="num">deck {b.get("deck_elevation", "—")} m</td>'
            f'<td class="num">—</td><td class="num">—</td>'
            f'<td class="num">{b.get("enquiry_gap", "—")}</td></tr>')
    if srows:
        sections.append(
            f'<section><h2>Structures — {n_struct} culverts, weirs &amp; bridges</h2>'
            f'<div class="tablewrap"><table><thead><tr><th>Structure</th><th>Type</th>'
            f'<th>Size</th><th>Manning</th><th>Apron</th><th>Enq. gap</th></tr></thead>'
            f'<tbody>{"".join(srows)}</tbody></table></div></section>')

    # ---- forcing ---------------------------------------------------------
    forcing = []
    if rain_entries:
        if rain:
            body = (
                f'{_hyetograph_svg(rain)}'
                f'<div class="rain-stats">'
                f'<div class="rs"><span class="rs-v">{rain["peak"]:.0f}</span><span class="rs-l">peak mm/hr</span></div>'
                f'<div class="rs"><span class="rs-v">{rain["peak_t_h"]:.1f} h</span><span class="rs-l">time of peak</span></div>'
                f'<div class="rs"><span class="rs-v">{rain["depth_mm"]:.0f} mm</span><span class="rs-l">total depth</span></div>'
                f'<div class="rs"><span class="rs-v">{rain["duration_h"]:.0f} h</span><span class="rs-l">wet duration</span></div>'
                f'<div class="rs"><span class="rs-v">{len(rain_entries)}</span><span class="rs-l">inputs</span></div>'
                f'</div>')
        else:
            body = ('<p class="note">Rainfall timeseries could not be read from disk; '
                    'showing input count only.</p>')
        forcing.append(
            f'<div class="card"><p class="eyebrow" style="margin-bottom:0">'
            f'Rainfall — catchment-mean hyetograph · {len(rain_entries)} input(s)</p>{body}</div>')
    if inlets:
        import numpy as np
        blocks = []
        for i in inlets:
            name = esc(str(i.get('name', 'inlet')).replace('_', ' '))
            ts = i.get('timeseries_file')
            series = _read_two_col(os.path.join(base_dir, ts)) if ts else None
            if series is None or len(series[0]) < 1:
                blocks.append(
                    f'<div class="kv"><span class="k">{name}</span>'
                    f'<span class="v">line source</span></div>')
            elif len(series[1]) == 1 or _is_constant(series[1]):
                q = float(series[1][0])
                blocks.append(
                    f'<div class="kv"><span class="k">{name}</span>'
                    f'<span class="v">constant {q:g} m³/s</span></div>')
            else:
                t, q = series
                blocks.append(
                    f'<div class="inlet-plot"><div class="kv" style="border:0;padding-bottom:.2rem">'
                    f'<span class="k">{name}</span>'
                    f'<span class="v">peak {float(q.max()):g} · mean {float(q.mean()):g} m³/s</span></div>'
                    f'{_line_chart_svg(t, q, "m³/s", accent="reed")}</div>')
        forcing.append(
            f'<div class="card" style="margin-top:1.1rem">'
            f'<p class="eyebrow" style="margin-bottom:.6rem">Inlets — line sources '
            f'(discharge)</p>{"".join(blocks)}</div>')
    if forcing:
        sections.append(f'<section><h2>Forcing</h2>{"".join(forcing)}</section>')

    title = f"{scenario.replace('_', ' ')} — ANUGA scenario"
    lede = (f"Shallow-water flood scenario declared in "
            f"<code>{esc(os.path.basename(config_path))}</code>. This is a dry-run summary "
            f"of the configuration — no mesh was built and no simulation was run.")
    body_html = (
        f'<div class="wrap"><p class="eyebrow">ANUGA scenario · anuga_toml_run --dry-run</p>'
        f'<h1>{esc(scenario.replace("_", " "))}</h1><p class="lede">{lede}</p>'
        f'<div class="badges">{"".join(badges)}</div>'
        f'<div class="stats">{"".join(stats)}</div>'
        f'{"".join(sections)}'
        f'<footer>Dry-run summary of <code>{esc(os.path.basename(config_path))}</code> '
        f'· no simulation was run</footer></div>')
    return f"<title>{esc(title)}</title>\n<style>{_CSS}</style>\n{body_html}"


def write_scenario_summary(config_path, output_html=None, base_dir=None,
                           open_browser=False):
    """Write the HTML summary and return its path.

    Parameters
    ----------
    config_path : str
        Path to the scenario ``.toml`` file.
    output_html : str, optional
        Output path. Defaults to ``<config_stem>_summary.html`` next to the TOML.
    base_dir : str, optional
        Directory that relative paths in the TOML resolve against
        (defaults to the TOML's directory).
    open_browser : bool
        If True, open the written file in the default web browser.
    """
    config_path = os.path.abspath(config_path)
    if output_html is None:
        stem = os.path.splitext(os.path.basename(config_path))[0]
        output_html = os.path.join(os.path.dirname(config_path), stem + '_summary.html')
    html_text = build_summary_html(config_path, base_dir=base_dir)
    with open(output_html, 'w') as fh:
        fh.write(html_text)
    if open_browser:
        _open_in_browser(output_html)
    return output_html


def _open_in_browser(path):
    """Open *path* in the default browser; fall back to printing on failure."""
    import webbrowser
    url = 'file://' + os.path.abspath(path)
    try:
        if webbrowser.open(url):
            return True
    except Exception:
        pass
    print(f'Could not open a browser automatically; open this file manually:\n  {path}')
    return False


# ---------------------------------------------------------------------------
# Styles (inline; hydrology palette; theme-aware)
# ---------------------------------------------------------------------------

_CSS = """
:root{
  --paper:#eef2f5;--surface:#fff;--surface-2:#e7edf1;--line:#d4dde4;
  --ink:#10222e;--ink-soft:#41586a;--ink-faint:#8496a3;
  --water:#0e6ba8;--silt:#c07f22;--reed:#4e8d5b;--stone:#6b7c8a;--focus:#0e6ba8;
  --sans:system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  --mono:ui-monospace,"SF Mono","Cascadia Code",Menlo,Consolas,monospace;
}
@media (prefers-color-scheme:dark){:root{
  --paper:#0b141b;--surface:#101d27;--surface-2:#0d1822;--line:#213240;
  --ink:#e7eef3;--ink-soft:#9db0bd;--ink-faint:#63788a;
  --water:#48b5e6;--silt:#e0a94b;--reed:#7ec48a;--stone:#8a9aa7;--focus:#48b5e6;}}
:root[data-theme="light"]{
  --paper:#eef2f5;--surface:#fff;--surface-2:#e7edf1;--line:#d4dde4;
  --ink:#10222e;--ink-soft:#41586a;--ink-faint:#8496a3;
  --water:#0e6ba8;--silt:#c07f22;--reed:#4e8d5b;--stone:#6b7c8a;--focus:#0e6ba8;}
:root[data-theme="dark"]{
  --paper:#0b141b;--surface:#101d27;--surface-2:#0d1822;--line:#213240;
  --ink:#e7eef3;--ink-soft:#9db0bd;--ink-faint:#63788a;
  --water:#48b5e6;--silt:#e0a94b;--reed:#7ec48a;--stone:#8a9aa7;--focus:#48b5e6;}
*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);font-family:var(--sans);
  line-height:1.55;-webkit-font-smoothing:antialiased}
.wrap{max-width:1000px;margin:0 auto;padding:clamp(1.2rem,4vw,3rem)}
.eyebrow{font-family:var(--mono);text-transform:uppercase;letter-spacing:.14em;
  font-size:.7rem;color:var(--ink-soft);margin:0 0 .5rem}
h1{font-size:clamp(1.9rem,5vw,2.9rem);line-height:1.05;margin:.1rem 0 .4rem;
  font-weight:680;letter-spacing:-.02em;text-wrap:balance}
.lede{color:var(--ink-soft);max-width:64ch;margin:0 0 1.3rem;font-size:1.02rem}
.badges{display:flex;flex-wrap:wrap;gap:.5rem;margin-bottom:2.4rem}
.badge{display:inline-flex;align-items:center;gap:.5rem;border:1px solid var(--line);
  border-radius:999px;padding:.32rem .8rem;background:var(--surface)}
.badge-k{font-family:var(--mono);font-size:.66rem;text-transform:uppercase;
  letter-spacing:.1em;color:var(--ink-faint)}
.badge-v{font-family:var(--mono);font-size:.82rem;color:var(--water);font-weight:600}
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(148px,1fr));gap:1px;
  background:var(--line);border:1px solid var(--line);border-radius:12px;
  overflow:hidden;margin-bottom:2.6rem}
.stat{background:var(--surface);padding:1.1rem 1.2rem}
.stat-v{font-family:var(--mono);font-size:1.75rem;font-weight:640;
  font-variant-numeric:tabular-nums;letter-spacing:-.02em}
.stat-l{font-size:.82rem;color:var(--ink-soft);margin-top:.15rem}
.stat-n{font-family:var(--mono);font-size:.68rem;color:var(--ink-faint);margin-top:.35rem}
section{margin-bottom:2.6rem}
h2{font-size:.78rem;font-family:var(--mono);text-transform:uppercase;letter-spacing:.14em;
  color:var(--ink-soft);font-weight:600;margin:0 0 1rem;padding-bottom:.55rem;
  border-bottom:1px solid var(--line)}
.card{background:var(--surface);border:1px solid var(--line);border-radius:12px;padding:1.2rem 1.35rem}
.grid2{display:grid;grid-template-columns:1fr 1fr;gap:1.1rem}
@media(max-width:640px){.grid2{grid-template-columns:1fr}}
.kv{display:flex;justify-content:space-between;gap:1rem;padding:.4rem 0;border-bottom:1px solid var(--line)}
.kv:last-child{border-bottom:0}
.kv .k{color:var(--ink-soft)}.kv .v{font-family:var(--mono);text-align:right}
.note{color:var(--ink-faint);font-size:.82rem;margin:.8rem 0 0}
ul.regions{list-style:none;margin:0;padding:0}
ul.regions li{display:flex;justify-content:space-between;align-items:center;
  padding:.4rem 0;border-bottom:1px solid var(--line)}
ul.regions li:last-child{border-bottom:0}
ul.regions code{font-family:var(--mono);font-size:.85rem}
.ir-res{font-family:var(--mono);color:var(--water);font-size:.85rem}
.bcs{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:.7rem}
.bc{border:1px solid var(--line);border-radius:10px;padding:.7rem .85rem;background:var(--surface-2)}
.bc-tag{font-family:var(--mono);font-weight:640;display:block}
.bc-type{font-size:.78rem;color:var(--ink-soft)}
.fbars{display:flex;flex-direction:column;gap:.55rem}
.fbar{display:grid;grid-template-columns:3rem 1fr 2.2rem auto;align-items:center;gap:.7rem}
.fbar-n{font-family:var(--mono);font-variant-numeric:tabular-nums;text-align:right;font-size:.9rem}
.fbar-track{background:var(--surface-2);border-radius:5px;height:14px;overflow:hidden}
.fbar-fill{display:block;height:100%;border-radius:5px}
.fbar-c{font-family:var(--mono);font-variant-numeric:tabular-nums;color:var(--ink-soft);font-size:.85rem}
.fbar-l{font-size:.85rem;color:var(--ink-soft)}
.c-silt{background:var(--silt)}.c-water{background:var(--water)}
.c-reed{background:var(--reed)}.c-stone{background:var(--stone)}
.tablewrap{overflow-x:auto;border:1px solid var(--line);border-radius:12px}
table{border-collapse:collapse;width:100%;font-size:.88rem}
th,td{text-align:left;padding:.6rem .9rem;border-bottom:1px solid var(--line);white-space:nowrap}
thead th{font-family:var(--mono);font-size:.66rem;text-transform:uppercase;letter-spacing:.1em;
  color:var(--ink-faint);font-weight:600;position:sticky;top:0;background:var(--surface)}
tbody tr:last-child td{border-bottom:0}
tbody tr:hover{background:var(--surface-2)}
.cv-label{font-weight:560}
td.num{font-family:var(--mono);font-variant-numeric:tabular-nums;color:var(--ink-soft)}
.chip{font-family:var(--mono);font-size:.68rem;text-transform:uppercase;letter-spacing:.06em;
  padding:.16rem .5rem;border-radius:6px;border:1px solid transparent}
.chip-box{background:color-mix(in srgb,var(--water) 16%,transparent);color:var(--water);
  border-color:color-mix(in srgb,var(--water) 35%,transparent)}
.chip-pipe{background:color-mix(in srgb,var(--reed) 16%,transparent);color:var(--reed);
  border-color:color-mix(in srgb,var(--reed) 35%,transparent)}
.chip-weir{background:color-mix(in srgb,var(--stone) 20%,transparent);color:var(--stone);
  border-color:color-mix(in srgb,var(--stone) 40%,transparent)}
.chip-bridge{background:color-mix(in srgb,var(--silt) 18%,transparent);color:var(--silt);
  border-color:color-mix(in srgb,var(--silt) 38%,transparent)}
.hyeto{width:100%;height:auto;display:block;margin:.5rem 0 .2rem;overflow:visible}
.hy-bar{fill:var(--water);opacity:.55}
.hy-cum{fill:none;stroke:var(--silt);stroke-width:1.8;stroke-linejoin:round}
.ln-water{fill:none;stroke:var(--water);stroke-width:1.8;stroke-linejoin:round}
.ln-reed{fill:none;stroke:var(--reed);stroke-width:1.8;stroke-linejoin:round}
.inlet-plot{padding:.4rem 0;border-bottom:1px solid var(--line)}
.inlet-plot:last-child{border-bottom:0}
.hy-axl{font-family:var(--mono);font-size:10px;font-variant-numeric:tabular-nums}
.hy-l{fill:var(--water)}.hy-r{fill:var(--silt)}.hy-reed{fill:var(--reed)}.hy-water{fill:var(--water)}
.hy-legend text{fill:var(--ink-soft);font-family:var(--mono);font-size:10px}
.grid{stroke:var(--line);stroke-width:1}
.ax{stroke:var(--ink-faint);stroke-width:1}
.axl{fill:var(--ink-faint);font-family:var(--mono);font-size:10px}
.rain-stats{display:flex;flex-wrap:wrap;gap:1.5rem;margin-top:.5rem;padding-top:.9rem;border-top:1px solid var(--line)}
.rs{display:flex;flex-direction:column;gap:.15rem}
.rs-v{font-family:var(--mono);font-size:1.2rem;font-weight:640;font-variant-numeric:tabular-nums}
.rs-l{font-size:.68rem;color:var(--ink-soft);font-family:var(--mono);text-transform:uppercase;letter-spacing:.09em}
footer{margin-top:2rem;padding-top:1.2rem;border-top:1px solid var(--line);
  color:var(--ink-faint);font-size:.8rem;font-family:var(--mono)}
footer code{color:var(--ink-soft)}
a{color:var(--water)}
"""
