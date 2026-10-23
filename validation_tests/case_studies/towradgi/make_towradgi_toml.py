#!/usr/bin/env python3
"""
Generate a faithful anuga_toml_run config (towradgi.toml) that reproduces
run_small_towradgi.py, plus the small helper CSVs it references.

Run from the towradgi case-study directory:

    python make_towradgi_toml.py

It writes:
  * towradgi.toml
  * bounding_polygon.csv, inlet_line.csv, inlet_hydrograph_<Q>.csv, east_stage.csv
  * culvert_exchange/<label>_{0,1}.csv        (exchange lines for each culvert)

The friction and rainfall lists are built from the actual data on disk
(Model/Mannings, Model/Creeks, Model/Buildings, Forcing/Rainfall/*), so this
stays in sync with the dataset. Culvert/bridge geometry and the Manning value
exceptions are transcribed from run_small_towradgi.py.

Requires the runner enhancements: rainfall .tms support and dict `losses`.

Three settings depart from run_small_towradgi.py, following the 2026 validation
sensitivity runs:
  * the creek inflow is 3 m^3/s, not 20. With 20 the Pioneer Rd gauge stays
    about 0.7 m above the record after the peak; with 3 it follows it.
  * the run starts at t = -3 h ([project] starttime), a burn-in that wets the
    creek with the inflow before the rain, which starts at t = 0. Model time
    then matches the rainfall clock; the Pioneer Rd logger runs 9 h ahead.
  * the flow algorithm is DE_ader2, not DE1. It gives the same flood-mark
    and gauge errors (peak levels within 0.05 m) and runs about 25% faster.
"""

import glob
import os

# ---- domain constants (project.py) -----------------------------------------
W, N, E, S = 303517, 6195670, 308570, 6193140
EPSG = 32756
MAX_TRI_AREA = 1000.0
BASE_FRICTION = 0.04
CHANNEL_MANNING = 0.03
SMOOTH_TS = 30.0
INLET_Q = 3.0            # creek inflow [m^3/s]; run_small_towradgi.py uses 20
# The inflow is in the file name, so a stale 20 m^3/s inlet_hydrograph.csv
# from an older run of this script is not picked up by mistake.
INLET_FILE = f'inlet_hydrograph_{INLET_Q:g}.csv'

# ---- Manning value rules (from ManningList) ---------------------------------
# Numbered Model/Mannings/<n>.csv are 0.15 unless listed here; named files and
# the Creeks/Buildings groups have fixed values.
MANNING_EXCEPTIONS = {   # basename (without .csv) -> n
    '1': 0.04, '9': 0.04, '13': 0.04, '18': 0.045, '18e': 0.08,
    '24': 0.05, '36': 0.05, '59': 0.08, '61': 0.08, '63': 0.08,
    '71': 0.05, '77': 0.07, 'Railway': 0.04, 'Escarpement': 0.15,
}
MANNING_DEFAULT = 0.15

# ---- Culverts & bridges (transcribed from run_small_towradgi.py) ------------
# el0/el1 are the two exchange-line coordinate pairs. Bridges are box culverts
# with pier losses. losses is a named-component dict (Boyd sums it).
_LC = {'inlet': 0.5, 'outlet': 1.0, 'bend': 0.0, 'grate': 0.0, 'pier': 0.0, 'other': 0.0}
_LP = {'inlet': 0.0, 'outlet': 0.0, 'bend': 0.0, 'grate': 0.0, 'pier': 1.0, 'other': 0.0}

CULVERTS = [
    dict(label='Branch_2_Brooker_St_Culvert', type='boyd_pipe', diameter=0.9,
         el0=[[305772.982, 6193988.557], [305772.378, 6193987.823]],
         el1=[[305794.592, 6193983.907], [305793.988, 6193983.173]],
         apron=3.0, enquiry_gap=10.0, manning=0.013, smoothing=0.0, losses=_LC),
    dict(label='Branch_2_Meadow_St_Culvert', type='boyd_box', width=5.4, height=0.6,
         el0=[[305886.333, 6193929.052], [305883.172, 6193922.986]],
         el1=[[305906.553, 6193910.461], [305903.393, 6193904.395]],
         apron=3.0, enquiry_gap=10.0, manning=0.013, smoothing=SMOOTH_TS, losses=_LC),
    dict(label='Branch_2_Williams_St_Culvert', type='boyd_pipe', diameter=1.2,
         el0=[[305945.955, 6193836.293], [305945.125, 6193835.387]],
         el1=[[306040.565, 6193827.573], [306039.735, 6193826.667]],
         apron=3.0, enquiry_gap=10.0, manning=0.013, smoothing=0.0, losses=_LC),
    dict(label='Branch_Towradgi_Meadow_St_Culvert', type='boyd_box', width=4.0, height=2.2,
         el0=[[305812.113, 6193591.972], [305809.390, 6193588.820]],
         el1=[[305834.913, 6193588.382], [305832.190, 6193585.230]],
         apron=3.0, enquiry_gap=10.0, manning=0.013, smoothing=SMOOTH_TS, losses=_LC),
    dict(label='Branch_5_Collins_St_Culverts', type='boyd_box', width=14.4, height=0.93,
         el0=[[306330.608, 6194817.116], [306320.768, 6194805.884]],
         el1=[[306369.483, 6194811.616], [306359.643, 6194800.384]],
         apron=3.0, enquiry_gap=10.0, manning=0.013, smoothing=SMOOTH_TS, losses=_LC),
    dict(label='Branch_5_Northern_Distributor_Culverts', type='boyd_box', width=9.09, height=0.85,
         el0=[[306956.242, 6194465.589], [306950.446, 6194457.411]],
         el1=[[307003.711, 6194446.089], [306997.916, 6194437.911]],
         apron=3.0, enquiry_gap=10.0, manning=0.013, smoothing=SMOOTH_TS, losses=_LC),
    dict(label='Branch_5_Coke_Works_Culverts', type='boyd_box', width=4.56, height=2.9,
         el0=[[307142.161, 6194181.3065], [307138.519, 6194174.394]],
         el1=[[307160.521, 6194164.8165], [307156.879, 6194157.904]],
         apron=3.1, enquiry_gap=10.0, manning=0.013, smoothing=SMOOTH_TS, losses=_LC),
    dict(label='Branch_6_Northern_Distributor_Culverts', type='boyd_box', width=3.6, height=1.2,
         el0=[[306950.758, 6193454.717], [306947.804, 6193453.283]],
         el1=[[306988.633, 6193474.217], [306985.679, 6193472.783]],
         apron=3.1, enquiry_gap=10.0, manning=0.013, smoothing=SMOOTH_TS, losses=_LC),
    dict(label='Branch_6_Railway_Culverts', type='boyd_box', width=1.0, height=3.5,
         el0=[[307139.134, 6193474.458], [307138.492, 6193473.542]],
         el1=[[307150.884, 6193469.458], [307150.242, 6193468.542]],
         apron=3.1, enquiry_gap=10.0, manning=0.013, smoothing=SMOOTH_TS, losses=_LC),
    dict(label='Branch_6_Colgong_St_Culverts', type='boyd_box', width=1.65, height=1.05,
         el0=[[307200.610, 6193476.765], [307199.140, 6193475.235]],
         el1=[[307224.610, 6193475.765], [307223.140, 6193474.235]],
         apron=3.1, enquiry_gap=10.0, manning=0.013, smoothing=SMOOTH_TS, losses=_LC),
    dict(label='Branch_3_Basin_Outlet_Culverts', type='boyd_box', width=6.0, height=0.86,
         el0=[[305629.639, 6194408.883], [305626.521, 6194400.457]],
         el1=[[305665.889, 6194347.183], [305662.771, 6194338.757]],
         apron=3.1, enquiry_gap=10.0, manning=0.013, smoothing=SMOOTH_TS, losses=_LC),
    dict(label='Branch_3_Bellambi_Rd_Culverts', type='boyd_box', width=1.65, height=1.05,
         el0=[[305777.182, 6194305.377], [305776.444, 6194304.623]],
         el1=[[305873.807, 6194303.377], [305873.069, 6194302.623]],
         apron=3.1, enquiry_gap=10.0, manning=0.013, smoothing=SMOOTH_TS, losses=_LC),
    dict(label='Branch_3_Meadow_St_Culverts', type='boyd_pipe', diameter=1.5,
         el0=[[305914.649, 6194322.375], [305913.477, 6194321.625]],
         el1=[[305950.711, 6194335.375], [305949.539, 6194334.625]],
         apron=3.1, enquiry_gap=10.0, manning=0.013, smoothing=0.0, losses=_LC),
    dict(label='Branch_3_13_Meadow_St_Culverts', type='boyd_pipe', diameter=1.5,
         el0=[[305911.280, 6194359.203], [305910.260, 6194358.017]],
         el1=[[305946.090, 6194353.573], [305945.070, 6194352.387]],
         apron=3.1, enquiry_gap=10.0, manning=0.013, smoothing=0.0, losses=_LC),
    dict(label='Branch_3_41_Angel_St_Culverts', type='boyd_box', width=10.0, height=0.35,
         el0=[[306196.779, 6194028.193], [306192.221, 6194010.807]],
         el1=[[306200.154, 6194018.693], [306195.596, 6194001.307]],
         apron=3.1, enquiry_gap=10.0, manning=0.013, smoothing=SMOOTH_TS, losses=_LC),
    dict(label='Branch_7_Carroll_St_Culverts', type='boyd_box', width=1.22, height=0.3,
         el0=[[308002.045, 6193820.163], [308001.215, 6193819.197]],
         el1=[[308021.965, 6193816.883], [308021.135, 6193815.917]],
         apron=3.1, enquiry_gap=10.0, manning=0.013, smoothing=SMOOTH_TS, losses=_LC),
    dict(label='Branch_7_Parker_Rd_Culverts', type='boyd_box', width=3.18, height=0.3,
         el0=[[308105.832, 6193803.622], [308103.648, 6193801.118]],
         el1=[[308126.782, 6193800.552], [308124.598, 6193798.048]],
         apron=3.1, enquiry_gap=10.0, manning=0.013, smoothing=SMOOTH_TS, losses=_LC),
    dict(label='Branch_7_Lake_Pde_Culverts', type='boyd_box', width=2.36, height=0.75,
         el0=[[308251.257, 6193614.658], [308248.343, 6193618.0]],
         el1=[[308232.0, 6193593.0], [308225.0, 6193596.0]],
         apron=3.1, enquiry_gap=10.0, manning=0.013, smoothing=SMOOTH_TS, losses=_LC),
    # Bridges (modelled as box culverts, pier losses, channel manning)
    dict(label='Branch_Towradgi_Princes_Hwy_Bridge', type='boyd_box', width=12.0, height=3.0,
         el0=[[306607.274, 6193707.421], [306602.635, 6193695.720]],
         el1=[[306626.205, 6193694.358], [306622.068, 6193683.138]],
         apron=0.0, enquiry_gap=10.0, manning=CHANNEL_MANNING, smoothing=SMOOTH_TS, losses=_LP),
    dict(label='Branch_Towradgi_Pioneer_Rd_Bridge', type='boyd_box', width=20.0, height=3.5,
         el0=[[307623.0, 6193610.0], [307622.0, 6193607.0]],
         el1=[[307610.0, 6193619.0], [307609.0, 6193616.0]],
         apron=0.0, enquiry_gap=10.0, manning=CHANNEL_MANNING, smoothing=SMOOTH_TS, losses=_LP),
    dict(label='Branch_Towradgi_Northern_Distributor_Bridge', type='boyd_box', width=45.0, height=6.0,
         el0=[[306985.0, 6193749.0], [306985.0, 6193736.0]],
         el1=[[306950.0, 6193745.0], [306950.0, 6193732.0]],
         apron=0.0, enquiry_gap=10.0, manning=CHANNEL_MANNING, smoothing=SMOOTH_TS, losses=_LP),
    dict(label='Branch_Towradgi_Railway_Bridge', type='boyd_box', width=20.0, height=8.0,
         el0=[[307236.0, 6193737.0], [307235.0, 6193733.0]],
         el1=[[307223.0, 6193738.0], [307222.0, 6193734.0]],
         apron=0.0, enquiry_gap=20.0, manning=CHANNEL_MANNING, smoothing=SMOOTH_TS, losses=_LP),
]


def _write_polyline(path, pts):
    with open(path, 'w') as fh:
        for x, y in pts:
            fh.write(f'{x},{y}\n')


def _manning_value(path):
    base = os.path.splitext(os.path.basename(path))[0]
    d = os.path.basename(os.path.dirname(path))
    if d == 'Creeks':
        return CHANNEL_MANNING
    if d == 'Buildings':
        return 10.0
    return MANNING_EXCEPTIONS.get(base, MANNING_DEFAULT)


def build_friction_entries():
    """Ordered [(polygon, value)] for TOML (earlier = higher priority).

    The script's Polygon_function gives LATER entries priority, with order
    Mannings -> Creeks -> Buildings. TOML is the reverse, so emit
    Buildings -> Creeks -> Mannings, then the base catch-all last.
    """
    buildings = sorted(glob.glob('Model/Buildings/*.csv'))
    creeks    = sorted(glob.glob('Model/Creeks/*.csv'))
    mannings  = sorted(glob.glob('Model/Mannings/*.csv'))
    entries = []
    for f in buildings + creeks + mannings:
        entries.append((f.replace('\\', '/'), _manning_value(f)))
    entries.append(('All', BASE_FRICTION))
    return entries


def build_rainfall_entries():
    """One [[rainfall]] per gauge: polygon = Gauge/<id>.csv, ts = Hort/<id>.tms."""
    out = []
    for g in sorted(glob.glob('Forcing/Rainfall/Gauge/*.csv')):
        gid = os.path.splitext(os.path.basename(g))[0]
        tms = f'Forcing/Rainfall/Hort/{gid}.tms'
        if os.path.exists(tms):
            out.append((tms, g.replace('\\', '/')))
    return out


def main():
    # ---- helper CSVs --------------------------------------------------------
    _write_polyline('bounding_polygon.csv',
                    [[W, S], [E, S], [E, N], [W, N]])
    _write_polyline('inlet_line.csv',
                    [[304000.0, 6194200.0], [304000.0, 6194600.0]])
    with open(INLET_FILE, 'w') as fh:
        fh.write(f'time,discharge\n0.0,{INLET_Q}\n1000000.0,{INLET_Q}\n')
    with open('east_stage.csv', 'w') as fh:
        fh.write('time,stage\n0.0,0.0\n1000000.0,0.0\n')

    os.makedirs('culvert_exchange', exist_ok=True)
    for c in CULVERTS:
        _write_polyline(f"culvert_exchange/{c['label']}_0.csv", c['el0'])
        _write_polyline(f"culvert_exchange/{c['label']}_1.csv", c['el1'])

    friction = build_friction_entries()
    rainfall = build_rainfall_entries()

    # ---- emit TOML ----------------------------------------------------------
    L = []
    A = L.append
    A('# =============================================================================')
    A('# Towradgi Creek storm — anuga_toml_run equivalent of run_small_towradgi.py')
    A('# GENERATED by make_towradgi_toml.py — edit that script, not this file.')
    A('# Run:  anuga_toml_run towradgi.toml   (needs the DEM_bridges/Model/Forcing data)')
    A('# =============================================================================')
    A('')
    A('[project]')
    A('scenario               = "Towradgi_historic_flood"')
    A('output_base_directory  = "OUTPUT/"')
    A('yieldstep              = "3*60"')
    A('starttime              = "-3*3600"   # 3 h burn-in; the rain starts at t = 0')
    A('finaltime              = "24*3600"')
    A('outputstep             = "15*60"')
    A(f'projection_information = "EPSG:{EPSG}"')
    A('flow_algorithm         = "DE_ader2"')   # run_small_towradgi.py uses DE1; see docstring
    A('compute_mode           = "unified"  # "legacy" (CPU OpenMP) or "unified" (CPU/GPU)')
    A('report_mass_conservation_statistics = true  # water balance each yieldstep')
    A('')
    A('[mesh]')
    A('bounding_polygon    = "bounding_polygon.csv"')
    A(f'default_res         = {MAX_TRI_AREA}')
    A('riverwall_csv_files = ["Model/Riverwalls/*.csv"]')
    for tag, edge in (('south', 0), ('east', 1), ('north', 2), ('west', 3)):
        A('[[mesh.boundary_tags]]')
        A(f'tag = "{tag}"')
        A(f'edges = [{edge}]')
    for poly, res in (('Model/Bdy/Catchment.csv', 100.0),
                      ('Model/Bdy/FineCatchment.csv', 36.0),
                      ('Model/Bdy/CreekBanks.csv', 8.0)):
        A('[[mesh.interior_regions]]')
        A(f'polygon    = "{poly}"')
        A(f'resolution = {res}')
    A('')
    A('[boundary_conditions]')
    for tag in ('west', 'south', 'north'):
        A('[[boundary_conditions.boundaries]]')
        A(f'tag  = "{tag}"')
        A('type = "Reflective"')
    A('[[boundary_conditions.boundaries]]')
    A('tag        = "east"')
    A('type       = "Stage"           # Transmissive_n_momentum_zero_t_momentum_set_stage')
    A('file       = "east_stage.csv"  # constant stage 0.0')
    A('start_time = 0.0')
    A('')
    A('[initial_conditions]')
    A('[[initial_conditions.elevation]]')
    A('polygon = "All"')
    A('value   = "DEM_bridges/towradgi.npy"   # x,y,z point cloud, nearest-neighbour')
    A('[[initial_conditions.stage]]')
    A('polygon = "All"')
    A('value   = 0.0')
    A('[[initial_conditions.xmomentum]]')
    A('polygon = "All"')
    A('value   = 0.0')
    A('[[initial_conditions.ymomentum]]')
    A('polygon = "All"')
    A('value   = 0.0')
    A(f'# Friction: {len(friction)} entries (Buildings 10.0, Creeks {CHANNEL_MANNING}, '
      f'Mannings mostly {MANNING_DEFAULT}), base {BASE_FRICTION} last.')
    for poly, val in friction:
        A('[[initial_conditions.friction]]')
        A(f'polygon = "{poly}"')
        A(f'value   = {val}')
    A('')
    A(f'# Inlet — constant {INLET_Q:g} m^3/s creek inflow (script creek_inlet uses 20)')
    A('[[inlets]]')
    A('name            = "Creek_Inlet_West"')
    A('line_file       = "inlet_line.csv"')
    A(f'timeseries_file = "{INLET_FILE}"')
    A('start_time      = 0.0')
    A('')
    A(f'# Culverts & bridges — {len(CULVERTS)} Boyd operators (exchange lines in culvert_exchange/).')
    for c in CULVERTS:
        A('[[culverts]]')
        A(f'type            = "{c["type"]}"')
        A(f'label           = "{c["label"]}"')
        if c['type'] == 'boyd_pipe':
            A(f'diameter        = {c["diameter"]}')
        else:
            A(f'width           = {c["width"]}')
            A(f'height          = {c["height"]}')
        A(f'exchange_line_0 = "culvert_exchange/{c["label"]}_0.csv"')
        A(f'exchange_line_1 = "culvert_exchange/{c["label"]}_1.csv"')
        A(f'apron           = {c["apron"]}')
        A(f'enquiry_gap     = {c["enquiry_gap"]}')
        A(f'manning         = {c["manning"]}')
        A(f'smoothing_timescale = {c["smoothing"]}')
        lo = c['losses']
        A('losses          = {' + ', '.join(f'{k} = {v}' for k, v in lo.items()) + '}')
    A('')
    A(f'# Rainfall — {len(rainfall)} gauges: polygon + ANUGA .tms rate timeseries.')
    A('# multiplier = 1e-3 is the Rate_operator unit factor (matches the script).')
    for tms, gauge in rainfall:
        A('[[rainfall]]')
        A(f'timeseries_file = "{tms}"')
        A(f'polygon         = "{gauge}"')
        A('start_time      = 0.0')
        A('multiplier      = 1.0e-3')
    L.append('')

    with open('towradgi.toml', 'w') as fh:
        fh.write('\n'.join(L))

    print(f'Wrote towradgi.toml: {len(friction)} friction, {len(CULVERTS)} culverts, '
          f'{len(rainfall)} rainfall gauges, + helper CSVs and '
          f'{2*len(CULVERTS)} exchange-line files.')


if __name__ == '__main__':
    main()
