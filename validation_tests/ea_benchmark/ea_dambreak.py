"""
EA 2D benchmark Tests 6A and 6B: dam break past a building.

The IMPACT isolated-building case (Soares-Frazao and Zech 2002): a channel
with a dam that has an opening, deep water behind it and shallow water in
front, and a rectangular building set at an angle downstream. All
boundaries are vertical walls. Test 6A is the laboratory flume; Test 6B
scales every length by 20.

The geometry is taken from the test's DEM, where walls stand 10 m above a
floor (with a shallow fillet along each side wall, kept here): the channel
spans the DEM's open cells, the dam's two blocks are cut out of the
bounding polygon and the building is a hole. Both DEMs share one grid, the
6A cell size being 1/20 of 6B's, and the same heights: in 6B units the
DEM's partly covered cells put the dam at -16 < x < 0 with its opening at
-10 < y < 10, and the building is 16 m x 8 m.

run() writes <name>_gauges.csv (level, depth and speed at G1-G6),
<name>_volume.txt (volumes and the largest speed in the domain),
<name>_geometry.json (for the report's set-up figure) and <name>.sww.

EA_<TEST>_MAX_AREA (for example EA_6B_MAX_AREA) sets the triangle area for
a mesh convergence study; the outputs are then named <test>_a<area>_*.
EA_<TEST>_BUILDING picks how the building is modelled: 'hole' (default, a
reflective hole in the mesh), 'block' (the fitted rectangle as a breakline
with the floor raised 10 m inside it) or 'dem' (the DEM's own cells,
sampled at the centroids, as the grid-based packages saw it); the outputs
gain a _<building> suffix. run(beta=...) overrides the slope limiter
(beta_w, beta_uh, beta_vh) for the limiter study of Test 6A; the outputs
gain a _b<beta> suffix. Finer meshes are best run in 'unified' compute
mode (ANUGA_DEFAULT_COMPUTE_MODE=unified) on a GPU build.
"""

import json
import os

import numpy as np
import anuga

import ea_common as ea

TESTS = {
    # scale: lengths relative to 6B; max_area: the specified grid squared
    '6A': dict(dem='Test6ADEM.asc', points='Test6Aoutput.csv', scale=0.05,
               manning=0.01, h_up=0.4, h_down=0.02, max_area=0.01,
               final_time=120.0, yield_step=0.1, output_step=1.0),
    '6B': dict(dem='Test6BDEM.asc', points='Test6Boutput.csv', scale=1.0,
               manning=0.05, h_up=8.0, h_down=0.4, max_area=4.0,
               final_time=1800.0, yield_step=1.0, output_step=60.0),
}
WALL = 5.0                 # DEM cells above this are walls (heights are not scaled)
LOW = 0.2                  # DEM cells above this are wall, dam or building


def run_name(test, beta=None):
    """The output stem for this test, the EA_<TEST>_* settings and beta."""
    p = TESTS[test]
    max_area = float(os.environ.get('EA_%s_MAX_AREA' % test, p['max_area']))
    building = os.environ.get('EA_%s_BUILDING' % test, 'hole')
    name = 'test%s' % test.lower()
    if max_area != p['max_area']:
        name += '_a%g' % max_area
    if building != 'hole':
        name += '_' + building
    if beta is not None:
        name += '_b%g' % beta
    return name, max_area, building


def geometry(test):
    """Bounding polygon (dam cut out) and building rectangle, from the DEM."""
    p = TESTS[test]
    s = p['scale']
    x_dem, y_dem, z_dem = ea.read_asc(ea.data_path('Test6', p['dem']))

    # Channel extent: the DEM's open rows and columns, to cell edges.
    open_cells = z_dem < WALL
    cs = x_dem[1] - x_dem[0]
    rows = np.where(open_cells.any(axis=1))[0]
    cols = np.where(open_cells.any(axis=0))[0]
    x0, x1 = x_dem[cols[0]] - cs / 2, x_dem[cols[-1]] + cs / 2
    y0, y1 = y_dem[rows[0]] - cs / 2, y_dem[rows[-1]] + cs / 2

    # Dam blocks at -16 < x < 0, leaving the opening -10 < y < 10 (6B units).
    dam_x0, dam_x1, gap = -16.0 * s, 0.0, 10.0 * s
    bounding_polygon = [[x0, y0], [dam_x0, y0], [dam_x0, -gap], [dam_x1, -gap],
                        [dam_x1, y0], [x1, y0], [x1, y1], [dam_x1, y1],
                        [dam_x1, gap], [dam_x0, gap], [dam_x0, y1], [x0, y1]]

    # The building as a rectangle fitted to its DEM cells.
    X, Y = np.meshgrid(x_dem, y_dem)
    in_building = (X > 10.0 * s) & (X < 200.0 * s) & (Y > y0) & (Y < y1)
    # Skip the low wall fillets: over the length of wall they would swamp
    # the covariance and turn the fit towards the channel axis.
    w = np.where(in_building & (z_dem > LOW), np.clip(z_dem / 10.0, 0.0, 1.0), 0.0)
    xc, yc = np.sum(w * X) / w.sum(), np.sum(w * Y) / w.sum()
    dx, dy = (X - xc)[w > 0], (Y - yc)[w > 0]
    cov = np.cov(np.vstack((dx, dy)), aweights=w[w > 0])
    evals, evecs = np.linalg.eigh(cov)
    major, minor = evecs[:, 1], evecs[:, 0]
    a, b = 8.0 * s, 4.0 * s                              # half sides
    corners = [(xc, yc) + sa * a * major + sb * b * minor
               for sa, sb in [(-1, -1), (1, -1), (1, 1), (-1, 1)]]
    # eigh differs in the last digits between OpenBLAS CPU kernels, and the mesh
    # generator turns that into a different mesh (55812 or 55601 triangles in 6A)
    # and gauge levels 6 cm apart. Round so every machine meshes the same corners.
    building = [[round(float(c[0]), 6), round(float(c[1]), 6)] for c in corners]
    return bounding_polygon, building, (x_dem, y_dem, z_dem)


def run(test, alg, verbose=False, beta=None, final_time=None):
    p = TESTS[test]
    s = p['scale']
    name, max_area, building_mode = run_name(test, beta)
    points = ea.read_points(ea.data_path('Test6', p['points']))
    bounding_polygon, building, (x_dem, y_dem, z_dem) = geometry(test)
    cs = x_dem[1] - x_dem[0]

    with open(name + '_geometry.json', 'w') as f:
        json.dump({'bounding_polygon': bounding_polygon, 'building': building}, f)

    domain = anuga.create_domain_from_regions(
        bounding_polygon,
        boundary_tags={'wall': list(range(len(bounding_polygon)))},
        maximum_triangle_area=max_area,
        interior_holes=[building] if building_mode == 'hole' else None,
        breaklines=[building + building[:1]] if building_mode == 'block' else None,
        use_cache=False, verbose=False)

    domain.set_name(name)
    domain.set_flow_algorithm(alg)
    if beta is not None:
        domain.beta_w = domain.beta_uh = domain.beta_vh = beta

    # Floor only: drop the walls and the partly covered cells beside them.
    z_floor = np.where(z_dem > LOW, 0.0, z_dem)
    domain.set_quantity('elevation',
                        ea.absolute(domain, ea.grid_function(x_dem, y_dem, z_floor)),
                        location='centroids')
    if building_mode != 'hole':
        xc, yc = domain.centroid_coordinates.T + [[domain.geo_reference.xllcorner],
                                                   [domain.geo_reference.yllcorner]]
        if building_mode == 'block':
            inside = anuga.inside_polygon(np.column_stack((xc, yc)), building)
        else:
            # The DEM cell under each centroid, as a grid model would see it.
            i = np.clip(np.round((yc - y_dem[0]) / cs).astype(int), 0, len(y_dem) - 1)
            j = np.clip(np.round((xc - x_dem[0]) / cs).astype(int), 0, len(x_dem) - 1)
            inside = np.where((xc > 10.0 * s) & (xc < 200.0 * s) & (z_dem[i, j] > WALL))[0]
        elev = domain.quantities['elevation'].centroid_values
        elev[inside] = 10.0
        print('building: %d centroids raised to 10 m' % len(inside))
    domain.set_quantity('friction', p['manning'], location='centroids')

    def initial_stage(x, y):
        return np.where(x < 0.0, p['h_up'], p['h_down'])

    domain.set_quantity('stage', ea.absolute(domain, initial_stage),
                        location='centroids')
    if building_mode != 'hole':   # dry on top of the building
        stage = domain.quantities['stage'].centroid_values
        stage[:] = np.maximum(stage, domain.quantities['elevation'].centroid_values)

    Br = anuga.Reflective_boundary(domain)
    domain.set_boundary({'wall': Br, 'interior': Br} if building_mode == 'hole'
                        else {'wall': Br})

    from anuga.validation_utilities import save_parameters_tex
    save_parameters_tex(domain)

    gauges = ea.Gauges(domain, points)
    volume_0 = ea.domain_volume(domain)
    max_speed = 0.0
    wet_depth = 0.05 * s
    for t in domain.evolve(yieldstep=p['yield_step'], outputstep=p['output_step'],
                           finaltime=final_time or p['final_time']):
        gauges.sample()
        h = (domain.quantities['stage'].centroid_values
             - domain.quantities['elevation'].centroid_values)
        wet = h > wet_depth
        if wet.any():
            u = np.hypot(domain.quantities['xmomentum'].centroid_values[wet],
                         domain.quantities['ymomentum'].centroid_values[wet]) / h[wet]
            max_speed = max(max_speed, float(u.max()))
        if verbose and abs(t % p['output_step']) < 1.0e-6:
            print(domain.timestepping_statistics())

    gauges.save(name + '_gauges.csv')
    with open(name + '_volume.txt', 'w') as f:
        f.write('initial %.6f\n' % volume_0)
        f.write('final %.6f\n' % ea.domain_volume(domain))
        f.write('max_speed %.6f\n' % max_speed)
        f.write('triangles %d\n' % domain.number_of_triangles)
    print('compute mode %s' % domain.compute_mode)


def plot_setup(test, filename):
    """The set-up figure: the whole channel, and the dam and building close up."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    p = TESTS[test]
    s = p['scale']
    with open('test%s_geometry.json' % test.lower()) as f:
        geo = json.load(f)
    points = ea.read_points(ea.data_path('Test6', p['points']))
    xs = [v[0] for v in geo['bounding_polygon']]
    ys = [v[1] for v in geo['bounding_polygon']]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 4.8), gridspec_kw={'height_ratios': [1, 3.2]})
    for ax in (ax1, ax2):
        ax.add_patch(Polygon(geo['bounding_polygon'], closed=True, fc='#dbe9f6', ec='k', lw=1))
        ax.add_patch(Polygon(geo['building'], closed=True, fc='0.4', ec='k', lw=1))
        ax.set_aspect('equal')
    ax1.set_xlim(min(xs) - 20 * s, max(xs) + 20 * s)
    ax1.set_ylim(min(ys) - 8 * s, max(ys) + 8 * s)
    ax1.set_title('The closed channel, %.4g m x %.3g m (dam at x = %g to 0 m)'
                  % (max(xs) - min(xs), max(ys) - min(ys), -16 * s), fontsize=10)
    ax2.set_xlim(-160 * s, 200 * s)
    ax2.set_ylim(-42 * s, 42 * s)
    ax2.fill_between([min(xs), -16 * s], min(ys), max(ys), color='C0', alpha=0.25,
                     label='%g m deep at t = 0' % p['h_up'])
    for n, (x, y) in points.items():
        ax2.plot(x, y, 'o', color='C3', ms=5)
        ax2.text(x + 3 * s, y + 3 * s, n, color='C3')
    ax2.text(75 * s, -30 * s, '%g m x %g m building' % (16 * s, 8 * s), ha='center', fontsize=8)
    ax2.set_xlabel('x (m)')
    ax2.set_ylabel('y (m)')
    ax2.legend(fontsize=8, loc='lower left')
    ax2.set_title('Reservoir, dam with its %g m opening, building and gauges' % (20 * s),
                  fontsize=10)
    fig.tight_layout()
    fig.savefig(filename, dpi=120)
    plt.close(fig)
