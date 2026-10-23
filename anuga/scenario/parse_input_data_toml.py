#!/usr/bin/python
"""
Read a TOML configuration file for an ANUGA scenario.

Drop-in replacement for the Excel-based AnugaXls / ProjectData interface in
parse_input_data.py.  Exposes exactly the same public attributes so that
prepare_data.py and all setup_*.py scripts work unchanged.

Requires Python 3.11+ (tomllib in stdlib) or the 'tomli' package for older
Python versions (pip install tomli).
"""

import ast
import operator
import os
import warnings
import glob as glob_module

from anuga.utilities import spatialInputUtil as su

# Must track Domain.set_flow_algorithm's accepted list.
_VALID_FLOW_ALGORITHMS = ('DE0', 'DE1', 'DE2', 'DE0_7', 'DE1_7', 'DE_ader2')

# ---------------------------------------------------------------------------
# Numeric values with optional arithmetic
# ---------------------------------------------------------------------------
# TOML has no arithmetic, so `finaltime = 5*60` is a parse error. As a
# convenience, numeric fields also accept a *quoted* arithmetic expression,
# e.g. finaltime = "5*60" -> 300.0. Only numeric literals and + - * / // % **
# with parentheses are allowed (no names, calls or attribute access), so it is
# safe to evaluate untrusted config.

_ARITH_BINOPS = {
    ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
    ast.Div: operator.truediv, ast.FloorDiv: operator.floordiv,
    ast.Mod: operator.mod, ast.Pow: operator.pow,
}
_ARITH_UNARYOPS = {ast.UAdd: operator.pos, ast.USub: operator.neg}


def _eval_arith(node):
    """Recursively evaluate a whitelisted arithmetic AST node."""
    if isinstance(node, ast.Expression):
        return _eval_arith(node.body)
    if isinstance(node, ast.Constant):
        # bool is an int subclass — reject so "true"-ish values don't sneak in.
        if isinstance(node.value, bool) or not isinstance(node.value, (int, float)):
            raise ValueError('only numbers are allowed')
        return node.value
    if isinstance(node, ast.BinOp) and type(node.op) in _ARITH_BINOPS:
        return _ARITH_BINOPS[type(node.op)](
            _eval_arith(node.left), _eval_arith(node.right))
    if isinstance(node, ast.UnaryOp) and type(node.op) in _ARITH_UNARYOPS:
        return _ARITH_UNARYOPS[type(node.op)](_eval_arith(node.operand))
    raise ValueError('unsupported expression')


def _num(value):
    """Coerce a TOML value to float, allowing a quoted arithmetic string.

    Numbers pass straight through; a string is evaluated as a simple arithmetic
    expression (``"5*60"`` -> ``300.0``). Raises ``ValueError``/``TypeError`` on
    anything else so callers can surface it as a configuration error.
    """
    if isinstance(value, bool):
        raise ValueError('expected a number, got a boolean')
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, str):
        try:
            tree = ast.parse(value, mode='eval')
        except SyntaxError as exc:
            raise ValueError('invalid arithmetic expression: %r' % value) from exc
        return float(_eval_arith(tree))
    raise TypeError('expected a number or arithmetic string, got %r' % (value,))

# ---------------------------------------------------------------------------
# Validation helper
# ---------------------------------------------------------------------------

class _Validator:
    """Accumulates TOML configuration errors for batch reporting.

    Each check method silently records an error entry rather than raising
    immediately, so that all problems in a config file are reported at once.
    Call :meth:`raise_if_errors` at the end of parsing to raise a single
    ``ValueError`` listing every problem found.
    """

    def __init__(self):
        self.errors = []

    def require(self, mapping, key, section):
        """Return ``mapping[key]``, recording an error if the key is absent."""
        if key not in mapping:
            self.errors.append(
                f'[{section}]: required field {key!r} is missing')
            return None
        return mapping[key]

    def to_float(self, mapping, key, section):
        """Return ``mapping[key]`` as a float, recording errors for absent or
        non-numeric values. Accepts a quoted arithmetic string (e.g. "5*60")."""
        val = self.require(mapping, key, section)
        if val is None:
            return None
        try:
            return _num(val)
        except (TypeError, ValueError):
            self.errors.append(
                f'[{section}] {key!r}: expected a number or arithmetic '
                f'string, got {val!r}')
            return None

    def positive(self, val, name, section):
        """Record an error if *val* is not ``> 0``."""
        if val is not None and val <= 0:
            self.errors.append(
                f'[{section}] {name!r}: must be > 0, got {val!r}')

    def non_negative(self, val, name, section):
        """Record an error if *val* is negative."""
        if val is not None and val < 0:
            self.errors.append(
                f'[{section}] {name!r}: must be >= 0, got {val!r}')

    def in_range(self, val, lo, hi, name, section):
        """Record an error if *val* is outside ``[lo, hi]``."""
        if val is not None and not (lo <= val <= hi):
            self.errors.append(
                f'[{section}] {name!r}: must be in [{lo}, {hi}], got {val!r}')

    def one_of(self, val, choices, name, section):
        """Record an error if *val* is not in *choices*."""
        if val is not None and val not in choices:
            self.errors.append(
                f'[{section}] {name!r}: must be one of {list(choices)}, '
                f'got {val!r}')

    def integer_multiple(self, val, base, name, section):
        """Record an error if *val* is not an integer multiple of *base*."""
        if val is not None and base is not None and base > 0:
            ratio = val / base
            if abs(ratio - round(ratio)) > 1e-6:
                self.errors.append(
                    f'[{section}] {name!r}: must be an integer multiple of '
                    f'yieldstep ({base}), got {val!r}')

    def raise_if_errors(self, filename):
        """Raise ``ValueError`` listing all accumulated errors, or do nothing."""
        if self.errors:
            bullets = '\n'.join(f'  • {e}' for e in self.errors)
            raise ValueError(
                f'Configuration errors in {filename!r}:\n{bullets}')


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _load_toml(filename):
    """Load a TOML file using tomllib (3.11+) or the tomli back-port."""
    try:
        import tomllib
    except ImportError:
        try:
            import tomli as tomllib
        except ImportError:
            raise ImportError(
                'TOML support requires Python 3.11+ or the tomli package '
                '(pip install tomli)')
    with open(filename, 'rb') as f:
        return tomllib.load(f)


def _normpath(path):
    """Normalise *path* to the OS-native directory separator.

    TOML files conventionally use forward slashes (``/``) for portability.
    On Windows this function converts them to backslashes so that downstream
    ANUGA code receives native paths.  On Linux/macOS the path is unchanged.

    The special quantity-polygon keywords ``'All'``, ``'None'``, and
    ``'Extent'`` are returned unchanged as they are not filesystem paths.
    """
    if path in ('All', 'None', 'Extent'):
        return path
    return os.path.normpath(path)


def _glob_files(patterns):
    """Expand a list of glob patterns into a flat, sorted list of file paths.

    Patterns may use either forward or backward slashes; both are normalised
    to the OS-native separator before globbing and in the returned paths.

    Raises FileNotFoundError if any pattern matches nothing.
    """
    files = []
    for pattern in patterns:
        native = os.path.normpath(pattern)
        matches = sorted(glob_module.glob(native))
        if not matches:
            raise FileNotFoundError(
                f'No files matched pattern: {pattern!r}')
        files.extend(os.path.normpath(m) for m in matches)
    return files


def _parse_losses(raw, sec, _v):
    """Normalise a culvert/weir ``losses`` value.

    Accepts a scalar, a list (Boyd operators sum the entries), or a dict of
    named components (e.g. ``{inlet, outlet, bend, grate, pier, other}``) — all
    forms that Boyd_box_operator / Boyd_pipe_operator accept directly. Every
    value is validated non-negative.
    """
    if isinstance(raw, dict):
        losses = {str(k): float(v) for k, v in raw.items()}
        for k, v in losses.items():
            _v.non_negative(v, f'losses[{k!r}]', sec)
    elif isinstance(raw, (list, tuple)):
        losses = [float(v) for v in raw]
        for j, v in enumerate(losses):
            _v.non_negative(v, f'losses[{j}]', sec)
    else:
        losses = float(raw)
        _v.non_negative(losses, 'losses', sec)
    return losses


def _parse_quantity_entries(entries, print_info):
    """Convert a list of TOML table entries for one quantity into the
    (data, clip_range) tuple that composite_quantity_setting_function expects.

    Each entry must have 'polygon' and 'value'.  clip_min / clip_max default
    to ±inf (no clipping).

    Handles the two-file wildcard polygon case: if 'polygon' is a glob pattern
    that matches exactly 2 line files, they are joined into a closed polygon
    (mirrors the Excel interface behaviour).
    """
    data = []
    clip_range = []

    for entry in entries:
        polygon  = _normpath(str(entry['polygon']))
        value    = entry['value']
        if isinstance(value, str):
            value = _normpath(value)
        clip_min = entry.get('clip_min', float('-inf'))
        clip_max = entry.get('clip_max', float('inf'))

        # Wildcard polygon: join 2 matching line files into a closed polygon
        if (isinstance(polygon, str) and
                polygon not in ('All', 'None', 'Extent') and
                ('*' in polygon or '?' in polygon)):
            matched = sorted(glob_module.glob(polygon))
            if len(matched) == 2:
                print_info.append(
                    f'Combining 2 files into polygon: {matched}')
                l0 = su.read_polygon(matched[0])
                l1 = su.read_polygon(matched[1])
                fake_bl = {matched[0]: l0, matched[1]: l1}
                prefix = polygon.split('*')[0]
                polygon = su.polygon_from_matching_breaklines(prefix, fake_bl)
            elif len(matched) > 2:
                raise ValueError(
                    f'Polygon pattern {polygon!r} matched {len(matched)} files; '
                    f'at most 2 are supported (they are joined into a polygon)')

        data.append([polygon, value])
        clip_range.append([clip_min, clip_max])

    return data, clip_range


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class ProjectDataTOML:
    """
    Parses a TOML scenario configuration file and exposes the same public
    attributes as the Excel-based ProjectData class so that prepare_data.py
    and the setup_*.py scripts work without modification.
    """

    def __init__(self, filename):
        self.config_filename = filename
        self.print_info = [
            '---------------------',
            'PARSING CONFIG FILE (TOML)',
            '--------------------',
            '',
        ]

        cfg = _load_toml(filename)
        _v = _Validator()

        self._parse_project(cfg.get('project', {}), _v)
        self._parse_mesh(cfg.get('mesh', {}), _v)
        # Initial conditions must be parsed before bridges / pumping stations
        # because those prepend entries to elevation_data / breakline_files.
        self._parse_initial_conditions(
            cfg.get('initial_conditions', {}),
            cfg.get('initial_condition_additions', {}))
        self._parse_boundary_conditions(cfg.get('boundary_conditions', {}))
        self._parse_inlets(cfg.get('inlets', []))
        self._parse_rainfall(cfg.get('rainfall', []))
        self._parse_bridges(cfg.get('bridges', []))
        self._parse_pumping_stations(cfg.get('pumping_stations', []))
        self._parse_culverts(cfg.get('culverts', []), _v)
        self._parse_weirs(cfg.get('weirs', []), _v)
        self._parse_erosion(cfg.get('erosion', []), _v)
        self._parse_sediment(cfg.get('sediment', None), _v)

        _v.raise_if_errors(filename)

    # -----------------------------------------------------------------------
    # Project settings
    # -----------------------------------------------------------------------

    def _parse_project(self, p, _v):
        raw_scenario = _v.require(p, 'scenario', 'project')
        self.scenario = str(raw_scenario) if raw_scenario is not None else ''

        raw_outdir = _v.require(p, 'output_base_directory', 'project')
        self.output_basedir = str(raw_outdir) if raw_outdir is not None else ''

        yieldstep = _v.to_float(p, 'yieldstep', 'project')
        _v.positive(yieldstep, 'yieldstep', 'project')
        self.yieldstep = yieldstep if yieldstep is not None else 1.0

        finaltime = _v.to_float(p, 'finaltime', 'project')
        _v.positive(finaltime, 'finaltime', 'project')
        self.finaltime = finaltime if finaltime is not None else float('inf')

        # Model time at the start of the run. A negative value runs a burn-in
        # before the forcing data (which start at t=0) begin.
        self.starttime = _num(p.get('starttime', 0.0))
        if self.starttime >= self.finaltime:
            _v.errors.append(
                '[project] starttime must be < finaltime'
                f' (got {self.starttime}, finaltime={self.finaltime})')

        raw_proj = _v.require(p, 'projection_information', 'project')
        if raw_proj is None:
            self.projection_information = None
        elif isinstance(raw_proj, float):
            self.projection_information = int(raw_proj)
        elif isinstance(raw_proj, int):
            self.projection_information = raw_proj
        else:
            self.projection_information = str(raw_proj)

        raw_alg = _v.require(p, 'flow_algorithm', 'project')
        _v.one_of(raw_alg, _VALID_FLOW_ALGORITHMS, 'flow_algorithm', 'project')
        self.flow_algorithm = str(raw_alg) if raw_alg is not None else 'DE0'

        self.output_tif_cellsize = _num(p.get('output_tif_cellsize', 50.0))
        _v.positive(self.output_tif_cellsize, 'output_tif_cellsize', 'project')

        otbp = p.get('output_tif_bounding_polygon', '')
        self.output_tif_bounding_polygon = str(otbp) if otbp else None

        self.max_quantity_update_frequency = int(
            p.get('max_quantity_update_frequency', 1))
        _v.positive(self.max_quantity_update_frequency,
                    'max_quantity_update_frequency', 'project')

        self.max_quantity_collection_start_time = _num(
            p.get('max_quantity_collection_starttime', 0.0))
        if finaltime is not None and \
                self.max_quantity_collection_start_time >= self.finaltime:
            _v.errors.append(
                '[project] max_quantity_collection_starttime must be < finaltime'
                f' (got {self.max_quantity_collection_start_time}, finaltime={self.finaltime})')

        self.store_vertices_uniquely        = bool(p.get('store_vertices_uniquely', False))
        # Track whether the user set this explicitly: with [[erosion]] present we
        # default it to true, but still warn if they deliberately chose false.
        self._store_elevation_explicit      = 'store_elevation_every_timestep' in p
        self.store_elevation_every_timestep = bool(p.get('store_elevation_every_timestep', False))
        self.spatial_text_output_dir        = str(p.get('spatial_text_output_dir', 'SPATIAL_TEXT'))

        self.report_mass_conservation_statistics      = bool(p.get('report_mass_conservation_statistics', False))
        self.report_peak_velocity_statistics          = bool(p.get('report_peak_velocity_statistics', False))
        self.report_smallest_edge_timestep_statistics = bool(p.get('report_smallest_edge_timestep_statistics', False))
        self.report_operator_statistics               = bool(p.get('report_operator_statistics', False))

        # Number of OpenMP threads (1 = single-threaded; None = use OMP_NUM_THREADS env var)
        raw_omp = p.get('omp_num_threads', None)
        if raw_omp is not None:
            self.omp_num_threads = int(raw_omp)
            _v.positive(self.omp_num_threads, 'omp_num_threads', 'project')
        else:
            self.omp_num_threads = None

        # Compute path:
        #   "legacy"  — the historical CPU OpenMP kernels (multiprocessor_mode 1)
        #   "unified" — the shared CPU/GPU (mode-2) kernels
        # The older integer `multiprocessor_mode` (1/2) is still accepted for
        # backward compatibility and mapped to the corresponding string.
        _INT_TO_MODE = {1: 'legacy', 2: 'unified'}
        _MODE_TO_INT = {'legacy': 1, 'unified': 2}
        if 'compute_mode' in p:
            self.compute_mode = str(p['compute_mode']).lower()
            _v.one_of(self.compute_mode, ('legacy', 'unified'),
                      'compute_mode', 'project')
        elif 'multiprocessor_mode' in p:
            mpm = int(p['multiprocessor_mode'])
            _v.one_of(mpm, (1, 2), 'multiprocessor_mode', 'project')
            self.compute_mode = _INT_TO_MODE.get(mpm, 'legacy')
        else:
            self.compute_mode = 'legacy'
        # Keep the integer form for any legacy consumers.
        self.multiprocessor_mode = _MODE_TO_INT.get(self.compute_mode, 1)

        # SWW output interval [seconds]; None means write every yieldstep.
        # Must be an integer multiple of yieldstep.
        raw_os = p.get('outputstep', None)
        if raw_os is not None:
            self.outputstep = _num(raw_os)
            _v.positive(self.outputstep, 'outputstep', 'project')
            _v.integer_multiple(self.outputstep, yieldstep, 'outputstep', 'project')
        else:
            self.outputstep = None

    # -----------------------------------------------------------------------
    # Mesh
    # -----------------------------------------------------------------------

    def _parse_mesh(self, m, _v):
        self.use_existing_mesh_pickle = bool(m.get('use_existing_mesh_pickle', False))

        raw_bp = _v.require(m, 'bounding_polygon', 'mesh')
        self.bounding_polygon_and_tags_file = (
            _normpath(str(raw_bp)) if raw_bp is not None else '')

        default_res = _v.to_float(m, 'default_res', 'mesh')
        _v.positive(default_res, 'default_res', 'mesh')
        self.default_res = default_res if default_res is not None else 1.0

        self.interior_regions_data = []
        for ir in m.get('interior_regions', []):
            poly = _normpath(str(ir['polygon']))
            res  = float(ir['resolution'])
            _v.positive(res, f'interior_regions resolution ({poly!r})', 'mesh')
            self.interior_regions_data.append([poly, res])

        # Interior holes — polygons excluded from the mesh entirely, leaving a
        # void rather than a refined region. Used for anything water should
        # neither enter nor flow through: building footprints, tank pads, solid
        # structures.
        #
        # create_pmesh_from_regions() has always accepted interior_holes and
        # hole_tags; there was simply no way to express them here.
        #
        # Distinct from interior_regions, which keeps the triangles and only
        # changes their size — so no resolution is taken. The optional 'tag'
        # maps to hole_tags, letting boundary conditions address the hole's
        # edges by name; omit it and the mesh generator uses its own default.
        self.interior_holes_data = []
        for ih in m.get('interior_holes', []):
            poly = _normpath(str(ih['polygon']))
            tag = ih.get('tag')
            self.interior_holes_data.append(
                [poly, str(tag) if tag is not None else None])

        # Explicit boundary tags — required when bounding_polygon is a CSV file,
        # ignored when it is a shapefile (tags come from the shapefile attributes)
        raw_tags = m.get('boundary_tags', [])
        if raw_tags:
            self.bounding_polygon_explicit_tags = [
                {'tag': str(t['tag']), 'edges': list(t['edges'])}
                for t in raw_tags
            ]
        else:
            self.bounding_polygon_explicit_tags = None

        self.breakline_files     = _glob_files(m.get('breakline_files', []))
        self.riverwall_csv_files = _glob_files(m.get('riverwall_csv_files', []))

        threshold = m.get('breakline_intersection_threshold', 'ignore')
        if isinstance(threshold, str):
            self.break_line_intersect_point_movement_threshold = threshold
        else:
            self.break_line_intersect_point_movement_threshold = float(threshold)

        region_file = m.get('region_areas_file', '')
        self.pt_areas = _normpath(str(region_file)) if region_file else None

        if self.pt_areas is not None:
            areas_type = str(m.get('region_areas_type', 'area'))
            if areas_type == 'length':
                self.region_resolutions_from_length = True
            elif areas_type == 'area':
                self.region_resolutions_from_length = False
            else:
                raise ValueError(
                    f'region_areas_type must be "area" or "length", '
                    f'got {areas_type!r}')

        # interior_regions and region_areas_file are two *alternative* ways to
        # specify mesh resolution, so they are mutually exclusive. Breaklines
        # and riverwalls are NOT — setup_mesh passes them to
        # create_pmesh_from_regions alongside interior_regions (e.g. a
        # catchment-refined model that also has riverwalls, as in the Towradgi
        # study), so they may be combined freely with either.
        if bool(self.interior_regions_data) and self.pt_areas is not None:
            raise ValueError(
                'Cannot specify both interior_regions and region_areas_file '
                '(they are alternative resolution-control methods)')

    # -----------------------------------------------------------------------
    # Initial conditions
    # -----------------------------------------------------------------------

    def _parse_initial_conditions(self, ic, ic_add):
        quantities = ['elevation', 'friction', 'stage', 'xmomentum', 'ymomentum']

        for q in quantities:
            entries  = ic.get(q, [])
            data, clip_range = _parse_quantity_entries(entries, self.print_info)

            # Optional per-quantity spatial averaging grid spacing [m]
            mean = ic.get(f'{q}_spatial_average', None)
            if mean is not None:
                mean = float(mean)

            setattr(self, f'{q}_data',       data)
            setattr(self, f'{q}_clip_range', clip_range)
            setattr(self, f'{q}_mean',       mean)

            # Additions — clip_range is not used downstream, so discard it
            add_entries = ic_add.get(q, [])
            additions, _ = _parse_quantity_entries(add_entries, self.print_info)
            setattr(self, f'{q}_additions', additions)

    # -----------------------------------------------------------------------
    # Boundary conditions
    # -----------------------------------------------------------------------

    def _parse_boundary_conditions(self, bc):
        self.boundary_tags_attribute_name = str(
            bc.get('boundary_tags_attribute_name', 'bnd_tag'))

        # Also accept the explicit ANUGA boundary class names (with or without
        # the trailing '_boundary') as aliases for the short type keywords.
        _BC_ALIASES = {
            'Transmissive_n_momentum_zero_t_momentum_set_stage':          'Stage',
            'Transmissive_n_momentum_zero_t_momentum_set_stage_boundary': 'Stage',
            'Flather_external_stage_zero_velocity':                        'Flather_Stage',
            'Flather_external_stage_zero_velocity_boundary':               'Flather_Stage',
        }
        self.boundary_data = []
        for b in bc.get('boundaries', []):
            tag   = str(b['tag'])
            raw   = str(b['type'])
            btype = _BC_ALIASES.get(raw, raw)
            if btype in ('Stage', 'Flather_Stage'):
                fpath = _normpath(str(b['file']))
                if not os.path.exists(fpath):
                    raise FileNotFoundError(
                        f"Boundary '{tag}': timeseries file not found: {fpath!r}")
                row = [tag, btype, fpath, float(b.get('start_time', 0.0))]
            elif btype == 'Reflective':
                row = [tag, btype]
            else:
                raise ValueError(
                    f'Unknown boundary type {raw!r} for tag {tag!r}. Valid types: '
                    f'"Reflective", "Stage" (alias '
                    f'Transmissive_n_momentum_zero_t_momentum_set_stage), '
                    f'"Flather_Stage" (alias Flather_external_stage_zero_velocity)')
            self.boundary_data.append(row)

    # -----------------------------------------------------------------------
    # Inlets
    # -----------------------------------------------------------------------

    def _parse_inlets(self, inlets):
        self.inlet_data = []
        for i in inlets:
            name  = str(i['name'])
            lfile = _normpath(str(i['line_file']))
            tfile = _normpath(str(i['timeseries_file']))
            if not os.path.exists(lfile):
                raise FileNotFoundError(
                    f"Inlet '{name}': line_file not found: {lfile!r}")
            if not os.path.exists(tfile):
                raise FileNotFoundError(
                    f"Inlet '{name}': timeseries_file not found: {tfile!r}")
            self.inlet_data.append(
                [name, lfile, tfile, float(i.get('start_time', 0.0))])

    # -----------------------------------------------------------------------
    # Rainfall
    # -----------------------------------------------------------------------

    # Erosion / scour operators. anuga/operators/ has provided these for years
    # (erosion_operators.py, sanddune_erosion_operator.py) but nothing in
    # anuga/scenario/ reached them, so any scenario involving bed erosion had to
    # drop out of the TOML workflow entirely.
    #
    # Five behaviours, not seven classes: Circular_erosion_operator and
    # Polygonal_erosion_operator are thin region-specification wrappers over the
    # base class and add no erosion physics, so they are reached here by giving
    # 'simple' a center/radius or a polygon rather than by naming them.
    EROSION_TYPES = ('simple', 'bed_shear', 'flat_slice', 'flat_fill', 'sand_dune')

    # Parameters accepted only by particular types. Validated rather than
    # ignored: silently dropping shear_factor from a flat_slice entry would let
    # a user believe they had configured something they had not.
    EROSION_TYPE_PARAMS = {
        'bed_shear': ('shear_factor',),
        'flat_slice': ('elevation',),
        'flat_fill': ('elevation',),
        'sand_dune': ('Ra',),
    }

    def _parse_erosion(self, erosion, _v):
        self.erosion_data = []
        for i, e in enumerate(erosion):
            sec = f'erosion[{i}]'
            etype = _v.require(e, 'type', sec)
            _v.one_of(etype, self.EROSION_TYPES, 'type', sec)

            # Region: polygon OR center+radius, exactly one. 'indices' is
            # deliberately not supported — raw triangle indices are a Python-API
            # convenience that cannot survive a re-mesh, so they have no stable
            # meaning in a declarative config.
            polygon = e.get('polygon')
            center, radius = e.get('center'), e.get('radius')
            has_poly = polygon is not None
            has_circle = center is not None or radius is not None
            if has_poly and has_circle:
                _v.errors.append(
                    f"[{sec}]: give either 'polygon' or 'center'+'radius', not both")
            elif not has_poly and not has_circle:
                _v.errors.append(
                    f"[{sec}]: needs a region — either 'polygon' or 'center'+'radius'")
            if has_circle:
                if center is None or radius is None:
                    _v.errors.append(
                        f"[{sec}]: 'center' and 'radius' must be given together")
                else:
                    if not (isinstance(center, (list, tuple)) and len(center) == 2):
                        _v.errors.append(
                            f"[{sec}] 'center': expected [x, y], got {center!r}")
                    _v.positive(float(radius), 'radius', sec)
            if has_poly:
                polygon = _normpath(str(polygon))

            row = {
                'type': etype,
                'polygon': polygon,
                'center': list(center) if has_circle and center is not None else None,
                'radius': float(radius) if has_circle and radius is not None else None,
                'threshold': float(e.get('threshold', 0.0)),
                'base': float(e.get('base', 0.0)),
                'description': e.get('description'),
                'label': e.get('label'),
                'logging': bool(e.get('logging', False)),
            }

            # Type-specific parameters: accept this type's own, reject others'.
            allowed = self.EROSION_TYPE_PARAMS.get(etype, ())
            for other_type, params in self.EROSION_TYPE_PARAMS.items():
                if other_type == etype:
                    continue
                for prm in params:
                    if prm in e and prm not in allowed:
                        _v.errors.append(
                            f"[{sec}] {prm!r}: only valid for type "
                            f"{other_type!r}, not {etype!r}")
            for prm in allowed:
                if prm in e:
                    row[prm] = float(e[prm])

            self.erosion_data.append(row)

        # Erosion changes elevation as the run proceeds, but ANUGA stores
        # elevation statically by default — so the .sww shows the initial
        # terrain forever and the erosion is invisible in every downstream
        # product. When erosion is present, default to storing elevation every
        # timestep so the changing bed is recorded. Only warn if the user
        # explicitly asked for static storage: don't silently override a choice
        # they wrote deliberately, but do tell them the eroded bed won't appear.
        if self.erosion_data and not self.store_elevation_every_timestep:
            if getattr(self, '_store_elevation_explicit', False):
                warnings.warn(
                    'Scenario defines [[erosion]] operators but [project] '
                    'store_elevation_every_timestep is explicitly false: '
                    'elevation will be written once at t=0, so the eroded bed '
                    'will not appear in the .sww or any raster derived from it.',
                    stacklevel=2)
            else:
                # Not set by the user: adopt the erosion-appropriate default.
                self.store_elevation_every_timestep = True

    # -----------------------------------------------------------------------
    # Sediment transport (suspended load with optional bedload and bed
    # evolution): the [sediment] table maps onto Domain.initialize_sediment_
    # operator and its setters, and each [[sediment.fractions]] onto
    # Domain.add_sediment_fraction. Every key is optional except a fraction's
    # name and diameter; an absent key leaves the domain's own default.
    # -----------------------------------------------------------------------
    SEDIMENT_SHEAR_CLOSURES = ('quadratic_drag', 'depth_slope', 'energy_slope')
    SEDIMENT_BED_MATERIALS = ('noncohesive', 'cohesive', 'partheniades')
    SEDIMENT_DEPOSITION_LAWS = ('d_star', 'threshold')
    SEDIMENT_NEAR_BED = ('constant', 'rouse')
    SEDIMENT_FRICTION_MODES = ('constant', 'larsen_lamb', 'wilson')
    SEDIMENT_WILSON_BEDS = ('sand', 'gravel', 'boulder')
    SEDIMENT_BEDLOAD = ('off', 'wong_parker_eq24', 'wong_parker_eq23',
                        'engelund_hansen', 'grass')
    SEDIMENT_KEYS = (
        'porosity', 'c_max', 'c_pack', 'bed_evolution', 'rho_w',
        'shear_closure', 'max_slope', 'freeze_slope', 'bed_material', 'tau_crit', 'K_e',
        'entrainment', 'de_leeuw_fit', 'skin_roughness',
        'deposition_law', 'tau_d', 'near_bed', 'reference_height_floor',
        'friction_mode', 'k_s', 'sigma_br', 'r_d', 'r_br', 'bed', 'grain_size',
        'bedload', 'bedload_K', 'bedload_m', 'bedload_tau_c_star',
        'bedload_open_boundaries',
        'angle_of_repose', 'repose_relax', 'repose_max_sweeps',
        'erodible_base_elevation', 'erodible_base_depth',
        'fractions', 'erodible_regions')
    SEDIMENT_FRACTION_KEYS = (
        'name', 'diameter', 'rho_s', 'tau_c_star', 'd_star',
        'initial_concentration', 'reference_height', 'nu', 'C1', 'C2',
        'boundary')

    def _parse_sediment(self, sed, _v):
        self.sediment_data = None
        if sed is None:
            return
        sec = 'sediment'
        if not isinstance(sed, dict):
            _v.errors.append(f'[{sec}]: expected a table, got {type(sed).__name__}')
            return

        for key in sed:
            if key not in self.SEDIMENT_KEYS:
                _v.errors.append(f'[{sec}] unknown key {key!r}')

        d = {}

        def opt_float(key, check=None):
            if key in sed:
                try:
                    val = _num(sed[key])
                except (TypeError, ValueError):
                    _v.errors.append(
                        f'[{sec}] {key!r}: expected a number, got {sed[key]!r}')
                    return
                if check is not None:
                    check(val, key, sec)
                d[key] = float(val)

        def opt_choice(key, choices):
            if key in sed:
                _v.one_of(sed[key], choices, key, sec)
                d[key] = sed[key]

        opt_float('porosity', lambda v, k, s: _v.in_range(v, 0.0, 0.99, k, s))
        opt_float('c_max', _v.positive)
        opt_float('c_pack', _v.positive)
        opt_float('rho_w', _v.positive)
        if 'bed_evolution' in sed:
            d['bed_evolution'] = bool(sed['bed_evolution'])

        opt_choice('shear_closure', self.SEDIMENT_SHEAR_CLOSURES)
        opt_float('max_slope', _v.non_negative)
        if 'freeze_slope' in sed:
            d['freeze_slope'] = bool(sed['freeze_slope'])
            if d['freeze_slope'] and sed.get('shear_closure') != 'depth_slope':
                _v.errors.append(
                    f"[{sec}] 'freeze_slope' applies to shear_closure "
                    f"'depth_slope' only")
        if ('max_slope' in sed) and sed.get('shear_closure', 'quadratic_drag') == 'quadratic_drag':
            _v.errors.append(
                f"[{sec}] 'max_slope' applies to shear_closure 'depth_slope' "
                f"or 'energy_slope'; quadratic_drag uses no slope")
        opt_choice('bed_material', self.SEDIMENT_BED_MATERIALS)
        opt_float('tau_crit', _v.positive)
        opt_float('K_e', _v.positive)
        if ('tau_crit' in sed or 'K_e' in sed) and \
                sed.get('bed_material', 'noncohesive') == 'noncohesive':
            _v.errors.append(
                f"[{sec}] 'tau_crit' / 'K_e' apply to bed_material 'cohesive' "
                f"or 'partheniades'; the noncohesive route uses each "
                f"fraction's 'tau_c_star'")

        opt_choice('entrainment', ('smith_mclean', 'de_leeuw'))
        opt_choice('de_leeuw_fit', ('de_leeuw_2020', 'nghiem_2022'))
        opt_float('skin_roughness', _v.positive)
        if sed.get('entrainment') == 'de_leeuw' and \
                sed.get('bed_material', 'noncohesive') != 'noncohesive':
            _v.errors.append(
                f"[{sec}] entrainment 'de_leeuw' is a non-cohesive relation; "
                f"use it with bed_material 'noncohesive'")
        if ('de_leeuw_fit' in sed or 'skin_roughness' in sed) and \
                sed.get('entrainment') != 'de_leeuw':
            _v.errors.append(
                f"[{sec}] 'de_leeuw_fit' / 'skin_roughness' apply to "
                f"entrainment 'de_leeuw'")

        opt_choice('deposition_law', self.SEDIMENT_DEPOSITION_LAWS)
        opt_float('tau_d', _v.non_negative)
        opt_choice('near_bed', self.SEDIMENT_NEAR_BED)
        opt_float('reference_height_floor', _v.non_negative)

        opt_choice('friction_mode', self.SEDIMENT_FRICTION_MODES)
        for key in ('k_s', 'sigma_br', 'r_d', 'r_br', 'grain_size'):
            opt_float(key, _v.positive)
        opt_choice('bed', self.SEDIMENT_WILSON_BEDS)
        mode = sed.get('friction_mode', 'constant')
        if mode == 'larsen_lamb' and 'k_s' not in sed and 'sigma_br' not in sed:
            _v.errors.append(
                f"[{sec}] friction_mode 'larsen_lamb' needs 'k_s' or 'sigma_br' "
                f"(site-measured; there is no default)")
        if mode == 'wilson' and 'grain_size' not in sed:
            _v.errors.append(
                f"[{sec}] friction_mode 'wilson' needs 'grain_size' "
                f"(D50 for sand, D84 for gravel/boulder)")

        opt_choice('bedload', self.SEDIMENT_BEDLOAD)
        opt_float('bedload_K', _v.positive)
        opt_float('bedload_m', _v.positive)
        opt_float('bedload_tau_c_star', _v.non_negative)
        if sed.get('bedload') == 'grass' and 'bedload_K' not in sed:
            _v.errors.append(
                f"[{sec}] bedload 'grass' needs 'bedload_K', the Grass "
                f"coefficient A_g (a calibration; there is no default)")
        if 'bedload_open_boundaries' in sed:
            tags = sed['bedload_open_boundaries']
            if isinstance(tags, str):
                tags = [tags]
            if not isinstance(tags, list) or \
                    not all(isinstance(t, str) for t in tags):
                _v.errors.append(
                    f"[{sec}] 'bedload_open_boundaries': expected a list of "
                    f"boundary tags, got {sed['bedload_open_boundaries']!r}")
            else:
                d['bedload_open_boundaries'] = list(tags)

        opt_float('angle_of_repose',
                  lambda v, k, s: _v.in_range(v, 0.0, 90.0, k, s))
        opt_float('repose_relax', lambda v, k, s: _v.in_range(v, 0.0, 1.0, k, s))
        if 'repose_max_sweeps' in sed:
            d['repose_max_sweeps'] = int(sed['repose_max_sweeps'])
            _v.positive(d['repose_max_sweeps'], 'repose_max_sweeps', sec)

        opt_float('erodible_base_elevation')
        opt_float('erodible_base_depth', _v.positive)
        if 'erodible_base_elevation' in sed and 'erodible_base_depth' in sed:
            _v.errors.append(
                f"[{sec}]: give 'erodible_base_elevation' or "
                f"'erodible_base_depth', not both")

        # Fractions: at least one, each with a name and a diameter
        fractions = sed.get('fractions', [])
        if not isinstance(fractions, list) or not fractions:
            _v.errors.append(
                f'[{sec}]: needs at least one [[sediment.fractions]] entry '
                f'(name and diameter)')
            fractions = []
        d['fractions'] = []
        names = set()
        for i, fr in enumerate(fractions):
            fsec = f'sediment.fractions[{i}]'
            for key in fr:
                if key not in self.SEDIMENT_FRACTION_KEYS:
                    _v.errors.append(f'[{fsec}] unknown key {key!r}')
            row = {}
            name = _v.require(fr, 'name', fsec)
            if name is not None:
                name = str(name)
                if name in names:
                    _v.errors.append(f'[{fsec}] duplicate fraction name {name!r}')
                names.add(name)
            row['name'] = name
            diameter = _v.to_float(fr, 'diameter', fsec)
            _v.positive(diameter, 'diameter', fsec)
            row['diameter'] = diameter
            for key, check in (('rho_s', _v.positive), ('tau_c_star', _v.non_negative),
                               ('d_star', _v.positive),
                               ('initial_concentration', _v.non_negative),
                               ('reference_height', _v.positive),
                               ('nu', _v.positive), ('C1', _v.positive),
                               ('C2', _v.positive)):
                if key in fr:
                    try:
                        val = _num(fr[key])
                    except (TypeError, ValueError):
                        _v.errors.append(
                            f'[{fsec}] {key!r}: expected a number, got {fr[key]!r}')
                        continue
                    check(val, key, fsec)
                    row[key] = float(val)
            boundary = fr.get('boundary')
            if boundary is not None:
                if not isinstance(boundary, dict):
                    _v.errors.append(
                        f"[{fsec}] 'boundary': expected a table of tag = "
                        f"concentration, got {boundary!r}")
                else:
                    row['boundary'] = {}
                    for tag, val in boundary.items():
                        try:
                            c = float(_num(val))
                        except (TypeError, ValueError):
                            _v.errors.append(
                                f"[{fsec}] boundary {tag!r}: expected a "
                                f"concentration, got {val!r}")
                            continue
                        _v.non_negative(c, f'boundary.{tag}', fsec)
                        row['boundary'][str(tag)] = c
            d['fractions'].append(row)

        # Erodible regions: polygon OR center+radius, as [[erosion]] takes them
        regions = sed.get('erodible_regions', [])
        d['erodible_regions'] = []
        for i, r in enumerate(regions if isinstance(regions, list) else []):
            rsec = f'sediment.erodible_regions[{i}]'
            polygon = r.get('polygon')
            center, radius = r.get('center'), r.get('radius')
            has_poly = polygon is not None
            has_circle = center is not None or radius is not None
            if has_poly and has_circle:
                _v.errors.append(
                    f"[{rsec}]: give either 'polygon' or 'center'+'radius', not both")
            elif not has_poly and not has_circle:
                _v.errors.append(
                    f"[{rsec}]: needs a region — either 'polygon' or 'center'+'radius'")
            if has_circle:
                if center is None or radius is None:
                    _v.errors.append(
                        f"[{rsec}]: 'center' and 'radius' must be given together")
                else:
                    if not (isinstance(center, (list, tuple)) and len(center) == 2):
                        _v.errors.append(
                            f"[{rsec}] 'center': expected [x, y], got {center!r}")
                    _v.positive(float(radius), 'radius', rsec)
            d['erodible_regions'].append({
                'polygon': _normpath(str(polygon)) if has_poly else None,
                'center': list(center) if has_circle and center is not None else None,
                'radius': float(radius) if has_circle and radius is not None else None,
                'erodible': bool(r.get('erodible', True)),
            })

        self.sediment_data = d

        # Sediment transport moves the bed unless bed_evolution is false; as
        # for [[erosion]], record the bed every timestep unless the user
        # explicitly chose static storage.
        if d.get('bed_evolution', True) and not self.store_elevation_every_timestep:
            if getattr(self, '_store_elevation_explicit', False):
                warnings.warn(
                    'Scenario defines [sediment] with bed evolution but '
                    '[project] store_elevation_every_timestep is explicitly '
                    'false: elevation will be written once at t=0, so the '
                    'evolving bed will not appear in the .sww.',
                    stacklevel=2)
            else:
                self.store_elevation_every_timestep = True

    def _parse_rainfall(self, rainfall):
        self.rain_data = []
        for r in rainfall:
            tfile = _normpath(str(r['timeseries_file']))
            if not os.path.exists(tfile):
                raise FileNotFoundError(
                    f"Rainfall: timeseries_file not found: {tfile!r}")
            polygon = _normpath(str(r.get('polygon', 'All')))
            if polygon not in ('All',) and not os.path.exists(polygon):
                raise FileNotFoundError(
                    f"Rainfall: polygon file not found: {polygon!r}")
            row = [
                tfile,
                float(r.get('start_time', 0.0)),
                str(r.get('interpolation_type', 'linear')),
                polygon,
                float(r.get('multiplier', 1.0)),
            ]
            self.rain_data.append(row)

    # -----------------------------------------------------------------------
    # Bridges
    # -----------------------------------------------------------------------

    def _parse_bridges(self, bridges):
        self.bridge_data = []
        for b in bridges:
            if not b.get('enabled', True):
                continue

            row = [
                str(b['label']),                             # [0]
                _normpath(str(b['deck_file'])),              # [1] — also used as breakline
                float(b['deck_elevation']),                  # [2] — also used as elevation value
                _normpath(str(b['exchange_line_0'])),        # [3]
                _normpath(str(b['exchange_line_1'])),        # [4]
                float(b['enquiry_gap']),                     # [5]
                _normpath(str(b['internal_boundary_curve_file'])),  # [6]
                float(b.get('vertical_datum_offset', 0.0)), # [7]
                float(b.get('smoothing_timescale', 0.0)),   # [8]
            ]

            for idx, desc in [(1, 'deck_file'),
                              (3, 'exchange_line_0'),
                              (4, 'exchange_line_1'),
                              (6, 'internal_boundary_curve_file')]:
                if not os.path.exists(row[idx]):
                    raise FileNotFoundError(
                        f"Bridge '{row[0]}': {desc} not found: {row[idx]!r}")

            # Prepend deck as a breakline and elevation source (highest priority)
            self.breakline_files.insert(0, row[1])
            self.elevation_data      = [[row[1], row[2]]] + self.elevation_data
            self.elevation_clip_range = [[float('-inf'), float('inf')]] + self.elevation_clip_range

            self.bridge_data.append(row)

    # -----------------------------------------------------------------------
    # Pumping stations
    # -----------------------------------------------------------------------

    def _parse_pumping_stations(self, stations):
        self.pumping_station_data = []
        for ps in stations:
            if not ps.get('enabled', True):
                continue

            row = [
                str(ps['label']),                           # [0]
                float(ps['pump_capacity']),                 # [1]
                float(ps['pump_rate_of_increase']),         # [2]
                float(ps['pump_rate_of_decrease']),         # [3]
                float(ps['hw_to_start_pumping']),           # [4]
                float(ps['hw_to_stop_pumping']),            # [5]
                _normpath(str(ps['basin_polygon_file'])),   # [6] — also breakline + elevation
                float(ps['basin_elevation']),               # [7] — elevation value
                _normpath(str(ps['exchange_line_0'])),      # [8]
                _normpath(str(ps['exchange_line_1'])),      # [9]
                float(ps.get('smoothing_timescale', 0.0)), # [10]
            ]

            for idx, desc in [(6, 'basin_polygon_file'),
                              (8, 'exchange_line_0'),
                              (9, 'exchange_line_1')]:
                if not os.path.exists(row[idx]):
                    raise FileNotFoundError(
                        f"Pumping station '{row[0]}': {desc} not found: {row[idx]!r}")

            # Prepend basin as a breakline and elevation source (highest priority)
            self.breakline_files.insert(0, row[6])
            self.elevation_data      = [[row[6], row[7]]] + self.elevation_data
            self.elevation_clip_range = [[float('-inf'), float('inf')]] + self.elevation_clip_range

            self.pumping_station_data.append(row)

    # -----------------------------------------------------------------------
    # Culverts  (Boyd box and Boyd pipe)
    # -----------------------------------------------------------------------

    def _parse_culverts(self, culverts, _v):
        """Parse ``[[culverts]]`` entries into ``self.culvert_data``.

        Each entry is stored as a dict so that ``setup_culverts.py`` can
        access parameters by name rather than by positional index.

        Supported ``type`` values: ``"boyd_box"`` (default), ``"boyd_pipe"``.

        Geometry is specified with either:
        * ``exchange_line_0`` / ``exchange_line_1`` — paths to CSV polyline
          files that define the upstream / downstream exchange zones, or
        * ``end_point_0`` / ``end_point_1`` — ``[x, y]`` arrays for the two
          culvert barrel ends (simpler; no polyline files required).
        """
        self.culvert_data = []
        for c in culverts:
            if not c.get('enabled', True):
                continue

            ctype = str(c.get('type', 'boyd_box'))
            if ctype not in ('boyd_box', 'boyd_pipe'):
                raise ValueError(
                    f"Culvert '{c.get('label', '?')}': unknown type {ctype!r}. "
                    f"Expected 'boyd_box' or 'boyd_pipe'.")

            label   = str(c['label'])
            sec     = f'culverts[{label!r}]'
            losses  = _parse_losses(c.get('losses', 0.0), sec, _v)
            barrels = float(c.get('barrels', 1.0))
            blockage= float(c.get('blockage', 0.0))
            manning = float(c.get('manning', 0.013))
            eq_gap  = float(c.get('enquiry_gap', 0.2))
            apron   = float(c.get('apron', 0.1))
            smoothing = float(c.get('smoothing_timescale', 0.0))

            _v.positive(barrels,  'barrels',  sec)
            _v.positive(manning,  'manning',  sec)
            _v.non_negative(eq_gap,   'enquiry_gap', sec)
            _v.non_negative(apron,    'apron',    sec)
            _v.non_negative(smoothing,'smoothing_timescale', sec)
            _v.in_range(blockage, 0.0, 1.0, 'blockage', sec)

            row = {
                'label':               label,
                'type':                ctype,
                'losses':              losses,
                'barrels':             barrels,
                'blockage':            blockage,
                'z1':                  float(c.get('z1', 0.0)),
                'z2':                  float(c.get('z2', 0.0)),
                'apron':               apron,
                'manning':             manning,
                'enquiry_gap':         eq_gap,
                'smoothing_timescale': smoothing,
                'use_momentum_jet':    bool(c.get('use_momentum_jet', True)),
                'use_velocity_head':   bool(c.get('use_velocity_head', True)),
            }

            # Type-specific geometry
            if ctype == 'boyd_box':
                row['width']    = float(c['width'])
                row['height']   = float(c['height']) if 'height' in c else None
                row['diameter'] = None
                _v.positive(row['width'], 'width', sec)
                if row['height'] is not None:
                    _v.positive(row['height'], 'height', sec)
            else:
                row['diameter'] = float(c['diameter'])
                row['width']    = None
                row['height']   = None
                _v.positive(row['diameter'], 'diameter', sec)

            # Exchange lines (file paths) or end points ([x, y] pairs)
            if 'exchange_line_0' in c:
                el0 = _normpath(str(c['exchange_line_0']))
                el1 = _normpath(str(c['exchange_line_1']))
                for fpath, key in [(el0, 'exchange_line_0'),
                                   (el1, 'exchange_line_1')]:
                    if not os.path.exists(fpath):
                        raise FileNotFoundError(
                            f"Culvert '{row['label']}': {key} not found: {fpath!r}")
                row['exchange_line_0'] = el0
                row['exchange_line_1'] = el1
                row['end_point_0']     = None
                row['end_point_1']     = None
            else:
                row['exchange_line_0'] = None
                row['exchange_line_1'] = None
                row['end_point_0']     = list(c['end_point_0'])
                row['end_point_1']     = list(c['end_point_1'])

            # Optional explicit invert elevations [upstream_m, downstream_m]
            if 'invert_elevations' in c:
                row['invert_elevations'] = [float(v) for v in c['invert_elevations']]
            else:
                row['invert_elevations'] = None

            self.culvert_data.append(row)

    # -----------------------------------------------------------------------
    # Weirs  (weir / orifice with trapezoidal cross-section)
    # -----------------------------------------------------------------------

    def _parse_weirs(self, weirs, _v):
        """Parse ``[[weirs]]`` entries into ``self.weir_data``.

        Uses the ``Weir_orifice_trapezoid_operator``.  Parameters and geometry
        specification follow the same conventions as ``[[culverts]]``.
        """
        self.weir_data = []
        for w in weirs:
            if not w.get('enabled', True):
                continue

            label    = str(w['label'])
            sec      = f'weirs[{label!r}]'
            width    = float(w['width'])
            height   = float(w['height']) if 'height' in w else None
            losses   = _parse_losses(w.get('losses', 0.0), sec, _v)
            barrels  = float(w.get('barrels', 1.0))
            blockage = float(w.get('blockage', 0.0))
            manning  = float(w.get('manning', 0.013))
            eq_gap   = float(w.get('enquiry_gap', 0.0))
            apron    = float(w.get('apron', 0.1))
            smoothing= float(w.get('smoothing_timescale', 0.0))

            _v.positive(width,   'width',   sec)
            _v.positive(barrels, 'barrels', sec)
            _v.positive(manning, 'manning', sec)
            _v.non_negative(eq_gap,    'enquiry_gap', sec)
            _v.non_negative(apron,     'apron',    sec)
            _v.non_negative(smoothing, 'smoothing_timescale', sec)
            _v.in_range(blockage, 0.0, 1.0, 'blockage', sec)
            if height is not None:
                _v.positive(height, 'height', sec)

            row = {
                'label':               label,
                'width':               width,
                'height':              height,
                'losses':              losses,
                'barrels':             barrels,
                'blockage':            blockage,
                'z1':                  float(w.get('z1', 0.0)),
                'z2':                  float(w.get('z2', 0.0)),
                'apron':               apron,
                'manning':             manning,
                'enquiry_gap':         eq_gap,
                'smoothing_timescale': smoothing,
                'use_momentum_jet':    bool(w.get('use_momentum_jet', True)),
                'use_velocity_head':   bool(w.get('use_velocity_head', True)),
            }

            # Exchange lines or end points
            if 'exchange_line_0' in w:
                el0 = _normpath(str(w['exchange_line_0']))
                el1 = _normpath(str(w['exchange_line_1']))
                for fpath, key in [(el0, 'exchange_line_0'),
                                   (el1, 'exchange_line_1')]:
                    if not os.path.exists(fpath):
                        raise FileNotFoundError(
                            f"Weir '{row['label']}': {key} not found: {fpath!r}")
                row['exchange_line_0'] = el0
                row['exchange_line_1'] = el1
                row['end_point_0']     = None
                row['end_point_1']     = None
            else:
                row['exchange_line_0'] = None
                row['exchange_line_1'] = None
                row['end_point_0']     = list(w['end_point_0'])
                row['end_point_1']     = list(w['end_point_1'])

            # Optional explicit invert elevations [upstream_m, downstream_m]
            if 'invert_elevations' in w:
                row['invert_elevations'] = [float(v) for v in w['invert_elevations']]
            else:
                row['invert_elevations'] = None

            self.weir_data.append(row)
