#!/usr/bin/env python3
"""GUI viewer for ANUGA SWW file animations.

Opens an SWW file, generates PNG frames for a chosen quantity using
SWW_plotter, then plays them as an animation.

Usage::

    anuga_sww_gui
    anuga_sww_gui --sww domain.sww
    anuga_sww_gui --sww domain.sww --qty depth
"""

import os
import glob
import sys

try:
    import tkinter as tk
    from tkinter import ttk, filedialog
except ImportError:
    sys.exit(
        "Error: tkinter is not available.\n"
        "On Debian/Ubuntu: sudo apt install python3-tk\n"
        "On conda: conda install tk"
    )

try:
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib import image as mpimage
except ImportError as _e:
    sys.exit(
        f"Error: matplotlib is required but could not be loaded: {_e}\n"
        "Install it with: conda install matplotlib  or  pip install matplotlib"
    )
except Exception as _e:
    sys.exit(
        f"Error: could not initialise the TkAgg matplotlib backend: {_e}\n"
        "Ensure a display is available and tkinter is installed."
    )


# ------------------------------------------------------------------ #
# Quantity configuration                                              #
# ------------------------------------------------------------------ #

_QUANTITIES = ('depth', 'stage', 'speed', 'speed_depth',
               'max_depth', 'max_stage', 'max_speed', 'max_speed_depth',
               'elev', 'elev_delta')

_QTY_DEFAULTS = {
    'depth':           dict(vmin=0.0,   vmax=20.0),
    'stage':           dict(vmin=-20.0, vmax=20.0),
    'speed':           dict(vmin=0.0,   vmax=10.0),
    'speed_depth':     dict(vmin=0.0,   vmax=20.0),
    'max_depth':       dict(vmin=0.0,   vmax=20.0),
    'max_stage':       dict(vmin=-20.0, vmax=20.0),
    'max_speed':       dict(vmin=0.0,   vmax=10.0),
    'max_speed_depth': dict(vmin=0.0,   vmax=20.0),
    'elev':            dict(vmin=-20.0, vmax=100.0),
    'elev_delta':      dict(vmin=-5.0,  vmax=5.0),
}

# Attribute on SWW_plotter used to compute auto vmin/vmax
_QTY_DATA_ATTR = {
    'depth':           'depth',
    'stage':           'stage',
    'speed':           'speed',
    'speed_depth':     'speed_depth',
    'max_depth':       'max_depth',
    'max_stage':       'max_stage',
    'max_speed':       'max_speed',
    'max_speed_depth': 'max_speed_depth',
    'elev':            'elev',
    'elev_delta':      'elev_delta',
}

# Method name on SWW_plotter to save a single frame
_QTY_SAVE_METHOD = {
    'depth':           'save_depth_frame',
    'stage':           'save_stage_frame',
    'speed':           'save_speed_frame',
    'speed_depth':     'save_speed_depth_frame',
    'max_depth':       'save_max_depth_frame',
    'max_stage':       'save_max_stage_frame',
    'max_speed':       'save_max_speed_frame',
    'max_speed_depth': 'save_max_speed_depth_frame',
    'elev':            'save_elev_frame',
    'elev_delta':      'save_elev_delta_frame',
}

# Colorbar label for each quantity (used by the clean-export renderer)
_QTY_CBAR_LABEL = {
    'depth':           'Depth (m)',
    'stage':           'Stage (m)',
    'speed':           'Speed (m/s)',
    'speed_depth':     'Flow depth (m2/s)',
    'max_depth':       'Max depth (m)',
    'max_stage':       'Max stage (m)',
    'max_speed':       'Max speed (m/s)',
    'max_speed_depth': 'Max flow depth (m2/s)',
    'elev':            'Elevation (m)',
    'elev_delta':      'Elevation change (m)',
}

# Human-readable quantity name for figure titles
_QTY_TITLE = {
    'depth':           'Water depth',
    'stage':           'Water stage',
    'speed':           'Flow speed',
    'speed_depth':     'Flow depth',
    'max_depth':       'Maximum water depth',
    'max_stage':       'Maximum water stage',
    'max_speed':       'Maximum flow speed',
    'max_speed_depth': 'Maximum flow depth',
    'elev':            'Bed elevation',
    'elev_delta':      'Bed elevation change',
}


# Tracers are user-defined, so unlike everything above they cannot be listed
# here: which ones exist is a property of the SWW file being opened. Loading a
# file appends its tracers (and their maxima, when the run stored them) to this
# list and to the five tables above, keyed 'tracer_<name>' / 'max_tracer_<name>'.
# Everything downstream then treats a tracer exactly like stage or depth.
_TRACER_QUANTITIES = []


def _all_quantities():
    """The fixed quantities plus whatever tracers the open file carries."""
    return tuple(_QUANTITIES) + tuple(_TRACER_QUANTITIES)


_CMAPS = ['viridis', 'plasma', 'inferno', 'magma', 'cividis',
          'Blues', 'Greens', 'Oranges', 'Reds', 'YlOrRd',
          'RdBu_r', 'coolwarm', 'seismic', 'jet', 'turbo',
          'terrain', 'gist_earth', 'gray']

# ------------------------------------------------------------------ #
# Frame helpers (shared with anuga_animate_gui convention)           #
# ------------------------------------------------------------------ #

def _find_frames(plot_dir, quantity, prefix):
    """Return sorted PNG paths matching prefix_quantity_*.png in plot_dir."""
    pattern = os.path.join(plot_dir, f'{prefix}_{quantity}_*.png')
    return sorted(glob.glob(pattern))


def _toml_value(v):
    """Format a Python value as an inline TOML literal."""
    if isinstance(v, bool):
        return 'true' if v else 'false'
    if isinstance(v, str):
        return '"' + v.replace('\\', '\\\\').replace('"', '\\"') + '"'
    return repr(v)


def _apply_config_to_gui(data, gui):
    """Apply a flat/sectioned TOML dict to a SWWAnimationGUI instance.

    Sections ``[render]`` and ``[generate]`` are merged with top-level keys
    so both flat and sectioned TOML layouts are accepted.
    """
    from anuga.utilities.animate import BASEMAP_PROVIDERS

    # Flatten render/generate sections into a single dict
    cfg = {}
    cfg.update(data.get('render',   {}))
    cfg.update(data.get('generate', {}))
    cfg.update(data.get('file',     {}))
    # Also accept top-level keys (flat layout)
    for k, v in data.items():
        if not isinstance(v, dict):
            cfg[k] = v

    if 'qty' in cfg and cfg['qty'] in _all_quantities():
        gui._qty_var.set(cfg['qty'])
    if 'vmin' in cfg:
        gui._vmin_var.set(str(cfg['vmin']))
        gui._auto_var.set(False)
    if 'vmax' in cfg:
        gui._vmax_var.set(str(cfg['vmax']))
        gui._auto_var.set(False)
    if 'cmap' in cfg and cfg['cmap'] in _CMAPS:
        gui._cmap_var.set(cfg['cmap'])
    if 'cmap_reverse' in cfg:
        gui._cmap_reverse_var.set(bool(cfg['cmap_reverse']))
    if 'mindepth' in cfg:
        gui._mindepth_var.set(str(cfg['mindepth']))
    if 'flat_view' in cfg:
        gui._show_edges_var.set(bool(cfg['flat_view']))
    if 'outdir' in cfg:
        gui._dir_var.set(str(cfg['outdir']))
    if 'dpi' in cfg:
        gui._dpi_var.set(str(cfg['dpi']))
    if 'stride' in cfg:
        gui._stride_var.set(str(cfg['stride']))
    if 'alpha' in cfg:
        gui._alpha_var.set(float(cfg['alpha']))
    if 'basemap_provider' in cfg:
        provider = cfg['basemap_provider']
        if provider in BASEMAP_PROVIDERS:
            gui._basemap_provider_var.set(provider)
    if 'epsg' in cfg and gui._splotter is not None:
        gui._epsg_var.set(str(cfg['epsg']))
        gui._on_set_epsg()
    if 'basemap' in cfg and gui._splotter is not None:
        gui._basemap_var.set(bool(cfg['basemap']))
        gui._on_basemap_toggle()


# ------------------------------------------------------------------ #
# GUI                                                                 #
# ------------------------------------------------------------------ #

class SWWAnimationGUI:
    """Tkinter + matplotlib GUI: open SWW -> generate frames -> animate."""

    def __init__(self, root, initial_sww=None, initial_qty=None,
                 initial_vmin=None, initial_vmax=None,
                 initial_cmap=None, initial_cmap_reverse=False,
                 initial_mindepth=None, initial_flat_view=False,
                 initial_outdir=None, initial_dpi=None, initial_stride=None,
                 initial_epsg=None, initial_basemap=None,
                 initial_basemap_provider=None, initial_alpha=None):
        self.root = root
        self.root.title('ANUGA SWW Animator')
        # Scale the minimum window size to the UI scale (set in main on HiDPI).
        ui_scale = getattr(root, 'ui_scale', 1.0) or 1.0
        self.root.minsize(int(860 * ui_scale), int(640 * ui_scale))
        self.root.protocol('WM_DELETE_WINDOW', self._on_close)

        # animation state
        self._frames = []
        self._current = 0
        self._playing = False
        self._after_id = None
        self._im = None

        # generation state
        self._splotter = None
        self._cancel_flag = False
        self._gen_after_id = None
        self._executor = None
        self._futures = []
        self._gen_plot_dir = ''
        self._gen_qty = 'depth'
        self._n_to_gen = 0

        # zoom state
        self._zoom_xlim = None
        self._zoom_ylim = None
        self._zoom_rect_patch = None
        self._zoom_selector = None
        self._zoom_mode = False

        # mesh overlay state
        self._mesh_overlay_lines = []

        # elevation contour overlay state
        self._elev_overlay_artists = []
        self._elev_contour_data = None    # cached (triang_px, elev_nodes, levels)

        # timeseries state
        self._ts_triangles = []       # list of picked triangle indices (in order)
        self._ts_vline = None
        self._pick_mode = False
        self._pick_cid = None
        self._pick_key_cid = None
        self._pick_overlays = []      # one artist per picked point
        self._pick_text = None
        self._plot_transform = None
        self._gen_used_basemap = False
        self._last_gen_dpi = 100
        self._last_gen_vmin = 0.0
        self._last_gen_vmax = 20.0
        self._last_gen_cmap = 'viridis'
        self._last_gen_qty = 'depth'
        self._last_gen_show_edges = False
        self._last_gen_show_elev = False
        self._last_gen_elev_levels = 10
        self._last_gen_show_mesh = False

        # cross-section state
        self._xs_mode = False
        self._xs_pts = []          # list of (x, y) in RELATIVE mesh coords, max 2
        self._xs_artists = []      # overlay artists on _ax
        self._xs_cid = None
        self._xs_key_cid = None
        self._xs_flux = None       # ndarray (n_times,) discharge m³/s
        self._xs_vline = None      # vertical cursor in xs plot

        # hover coordinate readout state
        self._hover_cid = None
        self._trifinder = None

        self._build_ui()

        # Apply render-tab params (independent of SWW file)
        if initial_qty and initial_qty in _all_quantities():
            self._qty_var.set(initial_qty)
        if initial_vmin is not None:
            self._vmin_var.set(str(initial_vmin))
            self._auto_var.set(False)
        if initial_vmax is not None:
            self._vmax_var.set(str(initial_vmax))
            self._auto_var.set(False)
        if initial_cmap and initial_cmap in _CMAPS:
            self._cmap_var.set(initial_cmap)
        if initial_cmap_reverse:
            self._cmap_reverse_var.set(True)
        if initial_mindepth is not None:
            self._mindepth_var.set(str(initial_mindepth))
        if initial_flat_view:
            self._show_edges_var.set(True)
        # Apply generate-tab params
        if initial_outdir:
            self._dir_var.set(initial_outdir)
        if initial_dpi is not None:
            self._dpi_var.set(str(initial_dpi))
        if initial_stride is not None:
            self._stride_var.set(str(initial_stride))
        if initial_alpha is not None:
            self._alpha_var.set(initial_alpha)
        if initial_basemap_provider:
            from anuga.utilities.animate import BASEMAP_PROVIDERS
            if initial_basemap_provider in BASEMAP_PROVIDERS:
                self._basemap_provider_var.set(initial_basemap_provider)

        if initial_sww:
            self._sww_var.set(os.path.abspath(initial_sww))
            self._on_sww_change()

        # EPSG and basemap override are applied after SWW load (needs _splotter)
        if initial_epsg is not None and self._splotter is not None:
            self._epsg_var.set(str(initial_epsg))
            self._on_set_epsg()
        if initial_basemap is not None and self._splotter is not None:
            self._basemap_var.set(initial_basemap)
            self._on_basemap_toggle()

    # -------------------------------------------------------------- #
    # UI construction                                                 #
    # -------------------------------------------------------------- #

    def _build_menubar(self):
        """File menu — Open/config/Quit, reusing the toolbar's handlers.

        Quit goes through _on_close so it does exactly what the window's close
        button does: cancel any running generation, stop playback and close the
        matplotlib figures.  Calling root.destroy() directly would leave a
        worker pool and open figures behind.
        """
        menubar = tk.Menu(self.root)

        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label='Open SWW...', accelerator='Ctrl+O',
                              command=self._browse_sww)
        file_menu.add_separator()
        file_menu.add_command(label='Load Config...', command=self._load_config)
        file_menu.add_command(label='Save Config...', command=self._save_config)
        file_menu.add_separator()
        file_menu.add_command(label='Quit', accelerator='Ctrl+Q',
                              command=self._on_close)
        menubar.add_cascade(label='File', menu=file_menu)

        self.root.config(menu=menubar)

        # Keyboard accelerators (the labels above are only cosmetic; tk needs
        # the bindings made explicitly).
        self.root.bind_all('<Control-o>', lambda e: self._browse_sww())
        self.root.bind_all('<Control-q>', lambda e: self._on_close())

    def _build_ui(self):
        self._build_menubar()

        # ---- status bar (packed first so it anchors to the bottom) ----
        status_frame = ttk.Frame(self.root)
        status_frame.pack(side=tk.BOTTOM, fill=tk.X)
        self._status_var = tk.StringVar(value='Open an SWW file to begin.')
        ttk.Label(status_frame, textvariable=self._status_var,
                  relief=tk.SUNKEN, anchor=tk.W,
                  padding=(4, 1)).pack(side=tk.LEFT, fill=tk.X, expand=True)
        self._coord_var = tk.StringVar()
        ttk.Label(status_frame, textvariable=self._coord_var,
                  relief=tk.SUNKEN, anchor=tk.E,
                  padding=(4, 1), width=44).pack(side=tk.RIGHT)

        # ================================================================
        # Top control panel
        # ================================================================
        ctrl = ttk.Frame(self.root, padding=(6, 4, 6, 2))
        ctrl.pack(side=tk.TOP, fill=tk.X)

        # ---- SWW file row (always visible) ----
        file_row = ttk.Frame(ctrl)
        file_row.pack(fill=tk.X, pady=(0, 4))
        ttk.Label(file_row, text='SWW file:').pack(side=tk.LEFT)
        self._sww_var = tk.StringVar()
        sww_entry = ttk.Entry(file_row, textvariable=self._sww_var)
        sww_entry.pack(side=tk.LEFT, padx=4, fill=tk.X, expand=True)
        sww_entry.bind('<Return>', lambda _e: self._on_sww_change())
        ttk.Button(file_row, text='Browse...',
                   command=self._browse_sww).pack(side=tk.LEFT)
        ttk.Button(file_row, text='Help',
                   command=self._show_help).pack(side=tk.LEFT, padx=(8, 0))
        self._sww_info_label = ttk.Label(file_row, text='No SWW file loaded.',
                                         foreground='grey')
        self._sww_info_label.pack(side=tk.LEFT, padx=(12, 0))

        # ================================================================
        # Notebook — Plot / Generate / Output tabs
        # ================================================================
        nb = ttk.Notebook(ctrl)
        nb.pack(fill=tk.X, pady=(0, 4))

        # ----------------------------------------------------------------
        # Tab 1 — Plot
        # ----------------------------------------------------------------
        tab_plot = ttk.Frame(nb, padding=(8, 6))
        nb.add(tab_plot, text='  Plot  ')

        # Row A: Quantity | vmin / vmax / Auto
        rA = ttk.Frame(tab_plot)
        rA.pack(fill=tk.X, pady=2)
        ttk.Label(rA, text='Quantity:').pack(side=tk.LEFT)
        self._qty_var = tk.StringVar(value='depth')
        qty_combo = ttk.Combobox(rA, textvariable=self._qty_var,
                                 values=list(_all_quantities()), width=13,
                                 state='readonly')
        # Kept so a newly loaded file can add its tracers to the list.
        self._qty_combo = qty_combo
        qty_combo.pack(side=tk.LEFT, padx=(2, 8))
        qty_combo.bind('<<ComboboxSelected>>', lambda _e: self._on_qty_change())

        ttk.Label(rA, text='vmin:').pack(side=tk.LEFT)
        self._vmin_var = tk.StringVar(value='0.0')
        ttk.Entry(rA, textvariable=self._vmin_var, width=8).pack(side=tk.LEFT, padx=2)
        ttk.Label(rA, text='vmax:').pack(side=tk.LEFT, padx=(6, 0))
        self._vmax_var = tk.StringVar(value='20.0')
        ttk.Entry(rA, textvariable=self._vmax_var, width=8).pack(side=tk.LEFT, padx=2)
        self._auto_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(rA, text='Auto from data',
                        variable=self._auto_var,
                        command=self._on_auto_toggle).pack(side=tk.LEFT, padx=(4, 0))

        # Row B: Colormap | Reverse | min depth | Flat View
        rB = ttk.Frame(tab_plot)
        rB.pack(fill=tk.X, pady=2)
        ttk.Label(rB, text='Colormap:').pack(side=tk.LEFT)
        self._cmap_var = tk.StringVar(value='viridis')
        ttk.Combobox(rB, textvariable=self._cmap_var,
                     values=_CMAPS, width=12).pack(side=tk.LEFT, padx=2)
        self._cmap_reverse_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(rB, text='Reverse',
                        variable=self._cmap_reverse_var).pack(side=tk.LEFT, padx=(2, 16))

        ttk.Label(rB, text='min depth:').pack(side=tk.LEFT)
        self._mindepth_var = tk.StringVar(value='0.001')
        ttk.Entry(rB, textvariable=self._mindepth_var, width=8).pack(side=tk.LEFT, padx=2)

        ttk.Separator(rB, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=12)
        self._show_edges_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(rB, text='Flat View',
                        variable=self._show_edges_var).pack(side=tk.LEFT, padx=(0, 12))

        ttk.Separator(rB, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=12)
        self._show_mesh_var = tk.BooleanVar(value=False)
        self._show_mesh_chk = ttk.Checkbutton(rB, text='Show Mesh',
                                               variable=self._show_mesh_var,
                                               command=self._on_show_mesh_toggle,
                                               state=tk.DISABLED)
        self._show_mesh_chk.pack(side=tk.LEFT, padx=(0, 8))

        self._show_elev_var = tk.BooleanVar(value=False)
        self._show_elev_chk = ttk.Checkbutton(rB, text='Show Elev',
                                               variable=self._show_elev_var,
                                               command=self._on_show_elev_toggle,
                                               state=tk.DISABLED)
        self._show_elev_chk.pack(side=tk.LEFT, padx=(0, 2))
        self._elev_levels_var = tk.StringVar(value='10')
        self._elev_levels_spin = ttk.Spinbox(rB, from_=3, to=50, increment=1,
                                              textvariable=self._elev_levels_var,
                                              width=3, state=tk.DISABLED,
                                              command=self._on_elev_levels_changed)
        self._elev_levels_spin.pack(side=tk.LEFT, padx=(0, 2))
        self._elev_levels_spin.bind('<Return>',   lambda _e: self._on_elev_levels_changed())
        self._elev_levels_spin.bind('<FocusOut>', lambda _e: self._on_elev_levels_changed())
        ttk.Label(rB, text='levels', foreground='grey').pack(side=tk.LEFT)

        # ----------------------------------------------------------------
        # Tab 2 — Generate
        # ----------------------------------------------------------------
        tab_gen = ttk.Frame(nb, padding=(8, 6))
        nb.add(tab_gen, text='  Generate  ')

        # Row A: Output dir | DPI | Every N
        gA = ttk.Frame(tab_gen)
        gA.pack(fill=tk.X, pady=2)
        ttk.Label(gA, text='Output dir:').pack(side=tk.LEFT)
        self._dir_var = tk.StringVar(value='_plot')
        ttk.Entry(gA, textvariable=self._dir_var, width=28).pack(side=tk.LEFT, padx=(2, 2))
        ttk.Button(gA, text='Browse...',
                   command=self._browse_dir).pack(side=tk.LEFT, padx=(0, 16))

        ttk.Label(gA, text='DPI:').pack(side=tk.LEFT)
        self._dpi_var = tk.StringVar(value='100')
        ttk.Entry(gA, textvariable=self._dpi_var, width=5).pack(side=tk.LEFT, padx=2)

        ttk.Separator(gA, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=12)
        ttk.Label(gA, text='Every N:').pack(side=tk.LEFT)
        self._stride_var = tk.StringVar(value='1')
        ttk.Spinbox(gA, from_=1, to=1000, increment=1,
                    textvariable=self._stride_var, width=5).pack(side=tk.LEFT, padx=2)
        ttk.Label(gA, text='frames', foreground='grey').pack(side=tk.LEFT, padx=(2, 0))

        # Row B: EPSG | Basemap | Alpha
        gB = ttk.Frame(tab_gen)
        gB.pack(fill=tk.X, pady=2)
        ttk.Label(gB, text='EPSG:').pack(side=tk.LEFT)
        self._epsg_var = tk.StringVar(value='')
        self._epsg_entry = ttk.Entry(gB, textvariable=self._epsg_var, width=8)
        self._epsg_entry.pack(side=tk.LEFT, padx=2)
        self._epsg_entry.bind('<Return>', lambda _e: self._on_set_epsg())
        ttk.Button(gB, text='Set', width=4,
                   command=self._on_set_epsg).pack(side=tk.LEFT, padx=(0, 12))

        from anuga.utilities.animate import BASEMAP_PROVIDERS
        self._basemap_var = tk.BooleanVar(value=False)
        self._basemap_chk = ttk.Checkbutton(gB, text='Basemap:',
                                             variable=self._basemap_var,
                                             command=self._on_basemap_toggle,
                                             state=tk.DISABLED)
        self._basemap_chk.pack(side=tk.LEFT, padx=(0, 2))
        self._basemap_provider_var = tk.StringVar(value='OpenStreetMap')
        self._basemap_provider_combo = ttk.Combobox(
            gB, textvariable=self._basemap_provider_var,
            values=list(BASEMAP_PROVIDERS.keys()),
            width=20, state=tk.DISABLED)
        self._basemap_provider_combo.pack(side=tk.LEFT, padx=(0, 12))

        ttk.Label(gB, text='Alpha:').pack(side=tk.LEFT)
        self._alpha_var = tk.DoubleVar(value=1.0)
        ttk.Spinbox(gB, from_=0.0, to=1.0, increment=0.05,
                    textvariable=self._alpha_var,
                    format='%.2f', width=5).pack(side=tk.LEFT, padx=2)

        # Row C: config file
        gC = ttk.Frame(tab_gen)
        gC.pack(fill=tk.X, pady=2)
        ttk.Button(gC, text='Save Config...', command=self._save_config).pack(side=tk.LEFT, padx=(0, 4))
        ttk.Button(gC, text='Load Config...', command=self._load_config).pack(side=tk.LEFT)

        # ----------------------------------------------------------------
        # Tab 3 — Output
        # ----------------------------------------------------------------
        tab_out = ttk.Frame(nb, padding=(8, 6))
        nb.add(tab_out, text='  Output  ')

        out_row = ttk.Frame(tab_out)
        out_row.pack(fill=tk.X, pady=4)
        self._save_frame_btn = ttk.Button(out_row, text='Save Frame...',
                                           command=self._save_frame,
                                           state=tk.DISABLED)
        self._save_frame_btn.pack(side=tk.LEFT, padx=(0, 4))
        self._export_frame_btn = ttk.Button(out_row, text='Export Frame...',
                                             command=self._export_frame,
                                             state=tk.DISABLED)
        self._export_frame_btn.pack(side=tk.LEFT, padx=4)
        self._save_anim_btn = ttk.Button(out_row, text='Save Animation...',
                                          command=self._save_animation,
                                          state=tk.DISABLED)
        self._save_anim_btn.pack(side=tk.LEFT, padx=4)

        ttk.Separator(out_row, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=12)
        self._mesh_btn = ttk.Button(out_row, text='View Mesh',
                                     command=self._show_mesh,
                                     state=tk.DISABLED)
        self._mesh_btn.pack(side=tk.LEFT, padx=4)
        self._save_mesh_btn = ttk.Button(out_row, text='Save Mesh...',
                                          command=self._save_mesh,
                                          state=tk.DISABLED)
        self._save_mesh_btn.pack(side=tk.LEFT, padx=4)

        # ================================================================
        # Generate Frames bar (always visible, below notebook)
        # ================================================================
        gen_bar = ttk.Frame(ctrl, padding=(0, 2, 0, 2))
        gen_bar.pack(fill=tk.X)

        self._gen_btn = ttk.Button(gen_bar, text='Generate Frames',
                                    command=self._start_generation,
                                    state=tk.DISABLED)
        self._gen_btn.pack(side=tk.LEFT, padx=(0, 4))
        self._cancel_btn = ttk.Button(gen_bar, text='Cancel',
                                       command=self._cancel_generation,
                                       state=tk.DISABLED)
        self._cancel_btn.pack(side=tk.LEFT, padx=(0, 8))

        self._progress_var = tk.IntVar(value=0)
        self._progress_bar = ttk.Progressbar(gen_bar, variable=self._progress_var,
                                              maximum=100)
        self._progress_bar.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 4))
        self._progress_label = ttk.Label(gen_bar, text='', width=14)
        self._progress_label.pack(side=tk.LEFT)

        ttk.Separator(ctrl, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=(4, 2))

        # ================================================================
        # Playback bar (always visible)
        # ================================================================
        pb = ttk.Frame(ctrl, padding=(0, 2))
        pb.pack(fill=tk.X)

        self._play_btn = ttk.Button(pb, text='Play', width=6,
                                     command=self._toggle_play)
        self._play_btn.pack(side=tk.LEFT, padx=(0, 2))
        ttk.Button(pb, text='|<', width=3, command=self._go_first).pack(side=tk.LEFT, padx=1)
        ttk.Button(pb, text='<',  width=3, command=self._step_back).pack(side=tk.LEFT, padx=1)
        ttk.Button(pb, text='>',  width=3, command=self._step_fwd).pack(side=tk.LEFT, padx=1)
        ttk.Button(pb, text='>|', width=3, command=self._go_last).pack(side=tk.LEFT, padx=1)
        ttk.Label(pb, text='FPS:').pack(side=tk.LEFT, padx=(8, 0))
        self._fps_var = tk.DoubleVar(value=5.0)
        ttk.Spinbox(pb, from_=0.5, to=30.0, increment=0.5,
                    textvariable=self._fps_var, width=5).pack(side=tk.LEFT, padx=2)

        ttk.Separator(pb, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)
        self._pick_btn = ttk.Button(pb, text='Pick timeseries',
                                     command=self._toggle_pick,
                                     state=tk.DISABLED)
        self._pick_btn.pack(side=tk.LEFT, padx=(0, 4))

        ttk.Separator(pb, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)
        self._xs_btn = ttk.Button(pb, text='Cross-section',
                                   command=self._toggle_xs_mode,
                                   state=tk.DISABLED)
        self._xs_btn.pack(side=tk.LEFT, padx=(0, 4))
        self._clear_xs_btn = ttk.Button(pb, text='Clear XS',
                                         command=self._clear_xs,
                                         state=tk.DISABLED)
        self._clear_xs_btn.pack(side=tk.LEFT, padx=(0, 4))

        ttk.Separator(pb, orient=tk.VERTICAL).pack(side=tk.LEFT, fill=tk.Y, padx=10)
        self._zoom_btn = ttk.Button(pb, text='Set Zoom',
                                     command=self._toggle_zoom_mode,
                                     state=tk.DISABLED)
        self._zoom_btn.pack(side=tk.LEFT, padx=(0, 4))
        self._reset_zoom_btn = ttk.Button(pb, text='Reset Zoom',
                                           command=self._reset_zoom,
                                           state=tk.DISABLED)
        self._reset_zoom_btn.pack(side=tk.LEFT, padx=(0, 4))

        # ---- Frame slider (always visible) ----
        slider_row = ttk.Frame(ctrl)
        slider_row.pack(fill=tk.X, pady=(4, 2))
        ttk.Label(slider_row, text='Frame:').pack(side=tk.LEFT)
        self._frame_label = ttk.Label(slider_row, text='-')
        self._frame_label.pack(side=tk.RIGHT, padx=10)
        self._slider_var = tk.IntVar(value=0)
        self._slider = ttk.Scale(slider_row, from_=0, to=0,
                                  variable=self._slider_var,
                                  orient=tk.HORIZONTAL,
                                  command=self._on_slider)
        self._slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=4)

        # ---- timeseries panel (hidden until a point is picked) ----
        self._ts_outer = ttk.Frame(self.root)
        # not packed yet

        ts_ctrl = ttk.Frame(self._ts_outer)
        ts_ctrl.pack(fill=tk.X)
        ttk.Label(ts_ctrl, text='Timeseries:').pack(side=tk.LEFT, padx=4)
        self._ts_qty_var = tk.StringVar(value='depth')
        ts_qty_combo = ttk.Combobox(ts_ctrl, textvariable=self._ts_qty_var,
                                    values=['depth', 'stage', 'speed',
                                            'speed_depth', 'elev'],
                                    width=12, state='readonly')
        ts_qty_combo.pack(side=tk.LEFT, padx=2)
        ts_qty_combo.bind('<<ComboboxSelected>>',
                          lambda _e: self._update_timeseries())
        self._ts_info_label = ttk.Label(ts_ctrl, text='', foreground='grey')
        self._ts_info_label.pack(side=tk.LEFT, padx=8)
        ttk.Button(ts_ctrl, text='Close',
                   command=self._close_timeseries).pack(side=tk.RIGHT, padx=4)
        ttk.Button(ts_ctrl, text='Clear picks',
                   command=self._close_timeseries).pack(side=tk.RIGHT, padx=4)
        ttk.Button(ts_ctrl, text='Export CSV',
                   command=self._export_timeseries).pack(side=tk.RIGHT, padx=4)

        # NB: do not set dpi here for the on-screen figures.  matplotlib's Tk
        # backend already scales the canvas by Tk's points-per-pixel (set from
        # the UI scale in tk_scaling), so passing a scaled dpi as well makes the
        # plots twice as large as intended.
        self._ts_fig, self._ts_ax = plt.subplots(figsize=(10, 1.8))
        self._ts_fig.tight_layout(pad=1.5)
        self._ts_canvas = FigureCanvasTkAgg(self._ts_fig, master=self._ts_outer)
        self._ts_canvas.draw()
        self._ts_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # ---- cross-section discharge panel (hidden until XS is defined) ----
        self._xs_outer = ttk.Frame(self.root)
        # not packed yet

        xs_ctrl = ttk.Frame(self._xs_outer)
        xs_ctrl.pack(fill=tk.X)
        ttk.Label(xs_ctrl, text='Cross-section discharge:').pack(side=tk.LEFT, padx=4)
        self._xs_info_lbl = ttk.Label(xs_ctrl, text='', foreground='grey')
        self._xs_info_lbl.pack(side=tk.LEFT, padx=8)
        ttk.Button(xs_ctrl, text='Close',
                   command=self._close_xs_panel).pack(side=tk.RIGHT, padx=4)
        ttk.Button(xs_ctrl, text='Export CSV',
                   command=self._export_xs).pack(side=tk.RIGHT, padx=4)

        self._xs_fig, self._xs_ax_plot = plt.subplots(figsize=(10, 1.8))
        self._xs_fig.tight_layout(pad=1.5)
        self._xs_fig_canvas = FigureCanvasTkAgg(self._xs_fig, master=self._xs_outer)
        self._xs_fig_canvas.draw()
        self._xs_fig_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # ---- matplotlib map canvas ----
        self._fig, self._ax = plt.subplots(figsize=(10, 6))
        self._ax.axis('off')
        self._fig.tight_layout(pad=0)

        self._canvas_frame = ttk.Frame(self.root)
        self._canvas_frame.pack(fill=tk.BOTH, expand=True)
        self._canvas = FigureCanvasTkAgg(self._fig, master=self._canvas_frame)
        self._canvas.draw()
        # ... and once more after the window has been laid out for real.
        self.root.after_idle(self._sync_canvas_sizes)
        self._canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self._hover_cid = self._canvas.mpl_connect(
            'motion_notify_event', self._on_hover)

    # -------------------------------------------------------------- #
    # SWW file selection                                              #
    # -------------------------------------------------------------- #

    def _browse_sww(self):
        path = filedialog.askopenfilename(
            title='Select SWW file',
            filetypes=[('SWW files', '*.sww'), ('All files', '*.*')],
            initialdir=os.path.dirname(self._sww_var.get()) or '.')
        if path:
            self._sww_var.set(path)
            self._on_sww_change()

    def _browse_dir(self):
        d = filedialog.askdirectory(
            title='Select output directory',
            initialdir=self._dir_var.get() or '.')
        if d:
            self._dir_var.set(d)

    def _on_sww_change(self):
        sww = self._sww_var.get().strip()
        if not sww or not os.path.isfile(sww):
            self._sww_info_label.config(text='File not found.', foreground='red')
            self._gen_btn.config(state=tk.DISABLED)
            self._pick_btn.config(state=tk.DISABLED)
            self._mesh_btn.config(state=tk.DISABLED)
            self._splotter = None
            return

        self._set_status(f'Loading {os.path.basename(sww)}...')
        self.root.update_idletasks()

        try:
            self._load_splotter(sww)
        except Exception as e:
            self._sww_info_label.config(text=f'Error: {e}', foreground='red')
            self._gen_btn.config(state=tk.DISABLED)
            self._set_status(f'Failed to load {sww}')
            return

        n = len(self._splotter.time)
        t0 = self._splotter.time[0]
        t1 = self._splotter.time[-1]
        self._sww_info_label.config(
            text=f'{n} timesteps  |  t={t0:.1f} .. {t1:.1f} s',
            foreground='grey')
        self._gen_btn.config(state=tk.NORMAL)
        self._pick_btn.config(state=tk.NORMAL)
        self._xs_btn.config(state=tk.NORMAL)
        self._mesh_btn.config(state=tk.NORMAL)
        self._set_status(f'Loaded {os.path.basename(sww)} -{n} timesteps')
        self._update_auto_limits()
        self._exit_pick_mode()
        self._close_timeseries()
        self._reset_zoom()
        self._clear_xs()

    def _load_splotter(self, sww):
        from anuga.utilities.animate import SWW_plotter

        self._trifinder = None   # invalidate cached trifinder for new mesh
        self._splotter = SWW_plotter(
            swwfile=sww,
            plot_dir=None,        # suppress automatic make_plot_dir
            min_depth=float(self._mindepth_var.get()))

        self._sww_prefix = self._splotter.name  # bare basename, no extension

        # Populate the EPSG entry from the file (may be empty for old SWW files)
        epsg_val = self._splotter.epsg
        self._epsg_var.set(str(epsg_val) if epsg_val is not None else '')

        self._register_tracer_quantities()

        self._refresh_basemap_state()

    def _register_tracer_quantities(self):
        """Add the open file's tracers to the quantity menu.

        Called after the plotter is built, which is where the tracer names and
        data become known. SWW_plotter has already bound a `tracer_<name>`
        attribute and a `save_tracer_<name>_frame` method for each, so the
        existing lookup tables are all that need filling in.
        """
        global _TRACER_QUANTITIES

        # Rebuild rather than append: opening a second file must not leave the
        # first file's tracers in the menu.
        for key in _TRACER_QUANTITIES:
            for table in (_QTY_DEFAULTS, _QTY_DATA_ATTR, _QTY_SAVE_METHOD,
                          _QTY_CBAR_LABEL, _QTY_TITLE):
                table.pop(key, None)
        _TRACER_QUANTITIES = []

        splotter = self._splotter
        for name in getattr(splotter, 'tracers', {}):
            for key, attr, label in (
                    ('tracer_' + name, 'tracer_' + name, name),
                    ('max_tracer_' + name, 'max_tracer_' + name, 'Max ' + name)):
                data = getattr(splotter, attr, None)
                if data is None:
                    continue          # no stored maximum for this tracer
                _TRACER_QUANTITIES.append(key)
                # A concentration's scale is the model's business -- sediment
                # runs at 1e-2, a normalised tracer at 1 -- so take the range
                # from the data rather than guessing a default.
                _QTY_DEFAULTS[key] = dict(vmin=float(data.min()),
                                          vmax=float(data.max()))
                _QTY_DATA_ATTR[key] = attr
                _QTY_SAVE_METHOD[key] = 'save_%s_frame' % key
                _QTY_CBAR_LABEL[key] = '%s (concentration)' % label
                _QTY_TITLE[key] = label

        if getattr(self, '_qty_combo', None) is not None:
            self._qty_combo.config(values=list(_all_quantities()))
            # A tracer selected for the previous file may not exist in this one.
            if self._qty_var.get() not in _all_quantities():
                self._qty_var.set('depth')

    def _refresh_basemap_state(self):
        """Enable/disable the basemap checkbox based on current EPSG + contextily."""
        if self._splotter is None:
            return
        has_epsg = self._splotter.epsg is not None
        has_contextily = False
        if has_epsg:
            try:
                import contextily  # noqa: F401
                has_contextily = True
            except ImportError:
                pass
        if has_epsg and has_contextily:
            self._basemap_chk.config(state=tk.NORMAL)
            self._basemap_var.set(True)
            self._on_basemap_toggle()
        elif has_epsg:
            self._basemap_chk.config(state=tk.DISABLED)
            self._basemap_var.set(False)
            self._basemap_provider_combo.config(state=tk.DISABLED)
        else:
            self._basemap_chk.config(state=tk.DISABLED)
            self._basemap_var.set(False)
            self._basemap_provider_combo.config(state=tk.DISABLED)

    def _on_set_epsg(self):
        """Apply the EPSG code typed in the entry field to the loaded SWW_plotter."""
        if self._splotter is None:
            self._set_status('Load an SWW file first.')
            return
        raw = self._epsg_var.get().strip()
        if raw == '':
            self._splotter.set_epsg(None)
            self._refresh_basemap_state()
            self._set_status('EPSG cleared - basemap disabled.')
            return
        try:
            code = int(raw)
        except ValueError:
            self._set_status(f'Invalid EPSG code "{raw}" - enter an integer, e.g. 32756.')
            return
        self._splotter.set_epsg(code)
        self._refresh_basemap_state()
        has_contextily = False
        try:
            import contextily  # noqa: F401
            has_contextily = True
        except ImportError:
            pass
        if has_contextily:
            self._set_status(f'EPSG set to {code} - basemap enabled.')
        else:
            self._set_status(
                f'EPSG set to {code} - install contextily for basemap support.')

    def _on_basemap_toggle(self):
        if self._basemap_var.get():
            self._alpha_var.set(0.6)
            self._basemap_provider_combo.config(state='readonly')
        else:
            self._alpha_var.set(1.0)
            self._basemap_provider_combo.config(state=tk.DISABLED)

    def _on_qty_change(self):
        qty = self._qty_var.get()
        defaults = _QTY_DEFAULTS[qty]
        self._vmin_var.set(str(defaults['vmin']))
        self._vmax_var.set(str(defaults['vmax']))
        if self._auto_var.get():
            self._update_auto_limits()

    def _on_auto_toggle(self):
        if self._auto_var.get():
            self._update_auto_limits()

    def _update_auto_limits(self):
        if self._splotter is None or not self._auto_var.get():
            return
        qty = self._qty_var.get()
        attr = _QTY_DATA_ATTR[qty]
        data = getattr(self._splotter, attr)
        if data is None:
            return
        if qty == 'elev_delta':
            # symmetric range around zero
            import numpy as np
            max_abs = float(max(abs(float(data.min())), abs(float(data.max()))))
            self._vmin_var.set(f'{-max_abs:.4g}')
            self._vmax_var.set(f'{max_abs:.4g}')
        else:
            self._vmin_var.set(f'{float(data.min()):.4g}')
            self._vmax_var.set(f'{float(data.max()):.4g}')

    def _is_static_qty(self, qty):
        """Return True if qty produces a single frame (not time-animated).

        max_* quantities always produce one frame.  elev and elev_delta produce
        one frame when elevation is static (1-D); they are animated when the SWW
        contains time-varying elevation (2-D, e.g. from an erosion run).
        """
        if qty.startswith('max_'):
            return True
        if qty in ('elev', 'elev_delta') and self._splotter is not None:
            return self._splotter.elev.ndim == 1
        return False

    # -------------------------------------------------------------- #
    # Frame generation                                                #
    # -------------------------------------------------------------- #

    def _start_generation(self):
        if self._splotter is None:
            return

        # Stop playback and discard frame list before deleting files so that
        # _advance cannot try to read a file that is about to be removed.
        self._stop_playback()
        self._frames = []

        # resolve output directory
        plot_dir = self._dir_var.get().strip() or '_plot'
        plot_dir = os.path.abspath(plot_dir)
        os.makedirs(plot_dir, exist_ok=True)

        qty   = self._qty_var.get()

        # Remove stale frames from a previous generation
        for old_file in _find_frames(plot_dir, qty, self._sww_prefix):
            try:
                os.remove(old_file)
            except OSError:
                pass

        try:
            vmin   = float(self._vmin_var.get())
            vmax   = float(self._vmax_var.get())
            dpi    = int(self._dpi_var.get())
            stride = max(1, int(self._stride_var.get()))
        except ValueError as e:
            self._set_status(f'Invalid parameter: {e}')
            return

        cmap = self._cmap_var.get().strip() or 'viridis'
        if self._cmap_reverse_var.get() and not cmap.endswith('_r'):
            cmap = cmap + '_r'
        basemap = self._basemap_var.get()
        alpha   = max(0.0, min(1.0, self._alpha_var.get()))

        from anuga.utilities.animate import BASEMAP_PROVIDERS, BASEMAP_DEFAULT
        provider_label = self._basemap_provider_var.get()
        basemap_provider = BASEMAP_PROVIDERS.get(provider_label, BASEMAP_DEFAULT)

        n_total = len(self._splotter.time)
        if self._is_static_qty(qty):
            # Static quantities (max_*, static elev) collapse to one image
            sww_frames = [0]
            stride_msg = ''
        else:
            sww_frames = list(range(0, n_total, stride))
            stride_msg = f' (every {stride})' if stride > 1 else ''
        n_to_gen = len(sww_frames)
        self._progress_bar.config(maximum=n_to_gen)
        self._progress_var.set(0)
        self._progress_label.config(text=f'0 / {n_to_gen}')
        self._gen_btn.config(state=tk.DISABLED, text='Generate Frames')
        self._cancel_btn.config(state=tk.NORMAL)
        self._cancel_flag = False
        self._set_status(
            f'Generating {n_to_gen} {qty} frames{stride_msg}...')
        self.root.update_idletasks()

        # Reset plotter frame counters and re-point it at the output dir
        self._splotter.plot_dir = plot_dir
        for attr in ('_depth_frame_count', '_stage_frame_count',
                     '_speed_frame_count', '_speed_depth_frame_count',
                     '_max_depth_frame_count', '_max_stage_frame_count',
                     '_max_speed_frame_count',
                     '_max_speed_depth_frame_count', '_elev_frame_count',
                     '_elev_delta_frame_count'):
            setattr(self._splotter, attr, 0)
        # Close any figures cached from a previous generation run
        self._splotter._clear_figure_cache()

        self._gen_used_basemap = basemap
        self._last_gen_dpi = dpi
        self._last_gen_vmin = vmin
        self._last_gen_vmax = vmax
        self._last_gen_cmap = cmap
        self._last_gen_qty = qty

        show_edges = self._show_edges_var.get()
        self._last_gen_show_edges = show_edges
        smooth = not show_edges

        show_elev = self._show_elev_var.get()
        try:
            elev_levels = max(3, int(self._elev_levels_var.get()))
        except (ValueError, AttributeError):
            elev_levels = 10
        show_mesh = self._show_mesh_var.get()

        self._last_gen_show_elev = show_elev
        self._last_gen_elev_levels = elev_levels
        self._last_gen_show_mesh = show_mesh

        # Single frame (max_* quantities) always runs sequentially.
        # For multi-frame quantities, attempt parallel generation and fall
        # back to sequential if multiprocessing is unavailable.
        if n_to_gen <= 1:
            save_method = getattr(self._splotter, _QTY_SAVE_METHOD[qty])
            self._generate_next_frame(0, sww_frames, save_method, dpi, vmin, vmax,
                                      plot_dir, qty, cmap, basemap, alpha,
                                      basemap_provider, smooth=smooth,
                                      show_elev=show_elev,
                                      elev_levels=elev_levels,
                                      show_mesh=show_mesh)
        else:
            try:
                self._start_parallel_generation(
                    sww_frames, plot_dir, qty, dpi, vmin, vmax,
                    cmap, basemap, alpha, basemap_provider,
                    smooth=smooth, show_elev=show_elev,
                    elev_levels=elev_levels, show_mesh=show_mesh)
            except Exception as e:
                import traceback
                print(f'[anuga_animate] parallel init failed: {e}',
                      file=sys.stderr, flush=True)
                traceback.print_exc(file=sys.stderr)
                self._set_status(
                    f'Parallel init failed ({e}); falling back to sequential.')
                save_method = getattr(self._splotter, _QTY_SAVE_METHOD[qty])
                self._generate_next_frame(0, sww_frames, save_method, dpi, vmin,
                                          vmax, plot_dir, qty, cmap, basemap,
                                          alpha, basemap_provider,
                                          smooth=smooth,
                                          show_elev=show_elev,
                                          elev_levels=elev_levels,
                                          show_mesh=show_mesh)

    def _generate_next_frame(self, pos, sww_frames, save_method,
                              dpi, vmin, vmax, plot_dir, qty,
                              cmap='viridis', basemap=False, alpha=1.0,
                              basemap_provider='OpenStreetMap.Mapnik',
                              smooth=False, show_elev=False, elev_levels=10,
                              show_mesh=False):
        n_to_gen  = len(sww_frames)
        sww_frame = sww_frames[pos]

        if self._cancel_flag:
            self._set_status(f'Generation cancelled after {pos} frames.')
            self._gen_btn.config(state=tk.NORMAL)
            self._cancel_btn.config(state=tk.DISABLED)
            return
        try:
            save_method(frame=sww_frame, dpi=dpi, vmin=vmin, vmax=vmax,
                        cmap=cmap, basemap=basemap, alpha=alpha,
                        basemap_provider=basemap_provider,
                        xlim=self._zoom_xlim, ylim=self._zoom_ylim,
                        smooth=smooth, show_elev=show_elev,
                        elev_levels=elev_levels, show_mesh=show_mesh)
        except Exception as e:
            self._set_status(f'Error generating frame {sww_frame}: {e}')
            self._gen_btn.config(state=tk.NORMAL)
            self._cancel_btn.config(state=tk.DISABLED)
            return

        self._progress_var.set(pos + 1)
        self._progress_label.config(text=f'{pos + 1} / {n_to_gen}')
        self.root.update_idletasks()

        if pos + 1 < n_to_gen:
            self._gen_after_id = self.root.after(
                1, lambda: self._generate_next_frame(
                    pos + 1, sww_frames, save_method,
                    dpi, vmin, vmax, plot_dir, qty, cmap, basemap, alpha,
                    basemap_provider, smooth=smooth,
                    show_elev=show_elev, elev_levels=elev_levels,
                    show_mesh=show_mesh))
        else:
            self._on_generation_done(plot_dir, qty, n_to_gen)

    def _start_parallel_generation(self, sww_frames, plot_dir, qty,
                                     dpi, vmin, vmax, cmap, basemap,
                                     alpha, basemap_provider,
                                     smooth=False, show_elev=False,
                                     elev_levels=10, show_mesh=False):
        """Spawn a ProcessPoolExecutor and submit all frames."""
        import multiprocessing as mp
        import platform
        from concurrent.futures import ProcessPoolExecutor
        from anuga.utilities._animate_worker import worker_init, worker_frame

        sww_path  = self._sww_var.get()
        min_depth = float(self._mindepth_var.get())
        epsg      = self._splotter.epsg
        n_workers = max(1, min(os.cpu_count() or 4, len(sww_frames), 4))

        # 'fork' on Linux: workers inherit parent memory, no reimport overhead.
        # 'spawn' on Windows/macOS: safer for GUI apps.
        ctx_name = 'fork' if platform.system() == 'Linux' else 'spawn'
        ctx = mp.get_context(ctx_name)
        self._executor = ProcessPoolExecutor(
            max_workers=n_workers,
            mp_context=ctx,
            initializer=worker_init,
            initargs=(sww_path, plot_dir, min_depth, epsg))

        save_name = _QTY_SAVE_METHOD[qty]
        self._futures = [
            self._executor.submit(
                worker_frame,
                frame, pos, save_name, qty,
                dpi, vmin, vmax, cmap, basemap, alpha, basemap_provider,
                self._zoom_xlim, self._zoom_ylim,
                smooth, show_elev, elev_levels, show_mesh)
            for pos, frame in enumerate(sww_frames)
        ]
        self._gen_plot_dir = plot_dir
        self._gen_qty      = qty
        self._n_to_gen     = len(sww_frames)
        self._set_status(
            f'Generating {self._n_to_gen} {qty} frames '
            f'({n_workers} workers)...')
        self._poll_generation()

    def _poll_generation(self):
        """Check parallel futures every 200 ms; update the progress bar."""
        if self._cancel_flag:
            # _cancel_generation already shut down the executor.
            return

        n_done = sum(1 for f in self._futures if f.done())
        self._progress_var.set(n_done)
        self._progress_label.config(text=f'{n_done} / {self._n_to_gen}')
        self.root.update_idletasks()

        # Surface first worker exception immediately.
        for fut in self._futures:
            if fut.done() and fut.exception() is not None:
                exc = fut.exception()
                for f in self._futures:
                    f.cancel()
                self._executor.shutdown(wait=False, cancel_futures=True)
                self._executor = None
                self._futures = []
                self._gen_btn.config(state=tk.NORMAL)
                self._cancel_btn.config(state=tk.DISABLED)
                self._set_status(f'Frame generation error: {exc}')
                return

        if n_done == self._n_to_gen:
            self._executor.shutdown(wait=False)
            self._executor = None
            self._futures = []
            self._on_generation_done(
                self._gen_plot_dir, self._gen_qty, self._n_to_gen)
        else:
            self._gen_after_id = self.root.after(200, self._poll_generation)

    def _cancel_generation(self):
        self._cancel_flag = True
        if self._gen_after_id is not None:
            self.root.after_cancel(self._gen_after_id)
            self._gen_after_id = None
        if self._executor is not None:
            for fut in self._futures:
                fut.cancel()
            self._executor.shutdown(wait=False, cancel_futures=True)
            self._executor = None
            self._futures = []
        self._gen_btn.config(state=tk.NORMAL)
        self._cancel_btn.config(state=tk.DISABLED)
        self._set_status('Generation cancelled.')

    def _on_generation_done(self, plot_dir, qty, n_frames):
        self._progress_var.set(n_frames)
        self._progress_label.config(text=f'{n_frames} / {n_frames}')
        self._gen_btn.config(state=tk.NORMAL)
        self._cancel_btn.config(state=tk.DISABLED)
        prefix = self._sww_prefix
        self._set_status(
            f'Done -{n_frames} frames saved to {plot_dir}')
        self._compute_plot_transform(self._last_gen_dpi)
        self._load_frames(plot_dir, qty, prefix)

    def _compute_plot_transform(self, dpi):
        """Render frame 0 off-screen (Agg) with the same params used to generate
        frames so that the resulting axes position/limits exactly match the PNGs."""
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        sp = self._splotter
        if sp is None:
            self._plot_transform = None
            return
        try:
            import numpy as np
            # For max quantities use the actual max; others use frame 0
            if self._last_gen_qty.startswith('max_'):
                depth = sp.max_depth
            else:
                depth = sp.depth[0, :]
            try:
                elev = sp.elev[0, :]
            except (IndexError, TypeError):
                elev = sp.elev

            vmin = self._last_gen_vmin
            vmax = self._last_gen_vmax
            cmap = self._last_gen_cmap
            use_basemap = self._gen_used_basemap and (sp.epsg is not None)
            triang = sp.triang_abs if use_basemap else sp.triang

            fig = Figure(figsize=(10, 6), dpi=dpi)
            FigureCanvasAgg(fig)
            ax = fig.add_subplot(111)

            # Pin axes to full mesh extent first — matches _animated_frame so
            # that pos/xlim/ylim here are identical to the saved PNGs.
            ax.set_xlim(triang.x.min(), triang.x.max())
            ax.set_ylim(triang.y.min(), triang.y.max())
            if not use_basemap:
                triang.set_mask(depth > sp.min_depth)
                ax.tripcolor(triang, facecolors=elev, cmap='Greys_r')
            triang.set_mask(depth < sp.min_depth)
            im = ax.tripcolor(triang, facecolors=depth, cmap=cmap,
                              vmin=vmin, vmax=vmax)
            triang.set_mask(None)
            ax.set_aspect('equal')
            ax.set_xlabel('Easting (m)')
            ax.set_ylabel('Northing (m)')
            if self._zoom_xlim is not None:
                ax.set_xlim(self._zoom_xlim)
            if self._zoom_ylim is not None:
                ax.set_ylim(self._zoom_ylim)
            fig.colorbar(im, ax=ax)

            if use_basemap:
                from anuga.utilities.animate import _add_basemap, BASEMAP_PROVIDERS
                provider_label = self._basemap_provider_var.get()
                provider_str = BASEMAP_PROVIDERS.get(
                    provider_label, 'OpenStreetMap.Mapnik')
                _add_basemap(ax, sp.epsg, provider_str, cache=sp._basemap_cache)

            fig.canvas.draw()

            self._plot_transform = dict(
                pos=ax.get_position(),
                xlim=ax.get_xlim(),
                ylim=ax.get_ylim(),
                W=int(fig.get_figwidth() * fig.get_dpi()),
                H=int(fig.get_figheight() * fig.get_dpi()),
                use_basemap=use_basemap,
            )
            self._elev_contour_data = None    # transform changed — invalidate cache
            self._update_mesh_overlay()
            self._update_elev_overlay()
            self._update_xs_overlay()
        except Exception:
            self._plot_transform = None

    # -------------------------------------------------------------- #
    # Frame loading and display                                       #
    # -------------------------------------------------------------- #

    def _load_frames(self, plot_dir, qty, prefix):
        self._stop_playback()
        self._frames = _find_frames(plot_dir, qty, prefix)
        n = len(self._frames)
        if n == 0:
            self._set_status(f'No frames found for {prefix}_{qty}_*.png in {plot_dir}')
            return
        self._slider.configure(to=n - 1)
        self._slider_var.set(0)
        self._current = 0
        # Remove any zoom rectangle — new frames already reflect the zoom
        self._remove_zoom_patch()
        self._show_frame(0)
        self._save_frame_btn.config(state=tk.NORMAL)
        self._export_frame_btn.config(state=tk.NORMAL)
        self._save_anim_btn.config(state=tk.NORMAL)
        self._zoom_btn.config(state=tk.NORMAL)
        self._save_mesh_btn.config(state=tk.NORMAL)
        self._show_mesh_chk.config(state=tk.NORMAL)
        self._show_elev_chk.config(state=tk.NORMAL)
        self._elev_levels_spin.config(state=tk.NORMAL)
        self._elev_contour_data = None    # invalidate on new generation
        self._update_xs_overlay()
        self._set_status(f'Loaded {n} frames  |  {plot_dir}')

    def _sync_canvas_sizes(self):
        """Match each embedded figure to its widget's real pixel size.

        matplotlib's Tk backend sizes a figure from the first <Configure> it
        receives, but on a HiDPI display that arrives before the canvas knows
        its device-pixel-ratio, so the figure comes out ui_scale times larger
        than the widget: the image then renders to an oversized canvas and only
        its corner is visible, un-centred and at the wrong scale.  Any later
        resize corrects it, which is why nudging the window "fixed" it.

        Re-assert the size once the widget geometry is real.  Safe to call
        repeatedly; a widget that is not laid out yet (1x1) is skipped.
        """
        for fig, canvas in ((self._fig, self._canvas),
                            (self._ts_fig, self._ts_canvas),
                            (self._xs_fig, self._xs_fig_canvas)):
            try:
                w = canvas.get_tk_widget()
                wpx, hpx = w.winfo_width(), w.winfo_height()
                if wpx <= 1 or hpx <= 1:
                    continue
                dpi = fig.get_dpi()
                if dpi <= 0:
                    continue
                want_w, want_h = wpx / dpi, hpx / dpi
                have_w, have_h = fig.get_size_inches()
                if (abs(want_w - have_w) > 0.01 or abs(want_h - have_h) > 0.01):
                    fig.set_size_inches(want_w, want_h, forward=False)
                    canvas.draw_idle()
            except (tk.TclError, AttributeError):
                continue

    def _show_frame(self, idx):
        if not self._frames:
            return
        idx = max(0, min(idx, len(self._frames) - 1))
        self._current = idx
        self._slider_var.set(idx)
        self._frame_label.config(text=f'Frame {idx + 1} / {len(self._frames)}')

        try:
            img = mpimage.imread(self._frames[idx])
        except FileNotFoundError:
            self._set_status('Frame file missing - click Generate Frames to rebuild.')
            self._stop_playback()
            return
        if self._im is None:
            self._im = self._ax.imshow(img, aspect='equal')
            self._im.set_extent([0, img.shape[1], img.shape[0], 0])
            # First image: the canvas may still be carrying the oversized
            # figure from its initial <Configure> (see _sync_canvas_sizes).
            self.root.after_idle(self._sync_canvas_sizes)
        else:
            self._im.set_data(img)
            self._im.set_extent([0, img.shape[1], img.shape[0], 0])
        self._update_pick_overlay()
        self._canvas.draw_idle()
        self._update_ts_cursor()

    def _on_slider(self, val):
        try:
            idx = int(float(val))
        except ValueError:
            return
        if idx != self._current:
            self._stop_playback()
            self._show_frame(idx)

    # -------------------------------------------------------------- #
    # Playback controls                                               #
    # -------------------------------------------------------------- #

    def _toggle_play(self):
        if self._playing:
            self._stop_playback()
        else:
            self._start_playback()

    def _start_playback(self):
        if not self._frames:
            return
        self._playing = True
        self._play_btn.config(text='Pause')
        self._schedule_next()

    def _stop_playback(self):
        self._playing = False
        self._play_btn.config(text='Play')
        if self._after_id is not None:
            self.root.after_cancel(self._after_id)
            self._after_id = None

    def _schedule_next(self):
        if not self._playing:
            return
        try:
            fps = max(0.1, self._fps_var.get())
        except tk.TclError:
            fps = 5.0
        self._after_id = self.root.after(int(1000 / fps), self._advance)

    def _advance(self):
        if not self._playing:
            return
        next_idx = (self._current + 1) % len(self._frames)
        self._show_frame(next_idx)
        self._schedule_next()

    def _step_fwd(self):
        self._stop_playback()
        self._show_frame(self._current + 1)

    def _step_back(self):
        self._stop_playback()
        self._show_frame(self._current - 1)

    def _go_first(self):
        self._stop_playback()
        self._show_frame(0)

    def _go_last(self):
        self._stop_playback()
        self._show_frame(len(self._frames) - 1)

    # -------------------------------------------------------------- #
    # Timeseries                                                      #
    # -------------------------------------------------------------- #

    def _toggle_pick(self):
        if self._pick_mode:
            self._exit_pick_mode()
        else:
            self._enter_pick_mode()

    def _enter_pick_mode(self):
        """Enter pick mode: keep the current animation frame as-is, overlay an
        instruction banner, and connect click events.  The frame image does not
        change size or appearance."""
        if self._splotter is None or not self._frames:
            return
        if self._plot_transform is None:
            self._set_status('Computing pick transform...')
            self.root.update_idletasks()
            self._compute_plot_transform(self._last_gen_dpi)
        if self._plot_transform is None:
            self._set_status('Pick mode unavailable: transform could not be computed')
            return

        # Overlay a semi-transparent instruction banner using axes-fraction coords
        # so it is independent of image size.
        self._pick_text = self._ax.text(
            0.5, 0.02,
            'Click to add timeseries points  |  Esc to finish',
            ha='center', va='bottom', fontsize=9,
            color='white', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.65),
            zorder=20,
            transform=self._ax.transAxes)

        self._canvas.draw_idle()
        self._canvas.get_tk_widget().config(cursor='crosshair')
        self._canvas.get_tk_widget().focus_set()  # needed for key_press_event (Esc)

        self._pick_cid = self._canvas.mpl_connect(
            'button_press_event', self._on_pick_click)
        self._pick_key_cid = self._canvas.mpl_connect(
            'key_press_event', self._on_pick_key)

        self._pick_mode = True
        self._pick_btn.config(text='Cancel pick')

    def _on_pick_click(self, event):
        if event.inaxes != self._ax or event.xdata is None:
            return
        import numpy as np
        sp = self._splotter

        xmesh, ymesh = self._imshow_to_mesh(event.xdata, event.ydata)

        # Centroids are stored in relative coords; for basemap the xlim/ylim
        # from the transform already include the xllcorner/yllcorner offset.
        use_basemap = self._gen_used_basemap and (sp.epsg is not None)
        if use_basemap:
            xc = sp.xc + sp.xllcorner
            yc = sp.yc + sp.yllcorner
        else:
            xc, yc = sp.xc, sp.yc

        dist = np.sqrt((xc - xmesh)**2 + (yc - ymesh)**2)
        tri = int(dist.argmin())
        if tri not in self._ts_triangles:
            self._ts_triangles.append(tri)
        # Stay in pick mode so the user can keep clicking to add more points.
        # _update_timeseries may pack the ts panel (resize), so force a
        # synchronous draw afterwards to ensure the markers are visible.
        self._update_timeseries()
        self._update_pick_overlay()
        self._canvas.draw()

    def _on_pick_key(self, event):
        if event.key == 'escape':
            self._exit_pick_mode()

    def _exit_pick_mode(self):
        """Disconnect pick events and remove the instruction banner."""
        if not self._pick_mode:
            return
        for cid in (self._pick_cid, self._pick_key_cid):
            if cid is not None:
                self._canvas.mpl_disconnect(cid)
        self._pick_cid = None
        self._pick_key_cid = None
        self._pick_mode = False
        self._pick_btn.config(text='Pick timeseries')
        self._canvas.get_tk_widget().config(cursor='')
        if self._pick_text is not None:
            try:
                self._pick_text.remove()
            except Exception:
                pass
            self._pick_text = None
        # Re-show the current frame: reapplies the correct imshow extent/limits
        # before adding the pick overlay (plain draw_idle is not enough because
        # ax.plot() in _update_pick_overlay can corrupt the axes limits).
        if self._frames:
            self._show_frame(self._current)
        else:
            self._canvas.draw_idle()

    # tab10 colour cycle (matches matplotlib default)
    _TS_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                  '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    @classmethod
    def _ts_color(cls, i):
        return cls._TS_COLORS[i % len(cls._TS_COLORS)]

    def _update_timeseries(self):
        """Plot quantity vs time for all picked centroids."""
        if not self._ts_triangles or self._splotter is None:
            return
        import numpy as _np

        sp  = self._splotter
        qty = self._ts_qty_var.get()
        t   = sp.time

        data_map = {
            'depth':       sp.depth,
            'stage':       sp.stage,
            'speed':       sp.speed,
            'speed_depth': sp.speed_depth,
        }
        ylabel_map = {
            'depth':       'Depth (m)',
            'stage':       'Stage (m)',
            'speed':       'Speed (m/s)',
            'speed_depth': 'Speed x Depth (m2/s)',
            'elev':        'Elevation (m)',
        }

        self._ts_ax.cla()
        for i, tri in enumerate(self._ts_triangles):
            if qty == 'elev':
                y = sp.elev[:, tri] if sp.elev.ndim == 2 else _np.full(len(t), sp.elev[tri])
            else:
                y = data_map[qty][:, tri]
            xc = sp.xc[tri] + sp.xllcorner
            yc = sp.yc[tri] + sp.yllcorner
            label = f'T{tri} ({xc:.0f}, {yc:.0f})'
            self._ts_ax.plot(t, y, color=self._ts_color(i),
                             linewidth=1.2, label=label)

        self._ts_ax.set_xlabel('Time (s)')
        self._ts_ax.set_ylabel(ylabel_map[qty])
        self._ts_ax.set_xlim(t[0], t[-1])
        self._ts_ax.grid(True, linestyle=':', alpha=0.5)
        if len(self._ts_triangles) > 1:
            self._ts_ax.legend(fontsize=7, loc='upper right', framealpha=0.7)

        # Vertical cursor at current animation time
        current_time = t[0]
        if self._frames and self._current < len(self._frames):
            frac = self._current / max(len(self._frames) - 1, 1)
            ts_idx = int(round(frac * (len(t) - 1)))
            current_time = t[ts_idx]

        self._ts_vline = self._ts_ax.axvline(current_time,
                                              color='red', linewidth=1.0,
                                              linestyle='--')
        self._ts_fig.tight_layout(pad=1.5)
        self._ts_canvas.draw_idle()

        n = len(self._ts_triangles)
        self._ts_info_label.config(
            text=f'{n} point{"s" if n != 1 else ""} picked')

        # Show the panel if not already visible.
        if not self._ts_outer.winfo_ismapped():
            self._ts_outer.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=False,
                                before=self._canvas_frame)

    def _export_timeseries(self):
        """Save all picked timeseries to a single CSV file."""
        if not self._ts_triangles or self._splotter is None:
            return
        import csv
        import numpy as _np2
        from tkinter import filedialog

        sp  = self._splotter
        qty = self._ts_qty_var.get()
        t   = sp.time
        tris = self._ts_triangles

        data_map = {
            'depth':       sp.depth,
            'stage':       sp.stage,
            'speed':       sp.speed,
            'speed_depth': sp.speed_depth,
        }
        qty_label = {
            'depth':       'depth_m',
            'stage':       'stage_m',
            'speed':       'speed_m_s',
            'speed_depth': 'speed_depth_m2_s',
            'elev':        'elevation_m',
        }[qty]

        # Build array of y-values: shape (n_times, n_tris)
        ys = []
        for tri in tris:
            if qty == 'elev':
                y = sp.elev[:, tri] if sp.elev.ndim == 2 else _np2.full(len(t), sp.elev[tri])
            else:
                y = data_map[qty][:, tri]
            ys.append(y)

        tri_str = '_'.join(str(tri) for tri in tris)
        default_name = f'{self._sww_prefix}_{qty}_tri{tri_str}.csv'
        path = filedialog.asksaveasfilename(
            title='Export timeseries',
            initialfile=default_name,
            defaultextension='.csv',
            filetypes=[('CSV files', '*.csv'), ('All files', '*.*')])
        if not path:
            return

        with open(path, 'w', newline='') as f:
            writer = csv.writer(f)
            # header comment rows — one per triangle
            for tri in tris:
                xc = sp.xc[tri] + sp.xllcorner
                yc = sp.yc[tri] + sp.yllcorner
                writer.writerow([f'# triangle: {tri}',
                                 f'x: {xc:.3f}', f'y: {yc:.3f}'])
            # column headers
            col_headers = ['time_s'] + [f'{qty_label}_tri{tri}' for tri in tris]
            writer.writerow(col_headers)
            for row_i, ti in enumerate(t):
                row = [f'{ti:.6g}'] + [f'{ys[j][row_i]:.6g}' for j in range(len(tris))]
                writer.writerow(row)

        self._set_status(f'Exported {len(t)} rows x {len(tris)} points to {path}')

    def _update_ts_cursor(self):
        """Move the vertical cursor line to the current animation time."""
        if self._splotter is None:
            return
        t = self._splotter.time
        if len(self._frames) > 1:
            frac = self._current / (len(self._frames) - 1)
        else:
            frac = 0.0
        ts_idx = int(round(frac * (len(t) - 1)))
        t_cur = t[ts_idx]
        if self._ts_vline is not None:
            self._ts_vline.set_xdata([t_cur, t_cur])
            self._ts_canvas.draw_idle()
        if self._xs_vline is not None:
            self._xs_vline.set_xdata([t_cur, t_cur])
            self._xs_fig_canvas.draw_idle()

    def _close_timeseries(self):
        """Hide the timeseries panel and clear all picked points."""
        self._ts_triangles = []
        self._ts_vline = None
        self._ts_ax.cla()
        self._ts_canvas.draw_idle()
        self._ts_info_label.config(text='')
        if self._ts_outer.winfo_ismapped():
            self._ts_outer.pack_forget()
        self._remove_pick_overlay()
        self._canvas.draw_idle()

    def _remove_pick_overlay(self):
        """Remove all centroid markers from the animation canvas."""
        for artist in self._pick_overlays:
            try:
                artist.remove()
            except Exception:
                pass
        self._pick_overlays = []

    def _update_pick_overlay(self):
        """Add/update one coloured star per picked centroid on the animation canvas."""
        self._remove_pick_overlay()
        if (not self._ts_triangles or self._plot_transform is None
                or not self._frames):
            return

        sp  = self._splotter
        pt  = self._plot_transform
        use_basemap = self._gen_used_basemap and (sp.epsg is not None)
        xlim, ylim = pt['xlim'], pt['ylim']
        pos = pt['pos']
        W, H = pt['W'], pt['H']

        for i, tri in enumerate(self._ts_triangles):
            if use_basemap:
                xd = sp.xc[tri] + sp.xllcorner
                yd = sp.yc[tri] + sp.yllcorner
            else:
                xd, yd = sp.xc[tri], sp.yc[tri]

            xfrac = (xd - xlim[0]) / (xlim[1] - xlim[0])
            yfrac = (yd - ylim[0]) / (ylim[1] - ylim[0])
            ax_x = (pos.x0 + pos.width  * xfrac) * W
            ax_y = (1.0 - (pos.y0 + pos.height * yfrac)) * H

            marker, = self._ax.plot(
                ax_x, ax_y, '*', markersize=9, zorder=10,
                color=self._ts_color(i),
                markeredgecolor='white', markeredgewidth=0.5,
                scalex=False, scaley=False)
            self._pick_overlays.append(marker)

    # -------------------------------------------------------------- #
    # Cross-section discharge                                         #
    # -------------------------------------------------------------- #

    def _toggle_xs_mode(self):
        if self._xs_mode:
            self._exit_xs_mode()
        else:
            self._enter_xs_mode()

    def _enter_xs_mode(self):
        if self._splotter is None:
            self._set_status('Load an SWW file first.')
            return
        self._exit_pick_mode()   # mutually exclusive with pick mode
        self._xs_pts = []
        self._remove_xs_overlay()
        self._xs_mode = True
        self._xs_btn.config(text='Cancel XS')
        self._canvas.get_tk_widget().config(cursor='crosshair')
        self._canvas.get_tk_widget().focus_set()
        self._xs_cid = self._canvas.mpl_connect(
            'button_press_event', self._on_xs_click)
        self._xs_key_cid = self._canvas.mpl_connect(
            'key_press_event', self._on_xs_key)
        self._set_status('Click to set cross-section start point.')

    def _exit_xs_mode(self):
        if not self._xs_mode:
            return
        for cid in (self._xs_cid, self._xs_key_cid):
            if cid is not None:
                self._canvas.mpl_disconnect(cid)
        self._xs_cid = None
        self._xs_key_cid = None
        self._xs_mode = False
        self._xs_btn.config(text='Cross-section')
        self._canvas.get_tk_widget().config(cursor='')

    def _on_xs_click(self, event):
        if event.inaxes != self._ax or event.xdata is None:
            return
        sp = self._splotter
        xm, ym = self._imshow_to_mesh(event.xdata, event.ydata)
        # Store in relative (mesh) coords
        use_basemap = self._gen_used_basemap and (sp.epsg is not None)
        if use_basemap:
            xm -= sp.xllcorner
            ym -= sp.yllcorner
        self._xs_pts.append((xm, ym))
        self._update_xs_overlay()
        self._canvas.draw()
        if len(self._xs_pts) == 1:
            self._set_status('Start point set.  Click to set end point.')
        elif len(self._xs_pts) >= 2:
            self._exit_xs_mode()
            self._compute_xs_flux()

    def _on_xs_key(self, event):
        if event.key == 'escape':
            self._xs_pts = []
            self._remove_xs_overlay()
            self._canvas.draw_idle()
            self._exit_xs_mode()

    def _compute_xs_flux(self):
        """Compute discharge Q(t) through the cross-section using exact mesh edge intersections."""
        import numpy as np
        sp = self._splotter
        P1, P2 = self._xs_pts[0], self._xs_pts[1]
        # xs_pts are in relative mesh coords; get_flow_through_cross_section
        # needs absolute coords (relative + geo_reference origin).
        ox, oy = sp.xllcorner, sp.yllcorner
        polyline = [[P1[0] + ox, P1[1] + oy], [P2[0] + ox, P2[1] + oy]]
        self._set_status('Computing cross-section discharge...')
        self.root.update_idletasks()
        try:
            _, Q = sp.get_flow_through_cross_section(polyline)
            self._xs_flux = np.asarray(Q)
        except Exception as e:
            self._set_status(f'Cross-section error: {e}')
            return
        dx = P2[0] - P1[0]
        dy = P2[1] - P1[1]
        L = float(np.sqrt(dx*dx + dy*dy))
        self._xs_info_lbl.config(
            text=f'length={L:.0f} m  |  '
                 f'P1=({P1[0]:.0f}, {P1[1]:.0f})  '
                 f'P2=({P2[0]:.0f}, {P2[1]:.0f})')
        self._update_xs_plot()
        if not self._xs_outer.winfo_ismapped():
            self._xs_outer.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=False,
                                before=self._canvas_frame)
        self._clear_xs_btn.config(state=tk.NORMAL)
        self._set_status(
            f'Cross-section discharge computed.  '
            f'Q range [{self._xs_flux.min():.3g}, {self._xs_flux.max():.3g}] m³/s')

    def _update_xs_plot(self):
        if self._xs_flux is None or self._splotter is None:
            return
        ax = self._xs_ax_plot
        ax.cla()
        t = self._splotter.time
        ax.plot(t, self._xs_flux, color='steelblue', linewidth=1.0)
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Q (m³/s)')
        ax.set_xlim(t[0], t[-1])
        ax.grid(True, linestyle=':', alpha=0.5)
        # Vertical cursor at current frame
        if len(self._frames) > 1:
            frac = self._current / (len(self._frames) - 1)
        else:
            frac = 0.0
        ts_idx = int(round(frac * (len(t) - 1)))
        self._xs_vline = ax.axvline(t[ts_idx], color='red',
                                    linestyle='--', linewidth=0.8)
        self._xs_fig.tight_layout(pad=1.5)
        self._xs_fig_canvas.draw()

    def _update_xs_overlay(self):
        """Draw cross-section endpoints and line on the animation canvas."""
        self._remove_xs_overlay()
        if not self._xs_pts or self._plot_transform is None:
            return
        sp = self._splotter
        use_basemap = self._gen_used_basemap and (sp is not None) and (sp.epsg is not None)
        for xr, yr in self._xs_pts:
            xa = xr + sp.xllcorner if use_basemap else xr
            ya = yr + sp.yllcorner if use_basemap else yr
            px, py = self._mesh_to_imshow(xa, ya)
            m, = self._ax.plot(px, py, 'o', color='cyan', markersize=8,
                               zorder=10, scalex=False, scaley=False)
            self._xs_artists.append(m)
        if len(self._xs_pts) == 2:
            (xr0, yr0), (xr1, yr1) = self._xs_pts
            xa0 = xr0 + sp.xllcorner if use_basemap else xr0
            ya0 = yr0 + sp.yllcorner if use_basemap else yr0
            xa1 = xr1 + sp.xllcorner if use_basemap else xr1
            ya1 = yr1 + sp.yllcorner if use_basemap else yr1
            px0, py0 = self._mesh_to_imshow(xa0, ya0)
            px1, py1 = self._mesh_to_imshow(xa1, ya1)
            ln, = self._ax.plot([px0, px1], [py0, py1], '-',
                                color='cyan', linewidth=1.5,
                                zorder=9, scalex=False, scaley=False)
            self._xs_artists.append(ln)
        self._canvas.draw_idle()

    def _remove_xs_overlay(self):
        for a in self._xs_artists:
            try:
                a.remove()
            except Exception:
                pass
        self._xs_artists = []

    def _close_xs_panel(self):
        if self._xs_outer.winfo_ismapped():
            self._xs_outer.pack_forget()

    def _clear_xs(self):
        self._exit_xs_mode()
        self._xs_pts = []
        self._xs_flux = None
        self._xs_vline = None
        self._remove_xs_overlay()
        self._xs_ax_plot.cla()
        self._xs_fig_canvas.draw()
        self._xs_info_lbl.config(text='')
        self._close_xs_panel()
        if hasattr(self, '_clear_xs_btn'):
            self._clear_xs_btn.config(state=tk.DISABLED)
        self._canvas.draw_idle()

    def _export_xs(self):
        if self._xs_flux is None or self._splotter is None:
            return
        path = filedialog.asksaveasfilename(
            title='Export cross-section discharge',
            defaultextension='.csv',
            filetypes=[('CSV', '*.csv'), ('All files', '*.*')],
            initialfile='cross_section_discharge.csv')
        if not path:
            return
        import csv
        t = self._splotter.time
        with open(path, 'w', newline='') as fh:
            writer = csv.writer(fh)
            writer.writerow(['time_s', 'discharge_m3_per_s'])
            for ti, qi in zip(t, self._xs_flux):
                writer.writerow([f'{float(ti):.4f}', f'{float(qi):.6f}'])
        self._set_status(f'Cross-section discharge exported to {os.path.basename(path)}')

    # -------------------------------------------------------------- #
    # Zoom                                                            #
    # -------------------------------------------------------------- #

    def _toggle_zoom_mode(self):
        if self._zoom_mode:
            self._exit_zoom_mode()
        else:
            self._enter_zoom_mode()

    def _enter_zoom_mode(self):
        """Activate rubber-band selection mode for choosing a zoom region."""
        if self._plot_transform is None:
            self._set_status('Generate frames first to enable zoom.')
            return
        from matplotlib.widgets import RectangleSelector
        self._zoom_selector = RectangleSelector(
            self._ax, self._on_zoom_select,
            useblit=False,
            button=[1],
            minspanx=5, minspany=5,
            spancoords='pixels',
            interactive=False)
        self._zoom_mode = True
        self._zoom_btn.config(text='Cancel Zoom')
        self._canvas.get_tk_widget().config(cursor='crosshair')
        self._set_status('Drag a rectangle on the image to set the zoom region.')

    def _exit_zoom_mode(self):
        if self._zoom_selector is not None:
            self._zoom_selector.set_active(False)
            self._zoom_selector = None
        self._zoom_mode = False
        self._zoom_btn.config(text='Set Zoom')
        self._canvas.get_tk_widget().config(cursor='')

    def _on_zoom_select(self, eclick, erelease):
        """Convert the rubber-band selection to mesh coordinates and store."""
        pt = self._plot_transform
        if pt is None:
            return
        x1, y1 = self._imshow_to_mesh(eclick.xdata,   eclick.ydata)
        x2, y2 = self._imshow_to_mesh(erelease.xdata, erelease.ydata)
        self._zoom_xlim = (min(x1, x2), max(x1, x2))
        self._zoom_ylim = (min(y1, y2), max(y1, y2))
        self._elev_contour_data = None   # recompute levels for the new zoom region
        self._exit_zoom_mode()
        self._draw_zoom_patch()
        self._reset_zoom_btn.config(state=tk.NORMAL)
        self._gen_btn.config(text='↻ Regenerate Frames')
        self._set_status(
            f'Zoom set - x: {self._zoom_xlim[0]:.1f} to {self._zoom_xlim[1]:.1f}  '
            f'y: {self._zoom_ylim[0]:.1f} to {self._zoom_ylim[1]:.1f}  '
            '- click ↻ Regenerate Frames to apply.')

    def _draw_zoom_patch(self):
        """Draw a yellow rectangle on the animation canvas showing the zoom region."""
        self._remove_zoom_patch()
        if self._zoom_xlim is None or self._plot_transform is None:
            return
        from matplotlib.patches import Rectangle
        px0, py0 = self._mesh_to_imshow(self._zoom_xlim[0], self._zoom_ylim[0])
        px1, py1 = self._mesh_to_imshow(self._zoom_xlim[1], self._zoom_ylim[1])
        x = min(px0, px1)
        y = min(py0, py1)
        w = abs(px1 - px0)
        h = abs(py1 - py0)
        self._zoom_rect_patch = Rectangle(
            (x, y), w, h,
            linewidth=2, edgecolor='yellow', facecolor='yellow',
            alpha=0.15, zorder=12)
        self._ax.add_patch(self._zoom_rect_patch)
        self._canvas.draw_idle()

    def _remove_zoom_patch(self):
        if self._zoom_rect_patch is not None:
            try:
                self._zoom_rect_patch.remove()
            except Exception:
                pass
            self._zoom_rect_patch = None
        self._canvas.draw_idle()

    def _reset_zoom(self):
        """Clear the zoom region and remove the overlay patch."""
        self._exit_zoom_mode()
        self._zoom_xlim = None
        self._zoom_ylim = None
        self._elev_contour_data = None   # recompute levels for full extent
        self._remove_zoom_patch()
        self._reset_zoom_btn.config(state=tk.DISABLED)
        if self._frames:
            self._gen_btn.config(text='↻ Regenerate Frames')
        self._set_status('Zoom reset - full extent will be used for generation.')

    # -------------------------------------------------------------- #
    # Mesh overlay                                                   #
    # -------------------------------------------------------------- #

    def _on_show_mesh_toggle(self):
        self._update_mesh_overlay()
        self._canvas.draw_idle()

    def _remove_mesh_overlay(self):
        for artist in self._mesh_overlay_lines:
            try:
                artist.remove()
            except Exception:
                pass
        self._mesh_overlay_lines = []

    def _update_mesh_overlay(self):
        """Build (or remove) a triplot overlay on the animation canvas.

        The overlay is drawn in image-pixel coordinates so it stays aligned
        with the imshow across all frames without per-frame rebuilding.
        """
        self._remove_mesh_overlay()
        if not self._show_mesh_var.get():
            return
        if self._plot_transform is None or self._splotter is None:
            return
        # Skip canvas overlay if mesh is already baked into frames
        if self._last_gen_show_mesh:
            return

        import numpy as np
        from matplotlib.tri import Triangulation

        sp  = self._splotter
        pt  = self._plot_transform
        use_basemap = self._gen_used_basemap and sp.epsg is not None
        src = sp.triang_abs if use_basemap else sp.triang

        xlim, ylim = pt['xlim'], pt['ylim']
        pos = pt['pos']
        W, H = pt['W'], pt['H']

        # Vectorised transform: mesh coords → image pixel coords
        xv = src.x
        yv = src.y
        px = (pos.x0 + pos.width  * (xv - xlim[0]) / (xlim[1] - xlim[0])) * W
        py = (1.0 - (pos.y0 + pos.height * (yv - ylim[0]) / (ylim[1] - ylim[0]))) * H

        # When zoomed, mask any triangle that has a vertex outside the visible
        # region so no edges extend beyond the plot boundary.
        mask = None
        if self._zoom_xlim is not None:
            tris = src.triangles
            mask = (
                (xv[tris[:, 0]] < xlim[0]) | (xv[tris[:, 0]] > xlim[1]) |
                (yv[tris[:, 0]] < ylim[0]) | (yv[tris[:, 0]] > ylim[1]) |
                (xv[tris[:, 1]] < xlim[0]) | (xv[tris[:, 1]] > xlim[1]) |
                (yv[tris[:, 1]] < ylim[0]) | (yv[tris[:, 1]] > ylim[1]) |
                (xv[tris[:, 2]] < xlim[0]) | (xv[tris[:, 2]] > xlim[1]) |
                (yv[tris[:, 2]] < ylim[0]) | (yv[tris[:, 2]] > ylim[1])
            )

        triang_px = Triangulation(px, py, src.triangles, mask=mask)
        self._mesh_overlay_lines = self._ax.triplot(
            triang_px, color='black', linewidth=0.25, alpha=0.45,
            scalex=False, scaley=False)

    # -------------------------------------------------------------- #
    # Elevation contour overlay                                       #
    # -------------------------------------------------------------- #

    def _on_show_elev_toggle(self):
        self._update_elev_overlay()
        self._canvas.draw_idle()

    def _on_elev_levels_changed(self):
        self._elev_contour_data = None    # force recompute at new level count
        if self._show_elev_var.get():
            self._update_elev_overlay()
            self._canvas.draw_idle()

    def _remove_elev_overlay(self):
        for obj in self._elev_overlay_artists:
            try:
                obj.remove()
            except Exception:
                # TriContourSet in older matplotlib needs per-collection removal
                try:
                    for coll in obj.collections:
                        coll.remove()
                except Exception:
                    pass
        self._elev_overlay_artists = []

    @staticmethod
    def _nice_contour_levels(vmin, vmax, n):
        """Return contour levels as multiples of a round step spanning vmin..vmax."""
        import numpy as _np
        span = vmax - vmin
        if span <= 0:
            return n
        raw_step = span / max(n, 1)
        magnitude = 10.0 ** _np.floor(_np.log10(raw_step))
        for factor in (1, 2, 5, 10, 20, 50, 100):
            step = magnitude * factor
            if span / step <= n * 1.5:
                break
        first = _np.ceil(vmin / step) * step
        levels = _np.arange(first, vmax + step * 0.01, step)
        return levels if len(levels) >= 2 else n

    def _build_elev_contour_data(self):
        """Build pixel-space triangulation + per-vertex elevation + nice levels.

        Cached until the plot transform or level count changes.
        Uses t=0 elevation for time-varying SWW files.
        """
        import numpy as _np
        from matplotlib.tri import Triangulation as _Tri
        from anuga.utilities.animate import _face_to_vertex

        sp  = self._splotter
        pt  = self._plot_transform
        use_basemap = self._gen_used_basemap and sp.epsg is not None
        src = sp.triang_abs if use_basemap else sp.triang

        elev_raw = sp.elev
        if elev_raw.ndim == 2:
            elev_face = elev_raw[0]
            self._set_status('Elev contours from t=0 (time-varying elevation)')
        else:
            elev_face = elev_raw

        elev_nodes = _face_to_vertex(src, elev_face)

        try:
            n_levels = max(3, int(self._elev_levels_var.get()))
        except ValueError:
            n_levels = 10

        # build pixel-space triangulation (same transform as mesh overlay)
        xlim, ylim = pt['xlim'], pt['ylim']
        pos = pt['pos']
        W, H = pt['W'], pt['H']
        xv = src.x
        yv = src.y
        px = (pos.x0 + pos.width  * (xv - xlim[0]) / (xlim[1] - xlim[0])) * W
        py = (1.0 - (pos.y0 + pos.height * (yv - ylim[0]) / (ylim[1] - ylim[0]))) * H

        # When zoomed, mask any triangle with a vertex outside the visible
        # region so contour lines don't extend beyond the plot boundary.
        mask = None
        if self._zoom_xlim is not None:
            tris = src.triangles
            mask = (
                (xv[tris[:, 0]] < xlim[0]) | (xv[tris[:, 0]] > xlim[1]) |
                (yv[tris[:, 0]] < ylim[0]) | (yv[tris[:, 0]] > ylim[1]) |
                (xv[tris[:, 1]] < xlim[0]) | (xv[tris[:, 1]] > xlim[1]) |
                (yv[tris[:, 1]] < ylim[0]) | (yv[tris[:, 1]] > ylim[1]) |
                (xv[tris[:, 2]] < xlim[0]) | (xv[tris[:, 2]] > xlim[1]) |
                (yv[tris[:, 2]] < ylim[0]) | (yv[tris[:, 2]] > ylim[1])
            )
        triang_px = _Tri(px, py, src.triangles, mask=mask)

        # Compute contour levels from the elevation within the visible region
        # only, so the requested number of levels actually appears in the view.
        if self._zoom_xlim is not None:
            in_zoom = ((xv >= xlim[0]) & (xv <= xlim[1]) &
                       (yv >= ylim[0]) & (yv <= ylim[1]))
            elev_vis = elev_nodes[in_zoom]
            if elev_vis.size >= 2:
                e_vmin = float(elev_vis.min())
                e_vmax = float(elev_vis.max())
            else:
                e_vmin = float(elev_nodes.min())
                e_vmax = float(elev_nodes.max())
        else:
            e_vmin = float(elev_nodes.min())
            e_vmax = float(elev_nodes.max())

        levels = self._nice_contour_levels(e_vmin, e_vmax, n_levels)

        self._elev_contour_data = (triang_px, elev_nodes, levels)

    def _update_elev_overlay(self):
        self._remove_elev_overlay()
        if not self._show_elev_var.get():
            return
        if self._plot_transform is None or self._splotter is None:
            return
        # Skip canvas overlay if elevation contours are already baked into frames
        if self._last_gen_show_elev:
            return

        if self._elev_contour_data is None:
            self._build_elev_contour_data()

        triang_px, elev_nodes, levels = self._elev_contour_data
        xlim, ylim = self._ax.get_xlim(), self._ax.get_ylim()
        cs = self._ax.tricontour(triang_px, elev_nodes, levels=levels,
                                 colors='dimgray', linewidths=0.6, alpha=0.7)
        labels = self._ax.clabel(cs, fmt='%g m', fontsize=6,
                                 inline=True, inline_spacing=2)
        self._ax.set_xlim(xlim)
        self._ax.set_ylim(ylim)
        self._elev_overlay_artists = [cs] + labels

    # -------------------------------------------------------------- #
    # Save mesh                                                      #
    # -------------------------------------------------------------- #

    def _save_mesh(self):
        """Open a modal dialog to configure and save the mesh render."""
        if self._splotter is None:
            return

        win = tk.Toplevel(self.root)
        win.title('Save mesh')
        win.resizable(False, False)
        win.grab_set()

        pad = dict(padx=8, pady=4)

        # --- format + DPI row ---
        fmt_frame = ttk.Frame(win, padding=6)
        fmt_frame.pack(fill=tk.X)
        ttk.Label(fmt_frame, text='Format:').pack(side=tk.LEFT, **pad)
        fmt_var = tk.StringVar(value='PDF')
        fmt_combo = ttk.Combobox(fmt_frame, textvariable=fmt_var, width=6,
                                 state='readonly',
                                 values=['PNG', 'PDF', 'SVG', 'PGF'])
        fmt_combo.pack(side=tk.LEFT, padx=(0, 16))

        ttk.Label(fmt_frame, text='DPI (PNG only):').pack(side=tk.LEFT, **pad)
        dpi_var = tk.IntVar(value=150)
        dpi_spin = ttk.Spinbox(fmt_frame, from_=72, to=600, increment=25,
                               textvariable=dpi_var, width=6)
        dpi_spin.pack(side=tk.LEFT)

        # --- PGF TeX engine row (shown only when PGF is selected) ---
        import shutil as _shutil
        _tex_engines = [e for e in ('pdflatex', 'xelatex', 'lualatex')
                        if _shutil.which(e)]
        pgf_row = ttk.Frame(win, padding=(6, 0, 6, 2))
        tex_var = tk.StringVar(value=_tex_engines[0] if _tex_engines else 'pdflatex')
        ttk.Label(pgf_row, text='TeX engine:').pack(side=tk.LEFT, **pad)
        tex_combo = ttk.Combobox(pgf_row, textvariable=tex_var, width=10,
                                 state='readonly',
                                 values=['pdflatex', 'xelatex', 'lualatex'])
        tex_combo.pack(side=tk.LEFT)
        if not _tex_engines:
            ttk.Label(pgf_row, text='[!] none found on PATH',
                      foreground='red').pack(side=tk.LEFT, padx=8)

        def _toggle_pgf_row(*_):
            dpi_spin.config(state='normal' if fmt_var.get() == 'PNG' else 'disabled')
            if fmt_var.get() == 'PGF':
                pgf_row.pack(fill=tk.X, after=fmt_frame)
            else:
                pgf_row.pack_forget()

        fmt_combo.bind('<<ComboboxSelected>>', _toggle_pgf_row)
        _toggle_pgf_row()

        ttk.Separator(win, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=6)

        # --- decoration checkboxes (no colorbar — mesh has none) ---
        chk_frame = ttk.Frame(win, padding=6)
        chk_frame.pack(fill=tk.X)
        labels_var  = tk.BooleanVar(value=True)
        title_var   = tk.BooleanVar(value=True)
        sp = self._splotter
        _bm_available = sp is not None and sp.epsg is not None
        basemap_var = tk.BooleanVar(
            value=self._gen_used_basemap and _bm_available)
        ttk.Checkbutton(chk_frame, text='Axis labels', variable=labels_var
                        ).pack(side=tk.LEFT, **pad)
        ttk.Checkbutton(chk_frame, text='Title',       variable=title_var
                        ).pack(side=tk.LEFT, **pad)
        bm_chk = ttk.Checkbutton(chk_frame, text='Basemap', variable=basemap_var)
        bm_chk.pack(side=tk.LEFT, **pad)
        if not _bm_available:
            bm_chk.config(state='disabled')

        ttk.Separator(win, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=6)

        # --- buttons ---
        btn_frame = ttk.Frame(win, padding=6)
        btn_frame.pack(fill=tk.X)
        ttk.Button(btn_frame, text='Cancel',
                   command=win.destroy).pack(side=tk.RIGHT, padx=4)

        def _do_save():
            fmt = fmt_var.get().lower()
            ext = f'.{fmt}'
            default = f'{self._sww_prefix}_mesh'
            filetypes = {
                'png': [('PNG image',    '*.png')],
                'pdf': [('PDF document', '*.pdf')],
                'svg': [('SVG image',    '*.svg')],
                'pgf': [('PGF/LaTeX',   '*.pgf')],
            }
            path = filedialog.asksaveasfilename(
                parent=win,
                title='Save mesh',
                initialfile=default,
                defaultextension=ext,
                filetypes=filetypes.get(fmt, [('All files', '*.*')]))
            if not path:
                return
            win.destroy()
            try:
                self._render_and_save_mesh(
                    path,
                    show_labels=labels_var.get(),
                    show_title=title_var.get(),
                    dpi=dpi_var.get(),
                    tex_engine=tex_var.get(),
                    use_basemap=basemap_var.get(),
                )
                self._set_status(f'Mesh saved to {path}')
            except Exception as exc:
                from tkinter import messagebox
                messagebox.showerror('Save mesh error', str(exc))

        ttk.Button(btn_frame, text='Save...',
                   command=_do_save).pack(side=tk.RIGHT, padx=4)

        win.wait_window()

    def _render_and_save_mesh(self, path, show_labels=True, show_title=True,
                               dpi=150, tex_engine='pdflatex', use_basemap=None):
        """Render the triangulation with triplot and save to *path*."""
        import os
        import matplotlib

        sp = self._splotter
        if use_basemap is None:
            use_basemap = self._gen_used_basemap and sp.epsg is not None
        triang = sp.triang_abs if use_basemap else sp.triang

        ext = os.path.splitext(path)[1].lower()
        if ext == '.pgf':
            old_backend = matplotlib.get_backend()
            old_texsystem = matplotlib.rcParams.get('pgf.texsystem', 'xelatex')
            matplotlib.rcParams['pgf.texsystem'] = tex_engine
            try:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(8, 6))
                ax.triplot(triang, color='steelblue', linewidth=0.4, alpha=0.7)
                ax.set_aspect('equal')
                if show_labels:
                    if use_basemap:
                        ax.set_xlabel('Easting (m)')
                        ax.set_ylabel('Northing (m)')
                    else:
                        ax.set_xlabel('x (m)')
                        ax.set_ylabel('y (m)')
                else:
                    ax.set_xlabel('')
                    ax.set_ylabel('')
                    ax.tick_params(labelbottom=False, labelleft=False)
                if show_title:
                    ax.set_title(f'{self._sww_prefix}  -  mesh  '
                                 f'({len(sp.triangles)} triangles)')
                if use_basemap:
                    from anuga.utilities.animate import _add_basemap, BASEMAP_PROVIDERS
                    provider_label = self._basemap_provider_var.get()
                    provider_str = BASEMAP_PROVIDERS.get(
                        provider_label, 'OpenStreetMap.Mapnik')
                    try:
                        _add_basemap(ax, sp.epsg, provider_str, cache=sp._basemap_cache)
                    except Exception:
                        pass
                fig.tight_layout()
                fig.savefig(path, bbox_inches='tight', backend='pgf')
                plt.close(fig)
            except Exception:
                matplotlib.rcParams['pgf.texsystem'] = old_texsystem
                raise
            matplotlib.rcParams['pgf.texsystem'] = old_texsystem
            return

        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg

        fig = Figure(figsize=(8, 6))
        FigureCanvasAgg(fig)
        ax = fig.add_subplot(111)
        ax.triplot(triang, color='steelblue', linewidth=0.4, alpha=0.7)
        ax.set_aspect('equal')

        if show_labels:
            if use_basemap:
                ax.set_xlabel('Easting (m)')
                ax.set_ylabel('Northing (m)')
            else:
                ax.set_xlabel('x (m)')
                ax.set_ylabel('y (m)')
        else:
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.tick_params(labelbottom=False, labelleft=False)

        if show_title:
            ax.set_title(f'{self._sww_prefix}  -  mesh  '
                         f'({len(sp.triangles)} triangles)')

        if use_basemap:
            from anuga.utilities.animate import _add_basemap, BASEMAP_PROVIDERS
            provider_label = self._basemap_provider_var.get()
            provider_str = BASEMAP_PROVIDERS.get(
                provider_label, 'OpenStreetMap.Mapnik')
            try:
                _add_basemap(ax, sp.epsg, provider_str, cache=sp._basemap_cache)
            except Exception:
                pass

        fig.tight_layout()
        fig.savefig(path, dpi=dpi, bbox_inches='tight')

    def _on_hover(self, event):
        """Update the coordinate readout as the cursor moves over the animation."""
        if (self._plot_transform is None
                or self._splotter is None
                or event.inaxes != self._ax
                or event.xdata is None):
            self._coord_var.set('')
            return

        xmesh, ymesh = self._imshow_to_mesh(event.xdata, event.ydata)

        if self._trifinder is None:
            import matplotlib.tri as _mtri
            sp = self._splotter
            use_basemap = self._gen_used_basemap and (sp.epsg is not None)
            src = sp.triang_abs if use_basemap else sp.triang
            # Build trifinder on a mask-free copy so every point in the domain
            # returns a valid triangle index (wet or dry).
            tf_triang = _mtri.Triangulation(src.x, src.y, src.triangles)
            self._trifinder = tf_triang.get_trifinder()

        tri = int(self._trifinder(xmesh, ymesh))
        if tri < 0:
            self._coord_var.set(f'x={xmesh:.1f}  y={ymesh:.1f}')
        else:
            self._coord_var.set(f'x={xmesh:.1f}  y={ymesh:.1f}  tri={tri}')

    def _imshow_to_mesh(self, px, py):
        """Convert imshow pixel coordinates to mesh/absolute coordinates."""
        pt = self._plot_transform
        pos = pt['pos']
        xlim, ylim = pt['xlim'], pt['ylim']
        W, H = pt['W'], pt['H']
        xfrac = (px / W - pos.x0) / pos.width
        yfrac = (1.0 - py / H - pos.y0) / pos.height
        mx = xlim[0] + xfrac * (xlim[1] - xlim[0])
        my = ylim[0] + yfrac * (ylim[1] - ylim[0])
        return mx, my

    def _mesh_to_imshow(self, mx, my):
        """Convert mesh/absolute coordinates to imshow pixel coordinates."""
        pt = self._plot_transform
        pos = pt['pos']
        xlim, ylim = pt['xlim'], pt['ylim']
        W, H = pt['W'], pt['H']
        xfrac = (mx - xlim[0]) / (xlim[1] - xlim[0])
        yfrac = (my - ylim[0]) / (ylim[1] - ylim[0])
        px = (pos.x0 + pos.width  * xfrac) * W
        py = (1.0 - (pos.y0 + pos.height * yfrac)) * H
        return px, py

    # -------------------------------------------------------------- #
    # Save frame / animation                                         #
    # -------------------------------------------------------------- #

    @staticmethod
    def _parse_times(text):
        """Parse a comma-separated string of floats; silently skip bad tokens."""
        times = []
        for part in text.split(','):
            part = part.strip()
            if part:
                try:
                    times.append(float(part))
                except ValueError:
                    pass
        return times

    def _time_to_sww_idx(self, t):
        """Return the SWW timestep index nearest to time *t*."""
        import numpy as _np
        return int(_np.argmin(_np.abs(self._splotter.time - t)))

    def _time_to_frame_idx(self, t):
        """Return the loaded-frame index nearest to SWW time *t*."""
        sww_idx = self._time_to_sww_idx(t)
        n_sww   = len(self._splotter.time)
        n_frames = len(self._frames)
        if n_frames <= 1 or n_sww <= 1:
            return 0
        frac = sww_idx / (n_sww - 1)
        return int(round(frac * (n_frames - 1)))

    def _save_frame(self):
        """Open a dialog to save the current frame, or a set of times, to PNG."""
        if not self._frames:
            return

        win = tk.Toplevel(self.root)
        win.title('Save Frame')
        win.resizable(False, False)
        win.grab_set()
        pad = dict(padx=8, pady=4)

        # --- frame selection ---
        sel_lf = ttk.LabelFrame(win, text='Frames', padding=6)
        sel_lf.pack(fill=tk.X, padx=8, pady=(8, 4))

        sel_var = tk.StringVar(value='current')
        ttk.Radiobutton(sel_lf, text='Current frame',
                        variable=sel_var, value='current').grid(
            row=0, column=0, sticky=tk.W)
        ttk.Radiobutton(sel_lf, text='Times (s):',
                        variable=sel_var, value='times').grid(
            row=1, column=0, sticky=tk.W, pady=(4, 0))

        sp = self._splotter
        t0, t1 = float(sp.time[0]), float(sp.time[-1])
        times_var = tk.StringVar(value=f'{t0:.1f}')
        times_entry = ttk.Entry(sel_lf, textvariable=times_var, width=40)
        times_entry.grid(row=1, column=1, sticky=tk.EW, padx=(4, 0), pady=(4, 0))
        ttk.Label(sel_lf, text=f'(comma-separated; range {t0:.1f} – {t1:.1f} s)',
                  foreground='grey').grid(row=2, column=1, sticky=tk.W, padx=(4, 0))
        sel_lf.columnconfigure(1, weight=1)

        ttk.Separator(win, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=6, pady=4)

        # --- buttons ---
        btn_frame = ttk.Frame(win, padding=6)
        btn_frame.pack(fill=tk.X)
        ttk.Button(btn_frame, text='Cancel',
                   command=win.destroy).pack(side=tk.RIGHT, padx=4)

        def _do_save():
            from tkinter import filedialog
            if sel_var.get() == 'current':
                default = (f'{self._sww_prefix}_{self._last_gen_qty}'
                           f'_frame{self._current + 1:04d}')
                path = filedialog.asksaveasfilename(
                    parent=win,
                    title='Save frame',
                    initialfile=default,
                    defaultextension='.png',
                    filetypes=[('PNG image', '*.png'),
                               ('PDF document', '*.pdf'),
                               ('SVG image', '*.svg'),
                               ('All files', '*.*')])
                if not path:
                    return
                win.destroy()
                self._fig.savefig(path, dpi=self._last_gen_dpi, bbox_inches='tight')
                self._set_status(f'Frame saved to {path}')
            else:
                times = self._parse_times(times_var.get())
                if not times:
                    from tkinter import messagebox
                    messagebox.showerror('No times', 'Enter at least one time value.',
                                         parent=win)
                    return
                if len(times) == 1:
                    # Single time → single file save dialog
                    fidx = self._time_to_frame_idx(times[0])
                    t_actual = sp.time[self._time_to_sww_idx(times[0])]
                    default = (f'{self._sww_prefix}_{self._last_gen_qty}'
                               f'_t{t_actual:.1f}')
                    path = filedialog.asksaveasfilename(
                        parent=win,
                        title='Save frame',
                        initialfile=default,
                        defaultextension='.png',
                        filetypes=[('PNG image', '*.png'), ('All files', '*.*')])
                    if not path:
                        return
                    win.destroy()
                    import shutil
                    shutil.copy2(self._frames[fidx], path)
                    self._set_status(f'Frame saved to {path}')
                else:
                    # Multiple times → choose output directory
                    out_dir = filedialog.askdirectory(
                        parent=win,
                        title='Choose output directory for frames')
                    if not out_dir:
                        return
                    win.destroy()
                    import shutil
                    saved = []
                    for t_req in times:
                        fidx = self._time_to_frame_idx(t_req)
                        t_actual = sp.time[self._time_to_sww_idx(t_req)]
                        fname = (f'{self._sww_prefix}_{self._last_gen_qty}'
                                 f'_t{t_actual:.1f}.png')
                        dest = os.path.join(out_dir, fname)
                        shutil.copy2(self._frames[fidx], dest)
                        saved.append(dest)
                    self._set_status(
                        f'Saved {len(saved)} frames to {out_dir}')

        ttk.Button(btn_frame, text='Save...',
                   command=_do_save).pack(side=tk.RIGHT, padx=4)

        win.wait_window()

    def _save_animation(self):
        """Save all loaded frames as a GIF or MP4 animation."""
        if not self._frames:
            return
        import shutil
        from tkinter import filedialog

        has_ffmpeg = shutil.which('ffmpeg') is not None
        filetypes = []
        if has_ffmpeg:
            filetypes.append(('MP4 video', '*.mp4'))
        filetypes.append(('GIF animation', '*.gif'))
        filetypes.append(('All files', '*.*'))

        default_ext = '.mp4' if has_ffmpeg else '.gif'
        default = f'{self._sww_prefix}_{self._last_gen_qty}'
        path = filedialog.asksaveasfilename(
            title='Save animation',
            initialfile=default,
            defaultextension=default_ext,
            filetypes=filetypes)
        if not path:
            return

        fps = max(0.5, self._fps_var.get())
        if path.lower().endswith('.mp4'):
            self._save_animation_mp4(path, fps)
        else:
            self._save_animation_gif(path, fps)

    def _save_animation_gif(self, path, fps):
        try:
            from PIL import Image
        except ImportError:
            from tkinter import messagebox
            messagebox.showerror(
                'Missing dependency',
                'Pillow is required to save GIF animations.\n'
                'Install with:  pip install Pillow')
            return
        self._set_status(f'Saving GIF ({len(self._frames)} frames)...')
        self.root.update_idletasks()
        duration_ms = max(20, int(1000 / fps))
        imgs = [Image.open(f).convert('RGBA') for f in self._frames]
        imgs[0].save(path, save_all=True, append_images=imgs[1:],
                     loop=0, duration=duration_ms, optimize=False)
        self._set_status(f'GIF saved ({len(imgs)} frames) to {path}')

    def _save_animation_mp4(self, path, fps):
        import shutil
        if not shutil.which('ffmpeg'):
            from tkinter import messagebox
            messagebox.showerror(
                'ffmpeg not found',
                'ffmpeg must be installed and on PATH to save MP4 videos.\n'
                'Download from https://ffmpeg.org or install via conda/apt/brew.')
            return
        import subprocess
        import tempfile
        import os
        self._set_status(f'Saving MP4 ({len(self._frames)} frames)...')
        self.root.update_idletasks()
        # Write a concat list so arbitrary frame sets work (stride, etc.)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt',
                                         delete=False, encoding='utf-8') as f:
            dur = 1.0 / fps
            for frame_path in self._frames:
                f.write(f"file '{frame_path}'\n")
                f.write(f"duration {dur:.6f}\n")
            # ffmpeg concat requires the last entry twice (no duration on final)
            f.write(f"file '{self._frames[-1]}'\n")
            concat_path = f.name
        try:
            result = subprocess.run(
                ['ffmpeg', '-y', '-f', 'concat', '-safe', '0',
                 '-i', concat_path,
                 '-c:v', 'libx264', '-pix_fmt', 'yuv420p',
                 '-movflags', '+faststart', path],
                capture_output=True, text=True)
            if result.returncode != 0:
                from tkinter import messagebox
                messagebox.showerror('ffmpeg error',
                                     result.stderr[-800:] or result.stdout[-800:])
                return
        finally:
            os.unlink(concat_path)
        self._set_status(f'MP4 saved ({len(self._frames)} frames) to {path}')

    # -------------------------------------------------------------- #
    # Export frame (re-rendered, with configurable decorations)       #
    # -------------------------------------------------------------- #

    def _export_frame(self):
        """Open a modal dialog to configure and export a re-rendered frame."""
        if not self._frames or self._splotter is None:
            return

        win = tk.Toplevel(self.root)
        win.title('Export frame')
        win.resizable(False, False)
        win.grab_set()

        pad = dict(padx=8, pady=4)

        # --- format + DPI row ---
        fmt_frame = ttk.Frame(win, padding=6)
        fmt_frame.pack(fill=tk.X)
        ttk.Label(fmt_frame, text='Format:').pack(side=tk.LEFT, **pad)
        fmt_var = tk.StringVar(value='PDF')
        fmt_combo = ttk.Combobox(fmt_frame, textvariable=fmt_var, width=6,
                                 state='readonly',
                                 values=['PNG', 'PDF', 'SVG', 'PGF'])
        fmt_combo.pack(side=tk.LEFT, padx=(0, 16))

        ttk.Label(fmt_frame, text='DPI (PNG only):').pack(side=tk.LEFT, **pad)
        dpi_var = tk.IntVar(value=self._last_gen_dpi)
        dpi_spin = ttk.Spinbox(fmt_frame, from_=72, to=600, increment=25,
                               textvariable=dpi_var, width=6)
        dpi_spin.pack(side=tk.LEFT)

        # --- PGF TeX engine row (shown only when PGF is selected) ---
        import shutil as _shutil
        _tex_engines = [e for e in ('pdflatex', 'xelatex', 'lualatex')
                        if _shutil.which(e)]
        pgf_row = ttk.Frame(win, padding=(6, 0, 6, 2))
        tex_var = tk.StringVar(value=_tex_engines[0] if _tex_engines else 'pdflatex')
        ttk.Label(pgf_row, text='TeX engine:').pack(side=tk.LEFT, **pad)
        tex_combo = ttk.Combobox(pgf_row, textvariable=tex_var, width=10,
                                 state='readonly',
                                 values=['pdflatex', 'xelatex', 'lualatex'])
        tex_combo.pack(side=tk.LEFT)
        if not _tex_engines:
            ttk.Label(pgf_row, text='[!] none found on PATH',
                      foreground='red').pack(side=tk.LEFT, padx=8)

        def _toggle_pgf_row(*_):
            dpi_spin.config(state='normal' if fmt_var.get() == 'PNG' else 'disabled')
            if fmt_var.get() == 'PGF':
                pgf_row.pack(fill=tk.X, after=fmt_frame)
            else:
                pgf_row.pack_forget()

        fmt_combo.bind('<<ComboboxSelected>>', _toggle_pgf_row)
        _toggle_pgf_row()

        ttk.Separator(win, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=6)

        # --- frame selection ---
        sel_lf = ttk.LabelFrame(win, text='Frames', padding=6)
        sel_lf.pack(fill=tk.X, padx=6, pady=4)

        exp_sel_var = tk.StringVar(value='current')
        ttk.Radiobutton(sel_lf, text='Current frame',
                        variable=exp_sel_var, value='current').grid(
            row=0, column=0, sticky=tk.W)
        ttk.Radiobutton(sel_lf, text='Times (s):',
                        variable=exp_sel_var, value='times').grid(
            row=1, column=0, sticky=tk.W, pady=(4, 0))

        sp = self._splotter
        t0e, t1e = float(sp.time[0]), float(sp.time[-1])
        exp_times_var = tk.StringVar(value=f'{t0e:.1f}')
        exp_times_entry = ttk.Entry(sel_lf, textvariable=exp_times_var, width=40)
        exp_times_entry.grid(row=1, column=1, sticky=tk.EW,
                             padx=(4, 0), pady=(4, 0))
        ttk.Label(sel_lf,
                  text=f'(comma-separated; range {t0e:.1f} – {t1e:.1f} s)',
                  foreground='grey').grid(row=2, column=1, sticky=tk.W, padx=(4, 0))
        sel_lf.columnconfigure(1, weight=1)

        ttk.Separator(win, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=6)

        # --- decoration checkboxes ---
        chk_frame = ttk.Frame(win, padding=6)
        chk_frame.pack(fill=tk.X)
        cbar_var   = tk.BooleanVar(value=True)
        labels_var = tk.BooleanVar(value=True)
        title_var  = tk.BooleanVar(value=True)
        edges_var  = tk.BooleanVar(value=self._last_gen_show_edges)
        mesh_var   = tk.BooleanVar(value=self._show_mesh_var.get())
        elev_var   = tk.BooleanVar(value=self._show_elev_var.get())
        ttk.Checkbutton(chk_frame, text='Colorbar',    variable=cbar_var
                        ).pack(side=tk.LEFT, **pad)
        ttk.Checkbutton(chk_frame, text='Axis labels', variable=labels_var
                        ).pack(side=tk.LEFT, **pad)
        ttk.Checkbutton(chk_frame, text='Title',       variable=title_var
                        ).pack(side=tk.LEFT, **pad)
        ttk.Checkbutton(chk_frame, text='Flat View',   variable=edges_var
                        ).pack(side=tk.LEFT, **pad)
        ttk.Checkbutton(chk_frame, text='Mesh',        variable=mesh_var
                        ).pack(side=tk.LEFT, **pad)
        ttk.Checkbutton(chk_frame, text='Elev contours', variable=elev_var
                        ).pack(side=tk.LEFT, **pad)

        ttk.Separator(win, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=6)

        # --- buttons ---
        btn_frame = ttk.Frame(win, padding=6)
        btn_frame.pack(fill=tk.X)
        ttk.Button(btn_frame, text='Cancel',
                   command=win.destroy).pack(side=tk.RIGHT, padx=4)

        def _do_export():
            from tkinter import messagebox
            fmt = fmt_var.get().lower()
            ext = f'.{fmt}'
            filetypes = {
                'png': [('PNG image', '*.png')],
                'pdf': [('PDF document', '*.pdf')],
                'svg': [('SVG image', '*.svg')],
                'pgf': [('PGF/LaTeX', '*.pgf')],
            }
            kw = dict(
                show_colorbar=cbar_var.get(),
                show_labels=labels_var.get(),
                show_title=title_var.get(),
                show_edges=edges_var.get(),
                show_mesh=mesh_var.get(),
                show_elev=elev_var.get(),
                elev_levels=max(3, int(self._elev_levels_var.get())),
                dpi=dpi_var.get(),
                tex_engine=tex_var.get(),
            )

            if exp_sel_var.get() == 'current':
                default = (f'{self._sww_prefix}_{self._last_gen_qty}'
                           f'_frame{self._current + 1:04d}')
                path = filedialog.asksaveasfilename(
                    parent=win,
                    title='Export frame',
                    initialfile=default,
                    defaultextension=ext,
                    filetypes=filetypes.get(fmt, [('All files', '*.*')]))
                if not path:
                    return
                win.destroy()
                try:
                    self._render_and_export_frame(path, **kw)
                    self._set_status(f'Frame exported to {path}')
                except Exception as exc:
                    messagebox.showerror('Export error', str(exc))
            else:
                times = self._parse_times(exp_times_var.get())
                if not times:
                    messagebox.showerror('No times',
                                         'Enter at least one time value.',
                                         parent=win)
                    return
                if len(times) == 1:
                    t_actual = sp.time[self._time_to_sww_idx(times[0])]
                    default = (f'{self._sww_prefix}_{self._last_gen_qty}'
                               f'_t{t_actual:.1f}')
                    path = filedialog.asksaveasfilename(
                        parent=win,
                        title='Export frame',
                        initialfile=default,
                        defaultextension=ext,
                        filetypes=filetypes.get(fmt, [('All files', '*.*')]))
                    if not path:
                        return
                    win.destroy()
                    try:
                        saved_current = self._current
                        self._current = self._time_to_frame_idx(times[0])
                        self._render_and_export_frame(path, **kw)
                        self._current = saved_current
                        self._set_status(f'Frame exported to {path}')
                    except Exception as exc:
                        messagebox.showerror('Export error', str(exc))
                else:
                    out_dir = filedialog.askdirectory(
                        parent=win,
                        title='Choose output directory for exported frames')
                    if not out_dir:
                        return
                    win.destroy()
                    saved_current = self._current
                    exported = []
                    errors = []
                    for t_req in times:
                        sww_idx = self._time_to_sww_idx(t_req)
                        t_actual = float(sp.time[sww_idx])
                        fname = (f'{self._sww_prefix}_{self._last_gen_qty}'
                                 f'_t{t_actual:.1f}{ext}')
                        path = os.path.join(out_dir, fname)
                        try:
                            self._current = self._time_to_frame_idx(t_req)
                            self._render_and_export_frame(path, **kw)
                            exported.append(path)
                        except Exception as exc:
                            errors.append(f't={t_actual:.1f}: {exc}')
                    self._current = saved_current
                    msg = f'Exported {len(exported)} frames to {out_dir}'
                    if errors:
                        msg += f' ({len(errors)} errors)'
                    self._set_status(msg)
                    if errors:
                        from tkinter import messagebox as _mb
                        _mb.showwarning('Export warnings',
                                        '\n'.join(errors[:10]))

        ttk.Button(btn_frame, text='Export...',
                   command=_do_export).pack(side=tk.RIGHT, padx=4)

        win.wait_window()

    def _render_and_export_frame(self, path, show_colorbar=True,
                                  show_labels=True, show_title=True,
                                  show_edges=False, show_mesh=False,
                                  show_elev=False, elev_levels=10,
                                  dpi=150, tex_engine='pdflatex'):
        """Re-render the current frame from raw SWW data and save to *path*.

        Unlike "Save Frame" (which screenshots the imshow canvas), this
        re-renders via tripcolor so the output is independent of the screen
        resolution and can be saved as a true vector PDF/SVG/PGF.

        Parameters
        ----------
        path : str
            Output file path.  Extension determines the format
            (png / pdf / svg / pgf).
        show_colorbar : bool
            Include the colorbar.
        show_labels : bool
            Include x/y axis labels and tick labels.
        show_title : bool
            Include a title showing quantity name and simulation time.
        dpi : int
            Raster resolution (only relevant for PNG output).
        tex_engine : str
            LaTeX engine for PGF output: 'pdflatex', 'xelatex', or 'lualatex'.
        """
        import numpy as np
        import matplotlib
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg

        sp   = self._splotter
        qty  = self._last_gen_qty
        vmin = self._last_gen_vmin
        vmax = self._last_gen_vmax
        cmap = self._last_gen_cmap
        use_basemap = self._gen_used_basemap and (sp.epsg is not None)
        triang = sp.triang_abs if use_basemap else sp.triang

        # --- resolve data array for the current frame ---
        is_static = qty.startswith('max_') or (
            qty == 'elev' and sp.elev.ndim == 1)

        if qty == 'elev':
            raw = sp.elev
            data = raw if raw.ndim == 1 else raw[self._frame_to_sww_idx(), :]
            elev_bg = None
        else:
            data_map = {
                'depth':           sp.depth,
                'stage':           sp.stage,
                'speed':           sp.speed,
                'speed_depth':     sp.speed_depth,
                'max_depth':       sp.depth,
                'max_speed':       sp.speed,
                'max_speed_depth': sp.speed_depth,
            }
            raw = data_map[qty]
            if is_static:
                data = np.max(raw, axis=0)
            else:
                data = raw[self._frame_to_sww_idx(), :]

            try:
                elev_bg = sp.elev[0, :] if sp.elev.ndim == 2 else sp.elev
            except (IndexError, TypeError):
                elev_bg = sp.elev

        # --- build figure ---
        is_pgf = path.lower().endswith('.pgf')
        fig = Figure(figsize=(10, 6))
        if is_pgf:
            matplotlib.rcParams['pgf.texsystem'] = tex_engine
            try:
                from matplotlib.backends.backend_pgf import FigureCanvasPgf
                FigureCanvasPgf(fig)
            except Exception:
                FigureCanvasAgg(fig)   # fall back; savefig will raise below
        else:
            FigureCanvasAgg(fig)

        ax = fig.add_subplot(111)

        from anuga.utilities.animate import _face_to_vertex
        import matplotlib.tri as _mtri
        import numpy as _np
        alpha = max(0.0, min(1.0, float(self._alpha_var.get())))
        smooth = not show_edges
        if smooth:
            wet_mask = (data < sp.min_depth) if elev_bg is not None else None

            if use_basemap and elev_bg is not None:
                # With a basemap, Gouraud gradient creases are visible against
                # the high-contrast imagery.  Interpolate to a regular grid and
                # use imshow — masked (dry) pixels are transparent.
                from matplotlib.tri import LinearTriInterpolator as _LTI
                triang_wet = _mtri.Triangulation(
                    triang.x, triang.y, triang.triangles, mask=wet_mask)
                v_data = _face_to_vertex(triang, data, face_mask=wet_mask)
                x0, x1 = float(triang.x.min()), float(triang.x.max())
                y0, y1 = float(triang.y.min()), float(triang.y.max())
                _span = max(x1 - x0, y1 - y0)
                nx = max(100, int(400 * (x1 - x0) / _span))
                ny = max(100, int(400 * (y1 - y0) / _span))
                xi = _np.linspace(x0, x1, nx)
                yi = _np.linspace(y0, y1, ny)
                Xi, Yi = _np.meshgrid(xi, yi)
                zi = _LTI(triang_wet, v_data)(Xi, Yi)
                im = ax.imshow(zi, extent=[x0, x1, y0, y1], origin='lower',
                               cmap=cmap, vmin=vmin, vmax=vmax, alpha=alpha,
                               interpolation='bilinear', zorder=1)
            else:
                # No basemap: Gouraud on separate Triangulations (works cleanly
                # against the flat grey elevation background).
                dry_mask = (data > sp.min_depth) if elev_bg is not None else None
                if elev_bg is not None:
                    triang_dry = _mtri.Triangulation(
                        triang.x, triang.y, triang.triangles, mask=dry_mask)
                    v_elev_bg = _face_to_vertex(triang, elev_bg, face_mask=dry_mask)
                    dry_elev = _np.asarray(elev_bg)[~_np.asarray(dry_mask, dtype=bool)]
                    e_vmin = float(dry_elev.min()) if dry_elev.size else 0.0
                    e_vmax = float(dry_elev.max()) if dry_elev.size else 1.0
                    if e_vmin >= e_vmax:
                        e_vmax = e_vmin + 1.0
                    ax.tripcolor(triang_dry, v_elev_bg, shading='gouraud', cmap='Greys_r',
                                 vmin=e_vmin, vmax=e_vmax).set_linewidths(0)
                if elev_bg is not None:
                    triang_wet = _mtri.Triangulation(
                        triang.x, triang.y, triang.triangles, mask=wet_mask)
                    v_data = _face_to_vertex(triang, data, face_mask=wet_mask)
                else:
                    triang_wet = _mtri.Triangulation(triang.x, triang.y, triang.triangles)
                    v_data = _face_to_vertex(triang, data)
                im = ax.tripcolor(triang_wet, v_data, shading='gouraud',
                                  cmap=cmap, vmin=vmin, vmax=vmax)
                im.set_linewidths(0)
        else:
            if elev_bg is not None and not use_basemap:
                mask_bg = data > sp.min_depth
                triang.set_mask(mask_bg)
                ax.tripcolor(triang, facecolors=elev_bg, cmap='Greys_r')
            if elev_bg is not None:
                triang.set_mask(data < sp.min_depth)
            else:
                triang.set_mask(None)
            im = ax.tripcolor(triang, facecolors=data, cmap=cmap,
                              vmin=vmin, vmax=vmax)
        triang.set_mask(None)
        ax.set_aspect('equal')

        if self._zoom_xlim is not None:
            ax.set_xlim(self._zoom_xlim)
        if self._zoom_ylim is not None:
            ax.set_ylim(self._zoom_ylim)

        if use_basemap:
            from anuga.utilities.animate import _add_basemap, BASEMAP_PROVIDERS
            provider_label = self._basemap_provider_var.get()
            provider_str = BASEMAP_PROVIDERS.get(
                provider_label, 'OpenStreetMap.Mapnik')
            _add_basemap(ax, sp.epsg, provider_str, cache=sp._basemap_cache)

        if show_mesh:
            ax.triplot(triang, color='black', linewidth=0.25, alpha=0.45)

        if show_elev:
            elev_raw = sp.elev
            elev_face = elev_raw[0] if elev_raw.ndim == 2 else elev_raw
            elev_v = _face_to_vertex(triang, elev_face)
            nice_levels = self._nice_contour_levels(
                float(elev_v.min()), float(elev_v.max()), elev_levels)
            cs_elev = ax.tricontour(triang, elev_v, levels=nice_levels,
                                    colors='dimgrey', linewidths=0.6, alpha=0.55)
            ax.clabel(cs_elev, inline=True, fontsize=7, fmt='%g')

        if show_colorbar:
            fig.colorbar(im, ax=ax,
                         label=_QTY_CBAR_LABEL.get(qty, qty))

        if show_labels:
            ax.set_xlabel('Easting (m)')
            ax.set_ylabel('Northing (m)')
        else:
            ax.set_xticklabels([])
            ax.set_yticklabels([])

        if show_title:
            if is_static:
                ax.set_title(f'{_QTY_TITLE.get(qty, qty)}  -  {self._sww_prefix}')
            else:
                t = sp.time[self._frame_to_sww_idx()]
                ax.set_title(
                    f'{_QTY_TITLE.get(qty, qty)}  -  '
                    f'{self._sww_prefix}  t = {t:.1f} s')

        try:
            fig.savefig(path, dpi=dpi, bbox_inches='tight')
        except FileNotFoundError as exc:
            if any(e in str(exc).lower()
                   for e in ('latex', 'xelatex', 'pdflatex', 'lualatex')):
                raise RuntimeError(
                    f'PGF export failed: {tex_engine!r} was not found.\n\n'
                    'Make sure a LaTeX distribution is installed '
                    '(TeX Live on Linux/macOS, MiKTeX on Windows) '
                    'and the binaries are on your PATH.\n\n'
                    'Then choose the matching engine in the "TeX engine" '
                    'selector and try again.'
                ) from exc
            raise

    def _frame_to_sww_idx(self):
        """Map self._current (frame index) to the SWW time-step index."""
        sp = self._splotter
        n_frames = len(self._frames)
        n_sww    = len(sp.time)
        if n_frames <= 1 or n_sww <= 1:
            return 0
        frac = self._current / (n_frames - 1)
        return int(round(frac * (n_sww - 1)))

    def _show_mesh(self):
        """Open a Toplevel window showing the triangulation."""
        if self._splotter is None:
            return
        sp = self._splotter
        _bm_available = sp.epsg is not None

        win = tk.Toplevel(self.root)
        win.title(f'Mesh - {len(sp.triangles)} triangles')
        win.resizable(True, True)

        from matplotlib.figure import Figure
        mesh_fig = Figure(figsize=(8, 6))
        ax = mesh_fig.add_subplot(111)
        mesh_canvas = FigureCanvasTkAgg(mesh_fig, master=win)
        mesh_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        bm_var = tk.BooleanVar(value=self._gen_used_basemap and _bm_available)

        def _redraw_mesh(*_):
            use_basemap = bm_var.get() and _bm_available
            triang = sp.triang_abs if use_basemap else sp.triang
            ax.cla()
            ax.triplot(triang, color='steelblue', linewidth=0.4, alpha=0.7)
            ax.set_aspect('equal')
            n_tri = len(sp.triangles)
            if use_basemap:
                ax.set_xlabel('Easting (m)')
                ax.set_ylabel('Northing (m)')
                ax.set_title(f'Mesh  ({n_tri} triangles)')
                from anuga.utilities.animate import BASEMAP_PROVIDERS, _add_basemap
                provider_label = self._basemap_provider_var.get()
                provider_str = BASEMAP_PROVIDERS.get(
                    provider_label, 'OpenStreetMap.Mapnik')
                try:
                    _add_basemap(ax, sp.epsg, provider_str,
                                 cache=sp._basemap_cache)
                except Exception as e:
                    ax.set_title(
                        f'Mesh  ({n_tri} triangles)  – basemap failed: {e}')
            else:
                ax.set_xlabel('x (m)')
                ax.set_ylabel('y (m)')
                ax.set_title(f'Mesh  ({n_tri} triangles)')
            mesh_fig.tight_layout()
            mesh_canvas.draw()

        _redraw_mesh()

        def _close_mesh():
            plt.close(mesh_fig)
            win.destroy()

        btn_row = ttk.Frame(win)
        btn_row.pack(pady=4)

        bm_chk = ttk.Checkbutton(btn_row, text='Basemap', variable=bm_var,
                                  command=_redraw_mesh)
        bm_chk.pack(side=tk.LEFT, padx=6)
        if not _bm_available:
            bm_chk.config(state='disabled')

        def _save_from_viewer():
            self._save_mesh()

        ttk.Button(btn_row, text='Save Mesh...',
                   command=_save_from_viewer).pack(side=tk.LEFT, padx=6)
        ttk.Button(btn_row, text='Close',
                   command=_close_mesh).pack(side=tk.LEFT, padx=6)
        win.protocol('WM_DELETE_WINDOW', _close_mesh)

    def _show_help(self):
        """Open a scrollable help window."""
        win = tk.Toplevel(self.root)
        win.title('ANUGA SWW Animation GUI - Help')
        win.resizable(True, True)

        text = tk.Text(win, wrap=tk.WORD, width=72, height=42,
                       padx=10, pady=8, relief=tk.FLAT)
        sb = ttk.Scrollbar(win, command=text.yview)
        text.configure(yscrollcommand=sb.set)
        ttk.Button(win, text='Close', command=win.destroy).pack(side=tk.BOTTOM, pady=6)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        text.tag_configure('h1', font=('TkDefaultFont', 13, 'bold'), spacing3=4)
        text.tag_configure('h2', font=('TkDefaultFont', 10, 'bold'), spacing1=6, spacing3=2)
        text.tag_configure('kw', font=('TkFixedFont', 9))

        def h1(s): text.insert(tk.END, s + '\n', 'h1')
        def h2(s): text.insert(tk.END, s + '\n', 'h2')
        def p(s):  text.insert(tk.END, s + '\n')
        def kw(s): text.insert(tk.END, s, 'kw')
        def nl():  text.insert(tk.END, '\n')

        h1('ANUGA SWW Animation GUI')
        p('Visualise and explore SWW simulation output as animated frames.')
        nl()

        h2('Quick-start workflow')
        p('1. Open an SWW file with Browse... or pass --sww on the command line.')
        p('2. On the Plot tab choose a Quantity; vmin/vmax default to Auto.')
        p('3. Click Generate Frames - PNGs are written to the Output dir.')
        p('4. Use Play / step buttons to animate, or drag the Frame slider.')
        nl()

        h2('Interface layout')
        p('The GUI uses three tabs plus always-visible bars:')
        p('  Plot tab     - quantity, colour range, colormap, overlays.')
        p('  Generate tab - output directory, DPI, stride, EPSG, basemap.')
        p('  Output tab   - Save/Export frame, Save Animation, View/Save Mesh.')
        p('Always visible: file row, Generate Frames bar, playback bar,')
        p('  frame slider, status bar.')
        nl()

        h2('Plot tab settings')
        h2('  Quantity')
        p('  depth, stage, speed, speed_depth - animated per timestep.')
        p('  max_depth, max_speed, max_speed_depth - single frame showing')
        p('  the maximum value at each triangle over all timesteps.')
        p('  elev - elevation (bed level).  Produces a single static frame')
        p('  when elevation is constant, or one frame per timestep when the')
        p('  SWW contains time-varying elevation (e.g. erosion simulations).')
        nl()
        h2('  vmin / vmax / Auto')
        p('  Colormap range.  "Auto from data" (ticked by default) sets the')
        p('  range automatically from the data in the current generation run.')
        nl()
        h2('  Colormap / Reverse')
        p('  Any matplotlib colormap name.  Tick Reverse to invert it.')
        nl()
        h2('  min depth')
        p('  Triangles with depth below this value are treated as dry and')
        p('  shown in grey (elevation shading).')
        nl()
        h2('  Flat View')
        p('  When unchecked (the default) the mesh is rendered with smooth')
        p('  Gouraud shading - colour is interpolated across each triangle')
        p('  so no triangle edges are visible.  Tick "Flat View" to revert')
        p('  to piecewise-constant (flat) rendering where each triangle is')
        p('  a single uniform colour.')
        nl()
        h2('  Show Mesh')
        p('  Tick to overlay the triangulation on every generated frame.')
        p('  The mesh is baked into the PNG so it appears at the correct')
        p('  z-order even over basemap tiles.  The canvas overlay toggle')
        p('  in the playback bar works independently for quick preview.')
        nl()
        h2('  Show Elev / levels')
        p('  Tick to overlay elevation contours on every generated frame.')
        p('  "levels" sets the target number of contour lines (rounded to')
        p('  a round step value).  Like Show Mesh, the contours are baked')
        p('  into the PNG during generation.')
        nl()

        h2('Generate tab settings')
        h2('  DPI')
        p('  Resolution of generated PNG frames (default 100).')
        nl()
        h2('  Every N frames')
        p('  Stride: generate one PNG every N SWW timesteps.')
        p('  Use a larger value for a quick preview of a long simulation.')
        nl()
        h2('  EPSG / Set')
        p('  Override or supply the coordinate-system code.  Older SWW')
        p('  files do not store an EPSG code - type the integer (e.g. 32756')
        p('  for UTM zone 56 S) and press Set or Enter to enable basemap.')
        p('  If the file already carries an EPSG it is pre-populated.')
        nl()
        h2('  Basemap / provider / Alpha')
        p('  Overlay an online tile basemap (OpenStreetMap, Esri Satellite,')
        p('  etc.).  Requires an EPSG code (from file or entered manually)')
        p('  and an internet connection.  Alpha controls mesh transparency.')
        p('  When the SWW file contains an EPSG code and contextily is')
        p('  installed, the basemap is enabled automatically on file load.')
        p('  Requires contextily:  pip install contextily')
        nl()

        h2('Generate Frames bar')
        p('  Generate Frames - start frame generation (parallelised).')
        p('  Cancel          - abort an in-progress generation.')
        p('  Progress label  - shows n / total frames completed.')
        nl()

        h2('Playback controls')
        p('  Play/Pause  - start or stop the animation.')
        p('  |< < > >|   - jump to first, step back, step forward, last frame.')
        p('  FPS         - playback speed in frames per second.')
        p('  Show Mesh   - canvas overlay toggle (quick preview only).')
        p('  Frame slider - drag to scrub through frames.')
        nl()

        h2('Pick timeseries')
        p('Click "Pick timeseries" to enter pick mode:')
        p('  * The cursor changes to a crosshair.')
        p('  * Click any point on the image to add the nearest triangle')
        p('    centroid to the timeseries panel.  Each picked point gets')
        p('    a distinct colour; its star marker appears on every frame.')
        p('  * Click further points to add more series to the same plot.')
        p('  * Press Esc or click "Cancel pick" to exit pick mode.')
        p('  * Click "Clear picks" to remove all picks and close the panel.')
        nl()

        h2('Timeseries panel')
        p('  Quantity dropdown - switch the plotted variable for all picks.')
        p('  Dashed lines     - mark the current animation frame time.')
        p('  Export CSV       - save all picked time series to one CSV file.')
        p('  Clear picks      - reset all picks and hide the panel.')
        p('  Close            - hide the panel (picks are also cleared).')
        nl()

        h2('Zoom region')
        p('Click "Set Zoom" to enter rubber-band selection mode:')
        p('  * Drag a rectangle on the animation frame to select a region.')
        p('  * A yellow highlight shows the selected area.')
        p('  * The status bar shows the mesh coordinate bounds of the selection.')
        p('  * Click "Generate Frames" to regenerate at full resolution for')
        p('    the selected region only.')
        p('  * Click "Reset Zoom" to clear the selection and return to')
        p('    full-extent generation.')
        p('  Note: Set Zoom is only available after frames have been generated.')
        nl()

        h2('Output tab')
        p('  Save Frame      - saves the current frame (with any pick-marker')
        p('  overlay).  Choose "Current frame" or enter a list of times (s)')
        p('  to save multiple frames at once.')
        nl()
        p('  Export frame... - like Save Frame but with explicit format choice.')
        nl()
        p('  Save Animation... - saves all loaded frames as:')
        p('    GIF - requires Pillow:  pip install Pillow')
        p('    MP4 - requires ffmpeg on PATH; produces smaller, higher-quality')
        p('          files.  MP4 is offered first when ffmpeg is detected.')
        p('          Install:  conda install ffmpeg  or  apt install ffmpeg')
        p('  Playback FPS is used as the animation frame rate.')
        nl()
        p('  View Mesh - open the mesh viewer.  The Basemap checkbox at the')
        p('  bottom of the viewer toggles the tile basemap on/off live.')
        nl()
        p('  Save Mesh... - save the mesh to PNG, PDF, SVG, or PGF.')
        p('  Options include axis labels, title, and a Basemap checkbox.')
        nl()

        h2('Performance')
        p('  Frame generation is automatically parallelised across up to 4')
        p('  CPU cores.  The progress bar advances in chunks as batches of')
        p('  frames complete.  Basemap tiles are cached in memory so the')
        p('  network is only hit once per generation run.')
        nl()

        h2('Optional dependencies')
        p('  Install all GUI extras with:  pip install anuga[gui]')
        p('  or individually:')
        p('    contextily - basemap tile overlay')
        p('    Pillow     - GIF animation export')
        p('    ffmpeg     - MP4 animation export (system package)')
        nl()

        h2('Command-line usage')
        kw('  anuga_sww_gui [--config FILE.toml] [--sww FILE] [--qty QUANTITY] [options]\n')
        nl()
        p('  --config FILE.toml   load settings from a TOML file (CLI args override)')
        nl()
        p('  Render options:')
        p('    --vmin / --vmax FLOAT    colour scale limits (disables auto-scale)')
        p('    --cmap NAME              colormap (viridis, jet, terrain, ...)')
        p('    --cmap-reverse           reverse the colormap')
        p('    --mindepth FLOAT         wet/dry depth threshold (m)')
        p('    --flat-view              flat per-triangle rendering')
        nl()
        p('  Generate options:')
        p('    --outdir PATH            output directory for frame PNGs')
        p('    --dpi INT                frame resolution')
        p('    --stride INT             generate every Nth frame')
        p('    --alpha FLOAT            quantity layer opacity (0–1)')
        p('    --epsg INT               EPSG code (enables basemap)')
        p('    --basemap / --no-basemap toggle basemap (requires --epsg)')
        p('    --basemap-provider NAME  tile provider (OpenStreetMap, Satellite (Esri), ...)')
        nl()

        text.configure(state=tk.DISABLED)

    # -------------------------------------------------------------- #
    # Config file (TOML)                                             #
    # -------------------------------------------------------------- #

    def _current_config(self):
        """Return a dict of all current GUI parameter values."""
        sww = self._sww_var.get().strip()
        epsg_raw = self._epsg_var.get().strip()
        try:
            epsg_val = int(epsg_raw)
        except (ValueError, TypeError):
            epsg_val = None
        cfg = {}
        if sww:
            cfg['sww'] = sww
        cfg['qty']              = self._qty_var.get()
        cfg['vmin']             = float(self._vmin_var.get())
        cfg['vmax']             = float(self._vmax_var.get())
        cfg['cmap']             = self._cmap_var.get()
        cfg['cmap_reverse']     = bool(self._cmap_reverse_var.get())
        cfg['mindepth']         = float(self._mindepth_var.get())
        cfg['flat_view']        = bool(self._show_edges_var.get())
        cfg['outdir']           = self._dir_var.get()
        cfg['dpi']              = int(self._dpi_var.get())
        cfg['stride']           = int(self._stride_var.get())
        cfg['alpha']            = float(self._alpha_var.get())
        if epsg_val is not None:
            cfg['epsg']         = epsg_val
        cfg['basemap']          = bool(self._basemap_var.get())
        cfg['basemap_provider'] = self._basemap_provider_var.get()
        return cfg

    def _save_config(self):
        """Write current settings to a TOML file chosen by the user."""
        path = filedialog.asksaveasfilename(
            title='Save config',
            defaultextension='.toml',
            filetypes=[('TOML config', '*.toml'), ('All files', '*.*')],
            initialfile='anuga_sww_gui.toml')
        if not path:
            return
        try:
            cfg = self._current_config()
            lines = ['# anuga_sww_gui configuration\n', '\n',
                     '[render]\n']
            for key in ('qty', 'vmin', 'vmax', 'cmap', 'cmap_reverse',
                        'mindepth', 'flat_view'):
                if key in cfg:
                    lines.append(f'{key} = {_toml_value(cfg[key])}\n')
            lines += ['\n', '[generate]\n']
            for key in ('outdir', 'dpi', 'stride', 'alpha', 'epsg',
                        'basemap', 'basemap_provider'):
                if key in cfg:
                    lines.append(f'{key} = {_toml_value(cfg[key])}\n')
            if 'sww' in cfg:
                lines += ['\n', '[file]\n',
                          f'sww = {_toml_value(cfg["sww"])}\n']
            with open(path, 'w') as fh:
                fh.writelines(lines)
            self._set_status(f'Config saved to {os.path.basename(path)}')
        except Exception as e:
            self._set_status(f'Failed to save config: {e}')

    def _load_config(self):
        """Load settings from a TOML file chosen by the user."""
        path = filedialog.askopenfilename(
            title='Load config',
            filetypes=[('TOML config', '*.toml'), ('All files', '*.*')])
        if not path:
            return
        try:
            import tomllib
        except ImportError:
            self._set_status('tomllib not available (requires Python 3.11+)')
            return
        try:
            with open(path, 'rb') as fh:
                data = tomllib.load(fh)
        except Exception as e:
            self._set_status(f'Failed to read config: {e}')
            return
        _apply_config_to_gui(data, self)
        self._set_status(f'Config loaded from {os.path.basename(path)}')

    def _set_status(self, msg):
        self._status_var.set(msg)

    def _on_close(self):
        self._cancel_generation()
        self._stop_playback()
        plt.close('all')
        self.root.quit()
        self.root.destroy()


# ------------------------------------------------------------------ #
# Entry point                                                         #
# ------------------------------------------------------------------ #

def main():
    import argparse
    from anuga.utilities.animate import BASEMAP_PROVIDERS

    parser = argparse.ArgumentParser(
        description='ANUGA SWW animation viewer',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Config file
    parser.add_argument('--config', default=None, metavar='FILE',
                        help='TOML config file; CLI args override config values')

    # File / quantity
    parser.add_argument('--sww', default=None,
                        help='SWW file to open on startup')
    parser.add_argument('--qty', default=None,
                        choices=list(_QUANTITIES),
                        help='Quantity to display')

    # Render tab
    parser.add_argument('--vmin', type=float, default=None,
                        help='Colour scale minimum (disables auto-scale)')
    parser.add_argument('--vmax', type=float, default=None,
                        help='Colour scale maximum (disables auto-scale)')
    parser.add_argument('--cmap', default=None,
                        choices=_CMAPS,
                        help='Colormap name')
    parser.add_argument('--cmap-reverse', action='store_true', default=False,
                        help='Reverse the colormap')
    parser.add_argument('--mindepth', type=float, default=None,
                        help='Minimum depth threshold for wet/dry (m)')
    parser.add_argument('--flat-view', action='store_true', default=False,
                        help='Use flat (per-triangle) rendering instead of smooth Gouraud')

    # Generate tab
    parser.add_argument('--outdir', default=None,
                        help='Output directory for generated frame PNGs')
    parser.add_argument('--dpi', type=int, default=None,
                        help='DPI for generated frames')
    parser.add_argument('--stride', type=int, default=None,
                        help='Generate every Nth frame (1 = all frames)')
    parser.add_argument('--alpha', type=float, default=None,
                        help='Quantity layer opacity (0–1)')
    parser.add_argument('--epsg', type=int, default=None,
                        help='EPSG code for coordinate reference system (enables basemap)')
    parser.add_argument('--basemap', dest='basemap', action='store_true',
                        default=None,
                        help='Enable basemap tile overlay (requires --epsg and contextily)')
    parser.add_argument('--no-basemap', dest='basemap', action='store_false',
                        help='Disable basemap tile overlay')
    parser.add_argument('--basemap-provider', default=None,
                        choices=list(BASEMAP_PROVIDERS.keys()),
                        metavar='PROVIDER',
                        help=('Basemap tile provider. Choices: '
                              + ', '.join(BASEMAP_PROVIDERS.keys())))
    parser.add_argument('--ui-scale', type=float, default=None, metavar='FACTOR',
                        help=('Scale the interface (fonts, spacing and on-screen '
                              'plots). Default: detected from the display. Use '
                              'e.g. 3 on a very high-DPI laptop panel, 1 to '
                              'disable. Equivalent to ANUGA_GUI_SCALE.'))

    parser.set_defaults(basemap=None)

    args = parser.parse_args()

    # Feed --ui-scale to the detector via the same env var users can export.
    if args.ui_scale is not None:
        os.environ['ANUGA_GUI_SCALE'] = str(args.ui_scale)

    # Load TOML config (if given); explicit CLI args take precedence.
    cfg = {}
    if args.config:
        try:
            import tomllib
        except ImportError:
            import sys
            sys.exit('--config requires Python 3.11+ (tomllib not available)')
        try:
            with open(args.config, 'rb') as fh:
                raw = tomllib.load(fh)
            # Flatten render/generate/file sections
            cfg.update(raw.get('render',   {}))
            cfg.update(raw.get('generate', {}))
            cfg.update(raw.get('file',     {}))
            for k, v in raw.items():
                if not isinstance(v, dict):
                    cfg[k] = v
        except Exception as e:
            import sys
            sys.exit(f'Failed to read config file: {e}')

    def _r(cli_val, key, default=None):
        """Return cli_val if set, else config value, else default."""
        return cli_val if cli_val is not None else cfg.get(key, default)

    def _rb(cli_val, key, default=False):
        """Boolean resolve: explicit True CLI flag wins; else config; else default."""
        if cli_val:
            return True
        return bool(cfg.get(key, default))

    root = tk.Tk()
    # Scale the UI up on HiDPI displays (Xft.dpi); no-op at normal 96 DPI.
    # Stash the factor so the GUI can scale its minsize to match.
    try:
        from anuga.utilities.tk_scaling import apply_hidpi_scaling
        root.ui_scale = apply_hidpi_scaling(root)
    except Exception:
        root.ui_scale = 1.0
    SWWAnimationGUI(
        root,
        initial_sww=             _r(args.sww,              'sww'),
        initial_qty=             _r(args.qty,              'qty'),
        initial_vmin=            _r(args.vmin,             'vmin'),
        initial_vmax=            _r(args.vmax,             'vmax'),
        initial_cmap=            _r(args.cmap,             'cmap'),
        initial_cmap_reverse=    _rb(args.cmap_reverse,    'cmap_reverse'),
        initial_mindepth=        _r(args.mindepth,         'mindepth'),
        initial_flat_view=       _rb(args.flat_view,       'flat_view'),
        initial_outdir=          _r(args.outdir,           'outdir'),
        initial_dpi=             _r(args.dpi,              'dpi'),
        initial_stride=          _r(args.stride,           'stride'),
        initial_alpha=           _r(args.alpha,            'alpha'),
        initial_epsg=            _r(args.epsg,             'epsg'),
        initial_basemap=         _r(args.basemap,          'basemap'),
        initial_basemap_provider=_r(args.basemap_provider, 'basemap_provider'),
    )
    root.mainloop()


if __name__ == '__main__':
    main()
