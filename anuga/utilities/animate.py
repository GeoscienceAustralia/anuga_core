"""
A module to allow interactive plotting in a Jupyter notebook of quantities and mesh
associated with an ANUGA domain and SWW file.
"""
import functools
import warnings
import numpy as np
import os

try:
    import matplotlib.pyplot as plt
except ImportError:
    raise ImportError(
        "matplotlib is required for anuga.utilities.animate. "
        "Install it with: conda install matplotlib  or  pip install matplotlib"
    )


# Curated basemap providers useful for hydrological/terrain work.
# Keys are human-readable labels; values are contextily provider strings
# in 'Provider' or 'Provider.Style' dot notation.
BASEMAP_PROVIDERS = {
    'OpenStreetMap':          'OpenStreetMap.Mapnik',
    'Satellite (Esri)':       'Esri.WorldImagery',
    'Shaded Relief (Esri)':   'Esri.WorldShadedRelief',
    'Topo Map (Esri)':        'Esri.WorldTopoMap',
    'Topo Map (OpenTopo)':    'OpenTopoMap',
    'Light (CartoDB)':        'CartoDB.Positron',
    'Dark (CartoDB)':         'CartoDB.DarkMatter',
}

# Default provider used when basemap=True but no provider is specified
BASEMAP_DEFAULT = 'OpenStreetMap.Mapnik'


def _resolve_provider(cx, provider_str):
    """Return the contextily provider dict for a dot-notation string."""
    parts = provider_str.split('.', 1)
    obj = cx.providers
    for part in parts:
        obj = obj[part]
    return obj


def _configure_tile_client(cx):
    """Identify ANUGA to tile servers and cache tiles between sessions.

    contextily defaults to ``USER_AGENT = "contextily-" + uuid4().hex`` — a
    fresh random agent every process.  OpenStreetMap's tile usage policy
    requires a stable User-Agent identifying the *application*, and treats
    rotating agents as an attempt to evade rate limits, so it blocks the
    request.  The block arrives as a normal HTTP 200 tile whose *image* reads
    "Access blocked ... osm.wiki/Blocked", so nothing raises and the notice is
    simply drawn onto the map.

    Sending a stable agent naming ANUGA and its project page fixes it, and a
    persistent on-disk cache keeps repeat views (every redraw, every frame of
    an animation) off the volunteer-run servers, which the policy also asks
    for.  Both are best-effort: contextily is a third-party package, so an
    unexpected version just leaves its defaults in place.
    """
    try:
        from anuga import __version__ as _v
    except Exception:
        _v = 'unknown'

    ua = ('ANUGA/%s (hydrodynamic modelling; '
          '+https://github.com/anuga-community/anuga_core)' % _v)

    try:
        import contextily.tile as _tile
        if getattr(_tile, 'USER_AGENT', '').startswith('contextily-'):
            _tile.USER_AGENT = ua
    except Exception:
        pass

    global _TILE_CACHE_SET
    if not _TILE_CACHE_SET:
        try:
            cache_dir = os.path.join(
                os.environ.get('XDG_CACHE_HOME',
                               os.path.join(os.path.expanduser('~'), '.cache')),
                'anuga', 'tiles')
            os.makedirs(cache_dir, exist_ok=True)
            cx.set_cache_dir(cache_dir)
        except Exception:
            pass
        _TILE_CACHE_SET = True


# set_cache_dir() only needs calling once per session
_TILE_CACHE_SET = False


def _add_basemap(ax, epsg, provider=BASEMAP_DEFAULT, cache=None):
    """Overlay a tile basemap on *ax* using contextily.

    Parameters
    ----------
    ax : matplotlib Axes
        Must already contain data plotted in the CRS given by *epsg*.
    epsg : int
        EPSG code of the data coordinate system.
    provider : str
        Dot-notation provider string, e.g. ``'Esri.WorldImagery'`` or
        ``'OpenTopoMap'``.  See :data:`BASEMAP_PROVIDERS` for the curated
        list, or use any key from ``contextily.providers``.
    cache : dict or None
        Optional in-memory cache dict.  On the first call the fetched tile
        image is stored; subsequent calls with the same bounds and provider
        reuse the cached image via ``ax.imshow``, avoiding repeated network
        and disk access.  Pass ``{}`` once per SWW_plotter instance.

    Warns and returns silently if contextily is not installed or tile
    fetching fails.
    """
    try:
        import contextily as cx
    except ImportError:
        warnings.warn(
            "contextily is not installed - basemap will be skipped. "
            "Install it with: conda install contextily  or  pip install contextily",
            stacklevel=3)
        return

    _configure_tile_client(cx)

    xl, xr = ax.get_xlim()
    yb, yt = ax.get_ylim()
    # Round to nearest metre — sub-metre differences never change the tile set
    cache_key = (epsg, provider, round(xl), round(xr), round(yb), round(yt))

    try:
        source = _resolve_provider(cx, provider)

        if cache is not None and cache_key in cache:
            img, extent = cache[cache_key]
            ax.imshow(img, extent=extent, interpolation='bilinear',
                      origin='upper', aspect='equal', zorder=0)
            ax.set_xlim(xl, xr)
            ax.set_ylim(yb, yt)
        else:
            # First call: fetch tiles via contextily and capture the result
            n_before = len(ax.get_images())
            cx.add_basemap(ax, crs=f'EPSG:{epsg}', source=source,
                           attribution_size=6)
            if cache is not None:
                images = ax.get_images()
                if len(images) > n_before:
                    img_artist = images[n_before]
                    arr = img_artist.get_array()
                    # get_array() may return a masked array; store plain ndarray
                    cache[cache_key] = (
                        arr.data if hasattr(arr, 'data') else arr,
                        img_artist.get_extent())
    except Exception as e:
        warnings.warn(f"Basemap '{provider}' could not be fetched: {e}",
                      stacklevel=3)


class Domain_plotter:
    """
    A class to wrap ANUGA domain centroid values for stage, height, elevation
    xmomentunm and ymomentum, and triangulation information.
    """


    def __init__(self, domain, plot_dir='_plot', min_depth=0.01, absolute=False):

        self.plot_dir = plot_dir
        self.make_plot_dir()

        self.zone = domain.geo_reference.zone

        self.min_depth = min_depth

        self.nodes = domain.nodes
        self.triangles = domain.triangles

        self.xllcorner = domain.geo_reference.xllcorner
        self.yllcorner = domain.geo_reference.yllcorner

        if absolute is False:
            self.x = domain.nodes[:, 0]
            self.y = domain.nodes[:, 1]

            self.xc = domain.centroid_coordinates[:, 0]
            self.yc = domain.centroid_coordinates[:, 1]
        else:
            self.x = domain.nodes[:, 0] + self.xllcorner
            self.y = domain.nodes[:, 1] + self.yllcorner

            self.xc = domain.centroid_coordinates[:, 0] + self.xllcorner
            self.yc = domain.centroid_coordinates[:, 1] + self.yllcorner


        import matplotlib.tri as tri
        self.triang = tri.Triangulation(self.x, self.y, self.triangles)

        self.epsg = domain.geo_reference.epsg

        if self.epsg is not None and not absolute:
            self.triang_abs = tri.Triangulation(
                self.x + self.xllcorner,
                self.y + self.yllcorner,
                self.triangles)
        else:
            self.triang_abs = self.triang

        self._basemap_cache = {}

        self.elev = domain.quantities['elevation'].centroid_values
        self.stage = domain.quantities['stage'].centroid_values

        self.xmom = domain.quantities['xmomentum'].centroid_values
        self.ymom = domain.quantities['ymomentum'].centroid_values

        self.friction = domain.quantities['friction'].centroid_values

        self.domain = domain
        self._depth_frame_count = 0
        self._stage_frame_count = 0
        self._speed_frame_count = 0

    # ------------------------------------------------------------------
    # Derived quantities.
    #
    # self.stage/elev/xmom/ymom are references into the domain's centroid
    # arrays, so they follow the evolution.  These are computed on access so
    # they follow it too -- caching them in __init__ froze them at t = 0, and
    # anything sampling e.g. plotter.speed inside an evolve loop silently
    # recorded the initial value at every yieldstep.
    # ------------------------------------------------------------------

    @property
    def depth(self):
        return self.stage - self.elev

    @property
    def xvel(self):
        depth = self.depth
        with np.errstate(invalid='ignore'):
            return np.where(depth > self.min_depth, self.xmom / depth, 0.0)

    @property
    def yvel(self):
        depth = self.depth
        with np.errstate(invalid='ignore'):
            return np.where(depth > self.min_depth, self.ymom / depth, 0.0)

    @property
    def speed(self):
        return np.sqrt(self.xvel**2 + self.yvel**2)

    @property
    def speed_depth(self):
        return self.speed * self.depth


    #------------------------------------------
    # General plots
    #------------------------------------------
    def plot_mesh(self, figsize=(10, 6), dpi=80):

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        ax.triplot(self.triang, linewidth=0.4)

        ax.set_title('Mesh')
        ax.set_aspect('equal')
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')

        #plt.show()

        return fig, ax

    def _depth_frame(self, figsize, dpi, vmin, vmax, cmap='viridis',
                     basemap=False, alpha=1.0,
                     basemap_provider=BASEMAP_DEFAULT):


        name = os.path.basename(self.domain.get_name())
        time = self.domain.get_time()

        md = self.min_depth

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        ax.set_title('Depth: Time {0:0>4}'.format(time))

        triang = self.triang_abs if (basemap and self.epsg) else self.triang
        if not (basemap and self.epsg):
            triang.set_mask(self.depth > md)
            ax.tripcolor(triang, facecolors=self.elev, cmap='Greys_r')

        triang.set_mask(self.depth <= md)
        im = ax.tripcolor(triang,
                      facecolors=self.depth,
                      cmap=cmap,
                      alpha=alpha,
                      vmin=vmin, vmax=vmax)

        triang.set_mask(None)

        ax.set_aspect('equal')
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        fig.colorbar(im, ax=ax)

        if basemap and self.epsg:
            _add_basemap(ax, self.epsg, basemap_provider, cache=self._basemap_cache)

        return fig, ax

    def save_depth_frame(self, figsize=(10, 6), dpi=80,
                         vmin=0.0, vmax=20, cmap='viridis',
                         basemap=False, alpha=1.0,
                     basemap_provider=BASEMAP_DEFAULT):



        plot_dir = self.plot_dir
        name = os.path.basename(self.domain.get_name())
        frame = self._depth_frame_count

        fig, ax = self._depth_frame(figsize, dpi, vmin, vmax, cmap, basemap, alpha, basemap_provider)

        if plot_dir is None:
            fig.savefig(name+'_depth_{0:0>10}.png'.format(frame))
        else:
            fig.savefig(os.path.join(plot_dir, name
                                     + '_depth_{0:0>10}.png'.format(frame)))

        self._depth_frame_count += 1
        fig.clf()
        plt.close()

    def plot_depth_frame(self, figsize=(10, 6), dpi=80,
                         vmin=0.0, vmax=20.0):

        import matplotlib.pyplot as plt

        fig, ax = self._depth_frame(figsize, dpi, vmin, vmax)

        #plt.show()

        return fig, ax

    def make_depth_animation(self):

        import numpy as np
        import glob
        from matplotlib import image, animation
        from matplotlib import pyplot as plt

        plot_dir = self.plot_dir
        name = os.path.basename(self.domain.get_name())
        time = self.domain.get_time()

        if plot_dir is None:
            expression = name+'_depth_*.png'
        else:
            expression = os.path.join(plot_dir, name+'_depth_*.png')
        img_files = sorted(glob.glob(expression))

        figsize = (10, 6)

        fig = plt.figure(figsize=figsize, dpi=80)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis('off')  # so there's not a second set of axes
        im = plt.imshow(image.imread(img_files[0]))

        def init():

            im.set_data(image.imread(img_files[0]))
            return im,

        def animate(i):

            image_i = image.imread(img_files[i])
            im.set_data(image_i)
            return im,

        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=len(img_files), interval=200, blit=True)

        plt.close()

        return anim

    def _stage_frame(self, figsize, dpi, vmin, vmax, cmap='viridis',
                     basemap=False, alpha=1.0,
                     basemap_provider=BASEMAP_DEFAULT):

        name = os.path.basename(self.domain.get_name())
        time = self.domain.get_time()

        md = self.min_depth

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        ax.set_title('Stage: Time {0:0>4}'.format(time))

        triang = self.triang_abs if (basemap and self.epsg) else self.triang
        if not (basemap and self.epsg):
            triang.set_mask(self.depth > md)
            ax.tripcolor(triang, facecolors=self.elev, cmap='Greys_r')

        triang.set_mask(self.depth <= md)
        im = ax.tripcolor(triang,
                      facecolors=self.stage,
                      cmap=cmap,
                      alpha=alpha,
                      vmin=vmin, vmax=vmax)

        triang.set_mask(None)

        ax.set_aspect('equal')
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        fig.colorbar(im, ax=ax)

        if basemap and self.epsg:
            _add_basemap(ax, self.epsg, basemap_provider, cache=self._basemap_cache)

        return fig, ax

    def save_stage_frame(self, figsize=(10, 6), dpi=80,
                         vmin=-20.0, vmax=20.0, cmap='viridis',
                         basemap=False, alpha=1.0,
                     basemap_provider=BASEMAP_DEFAULT):

        import matplotlib.pyplot as plt

        plot_dir = self.plot_dir
        name = os.path.basename(self.domain.get_name())
        frame = self._stage_frame_count

        fig, ax = self._stage_frame(figsize, dpi, vmin, vmax, cmap, basemap, alpha, basemap_provider)

        if plot_dir is None:
            fig.savefig(name+'_stage_{0:0>10}.png'.format(frame))
        else:
            fig.savefig(os.path.join(plot_dir, name
                                     + '_stage_{0:0>10}.png'.format(frame)))

        self._stage_frame_count += 1
        fig.clf()
        plt.close()

    def plot_stage_frame(self, figsize=(10, 6), dpi=80,
                         vmin=-20.0, vmax=20.0):

        import matplotlib.pyplot as plt

        fig, ax = self._stage_frame(figsize, dpi, vmin, vmax)

        #plt.show()

        return fig, ax

    def make_stage_animation(self):

        import numpy as np
        import glob
        from matplotlib import image, animation
        from matplotlib import pyplot as plt

        plot_dir = self.plot_dir
        name = os.path.basename(self.domain.get_name())
        time = self.domain.get_time()

        if plot_dir is None:
            expression = name+'_stage_*.png'
        else:
            expression = os.path.join(plot_dir, name+'_stage_*.png')

        img_files = sorted(glob.glob(expression))

        figsize = (10, 6)

        fig = plt.figure(figsize=figsize, dpi=80)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis('off')  # so there's not a second set of axes
        im = plt.imshow(image.imread(img_files[0]))

        def init():

            im.set_data(image.imread(img_files[0]))
            return im,

        def animate(i):

            image_i = image.imread(img_files[i])
            im.set_data(image_i)
            return im,

        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=len(img_files), interval=200, blit=True)

        plt.close()

        return anim

    def _speed_frame(self, figsize, dpi, vmin, vmax, cmap='viridis',
                     basemap=False, alpha=1.0,
                     basemap_provider=BASEMAP_DEFAULT):


        name = os.path.basename(self.domain.get_name())
        time = self.domain.get_time()

        md = self.min_depth

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        ax.set_title('Speed: Time {0:0>4}'.format(time))

        triang = self.triang_abs if (basemap and self.epsg) else self.triang
        if not (basemap and self.epsg):
            triang.set_mask(self.depth > md)
            ax.tripcolor(triang, facecolors=self.elev, cmap='Greys_r')

        triang.set_mask(self.depth <= md)
        im = ax.tripcolor(triang,
                      facecolors=self.speed,
                      cmap=cmap,
                      alpha=alpha,
                      vmin=vmin, vmax=vmax)

        triang.set_mask(None)

        ax.set_aspect('equal')
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        fig.colorbar(im, ax=ax)

        if basemap and self.epsg:
            _add_basemap(ax, self.epsg, basemap_provider, cache=self._basemap_cache)

        return fig, ax

    def save_speed_frame(self, figsize=(10, 6), dpi=80,
                         vmin=-20.0, vmax=20.0, cmap='viridis',
                         basemap=False, alpha=1.0,
                     basemap_provider=BASEMAP_DEFAULT):

        import matplotlib.pyplot as plt

        plot_dir = self.plot_dir
        name = os.path.basename(self.domain.get_name())
        frame = self._speed_frame_count

        fig, ax = self._speed_frame(figsize, dpi, vmin, vmax, cmap, basemap, alpha, basemap_provider)

        if plot_dir is None:
            fig.savefig(name+'_speed_{0:0>10}.png'.format(frame))
        else:
            fig.savefig(os.path.join(plot_dir, name
                                     + '_speed_{0:0>10}.png'.format(frame)))

        self._speed_frame_count += 1
        fig.clf()
        plt.close()

    def plot_speed_frame(self, figsize=(5, 3), dpi=80,
                         vmin=-20.0, vmax=20.0):

        import matplotlib.pyplot as plt

        fig, ax = self._speed_frame(figsize, dpi, vmin, vmax)

        #plt.show()

        return fig, ax

    def make_speed_animation(self):

        import numpy as np
        import glob
        from matplotlib import image, animation
        from matplotlib import pyplot as plt

        plot_dir = self.plot_dir
        name = os.path.basename(self.domain.get_name())
        time = self.domain.get_time()

        if plot_dir is None:
            expression = name+'_speed_*.png'
        else:
            expression = os.path.join(plot_dir, name+'_speed_*.png')

        img_files = sorted(glob.glob(expression))

        figsize = (10, 6)

        fig = plt.figure(figsize=figsize, dpi=80)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis('off')  # so there's not a second set of axes
        im = plt.imshow(image.imread(img_files[0]))

        def init():

            im.set_data(image.imread(img_files[0]))
            return im,

        def animate(i):

            image_i = image.imread(img_files[i])
            im.set_data(image_i)
            return im,

        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=len(img_files), interval=200, blit=True)

        plt.close()

        return anim

    def make_plot_dir(self, clobber=True):
        """
        Utility function to create a directory for storing a sequence of plot
        files, or if the directory already exists, clear out any old plots.
        If clobber==False then it will abort instead of deleting existing files.
        """

        plot_dir = self.plot_dir
        if plot_dir is None:
            return
        else:
            if os.path.isdir(plot_dir):
                if clobber:
                    try:
                        os.remove("%s/*" % plot_dir)
                    except OSError:
                        pass
                else:
                    raise OSError(
                        '*** Cannot clobber existing directory %s' % plot_dir)
            else:
                os.mkdir("%s" % plot_dir)
                print("Figure files for each frame will be stored in " + plot_dir)


def _face_to_vertex(triang, face_values, face_mask=None):
    """Average per-triangle (centroid) scalar values onto mesh vertices.

    Used to convert flat-shaded centroid data to per-vertex values for
    Gouraud (smooth) shading in tripcolor.

    Parameters
    ----------
    triang : matplotlib.tri.Triangulation
    face_values : array-like, shape (n_tris,)
    face_mask : boolean array, shape (n_tris,), optional
        If given, only triangles where face_mask is *False* contribute to
        the vertex averages.  (Convention matches Triangulation.set_mask:
        True = masked = excluded.)
    """
    import numpy as np
    n_v = len(triang.x)
    tris = triang.triangles          # shape (n_tris, 3)
    fv   = np.asarray(face_values, dtype=float)
    if face_mask is not None:
        active = ~np.asarray(face_mask, dtype=bool)
        tris = tris[active]
        fv   = fv[active]
    if tris.shape[0] == 0:
        return np.zeros(n_v)
    tris_flat = tris.ravel()
    vals_flat = np.repeat(fv, 3)
    v_sum = np.bincount(tris_flat, weights=vals_flat, minlength=n_v)
    v_cnt = np.bincount(tris_flat, minlength=n_v).astype(float)
    v_cnt[v_cnt == 0] = 1            # guard isolated vertices
    return v_sum / v_cnt


def _nice_contour_levels(vmin, vmax, n):
    """Return elevation contour levels as round multiples spanning [vmin, vmax].

    Picks the smallest step from 1/2/5/10/20/50/100 × order-of-magnitude so
    that levels fall on clean values (e.g. 5 m, 50 m, 500 m).
    """
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


def _draw_elev_contours(ax, triang, elev_face, n_levels):
    """Draw elevation contours + labels onto *ax* using data-space *triang*.

    Saves/restores axis limits to prevent tricontour autoscale side-effects.
    """
    v_elev = _face_to_vertex(triang, elev_face)
    levels = _nice_contour_levels(float(v_elev.min()), float(v_elev.max()), n_levels)
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    cs = ax.tricontour(triang, v_elev, levels=levels,
                       colors='dimgray', linewidths=0.6, alpha=0.7, zorder=3)
    ax.clabel(cs, fmt='%g m', fontsize=6, inline=True, inline_spacing=2)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)


class SWW_plotter:
    """
    A class to wrap ANUGA swwfile centroid values for stage, height, elevation
    xmomentum and ymomentum, and triangulation information.

    Plotting functions are provided to plot depth, speed, speed_depth, stage and mesh,
    and to create animations of depth, speed, speed_depth (momentum) and stage.

    Plotting based on matplotlib's tripcolor and triplot functions.

    Example:

        Open an SWW file and plot various quantities.

    >>> import anuga
    >>> import matplotlib.pyplot as plt
    >>>
    >>> # Enable interactive mode
    >>> plot.ion()
    >>>
    >>> # Create SWW plotter object by opening an SWW file
    >>> splotter = anuga.SWW_plotter('domain.sww')
    >>>
    >>> # Find min and max depth for plotting
    >>> vmin = splotter.depth.min()
    >>> vmax = splotter.depth.max()
    >>>
    >>> # Plot depth at last frame.
    >>> fig, ax = splotter.plot_depth_frame(-1, vmin=vmin, vmax=vmax)
    >>> ax.set_title('Depth at final time')
    >>>
    >>> # Plot speed at second frame
    >>> fig, ax = splotter.plot_speed_frame(1)
    >>> ax.set_title('Speed at second frame')
    >>>
    >>> # Plot Mesh
    >>> fig, ax, im = splotter.plot_mesh()
    >>> ax.set_title('Mesh')
    >>>
    >>> # Animate depth
    >>> for frame in range(len(splotter.time)):
    ...     splotter.save_depth_frame(frame, vmin=vmin, vmax=vmax)
    >>> anim = splotter.make_depth_animation()
    """

    def __init__(self, swwfile='domain.sww', plot_dir='_plot',
                 min_depth = 0.001,
                 absolute=False):

        self.filename = swwfile

        self.plot_dir = plot_dir
        self.make_plot_dir()

        self.min_depth = min_depth

        import matplotlib.tri as tri
        import numpy as np

        self.name = os.path.basename(os.path.splitext(swwfile)[0])

        from anuga.file.netcdf import NetCDFFile
        p = NetCDFFile(swwfile)

        self.x = p.variables['x'][:]
        self.y = p.variables['y'][:]
        self.triangles = p.variables['volumes'][:]

        vols0 = self.triangles[:, 0]
        vols1 = self.triangles[:, 1]
        vols2 = self.triangles[:, 2]

        self.xc = (self.x[vols0]+self.x[vols1]+self.x[vols2])/3.0
        self.yc = (self.y[vols0]+self.y[vols1]+self.y[vols2])/3.0

        self.xllcorner = p.xllcorner
        self.yllcorner = p.yllcorner
        self.zone = p.zone
        self.timezone = p.timezone
        self.starttime = p.starttime

        try:
            self.epsg = int(p.epsg)
        except AttributeError:
            self.epsg = None

        if absolute is True:
            self.x[:] = self.x + self.xllcorner
            self.y[:] = self.y + self.yllcorner

            self.xc[:] = self.xc + self.xllcorner
            self.yc[:] = self.yc + self.yllcorner

        # Built after the absolute shift above: Triangulation copies the
        # coordinate arrays it is given, so constructing it first would leave
        # self.triang in relative coordinates while self.x/self.xc are absolute.
        self.triang = tri.Triangulation(self.x, self.y, self.triangles)

        # Absolute-coordinate triangulation used when drawing a basemap
        if self.epsg is not None and not absolute:
            self.triang_abs = tri.Triangulation(
                self.x + self.xllcorner,
                self.y + self.yllcorner,
                self.triangles)
        else:
            self.triang_abs = self.triang


        self.time = p.variables['time'][:]

        # The (time, triangle) quantities are read and derived a timestep at
        # a time as they are used, rather than all held in memory: see
        # _sww_frames.  A static elevation is small, so it is read now.
        from anuga.utilities._sww_frames import FrameArray, SWWSource
        self._source = source = SWWSource(swwfile)
        shape = (len(self.time), len(self.triangles))

        def quantity(compute):
            return FrameArray(source, compute, *shape)

        min_depth = self.min_depth

        def depth(get):
            return get('stage_c') - get('elevation_c')

        def velocity(get, momentum):
            d = depth(get)
            with np.errstate(invalid='ignore', divide='ignore'):
                return np.where(d > min_depth, get(momentum) / d, 0.0)

        def speed(get):
            return np.sqrt(velocity(get, 'xmomentum_c')**2
                           + velocity(get, 'ymomentum_c')**2)

        if p.variables['elevation_c'].ndim == 2:
            self.elev = quantity(lambda get: get('elevation_c'))
        else:
            self.elev = p.variables['elevation_c'][:]
        self.stage = quantity(lambda get: get('stage_c'))
        self.xmom = quantity(lambda get: get('xmomentum_c'))
        self.ymom = quantity(lambda get: get('ymomentum_c'))
        self.depth = quantity(depth)
        self.xvel = quantity(lambda get: velocity(get, 'xmomentum_c'))
        self.yvel = quantity(lambda get: velocity(get, 'ymomentum_c'))
        self.speed = quantity(speed)
        self.speed_depth = quantity(lambda get: speed(get) * depth(get))

        # Tracers -- and suspended sediment classes, which are tracers -- are
        # stored one dynamic centroid variable each. The file names them in a
        # tracer_names attribute, since <name>_c is indistinguishable from
        # stage_c otherwise.
        from anuga.utilities.plot_utils import _tracer_names_from
        self.tracers = {}
        self.max_tracers = {}
        self._tracer_frame_count = {}
        for _tname in _tracer_names_from(p):
            if _tname + '_c' not in p.variables:
                continue
            self.tracers[_tname] = quantity(
                lambda get, name=_tname + '_c': get(name))
            self._tracer_frame_count[_tname] = 0
            # Named attribute and bound save method per tracer, so a caller
            # that drives this class by attribute name -- anuga_sww_gui does
            # -- reaches a tracer the same way it reaches stage.
            setattr(self, 'tracer_' + _tname, self.tracers[_tname])
            setattr(self, 'save_tracer_%s_frame' % _tname,
                    functools.partial(self._save_tracer_frame, _tname))
            if 'max_' + _tname + '_c' in p.variables:
                self.max_tracers[_tname] = p.variables['max_' + _tname + '_c'][:]
                setattr(self, 'max_tracer_' + _tname, self.max_tracers[_tname])
                setattr(self, 'save_max_tracer_%s_frame' % _tname,
                        functools.partial(self._save_max_tracer_frame, _tname))

        # Load precomputed max quantities if stored by Collect_max_quantities_operator
        pvars = p.variables
        self._max_depth_c      = pvars['max_depth_c'][:] if 'max_depth_c' in pvars else None
        self._max_speed_c      = pvars['max_speed_c'][:] if 'max_speed_c' in pvars else None
        self._max_uh_c         = pvars['max_uh_c'][:]    if 'max_uh_c'    in pvars else None
        self._max_stage_c      = pvars['max_stage_c'][:] if 'max_stage_c' in pvars else None

        self._depth_frame_count = 0
        self._stage_frame_count = 0
        self._speed_depth_frame_count = 0
        self._speed_frame_count = 0
        self._max_depth_frame_count = 0
        self._max_stage_frame_count = 0
        self._max_speed_frame_count = 0
        self._max_speed_depth_frame_count = 0
        self._elev_frame_count = 0
        self._elev_delta_frame_count = 0

        self._basemap_cache = {}
        self._figure_cache = {}

        p.close()

    def close(self):
        """Close the SWW file.  It is opened again if a quantity is read."""
        self._source.close()

    #------------------------------------------
    # General plots
    #------------------------------------------
    def plot_mesh(self, figsize=(10, 8), dpi=80, **kwargs):

        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        im = ax.triplot(self.triang, linewidth=0.4, **kwargs)

        ax.set_title('Mesh')
        ax.set_aspect('equal')
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')

        #plt.show()

        return fig, ax, im

    #------------------------------------------
    # Frame rendering helpers (figure reuse)
    #------------------------------------------


    #------------------------------------------
    # Tracer procedures
    #
    # One pair of methods for every tracer rather than a hand-written pair
    # per name: a tracer is user-defined, so there is no fixed list to write
    # out. __init__ binds `save_tracer_<name>_frame` for each.
    #------------------------------------------
    def _tracer_frame(self, tracer, frame, figsize, dpi, vmin, vmax,
                      cmap='viridis', basemap=False, alpha=1.0,
                      basemap_provider=BASEMAP_DEFAULT,
                      xlim=None, ylim=None, smooth=False,
                      show_elev=False, elev_levels=10, show_mesh=False):
        return self._animated_frame(
            frame, 'tracer_' + tracer, self.tracers[tracer][frame, :],
            tracer, figsize, dpi, vmin, vmax, cmap, basemap, alpha,
            basemap_provider, xlim=xlim, ylim=ylim, smooth=smooth,
            show_elev=show_elev, elev_levels=elev_levels, show_mesh=show_mesh)

    def _save_tracer_frame(self, tracer, frame=-1, figsize=(10, 6), dpi=160,
                           vmin=0.0, vmax=1.0, cmap='viridis', basemap=False,
                           alpha=1.0, basemap_provider=BASEMAP_DEFAULT,
                           xlim=None, ylim=None, smooth=False,
                           show_elev=False, elev_levels=10, show_mesh=False):
        frame_num = self._tracer_frame_count[tracer]
        fig, ax = self._tracer_frame(tracer, frame, figsize, dpi, vmin, vmax,
                                     cmap, basemap, alpha, basemap_provider,
                                     xlim=xlim, ylim=ylim, smooth=smooth,
                                     show_elev=show_elev,
                                     elev_levels=elev_levels,
                                     show_mesh=show_mesh)
        fname = '%s_tracer_%s_%010d.png' % (self.name, tracer, frame_num)
        if self.plot_dir is None:
            fig.savefig(fname)
        else:
            fig.savefig(os.path.join(self.plot_dir, fname))
        plt.close(fig)
        self._tracer_frame_count[tracer] += 1

    def _save_max_tracer_frame(self, tracer, frame=None, figsize=(10, 6),
                               dpi=160, vmin=0.0, vmax=1.0, cmap='viridis',
                               basemap=False, alpha=1.0,
                               basemap_provider=BASEMAP_DEFAULT,
                               xlim=None, ylim=None, smooth=False,
                               show_elev=False, elev_levels=10,
                               show_mesh=False):
        """The running maximum, if the run stored one. Time-independent, so
        `frame` is accepted and ignored, as for the other max quantities."""
        data = self.max_tracers[tracer]
        # Stored flag-4, i.e. overwritten each yield step: the last slice is
        # the maximum over the whole run.
        values = data[-1, :] if data.ndim == 2 else data
        fig, ax = self._animated_frame(
            -1, 'max_tracer_' + tracer, values, 'Max ' + tracer,
            figsize, dpi, vmin, vmax, cmap, basemap, alpha, basemap_provider,
            xlim=xlim, ylim=ylim, smooth=smooth, show_elev=show_elev,
            elev_levels=elev_levels, show_mesh=show_mesh)
        fname = self.name + '_max_tracer_%s_0000000000.png' % tracer
        if self.plot_dir is None:
            fig.savefig(fname)
        else:
            fig.savefig(os.path.join(self.plot_dir, fname))
        plt.close(fig)

    def _animated_frame(self, frame, qty_name, qty_data, qty_label,
                        figsize, dpi, vmin, vmax, cmap, basemap, alpha,
                        basemap_provider, xlim=None, ylim=None,
                        smooth=False, show_elev=False, elev_levels=10,
                        show_mesh=False):
        """Shared core for depth/stage/speed/speed_depth frame rendering.

        smooth=True uses Gouraud shading (per-vertex interpolation) which
        eliminates visible triangle-edge artefacts at the cost of slightly
        blurring sharp gradients.
        """
        import matplotlib.pyplot as plt

        time = self.time[frame]
        depth = self.depth[frame, :]
        try:
            elev = self.elev[frame, :]
        except IndexError:
            elev = self.elev
        md = self.min_depth

        cache_key = (qty_name, figsize, dpi, basemap, basemap_provider,
                     cmap, vmin, vmax, alpha,
                     tuple(xlim) if xlim is not None else None,
                     tuple(ylim) if ylim is not None else None,
                     show_elev, elev_levels, show_mesh)

        if cache_key not in self._figure_cache:
            self._figure_cache[cache_key] = plt.figure(figsize=figsize, dpi=dpi)

        fig = self._figure_cache[cache_key]
        fig.clf()
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        ax.set_title(f'{qty_label}: Time {time:.2f}')

        triang = self.triang_abs if (basemap and self.epsg) else self.triang

        # Pin axes to full mesh extent before any plotting so that flat
        # (PolyCollection) and smooth (TriMesh / imshow) modes produce
        # identical axis limits.  TriMesh reports limits from ALL vertex
        # coordinates while PolyCollection only uses unmasked triangle
        # vertices, causing a small but visible jump when toggling show_edges.
        ax.set_xlim(triang.x.min(), triang.x.max())
        ax.set_ylim(triang.y.min(), triang.y.max())

        if basemap and self.epsg:
            # When a basemap is active, always render via imshow regardless of
            # the smooth flag.  Flat tripcolor leaves blocky triangle-edge
            # artefacts at the wet/dry interface because the basemap imagery is
            # visible between adjacent triangles of opposite wetness.  Imshow
            # with a LinearTriInterpolator gives a smooth boundary: dry pixels
            # are masked (alpha=0) so the basemap shows through cleanly.
            #
            # Resolution: compute the interpolation grid over the *visible*
            # region (zoom window if set, otherwise the full mesh extent).
            # This gives proper pixel density even when zoomed into a small
            # coastal area of a large domain.
            import matplotlib.tri as _mtri
            import numpy as _np
            from matplotlib.tri import LinearTriInterpolator as _LTI

            wet_mask = depth < md
            triang_wet = _mtri.Triangulation(
                triang.x, triang.y, triang.triangles, mask=wet_mask)
            # Use face-to-vertex averaging so the interpolator gets smooth data
            v_qty = _face_to_vertex(triang, qty_data, face_mask=wet_mask)

            # Visible region: clip zoom limits to mesh bounds
            mx0 = float(triang.x.min()); mx1 = float(triang.x.max())
            my0 = float(triang.y.min()); my1 = float(triang.y.max())
            x0 = max(mx0, xlim[0]) if xlim is not None else mx0
            x1 = min(mx1, xlim[1]) if xlim is not None else mx1
            y0 = max(my0, ylim[0]) if ylim is not None else my0
            y1 = min(my1, ylim[1]) if ylim is not None else my1
            if x0 >= x1:
                x0, x1 = mx0, mx1
            if y0 >= y1:
                y0, y1 = my0, my1

            # Target ~600 pixels across the longer visible axis
            _span = max(x1 - x0, y1 - y0)
            nx = max(200, int(600 * (x1 - x0) / _span))
            ny = max(200, int(600 * (y1 - y0) / _span))
            xi = _np.linspace(x0, x1, nx)
            yi = _np.linspace(y0, y1, ny)
            Xi, Yi = _np.meshgrid(xi, yi)
            zi = _LTI(triang_wet, v_qty)(Xi, Yi)   # masked outside wet region
            im = ax.imshow(zi, extent=[x0, x1, y0, y1], origin='lower',
                           cmap=cmap, vmin=vmin, vmax=vmax, alpha=alpha,
                           interpolation='bilinear', zorder=1)

        elif smooth:
            import matplotlib.tri as _mtri
            import numpy as _np
            wet_mask = depth < md   # True = dry  → quantity hides dry cells

            # No basemap: Gouraud shading on separate Triangulations.
            # Gouraud needs separate Triangulation objects per collection
            # because TriMesh reads the mask dynamically at draw time.
            dry_mask = depth > md
            triang_dry = _mtri.Triangulation(
                triang.x, triang.y, triang.triangles, mask=dry_mask)
            v_elev = _face_to_vertex(triang, elev, face_mask=dry_mask)
            # Normalise colormap to dry-cell elevation range; wet-only
            # vertices get v_elev=0 from _face_to_vertex which would
            # otherwise skew the auto-scale and make terrain appear white.
            dry_elev = _np.asarray(elev)[~_np.asarray(dry_mask, dtype=bool)]
            e_vmin = float(dry_elev.min()) if dry_elev.size else 0.0
            e_vmax = float(dry_elev.max()) if dry_elev.size else 1.0
            if e_vmin >= e_vmax:
                e_vmax = e_vmin + 1.0
            ax.tripcolor(triang_dry, v_elev, shading='gouraud', cmap='Greys_r',
                         vmin=e_vmin, vmax=e_vmax).set_linewidths(0)
            triang_wet = _mtri.Triangulation(
                triang.x, triang.y, triang.triangles, mask=wet_mask)
            v_qty = _face_to_vertex(triang, qty_data, face_mask=wet_mask)
            im = ax.tripcolor(triang_wet, v_qty, shading='gouraud',
                              cmap=cmap, alpha=alpha, vmin=vmin, vmax=vmax)
            im.set_linewidths(0)
        else:
            # Flat tripcolor, no basemap
            triang.set_mask(depth > md)
            ax.tripcolor(triang, facecolors=elev, cmap='Greys_r')
            triang.set_mask(depth < md)
            im = ax.tripcolor(triang, facecolors=qty_data,
                              cmap=cmap, alpha=alpha, vmin=vmin, vmax=vmax)
        triang.set_mask(None)

        # Apply zoom limits before colorbar/basemap so tiles are fetched
        # for the correct region.
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)

        if show_elev:
            _draw_elev_contours(ax, triang, elev, elev_levels)
        if show_mesh:
            ax.triplot(triang, color='black', linewidth=0.25, alpha=0.45,
                       zorder=4)

        fig.colorbar(im, ax=ax)

        if basemap and self.epsg:
            _add_basemap(ax, self.epsg, basemap_provider,
                         cache=self._basemap_cache)

        return fig, ax

    def _clear_figure_cache(self):
        """Close all cached figures and release their memory."""
        import matplotlib.pyplot as plt
        for fig in self._figure_cache.values():
            plt.close(fig)
        self._figure_cache = {}

    #------------------------------------------
    # Depth procedures
    #------------------------------------------
    def _depth_frame(self, frame, figsize, dpi, vmin, vmax, cmap='viridis',
                     basemap=False, alpha=1.0,
                     basemap_provider=BASEMAP_DEFAULT, xlim=None, ylim=None,
                     smooth=False, show_elev=False, elev_levels=10,
                     show_mesh=False):
        return self._animated_frame(
            frame, 'depth', self.depth[frame, :], 'Depth',
            figsize, dpi, vmin, vmax, cmap, basemap, alpha, basemap_provider,
            xlim=xlim, ylim=ylim, smooth=smooth,
            show_elev=show_elev, elev_levels=elev_levels, show_mesh=show_mesh)

    def save_depth_frame(self, frame=-1, figsize=(10, 6), dpi=160,
                         vmin=0.0, vmax=20.0, cmap='viridis', basemap=False,
                         alpha=1.0, basemap_provider=BASEMAP_DEFAULT,
                         xlim=None, ylim=None, smooth=False,
                         show_elev=False, elev_levels=10, show_mesh=False):

        name = self.name
        plot_dir = self.plot_dir
        frame_num = self._depth_frame_count

        fig, ax = self._depth_frame(frame, figsize, dpi, vmin, vmax, cmap,
                                    basemap, alpha, basemap_provider,
                                    xlim=xlim, ylim=ylim, smooth=smooth,
                                    show_elev=show_elev,
                                    elev_levels=elev_levels,
                                    show_mesh=show_mesh)

        if plot_dir is None:
            fig.savefig(name+'_depth_{0:0>10}.png'.format(frame_num))
        else:
            fig.savefig(os.path.join(plot_dir, name
                                     + '_depth_{0:0>10}.png'.format(frame_num)))
        plt.close(fig)
        self._depth_frame_count += 1

    def plot_depth_frame(self, frame=-1, figsize=(10, 6), dpi = 80,
                         vmin=0.0, vmax=20.0):

        import matplotlib.pyplot as plt

        fig, ax = self._depth_frame(frame, figsize, dpi, vmin, vmax)

        #plt.show()

        return fig, ax


    #------------------------------------------
    # Stage procedures
    #------------------------------------------
    def _stage_frame(self, frame, figsize, dpi, vmin, vmax, cmap='viridis',
                     basemap=False, alpha=1.0,
                     basemap_provider=BASEMAP_DEFAULT, xlim=None, ylim=None,
                     smooth=False, show_elev=False, elev_levels=10,
                     show_mesh=False):
        return self._animated_frame(
            frame, 'stage', self.stage[frame, :], 'Stage',
            figsize, dpi, vmin, vmax, cmap, basemap, alpha, basemap_provider,
            xlim=xlim, ylim=ylim, smooth=smooth,
            show_elev=show_elev, elev_levels=elev_levels, show_mesh=show_mesh)

    def save_stage_frame(self, frame=-1, figsize=(10, 6), dpi=160,
                         vmin=-20.0, vmax=20.0, cmap='viridis', basemap=False,
                         alpha=1.0, basemap_provider=BASEMAP_DEFAULT,
                         xlim=None, ylim=None, smooth=False,
                         show_elev=False, elev_levels=10, show_mesh=False):

        name = self.name
        plot_dir = self.plot_dir
        frame_num = self._stage_frame_count

        fig, ax = self._stage_frame(frame, figsize, dpi, vmin, vmax, cmap,
                                    basemap, alpha, basemap_provider,
                                    xlim=xlim, ylim=ylim, smooth=smooth,
                                    show_elev=show_elev,
                                    elev_levels=elev_levels,
                                    show_mesh=show_mesh)

        if plot_dir is None:
            fig.savefig(name+'_stage_{0:0>10}.png'.format(frame_num))
        else:
            fig.savefig(os.path.join(plot_dir, name
                                     + '_stage_{0:0>10}.png'.format(frame_num)))
        plt.close(fig)
        self._stage_frame_count += 1

    def plot_stage_frame(self, frame=-1, figsize=(5, 3), dpi=80,
                         vmin=-20, vmax=20.0):

        import matplotlib.pyplot as plt

        fig, ax = self._stage_frame(frame, figsize, dpi, vmin, vmax)

        #plt.show()

        return fig, ax

    #------------------------------------------
    # Depth Speed procedures
    #------------------------------------------
    def _speed_depth_frame(self, frame, figsize, dpi, vmin, vmax,
                           cmap='viridis', basemap=False, alpha=1.0,
                           basemap_provider=BASEMAP_DEFAULT,
                           xlim=None, ylim=None, smooth=False,
                           show_elev=False, elev_levels=10, show_mesh=False):
        return self._animated_frame(
            frame, 'speed_depth', self.speed_depth[frame, :], 'Speed x Depth',
            figsize, dpi, vmin, vmax, cmap, basemap, alpha, basemap_provider,
            xlim=xlim, ylim=ylim, smooth=smooth,
            show_elev=show_elev, elev_levels=elev_levels, show_mesh=show_mesh)

    def save_speed_depth_frame(self, frame=-1, figsize=(10, 6), dpi=160,
                               vmin=0.0, vmax=20.0, cmap='viridis',
                               basemap=False, alpha=1.0,
                               basemap_provider=BASEMAP_DEFAULT,
                               xlim=None, ylim=None, smooth=False,
                               show_elev=False, elev_levels=10,
                               show_mesh=False):

        name = self.name
        plot_dir = self.plot_dir
        frame_num = self._speed_depth_frame_count

        fig, ax = self._speed_depth_frame(frame, figsize, dpi, vmin, vmax,
                                          cmap, basemap, alpha, basemap_provider,
                                          xlim=xlim, ylim=ylim,
                                          smooth=smooth,
                                          show_elev=show_elev,
                                          elev_levels=elev_levels,
                                          show_mesh=show_mesh)

        if plot_dir is None:
            fig.savefig(name+'_speed_depth_{0:0>10}.png'.format(frame_num))
        else:
            fig.savefig(os.path.join(plot_dir, name
                                     + '_speed_depth_{0:0>10}.png'.format(frame_num)))
        plt.close(fig)
        self._speed_depth_frame_count += 1

    def plot_speed_depth_frame(self, frame=-1, figsize=(5, 3), dpi=80,
                         vmin=-20, vmax=20.0):

        import matplotlib.pyplot as plt

        fig, ax = self._speed_depth_frame(frame, figsize, dpi, vmin, vmax)

        #plt.show()

        return fig, ax

    #------------------------------------------
    # Speed procedures
    #------------------------------------------
    def _speed_frame(self, frame, figsize, dpi, vmin, vmax, cmap='viridis',
                     basemap=False, alpha=1.0,
                     basemap_provider=BASEMAP_DEFAULT, xlim=None, ylim=None,
                     smooth=False, show_elev=False, elev_levels=10,
                     show_mesh=False):
        return self._animated_frame(
            frame, 'speed', self.speed[frame, :], 'Speed',
            figsize, dpi, vmin, vmax, cmap, basemap, alpha, basemap_provider,
            xlim=xlim, ylim=ylim, smooth=smooth,
            show_elev=show_elev, elev_levels=elev_levels, show_mesh=show_mesh)

    def save_speed_frame(self, frame=-1, figsize=(10, 6), dpi=160,
                         vmin=0.0, vmax=10.0, cmap='viridis', basemap=False,
                         alpha=1.0, basemap_provider=BASEMAP_DEFAULT,
                         xlim=None, ylim=None, smooth=False,
                         show_elev=False, elev_levels=10, show_mesh=False):

        name = self.name
        plot_dir = self.plot_dir
        frame_num = self._speed_frame_count

        fig, ax = self._speed_frame(frame, figsize, dpi, vmin, vmax, cmap,
                                    basemap, alpha, basemap_provider,
                                    xlim=xlim, ylim=ylim, smooth=smooth,
                                    show_elev=show_elev,
                                    elev_levels=elev_levels,
                                    show_mesh=show_mesh)

        if plot_dir is None:
            fig.savefig(name+'_speed_{0:0>10}.png'.format(frame_num))
        else:
            fig.savefig(os.path.join(plot_dir, name
                                     + '_speed_{0:0>10}.png'.format(frame_num)))
        plt.close(fig)
        self._speed_frame_count += 1

    def plot_speed_frame(self, frame=-1, figsize=(10, 6), dpi=80,
                         vmin=0.0, vmax=10.0):

        import matplotlib.pyplot as plt

        fig, ax = self._speed_frame(frame, figsize, dpi, vmin, vmax)

        #plt.show()

        return fig, ax

    #------------------------------------------
    # Elevation procedures (static or time-varying)
    #------------------------------------------

    def _elev_frame(self, frame, figsize, dpi, vmin, vmax, cmap='terrain',
                    basemap=False, alpha=1.0,
                    basemap_provider=BASEMAP_DEFAULT, xlim=None, ylim=None,
                    smooth=False, show_mesh=False):
        """Render one elevation frame.

        Works for both static (1-D) and time-varying (2-D) elevation arrays.
        No wet/dry masking is applied — elevation is shown everywhere.
        """
        import matplotlib.pyplot as plt

        if self.elev.ndim == 2:
            elev_data = self.elev[frame, :]
            title = f'Elevation: Time {self.time[frame]:.2f}'
        else:
            elev_data = self.elev
            title = 'Elevation'

        cache_key = ('elev', figsize, dpi, basemap, basemap_provider,
                     cmap, vmin, vmax, alpha,
                     tuple(xlim) if xlim is not None else None,
                     tuple(ylim) if ylim is not None else None,
                     show_mesh)

        if cache_key not in self._figure_cache:
            self._figure_cache[cache_key] = plt.figure(figsize=figsize, dpi=dpi)

        fig = self._figure_cache[cache_key]
        fig.clf()
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        ax.set_title(title)

        triang = self.triang_abs if (basemap and self.epsg) else self.triang
        triang.set_mask(None)
        if smooth:
            import matplotlib.tri as _mtri
            triang_all = _mtri.Triangulation(triang.x, triang.y, triang.triangles)
            v_elev = _face_to_vertex(triang, elev_data)
            im = ax.tripcolor(triang_all, v_elev, shading='gouraud',
                              cmap=cmap, alpha=alpha, vmin=vmin, vmax=vmax)
            im.set_linewidths(0)
        else:
            im = ax.tripcolor(triang, facecolors=elev_data,
                              cmap=cmap, alpha=alpha, vmin=vmin, vmax=vmax)

        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)

        if show_mesh:
            ax.triplot(triang, color='black', linewidth=0.25, alpha=0.45,
                       zorder=4)

        fig.colorbar(im, ax=ax, label='Elevation (m)')

        if basemap and self.epsg:
            _add_basemap(ax, self.epsg, basemap_provider,
                         cache=self._basemap_cache)

        return fig, ax

    def save_elev_frame(self, frame=-1, figsize=(10, 6), dpi=160,
                        vmin=-20.0, vmax=100.0, cmap='terrain', basemap=False,
                        alpha=1.0, basemap_provider=BASEMAP_DEFAULT,
                        xlim=None, ylim=None, smooth=False,
                        show_elev=False, elev_levels=10, show_mesh=False):
        """Save one elevation frame to disk.

        For static elevation (1-D) the *frame* argument is ignored.
        """
        name = self.name
        plot_dir = self.plot_dir
        frame_num = self._elev_frame_count

        # For static elevation there is only one frame; always render it.
        render_frame = 0 if self.elev.ndim == 1 else frame
        fig, ax = self._elev_frame(render_frame, figsize, dpi, vmin, vmax, cmap,
                                   basemap, alpha, basemap_provider,
                                   xlim=xlim, ylim=ylim, smooth=smooth,
                                   show_mesh=show_mesh)

        fname = name + '_elev_{0:0>10}.png'.format(frame_num)
        if plot_dir is None:
            fig.savefig(fname)
        else:
            fig.savefig(os.path.join(plot_dir, fname))
        plt.close(fig)
        self._elev_frame_count += 1

    def plot_elev_frame(self, frame=-1, figsize=(10, 6), dpi=80,
                        vmin=-20.0, vmax=100.0):

        import matplotlib.pyplot as plt

        render_frame = 0 if self.elev.ndim == 1 else frame
        fig, ax = self._elev_frame(render_frame, figsize, dpi, vmin, vmax)
        return fig, ax

    #------------------------------------------
    # Max-quantity accessors (precomputed or derived)
    #------------------------------------------

    @property
    def max_depth(self):
        """Max depth per triangle: precomputed from SWW if available, else max over time."""
        if self._max_depth_c is not None:
            return self._max_depth_c
        import numpy as np
        return np.max(self.depth, axis=0)

    @property
    def max_speed(self):
        """Max speed per triangle: precomputed from SWW if available, else max over time."""
        if self._max_speed_c is not None:
            return self._max_speed_c
        import numpy as np
        return np.max(self.speed, axis=0)

    @property
    def max_speed_depth(self):
        """Max momentum magnitude per triangle: precomputed (max_uh_c) or max over time."""
        if self._max_uh_c is not None:
            return self._max_uh_c
        import numpy as np
        return np.max(self.speed_depth, axis=0)

    @property
    def max_stage(self):
        """Max stage per triangle: precomputed from SWW if available, else max over time."""
        if self._max_stage_c is not None:
            return self._max_stage_c
        import numpy as np
        return np.max(self.stage, axis=0)

    #------------------------------------------
    # Elevation-change (delta) procedures
    #------------------------------------------

    @property
    def elev_delta(self):
        """Bed elevation change from t=0: shape (n_time, n_tri), or None for static."""
        if self.elev.ndim != 2:
            return None
        if getattr(self, '_elev_delta', None) is None:
            from anuga.utilities._sww_frames import FrameArray
            first = self.elev[0, :]
            self._elev_delta = FrameArray(
                self._source, lambda get: get('elevation_c') - first,
                *self.elev.shape)
        return self._elev_delta

    def _elev_delta_frame(self, frame, figsize, dpi, vmin, vmax, cmap='RdBu_r',
                          basemap=False, alpha=1.0,
                          basemap_provider=BASEMAP_DEFAULT, xlim=None, ylim=None,
                          smooth=False, show_mesh=False):
        """Render one elevation-change frame: elev[t] - elev[0]."""
        import matplotlib.pyplot as plt

        if self.elev.ndim == 2:
            delta_data = self.elev[frame, :] - self.elev[0, :]
            title = f'Elevation change: Time {self.time[frame]:.2f}'
        else:
            import numpy as np
            delta_data = np.zeros(self.elev.shape)
            title = 'Elevation change (static)'

        cache_key = ('elev_delta', figsize, dpi, basemap, basemap_provider,
                     cmap, vmin, vmax, alpha,
                     tuple(xlim) if xlim is not None else None,
                     tuple(ylim) if ylim is not None else None,
                     show_mesh)

        if cache_key not in self._figure_cache:
            self._figure_cache[cache_key] = plt.figure(figsize=figsize, dpi=dpi)

        fig = self._figure_cache[cache_key]
        fig.clf()
        ax = fig.add_subplot(111)
        ax.set_aspect('equal')
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        ax.set_title(title)

        triang = self.triang_abs if (basemap and self.epsg) else self.triang
        triang.set_mask(None)
        if smooth:
            import matplotlib.tri as _mtri
            triang_all = _mtri.Triangulation(triang.x, triang.y, triang.triangles)
            v_delta = _face_to_vertex(triang, delta_data)
            im = ax.tripcolor(triang_all, v_delta, shading='gouraud',
                              cmap=cmap, alpha=alpha, vmin=vmin, vmax=vmax)
            im.set_linewidths(0)
        else:
            im = ax.tripcolor(triang, facecolors=delta_data,
                              cmap=cmap, alpha=alpha, vmin=vmin, vmax=vmax)

        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)

        if show_mesh:
            ax.triplot(triang, color='black', linewidth=0.25, alpha=0.45, zorder=4)

        fig.colorbar(im, ax=ax, label='Elevation change (m)')

        if basemap and self.epsg:
            _add_basemap(ax, self.epsg, basemap_provider, cache=self._basemap_cache)

        return fig, ax

    def save_elev_delta_frame(self, frame=-1, figsize=(10, 6), dpi=160,
                              vmin=-5.0, vmax=5.0, cmap='RdBu_r', basemap=False,
                              alpha=1.0, basemap_provider=BASEMAP_DEFAULT,
                              xlim=None, ylim=None, smooth=False,
                              show_elev=False, elev_levels=10, show_mesh=False):
        """Save one elevation-change frame (elev[t] - elev[0]) to disk."""
        name = self.name
        plot_dir = self.plot_dir
        frame_num = self._elev_delta_frame_count

        render_frame = 0 if self.elev.ndim == 1 else frame
        fig, ax = self._elev_delta_frame(render_frame, figsize, dpi, vmin, vmax, cmap,
                                         basemap, alpha, basemap_provider,
                                         xlim=xlim, ylim=ylim, smooth=smooth,
                                         show_mesh=show_mesh)

        fname = name + '_elev_delta_{0:0>10}.png'.format(frame_num)
        if plot_dir is None:
            fig.savefig(fname)
        else:
            fig.savefig(os.path.join(plot_dir, fname))
        plt.close(fig)
        self._elev_delta_frame_count += 1

    #------------------------------------------
    # Maximum-over-time procedures
    #------------------------------------------

    def _render_max_qty(self, ax, triang, max_depth, qty_data, elev, md,
                        basemap, cmap, alpha, vmin, vmax, smooth,
                        show_elev=False, elev_levels=10, show_mesh=False,
                        xlim=None, ylim=None):
        """Shared renderer for save_max_*_frame; returns the ScalarMappable."""
        import numpy as np
        ax.set_xlim(triang.x.min(), triang.x.max())
        ax.set_ylim(triang.y.min(), triang.y.max())

        wet_mask = max_depth <= md   # True = never-wet → exclude from quantity

        if basemap and self.epsg:
            # Always use imshow for basemap regardless of smooth flag.
            # Flat tripcolor leaves blocky triangle artefacts over the map tiles;
            # imshow + LinearTriInterpolator gives a smooth anti-aliased boundary.
            # Resolution targets the visible (zoomed) region, matching _animated_frame.
            import matplotlib.tri as _mtri
            from matplotlib.tri import LinearTriInterpolator as _LTI

            triang_wet = _mtri.Triangulation(
                triang.x, triang.y, triang.triangles, mask=wet_mask)
            v_qty = _face_to_vertex(triang, qty_data, face_mask=wet_mask)

            mx0 = float(triang.x.min()); mx1 = float(triang.x.max())
            my0 = float(triang.y.min()); my1 = float(triang.y.max())
            x0 = max(mx0, xlim[0]) if xlim is not None else mx0
            x1 = min(mx1, xlim[1]) if xlim is not None else mx1
            y0 = max(my0, ylim[0]) if ylim is not None else my0
            y1 = min(my1, ylim[1]) if ylim is not None else my1
            if x0 >= x1:
                x0, x1 = mx0, mx1
            if y0 >= y1:
                y0, y1 = my0, my1

            _span = max(x1 - x0, y1 - y0)
            nx = max(200, int(600 * (x1 - x0) / _span))
            ny = max(200, int(600 * (y1 - y0) / _span))
            xi = np.linspace(x0, x1, nx)
            yi = np.linspace(y0, y1, ny)
            Xi, Yi = np.meshgrid(xi, yi)
            zi = _LTI(triang_wet, v_qty)(Xi, Yi)
            im = ax.imshow(zi, extent=[x0, x1, y0, y1], origin='lower',
                           cmap=cmap, vmin=vmin, vmax=vmax, alpha=alpha,
                           interpolation='bilinear', zorder=1)

        elif smooth:
            import matplotlib.tri as _mtri
            dry_mask = max_depth > md
            triang_dry = _mtri.Triangulation(
                triang.x, triang.y, triang.triangles, mask=dry_mask)
            v_elev = _face_to_vertex(triang, elev, face_mask=dry_mask)
            dry_elev = np.asarray(elev)[~np.asarray(dry_mask, dtype=bool)]
            e_vmin = float(dry_elev.min()) if dry_elev.size else 0.0
            e_vmax = float(dry_elev.max()) if dry_elev.size else 1.0
            if e_vmin >= e_vmax:
                e_vmax = e_vmin + 1.0
            ax.tripcolor(triang_dry, v_elev, shading='gouraud', cmap='Greys_r',
                         vmin=e_vmin, vmax=e_vmax).set_linewidths(0)
            triang_wet = _mtri.Triangulation(
                triang.x, triang.y, triang.triangles, mask=wet_mask)
            v_qty = _face_to_vertex(triang, qty_data, face_mask=wet_mask)
            im = ax.tripcolor(triang_wet, v_qty, shading='gouraud',
                              cmap=cmap, alpha=alpha, vmin=vmin, vmax=vmax)
            im.set_linewidths(0)

        else:
            triang.set_mask(max_depth > md)
            ax.tripcolor(triang, facecolors=elev, cmap='Greys_r')
            triang.set_mask(max_depth <= md)
            im = ax.tripcolor(triang, facecolors=qty_data, cmap=cmap,
                              alpha=alpha, vmin=vmin, vmax=vmax)
            triang.set_mask(None)
        if show_elev:
            _draw_elev_contours(ax, triang, elev, elev_levels)
        if show_mesh:
            ax.triplot(triang, color='black', linewidth=0.25, alpha=0.45,
                       zorder=4)

        # Use offset notation on both axes (e.g. "1e6" header + short ticks)
        # and limit tick count so large UTM values don't collide.
        from matplotlib.ticker import ScalarFormatter, MaxNLocator
        for _axis in (ax.xaxis, ax.yaxis):
            _fmt = ScalarFormatter(useOffset=True, useMathText=False)
            _fmt.set_powerlimits((-3, 4))
            _axis.set_major_formatter(_fmt)
            _axis.set_major_locator(MaxNLocator(nbins=5, integer=True))

        return im

    def save_max_depth_frame(self, frame=None, figsize=(10, 6), dpi=160,
                             vmin=0.0, vmax=20.0, cmap='viridis', basemap=False,
                             alpha=1.0, basemap_provider=BASEMAP_DEFAULT,
                             xlim=None, ylim=None, smooth=False,
                             show_elev=False, elev_levels=10, show_mesh=False):
        """Save a single frame showing the maximum depth at each triangle."""
        import matplotlib.pyplot as plt

        max_depth = self.max_depth
        md = self.min_depth
        try:
            elev = self.elev[0, :]
        except (IndexError, TypeError):
            elev = self.elev

        triang = self.triang_abs if (basemap and self.epsg) else self.triang
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        plt.title('Max Depth')
        im = self._render_max_qty(ax, triang, max_depth, max_depth, elev, md,
                                   basemap, cmap, alpha, vmin, vmax, smooth,
                                   show_elev=show_elev, elev_levels=elev_levels,
                                   show_mesh=show_mesh, xlim=xlim, ylim=ylim)
        triang.set_mask(None)
        ax.set_aspect('equal')
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        fig.colorbar(im, ax=ax)
        if basemap and self.epsg:
            _add_basemap(ax, self.epsg, basemap_provider, cache=self._basemap_cache)
        fname = '{}_max_depth_{:0>10}.png'.format(self.name, self._max_depth_frame_count)
        path = fname if self.plot_dir is None else os.path.join(self.plot_dir, fname)
        fig.savefig(path)
        self._max_depth_frame_count += 1
        plt.close()
        fig.clf()

    def save_max_stage_frame(self, frame=None, figsize=(10, 6), dpi=160,
                             vmin=-20.0, vmax=20.0, cmap='viridis', basemap=False,
                             alpha=1.0, basemap_provider=BASEMAP_DEFAULT,
                             xlim=None, ylim=None, smooth=False,
                             show_elev=False, elev_levels=10, show_mesh=False):
        """Save a single frame showing the maximum stage at each triangle."""
        import matplotlib.pyplot as plt

        max_depth = self.max_depth
        max_stage = self.max_stage
        md = self.min_depth
        try:
            elev = self.elev[0, :]
        except (IndexError, TypeError):
            elev = self.elev

        triang = self.triang_abs if (basemap and self.epsg) else self.triang
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        plt.title('Max Stage')
        im = self._render_max_qty(ax, triang, max_depth, max_stage, elev, md,
                                   basemap, cmap, alpha, vmin, vmax, smooth,
                                   show_elev=show_elev, elev_levels=elev_levels,
                                   show_mesh=show_mesh, xlim=xlim, ylim=ylim)
        triang.set_mask(None)
        ax.set_aspect('equal')
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        fig.colorbar(im, ax=ax)
        if basemap and self.epsg:
            _add_basemap(ax, self.epsg, basemap_provider, cache=self._basemap_cache)
        fname = '{}_max_stage_{:0>10}.png'.format(self.name, self._max_stage_frame_count)
        path = fname if self.plot_dir is None else os.path.join(self.plot_dir, fname)
        fig.savefig(path)
        self._max_stage_frame_count += 1
        plt.close()
        fig.clf()

    def save_max_speed_frame(self, frame=None, figsize=(10, 6), dpi=160,
                             vmin=0.0, vmax=10.0, cmap='viridis', basemap=False,
                             alpha=1.0, basemap_provider=BASEMAP_DEFAULT,
                             xlim=None, ylim=None, smooth=False,
                             show_elev=False, elev_levels=10, show_mesh=False):
        """Save a single frame showing the maximum speed at each triangle."""
        import matplotlib.pyplot as plt

        max_depth = self.max_depth
        max_speed = self.max_speed
        md = self.min_depth
        try:
            elev = self.elev[0, :]
        except (IndexError, TypeError):
            elev = self.elev

        triang = self.triang_abs if (basemap and self.epsg) else self.triang
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        plt.title('Max Speed')
        im = self._render_max_qty(ax, triang, max_depth, max_speed, elev, md,
                                   basemap, cmap, alpha, vmin, vmax, smooth,
                                   show_elev=show_elev, elev_levels=elev_levels,
                                   show_mesh=show_mesh, xlim=xlim, ylim=ylim)
        triang.set_mask(None)
        ax.set_aspect('equal')
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        fig.colorbar(im, ax=ax)
        if basemap and self.epsg:
            _add_basemap(ax, self.epsg, basemap_provider, cache=self._basemap_cache)
        fname = '{}_max_speed_{:0>10}.png'.format(self.name, self._max_speed_frame_count)
        path = fname if self.plot_dir is None else os.path.join(self.plot_dir, fname)
        fig.savefig(path)
        self._max_speed_frame_count += 1
        plt.close()
        fig.clf()

    def save_max_speed_depth_frame(self, frame=None, figsize=(10, 6), dpi=160,
                                   vmin=0.0, vmax=20.0, cmap='viridis', basemap=False,
                                   alpha=1.0, basemap_provider=BASEMAP_DEFAULT,
                                   xlim=None, ylim=None, smooth=False,
                                   show_elev=False, elev_levels=10,
                                   show_mesh=False):
        """Save a single frame showing the maximum speed×depth at each triangle."""
        import matplotlib.pyplot as plt

        max_depth = self.max_depth
        max_speed_depth = self.max_speed_depth
        md = self.min_depth
        try:
            elev = self.elev[0, :]
        except (IndexError, TypeError):
            elev = self.elev

        triang = self.triang_abs if (basemap and self.epsg) else self.triang
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        plt.title('Max Speed x Depth')
        im = self._render_max_qty(ax, triang, max_depth, max_speed_depth, elev, md,
                                   basemap, cmap, alpha, vmin, vmax, smooth,
                                   show_elev=show_elev, elev_levels=elev_levels,
                                   show_mesh=show_mesh, xlim=xlim, ylim=ylim)
        triang.set_mask(None)
        ax.set_aspect('equal')
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        fig.colorbar(im, ax=ax)
        if basemap and self.epsg:
            _add_basemap(ax, self.epsg, basemap_provider, cache=self._basemap_cache)
        fname = '{}_max_speed_depth_{:0>10}.png'.format(
            self.name, self._max_speed_depth_frame_count)
        path = fname if self.plot_dir is None else os.path.join(self.plot_dir, fname)
        fig.savefig(path)
        self._max_speed_depth_frame_count += 1
        plt.close()
        fig.clf()

    #------------------------------------------
    # Animation procedures
    #------------------------------------------
    def make_depth_animation(self):

        return self._make_quantity_animation(quantity='depth')

    def make_speed_animation(self):

        return self._make_quantity_animation(quantity='speed')

    def make_stage_animation(self):

        return self._make_quantity_animation(quantity='stage')

    def make_speed_depth_animation(self):

        return self._make_quantity_animation(quantity='speed_depth')


    def _make_quantity_animation(self, quantity='depth'):

        import numpy as np
        import glob
        from matplotlib import image, animation
        from matplotlib import pyplot as plt

        plot_dir = self.plot_dir
        name = self.name

        if plot_dir is None:
            expression = name+'_'+quantity+'_*.png'
        else:
            expression = os.path.join(plot_dir, name+'_'+quantity+'_*.png')
        img_files = sorted(glob.glob(expression))

        figsize = (10, 6)

        fig = plt.figure(figsize=figsize, dpi=80)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.axis('off')  # so there's not a second set of axes
        im = plt.imshow(image.imread(img_files[0]))

        def init():
            im.set_data(image.imread(img_files[0]))
            return im,

        def animate(i):
            image_i = image.imread(img_files[i])
            im.set_data(image_i)
            return im,

        anim = animation.FuncAnimation(fig, animate, init_func=init,
                                       frames=len(img_files),
                                       interval=200, blit=True)

        plt.close()

        return anim

    def make_plot_dir(self, clobber=True):
        """
        Utility function to create a directory for storing a sequence of plot
        files, or if the directory already exists, clear out any old plots.
        If clobber==False then it will abort instead of deleting existing files.
        """

        plot_dir = self.plot_dir
        if plot_dir is None:
            return
        else:
            if os.path.isdir(plot_dir):
                if clobber:
                    try:
                        os.remove("%s/*" % plot_dir)
                    except OSError:
                        pass
                else:
                    raise OSError(
                      '*** Cannot clobber existing directory %s' % plot_dir)
            else:
                os.mkdir("%s" % plot_dir)
                print("Figure files for each frame will be stored in " + plot_dir)

    def set_epsg(self, epsg):
        """Set (or override) the EPSG code and rebuild the absolute-coordinate
        triangulation used for basemap overlays.

        Useful when the SWW file pre-dates EPSG storage (older files lack the
        attribute).  Call this before generating frames with ``basemap=True``.

        Parameters
        ----------
        epsg : int or None
            EPSG code of the coordinate system, e.g. 32756 for UTM zone 56 S.
            Pass ``None`` to clear the code and disable basemap support.
        """
        import matplotlib.tri as tri
        self.epsg = int(epsg) if epsg is not None else None
        self._basemap_cache = {}  # tiles fetched under old EPSG are invalid
        if self.epsg is not None:
            self.triang_abs = tri.Triangulation(
                self.x + self.xllcorner,
                self.y + self.yllcorner,
                self.triangles)
        else:
            self.triang_abs = self.triang

    def triplot(self, figsize = (10, 6), dpi=80, **kwargs):
        """
        Create a triplot of the mesh.

        Args:
            figsize: Figure size
            dpi: Dots per inch
            **kwargs: Additional arguments to pass to triplot
        Returns:
            fig: Matplotlib figure object
            ax: Matplotlib axes object
            lines: The lines created by triplot
        """

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        lines = ax.triplot(self.triang, **kwargs)
        return fig, ax, lines


    def tripcolor(self, figsize = (10, 6), dpi=80, **kwargs):
        """
        Create a tripcolor plot of the mesh.

        Args:
            figsize: Figure size
            dpi: Dots per inch
            **kwargs: Additional arguments to pass to tripcolor
        Returns:
            fig: Matplotlib figure object
            ax: Matplotlib axes object
            im: The image created by tripcolor
        """

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

        im = ax.tripcolor(self.triang, **kwargs)
        return fig, ax, im

    def get_flow_through_cross_section(self, polyline: list, verbose: bool = False) -> tuple[np.ndarray, list]:
        """
        Calculate flow through a cross-section defined by a polyline.

        Args:
            polyline: List of [x, y] coordinates defining the cross-section
            verbose: Whether to print debug information

        Returns:
            tuple containing:
                - time: numpy array of time values
                - Q: list of flow values at each timestep
        """

        if not hasattr(self, 'mesh'):
            from anuga.file.sww import get_mesh_and_quantities_from_file
            self.mesh, _, __ = get_mesh_and_quantities_from_file(self.filename, verbose=verbose)

        segments = self.mesh.get_intersecting_segments(polyline, verbose=verbose)

        if verbose:
            print(f'Cross-section intersects {len(segments)} triangles')
            print('Computing hydrograph')

        # Pre-extract per-segment geometry as arrays for vectorised computation
        tri_ids = np.array([seg.triangle_id for seg in segments])  # (S,)
        normals = np.array([seg.normal for seg in segments])        # (S, 2)
        lengths = np.array([seg.length for seg in segments])        # (S,)

        # Slice momentum arrays for the intersected triangles: (T, S)
        uh = self.xmom[:, tri_ids]
        vh = self.ymom[:, tri_ids]

        # Normal momentum at each segment and timestep, summed to hydrograph
        normal_mom = uh * normals[:, 0] + vh * normals[:, 1]  # (T, S)
        Q = (normal_mom * lengths).sum(axis=1)                 # (T,)

        if verbose:
            print(f'Done: {len(self.time)} timesteps, Q range [{Q.min():.3g}, {Q.max():.3g}] m^3/s')

        return self.time, Q

    def get_triangles_inside_polygon(self, polygon: list | np.ndarray, verbose: bool = False) -> list | np.ndarray:
        """
        Get list of triangle IDs whose centroids lie within a given polygon.

        Args:
            polygon: List of [x, y] coordinates defining the polygon
            verbose: Whether to print debug information

        Returns:
            List of triangle IDs inside the polygon
        """

        try:
            _ = self.mesh
        except AttributeError:
            from anuga.file.sww import get_mesh_and_quantities_from_file
            self.mesh, _, __ = get_mesh_and_quantities_from_file(self.filename, verbose=verbose)

        triangle_ids = self.mesh.get_triangles_inside_polygon(polygon, verbose=verbose)

        return triangle_ids

    def water_volume(self, per_unit_area=False, triangle_ids=None, polygon=None, verbose=False) -> np.ndarray:
        """
        Compute the water volume associated within a subset of triangles or within a polygon.

        Args:
            per_unit_area: If True, return volume per unit area
            triangle_ids: List of triangle IDs to include
            polygon: Polygon defining area of interest
            verbose: Whether to print debug information

        Returns:
            Numpy array of water volume at each timestep
        """

        if not hasattr(self, "mesh"):
            from anuga.file.sww import get_mesh_and_quantities_from_file
            self.mesh, _, __ = get_mesh_and_quantities_from_file(self.filename, verbose=verbose)

        self.areas = self.mesh.areas



        if(triangle_ids is None and polygon is None):
            triangle_ids=list(range(len(self.xc)))
        elif(triangle_ids is not None):
            pass
        else:
            triangle_ids = self.get_triangles_inside_polygon(polygon, verbose=verbose)


        l=len(self.time)
        areas=self.areas[triangle_ids]

        total_area=areas.sum()
        volume=self.time*0.

        if self.elev.ndim ==1:
            for i in range(l):
                volume[i]=((self.stage[i,triangle_ids]-self.elev[triangle_ids])*(self.stage[i,triangle_ids]>self.elev[triangle_ids])*areas).sum()
        else:
            for i in range(l):
                volume[i]=((self.stage[i,triangle_ids]-self.elev[i,triangle_ids])*(self.stage[i,triangle_ids]>self.elev[i,triangle_ids])*areas).sum()

        if(per_unit_area):
            volume = volume / total_area

        return volume



