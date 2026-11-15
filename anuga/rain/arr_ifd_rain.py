"""
ARR 2016 IFD (Intensity-Frequency-Duration) rainfall grid data.

Classes
-------
Arr_ifd_rain
    Scans a directory of IFD zip files and opens individual grids.
Arr_grd
    One IFD rainfall depth grid (one duration/AEP combination).
"""

import re
import warnings
import zipfile
import glob
import os
import numpy as np


# catchment_depth_<duration in minutes>min_<AEP in percent * 1000>aep.txt.asc
_GRID_NAME = re.compile(r'catchment_depth_(\d+)min_(\d+)aep\.txt\.asc')


def _decode(line):
    """Decode a bytes line (or list of bytes lines) from a zip file."""
    if isinstance(line, list):
        return [_decode(x) for x in line]
    return line.decode('utf-8').strip('\r\n').split(',')


class Arr_ifd_rain:
    """Scan a directory of ARR 2016 IFD zip files and expose grid data.

    Parameters
    ----------
    IFD_DIR : str
        Directory containing the zip files from the Bureau of Meteorology's
        design rainfall system, under any name. Zip files without grids in them
        are ignored.
    Lat : float
        Latitude of the point of interest (degrees, negative = south).
    Lon : float
        Longitude of the point of interest (degrees).
    debug : int, optional
        Verbosity level (0 = quiet).
    """

    def __init__(self, IFD_DIR, Lat, Lon, debug=0):
        self.IFD_DIR = IFD_DIR
        self.Lat     = Lat
        self.Lon     = Lon
        self.debug   = debug

        dur_list  = []
        frq_list  = []
        file_list = []
        zip_list  = []

        for zip_file in sorted(glob.glob(os.path.join(IFD_DIR, '*.zip'))):
            try:
                files_in_zip = zipfile.ZipFile(zip_file).namelist()
            except zipfile.BadZipFile:
                warnings.warn(f'Skipping {zip_file}: not a valid zip file.',
                              stacklevel=2)
                continue
            if debug > 0:
                print(files_in_zip)

            for file_in_zip in files_in_zip:
                match = _GRID_NAME.fullmatch(os.path.basename(file_in_zip))
                if match is None:
                    continue

                if debug > 0:
                    print(file_in_zip)

                dur_list.append(match.group(1))
                frq_list.append(match.group(2))
                file_list.append(file_in_zip)
                zip_list.append(zip_file)

        if not file_list:
            warnings.warn(
                f'No IFD grids (catchment_depth_*.asc in a zip file) found in '
                f'{IFD_DIR}.', stacklevel=2)

        self.dur_list  = np.asarray(dur_list, dtype=int)
        self.frq_list  = np.asarray(frq_list, dtype=int)
        self.file_list = file_list
        self.zip_list  = zip_list

    def open_grd(self, Dur, Frq=None, aep_percent=None):
        """Return an :class:`Arr_grd` for the given duration and AEP.

        Give the probability as either ``aep_percent`` or ``Frq``, not both.

        Parameters
        ----------
        Dur : int
            Duration in minutes.
        Frq : int, optional
            The AEP label in the Bureau of Meteorology's file name, which is
            the AEP in percent times 1000: 1000 is 1% AEP (``01000aep``), and
            100 is 0.1% AEP, or 1 in 1000 (``00100aep``).
        aep_percent : float, optional
            Annual exceedance probability in percent, for example 1.0 for 1%
            AEP or 0.05 for 1 in 2000. The exact labels for 18.127, 39.347 and
            63.212 can also be given as they are rounded in tables, such as
            39.35.

        Returns
        -------
        Arr_grd

        Raises
        ------
        ValueError
            If the probability is not given exactly once, or no grid matches.
        """
        if (Frq is None) == (aep_percent is None):
            raise ValueError('Give exactly one of Frq and aep_percent.')

        if aep_percent is not None:
            # 0.1% of the value separates the rounded labels from the exact
            # ones without reaching a neighbouring AEP.
            labels = [int(f) for f in set(self.frq_list)
                      if abs(f / 1000.0 - aep_percent) <= 1.0e-3 * aep_percent]
            label = labels[0] if labels else None
        else:
            label = int(Frq)

        matches = [(dur_r, frq_r, file_in_zip, zip_name)
                   for dur_r, frq_r, file_in_zip, zip_name
                   in zip(self.dur_list, self.frq_list, self.file_list, self.zip_list)
                   if Dur == int(dur_r) and label == int(frq_r)]
        if not matches:
            asked = f'{aep_percent}% AEP' if aep_percent is not None else f'Frq {Frq}'
            raise ValueError(
                f'No IFD grid in {self.IFD_DIR} for duration {Dur} min and {asked}. '
                f'Durations (min): {sorted(set(int(d) for d in self.dur_list))}. '
                f'AEPs (%): {sorted(set(int(f) / 1000.0 for f in self.frq_list))}.'
            )

        # A grid in more than one zip (1% AEP, say, is in two) is the same grid.
        dur_r, frq_r, file_in_zip, zip_name = matches[-1]
        return Arr_grd(dur_r, frq_r, file_in_zip, zip_name,
                       self.IFD_DIR, debug=self.debug)


class Arr_grd:
    """One ARR 2016 IFD rainfall depth grid (ASCII raster inside a zip).

    Parameters
    ----------
    Dur : str or int
        Duration label (minutes).
    Frq : str or int
        AEP label, the AEP in percent times 1000. See
        :meth:`Arr_ifd_rain.open_grd`.
    file_in_zip : str
        Path of the ``.asc`` file inside the zip archive.
    zip_name : str
        Path to the outer zip archive.
    IFD_DIR : str
        Parent directory (kept for reference).
    debug : int, optional
        Verbosity level.
    """

    def __init__(self, Dur, Frq, file_in_zip, zip_name, IFD_DIR, debug=0):
        self.debug   = debug
        self.dur     = Dur
        self.frq     = Frq
        self.IFD_DIR = IFD_DIR

        zf = zipfile.ZipFile(zip_name)
        with zf.open(file_in_zip) as f:
            lines = f.readlines()

        # Parse ASCII raster header
        cols = rows = None
        xllcorner = yllcorner = cellsize = None
        NODATA_value = None

        for line in lines:
            decoded = _decode(line)
            text = ' '.join(decoded[0].split())
            if text.startswith('ncols'):
                cols = int(text.split()[1])
            elif text.startswith('nrows'):
                rows = int(text.split()[1])
            elif text.startswith('xllcorner'):
                xllcorner = float(text.split()[1])
            elif text.startswith('yllcorner'):
                yllcorner = float(text.split()[1])
            elif text.startswith('cellsize'):
                cellsize = float(text.split()[1])
            elif text.startswith('NODATA_value'):
                NODATA_value = float(text.split()[1])
                xurcorner = xllcorner + cols * cellsize
                yurcorner = yllcorner + rows * cellsize
                if debug > 2:
                    print('xurcorner', xurcorner, 'yurcorner', yurcorner)

        decoded_lines = _decode(lines)
        IFD_data = []
        for l in decoded_lines[6:]:
            row_text = l[0].rstrip(' ')
            IFD_data.append([float(s) for s in row_text.split(' ')])

        self.IFD_Data    = np.asarray(IFD_data, dtype=float)
        if NODATA_value is not None:
            self.IFD_Data[self.IFD_Data == NODATA_value] = np.nan
        self.lons        = xllcorner + np.arange(cols) * cellsize
        lats_r           = yllcorner + np.arange(rows) * cellsize
        self.lats        = lats_r[::-1]  # flip so lats[0] is northernmost
        self.cols        = cols
        self.rows        = rows
        self.xllcorner   = xllcorner
        self.yllcorner   = yllcorner
        self.cellsize    = cellsize
        self.NODATA_value = NODATA_value

    def __repr__(self):
        return (
            f'Arr_grd(\n'
            f'  cols={self.cols}, rows={self.rows}\n'
            f'  xllcorner={self.xllcorner}, yllcorner={self.yllcorner}\n'
            f'  cellsize={self.cellsize}, NODATA={self.NODATA_value}\n'
            f')'
        )

    def _lookup(self, Lons, Lats):
        """Return the value of the cell holding each point.

        Cell *i* spans ``xllcorner + i*cellsize`` to ``xllcorner + (i+1)*cellsize``;
        a point on the outer edge of the grid belongs to the last cell.
        """
        eps = 1.0e-9  # of a cell, for round-off on the edges of cells
        cols = (np.asarray(Lons, dtype=float) - self.xllcorner) / self.cellsize
        rows = (np.asarray(Lats, dtype=float) - self.yllcorner) / self.cellsize

        inside = ((cols >= -eps) & (cols <= self.cols + eps) &
                  (rows >= -eps) & (rows <= self.rows + eps))
        if not inside.all():
            raise ValueError(
                f'{int((~inside).sum())} of {inside.size} points are outside the '
                f'IFD grid (Lon {self.xllcorner} to '
                f'{self.xllcorner + self.cols * self.cellsize}, Lat '
                f'{self.yllcorner} to {self.yllcorner + self.rows * self.cellsize}).'
            )

        cols = np.clip(np.floor(cols + eps).astype(int), 0, self.cols - 1)
        rows = np.clip(np.floor(rows + eps).astype(int), 0, self.rows - 1)
        values = self.IFD_Data[self.rows - 1 - rows, cols]  # row 0 is northernmost

        zero = values == 0.0
        if zero.any():
            warnings.warn(
                f'{int(zero.sum())} of {zero.size} points read a depth of 0 from '
                f'the IFD grid. These grids hold 0, not NODATA, over the ocean, so '
                f'those points get no rain.', stacklevel=3)

        nodata = np.isnan(values)
        if nodata.any():
            warnings.warn(
                f'{int(nodata.sum())} of {nodata.size} points fall on NODATA cells '
                f'of the IFD grid and return NaN.', stacklevel=3)
        return values

    def get_rain_at_point(self, Lon=151.0, Lat=-33.5):
        """Return total event depth (mm) of the grid cell holding the point.

        Parameters
        ----------
        Lon : float
            Longitude (degrees).
        Lat : float
            Latitude (degrees, negative = south).

        Raises
        ------
        ValueError
            If the point is outside the grid.
        """
        return self._lookup(np.array([Lon]), np.array([Lat]))[0]

    def get_rain_at_points(self, Lons, Lats):
        """Return total event depths (mm) of the grid cells holding the points.

        Parameters
        ----------
        Lons : array-like
            Longitudes (degrees).
        Lats : array-like
            Latitudes (degrees, negative = south).

        Returns
        -------
        numpy.ndarray, shape (N,)
            NaN (with a warning) for points on NODATA cells.

        Raises
        ------
        ValueError
            If any point is outside the grid.
        """
        return self._lookup(Lons, Lats)

    def plot(self, Lon=None, Lat=None, SiteLabel='SiteLabel',
             ax=None, close_plot=False):
        """Plot the IFD depth grid with optional site marker."""
        import matplotlib.pyplot as plt

        try:
            import cmaps as nclcmaps
            cmap = nclcmaps.MPL_plasma_r
        except ImportError:
            cmap = 'plasma_r'

        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 6))
        else:
            fig = ax.get_figure()

        IFD_Data = self.IFD_Data
        plot_marker = Lon is not None and Lat is not None
        if plot_marker:
            site_rain = self.get_rain_at_point(Lon=Lon, Lat=Lat)

        mask = IFD_Data > 0.0
        vmax = np.nanmax(IFD_Data)
        vmin = np.min(IFD_Data[mask])
        # lons and lats are the cells' lower-left corners; plot at the centres
        half = 0.5 * self.cellsize
        cs = ax.pcolormesh(self.lons + half, self.lats + half, IFD_Data,
                           cmap=cmap, vmin=vmin, vmax=vmax)
        fig.colorbar(cs)

        for poly_file in glob.glob(os.path.join(self.IFD_DIR, '*.ply')):
            xp, yp = [], []
            with open(poly_file) as f:
                for ln in f.readlines():
                    xp.append(float(ln.split(',')[0]))
                    yp.append(float(ln.split(',')[1]))
            plt.plot(xp, yp, c='white', linestyle='--', linewidth=1.5)

        if plot_marker:
            plt.scatter(Lon, Lat, marker='+', color='white', s=10, zorder=10)

        fig.suptitle(f'2016 IFD Grid  AEP {float(self.frq) / 1000.0:g}%  {self.dur} min',
                     color='red')
        if plot_marker:
            ax.set_title(
                f'{SiteLabel}  Lat: {Lat:.4f}  Lon: {Lon:.4f}  Depth: {site_rain:.1f} mm',
                color='blue')
        else:
            ax.set_title(SiteLabel, color='blue')

        plt.show()
        if close_plot:
            plt.close(fig=fig)

        return ax
