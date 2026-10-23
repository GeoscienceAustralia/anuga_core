"""
ARR 2016 IFD (Intensity-Frequency-Duration) rainfall grid data.

Classes
-------
Arr_ifd_rain
    Scans a directory of IFD zip files and opens individual grids.
Arr_grd
    One IFD rainfall depth grid (one duration/AEP combination).
"""

import zipfile
import glob
import os
import numpy as np


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
        Directory containing ``*IFD*.zip`` files.
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

        for zip_file in glob.glob(os.path.join(IFD_DIR, '*IFD*.zip')):
            zf = zipfile.ZipFile(zip_file)
            files_in_zip = zf.namelist()
            if debug > 0:
                print(files_in_zip)

            for file_in_zip in files_in_zip:
                if 'epsg' in file_in_zip:
                    continue

                if debug > 0:
                    print(file_in_zip)

                # Filename pattern: catchment_depth_<Dur>min_<Frq>aep.txt.asc
                try:
                    dur_r = file_in_zip.split('min')[0].split('_')[2]
                    frq_r = file_in_zip.split('aep')[0].split('_')[3]
                except IndexError:
                    continue

                if not dur_r.isdigit() or not frq_r.isdigit():
                    continue

                dur_list.append(dur_r)
                frq_list.append(frq_r)
                file_list.append(file_in_zip)
                zip_list.append(zip_file)

        self.dur_list  = np.asarray(dur_list, dtype=int)
        self.frq_list  = np.asarray(frq_list, dtype=int)
        self.file_list = file_list
        self.zip_list  = zip_list

    def open_grd(self, Dur, Frq):
        """Return an :class:`Arr_grd` for the given duration and AEP.

        Parameters
        ----------
        Dur : int
            Duration in minutes.
        Frq : int
            Annual exceedance probability (AEP), e.g. 100 for 1% AEP.

        Returns
        -------
        Arr_grd or None
            ``None`` when no matching file is found (a warning is printed).
        """
        found = False
        arr_grd = None

        for dur_r, frq_r, file_in_zip, zip_name in zip(
                self.dur_list, self.frq_list, self.file_list, self.zip_list):
            if Dur == int(dur_r) and Frq == int(frq_r):
                found   = True
                arr_grd = Arr_grd(dur_r, frq_r, file_in_zip, zip_name,
                                  self.IFD_DIR, debug=self.debug)

        if not found:
            print(
                f'WARNING: No file for duration {Dur} and frequency {Frq}.\n'
                f'Check dur_list and frq_list for available pairs.'
            )

        return arr_grd


class Arr_grd:
    """One ARR 2016 IFD rainfall depth grid (ASCII raster inside a zip).

    Parameters
    ----------
    Dur : str or int
        Duration label (minutes).
    Frq : str or int
        AEP label.
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
                NODATA_value = text.split()[1]
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

    def get_rain_at_point(self, Lon=151.0, Lat=-33.5):
        """Return total event depth (mm) at the nearest grid cell.

        Parameters
        ----------
        Lon : float
            Longitude (degrees).
        Lat : float
            Latitude (degrees, negative = south).
        """
        col = round(float((Lon - self.xllcorner) / self.cellsize))
        row = round(float((Lat - self.yllcorner) / self.cellsize))
        row = int(self.rows) - row
        return self.IFD_Data[row, col]

    def get_rain_at_points(self, Lons, Lats):
        """Return total event depths (mm) at an array of grid points (nearest-neighbour).

        Parameters
        ----------
        Lons : array-like
            Longitudes (degrees).
        Lats : array-like
            Latitudes (degrees, negative = south).

        Returns
        -------
        numpy.ndarray, shape (N,)
        """
        Lons = np.asarray(Lons, dtype=float)
        Lats = np.asarray(Lats, dtype=float)
        cols = np.clip(
            np.round((Lons - self.xllcorner) / self.cellsize).astype(int),
            0, self.cols - 1)
        rows_r = np.clip(
            np.round((Lats - self.yllcorner) / self.cellsize).astype(int),
            0, self.rows - 1)
        rows = self.rows - rows_r
        rows = np.clip(rows, 0, self.rows - 1)
        return self.IFD_Data[rows, cols]

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
        vmax = np.max(IFD_Data)
        vmin = np.min(IFD_Data[mask])
        cs = ax.pcolormesh(self.lons, self.lats, IFD_Data,
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

        fig.suptitle(f'2016 IFD Grid  AEP {self.frq}  {self.dur} min', color='red')
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
