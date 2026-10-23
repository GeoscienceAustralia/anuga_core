"""
ARR 2016 hub rainfall data parser.

Classes
-------
Arr_hub_rain
    Reads metadata from an ARR hub export text file.
ARR_point_rainfall_patterns
    Loads temporal rainfall patterns from an ARR 2016 zip archive.
Single_pattern
    One temporal pattern (hyetograph) extracted from ARR pattern data.
"""

import zipfile
import numpy as np


def _decode(line):
    """Decode a bytes line (or list of bytes lines) from a zip file."""
    if isinstance(line, list):
        return [_decode(x) for x in line]
    return line.decode('utf-8').strip('\r\n').split(',')


class Arr_hub_rain:
    """Parse metadata from an ARR 2016 hub export text file.

    Parameters
    ----------
    hub_filename : str
        Path to the ARR hub ``.txt`` export file.
    """

    def __init__(self, hub_filename):
        self.hub_filename = hub_filename
        self._open_hubfile()

    def _open_hubfile(self):
        with open(self.hub_filename) as fh:
            lines = fh.readlines()
        self.lines = lines
        for lcount, line in enumerate(lines):
            if line.startswith('[INPUTDATA]'):
                self.Loc_Lat = lines[lcount+1].split(',')[1].strip()
                self.Loc_Lon = lines[lcount+2].split(',')[1].strip()
            elif line.startswith('[RIVREG]'):
                self.Divis   = lines[lcount+1].split(',')[1].strip()
                self.RivNum  = lines[lcount+2].split(',')[1].strip()
                self.RivName = lines[lcount+3].split(',')[1].strip()
            elif line.startswith('[RIVREG_META]'):
                self.TimeAccessed = lines[lcount+1].split(',')[1].strip()
                self.Version      = lines[lcount+2].split(',')[1].strip()
            elif line.startswith('[LONGARF]'):
                self.ARF_a = lines[lcount+2].split(',')[1].strip()
                self.ARF_b = lines[lcount+3].split(',')[1].strip()
                self.ARF_c = lines[lcount+4].split(',')[1].strip()
                self.ARF_d = lines[lcount+5].split(',')[1].strip()
                self.ARF_e = lines[lcount+6].split(',')[1].strip()
                self.ARF_f = lines[lcount+7].split(',')[1].strip()
                self.ARF_g = lines[lcount+8].split(',')[1].strip()
                self.ARF_h = lines[lcount+9].split(',')[1].strip()
                self.ARF_i = lines[lcount+10].split(',')[1].strip()
            elif line.startswith('[LOSSES]'):
                self.ARR_IL = lines[lcount+2].split(',')[1].strip()
                self.ARR_CL = lines[lcount+3].split(',')[1].strip()
            elif line.startswith('[TP]'):
                self.Tpat_code  = lines[lcount+1].split(',')[1].strip()
                self.Tpatlabel  = lines[lcount+2].split(',')[1].strip()
            elif line.startswith('[ATP]'):
                self.ATpat_code = lines[lcount+1].split(',')[1].strip()
                self.ATpatlabel = lines[lcount+2].split(',')[1].strip()


class ARR_point_rainfall_patterns:
    """Load ARR 2016 temporal rainfall patterns from a zip archive.

    Parameters
    ----------
    pattern_zip_file : str
        Path to the zip file containing ``<Tpat_code>_AllStats.csv`` and
        ``<Tpat_code>_Increments.csv``.
    Tpat_code : str
        Temporal pattern code (from :class:`Arr_hub_rain`).
    """

    def __init__(self, pattern_zip_file, Tpat_code, debug=0):
        try:
            zf = zipfile.ZipFile(pattern_zip_file)
            files_in_zip = zf.namelist()
        except Exception as exc:
            raise OSError(f'Cannot open pattern zip file: {pattern_zip_file}') from exc

        allstats_name   = Tpat_code + '_AllStats.csv'
        increments_name = Tpat_code + '_Increments.csv'

        if allstats_name not in files_in_zip:
            raise OSError(f'{allstats_name} not found in {pattern_zip_file}')
        if increments_name not in files_in_zip:
            raise OSError(f'{increments_name} not found in {pattern_zip_file}')

        with zf.open(allstats_name) as f:
            lines_astat = f.readlines()
        with zf.open(increments_name) as f:
            lines_inc = f.readlines()

        self.linesAStat     = _decode(lines_astat)
        self.linesInc       = _decode(lines_inc)
        self.STATS_Labels   = _decode(lines_astat[0])
        self.INCS_Labels    = _decode(lines_inc[0])


class Single_pattern:
    """One temporal rainfall pattern (hyetograph) from ARR 2016 data.

    Parameters
    ----------
    prp : ARR_point_rainfall_patterns
        Pattern data container.
    index : int
        Row index in the patterns data (1–720).
    Ev_dep : float
        Total event depth (mm) used to scale the pattern increments.
    """

    def __init__(self, prp, index=1, Ev_dep=1.0):
        assert 1 <= index <= 720, f'index = {index} must be in [1, 720]'

        self.index  = index
        self.prp    = prp

        bits_astat = prp.linesAStat[index]
        bits_inc   = prp.linesInc[index]

        self.bitsAStat = bits_astat
        self.bitsInc   = bits_inc

        self.Ev_ID  = bits_inc[0]
        self.Ev_dur = int(bits_inc[1])
        self.Tstep  = int(bits_inc[2])
        self.Zone   = bits_inc[3]
        self.Ev_Frq = bits_inc[4]
        self.Ev_dep = Ev_dep
        self.Tstps  = int(self.Ev_dur / self.Tstep)

        Rplot = [0.0]
        Tplot = [0.0]
        for tcount in range(self.Tstps):
            R = float(bits_inc[5 + tcount]) * Ev_dep / 100.0
            Rplot.append(R)
            Tplot.append(self.Tstep * (tcount + 1))

        self.Tplot = np.asarray(Tplot)  # minutes
        self.Rplot = np.asarray(Rplot)  # mm incremental depth per Tstep interval

    def plot(self, title=None, ax=None):
        """Plot the hyetograph (requires matplotlib)."""
        plot_single_pattern(self.bitsInc, self.Tstep, self.Tstps,
                            Ev_dep=self.Ev_dep, title=title, ax=ax)


def plot_single_pattern(bitsInc, Tstep, Tstps, Ev_dep=100.0, title=None, ax=None):
    """Plot a single rainfall pattern as a bar chart with cumulative line."""
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()

    Rplot = [0.0]
    Tplot = [0.0]
    Ev_ID  = bitsInc[0]
    Ev_Frq = bitsInc[4]

    if title is None:
        title = Ev_ID + ':' + Ev_Frq

    for tcount in range(Tstps):
        R = float(bitsInc[5 + tcount]) * Ev_dep / 100.0
        Rplot.append(R)
        Tplot.append(Tstep * (tcount + 1))

    inv_Rplot = np.array(Rplot)
    cum_pmp   = np.cumsum(inv_Rplot)

    ax.set_title(title, fontsize=10)
    ax.bar(Tplot[1:], Rplot[1:], -Tstep, edgecolor='blue', color='None', align='edge')
    ax2 = ax.twinx()
    ax2.plot(Tplot, cum_pmp, color='red')
    ax2.yaxis.label.set_color('red')
    ax2.set_ylabel('Cumulative rain (mm)', fontsize=8)
    ax2.tick_params(axis='y', colors='red')
    ax.set_xlabel('Time Steps (min)', fontsize=8)
    ax.tick_params(axis='y', colors='green')
    ax.set_ylabel('Incremental rain (mm)', fontsize=8)
    ax.yaxis.label.set_color('green')

    plt.show()
    plt.close(fig=fig)


def plot_frq_patterns(PatternType, Dur_Incs, Tstep, Tstps, Zone, Ev_dur):
    """Plot up to 30 patterns for one duration on a single figure."""
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(12, 8))
    Ev_dep = 100.0

    n = len(Dur_Incs)
    if n == 30:
        v, h = 6, 5
    else:
        v, h = 2, 5

    for i, dur_inc in enumerate(Dur_Incs, start=1):
        ax = fig.add_subplot(v, 5, i)
        plot_single_pattern(dur_inc, Tstep, Tstps, Ev_dep, ax=ax)

    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    stitle = f'{PatternType}, {Zone}  Rain Proportion Patterns for {Ev_dur}-minute Duration'
    fig.suptitle(stitle, fontsize=16, color='red')
    plt.show()


def plot_all_patterns_for_duration(PatternType, STATS_Labels, INCS_Labels,
                                   Dur_Astat, Dur_Incs, debug=0):
    """Plot all 30 patterns (10 frequent, 10 intermediate, 10 rare) for one duration.

    Parameters
    ----------
    PatternType : str
    STATS_Labels, INCS_Labels : list
        Column headers from the patterns CSVs.
    Dur_Astat, Dur_Incs : list
        30-row slices of ``linesAStat`` / ``linesInc`` for one duration.
    debug : int, optional
    """
    Ev_dur = int(Dur_Incs[0][1])
    Tstep  = int(Dur_Incs[0][2])
    Tstps  = int(Ev_dur / Tstep)
    Zone   = Dur_Incs[0][3]

    if debug > 0:
        print(f'E_Dur: {Ev_dur}, TStep: {Tstep}, TSteps: {Tstps}, Zone: {Zone}')

    starts = [0, 10, 20]
    ends   = [10, 20, 30]
    for st, en, frq_label in zip(starts, ends, ['Frequent', 'Intermediate', 'Rare']):
        if debug > 0:
            print(frq_label)
        plot_frq_patterns(PatternType, Dur_Incs[st:en], Tstep, Tstps, Zone, Ev_dur)
