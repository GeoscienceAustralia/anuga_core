#!/usr/bin/env python3
"""GUI viewer for ANUGA frame animations.

Loads PNG frames saved by Domain_plotter / SWW_plotter (save_depth_frame,
save_stage_frame, save_speed_frame) and plays them as an animation.

Usage::

    anuga_animate_gui
    anuga_animate_gui --dir _plot
    anuga_animate_gui --dir _plot --qty depth
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
# Helpers                                                             #
# ------------------------------------------------------------------ #

_QUANTITIES = ('depth', 'stage', 'speed')


def _find_frames(plot_dir, quantity, prefix):
    """Return sorted list of PNG paths matching prefix_quantity_*.png."""
    pattern = os.path.join(plot_dir, f'{prefix}_{quantity}_*.png')
    return sorted(glob.glob(pattern))


def _detect_prefixes(plot_dir):
    """Scan *plot_dir* and return sorted list of unique domain name prefixes."""
    prefixes = set()
    for qty in _QUANTITIES:
        for path in glob.glob(os.path.join(plot_dir, f'*_{qty}_*.png')):
            fname = os.path.basename(path)
            suffix = f'_{qty}_'
            idx = fname.find(suffix)
            if idx > 0:
                prefixes.add(fname[:idx])
    return sorted(prefixes)


def _detect_quantities(plot_dir, prefix):
    """Return list of quantities that have at least one frame for *prefix*."""
    return [q for q in _QUANTITIES if _find_frames(plot_dir, q, prefix)]


# ------------------------------------------------------------------ #
# GUI                                                                 #
# ------------------------------------------------------------------ #

class AnimationGUI:
    """Tkinter + matplotlib GUI for browsing and playing ANUGA frame PNGs."""

    def __init__(self, root, initial_dir=None, initial_qty=None):
        self.root = root
        self.root.title('ANUGA Frame Animator')
        self.root.minsize(800, 560)
        self.root.protocol('WM_DELETE_WINDOW', self._on_close)

        self._frames = []
        self._current = 0
        self._playing = False
        self._after_id = None
        self._im = None

        self._build_ui()

        start_dir = initial_dir or '_plot'
        self._dir_var.set(start_dir)
        if initial_qty and initial_qty in _QUANTITIES:
            self._qty_var.set(initial_qty)
        self._on_dir_change()

    # -------------------------------------------------------------- #
    # UI construction                                                 #
    # -------------------------------------------------------------- #

    def _build_ui(self):
        # ---------- control panel (top) ----------
        ctrl = ttk.Frame(self.root, padding=6)
        ctrl.pack(side=tk.TOP, fill=tk.X)

        # Row 1: directory
        row1 = ttk.Frame(ctrl)
        row1.pack(fill=tk.X, pady=2)
        ttk.Label(row1, text='Directory:').pack(side=tk.LEFT)
        self._dir_var = tk.StringVar()
        dir_entry = ttk.Entry(row1, textvariable=self._dir_var, width=45)
        dir_entry.pack(side=tk.LEFT, padx=4, fill=tk.X, expand=True)
        dir_entry.bind('<Return>', lambda _e: self._on_dir_change())
        ttk.Button(row1, text='Browse...', command=self._browse_dir).pack(side=tk.LEFT)
        ttk.Button(row1, text='Reload', command=self._on_dir_change).pack(side=tk.LEFT, padx=4)

        # Row 2: prefix + quantity
        row2 = ttk.Frame(ctrl)
        row2.pack(fill=tk.X, pady=2)
        ttk.Label(row2, text='Name prefix:').pack(side=tk.LEFT)
        self._prefix_var = tk.StringVar()
        self._prefix_combo = ttk.Combobox(row2, textvariable=self._prefix_var,
                                          width=22, state='readonly')
        self._prefix_combo.pack(side=tk.LEFT, padx=4)
        self._prefix_combo.bind('<<ComboboxSelected>>', lambda _e: self._on_selection_change())

        ttk.Label(row2, text='   Quantity:').pack(side=tk.LEFT)
        self._qty_var = tk.StringVar(value='depth')
        self._qty_combo = ttk.Combobox(row2, textvariable=self._qty_var,
                                        values=list(_QUANTITIES), width=10,
                                        state='readonly')
        self._qty_combo.pack(side=tk.LEFT, padx=4)
        self._qty_combo.bind('<<ComboboxSelected>>', lambda _e: self._on_selection_change())

        self._info_label = ttk.Label(row2, text='', foreground='grey')
        self._info_label.pack(side=tk.LEFT, padx=12)

        # Row 3: playback controls
        row3 = ttk.Frame(ctrl)
        row3.pack(fill=tk.X, pady=4)

        self._play_btn = ttk.Button(row3, text='Play', width=6,
                                     command=self._toggle_play)
        self._play_btn.pack(side=tk.LEFT, padx=2)

        ttk.Button(row3, text='|<',  width=3, command=self._go_first).pack(side=tk.LEFT, padx=1)
        ttk.Button(row3, text='<',   width=3, command=self._step_back).pack(side=tk.LEFT, padx=1)
        ttk.Button(row3, text='>',   width=3, command=self._step_fwd).pack(side=tk.LEFT, padx=1)
        ttk.Button(row3, text='>|',  width=3, command=self._go_last).pack(side=tk.LEFT, padx=1)

        ttk.Label(row3, text='   FPS:').pack(side=tk.LEFT)
        self._fps_var = tk.DoubleVar(value=5.0)
        ttk.Spinbox(row3, from_=0.5, to=30.0, increment=0.5,
                    textvariable=self._fps_var, width=6).pack(side=tk.LEFT, padx=4)

        self._frame_label = ttk.Label(row3, text='-')
        self._frame_label.pack(side=tk.RIGHT, padx=10)

        # Row 4: frame slider
        row4 = ttk.Frame(ctrl)
        row4.pack(fill=tk.X, pady=2)
        ttk.Label(row4, text='Frame:').pack(side=tk.LEFT)
        self._slider_var = tk.IntVar(value=0)
        self._slider = ttk.Scale(row4, from_=0, to=0,
                                  variable=self._slider_var,
                                  orient=tk.HORIZONTAL,
                                  command=self._on_slider)
        self._slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=4)

        # ---------- matplotlib canvas ----------
        self._fig, self._ax = plt.subplots(figsize=(10, 6))
        self._ax.axis('off')
        self._fig.tight_layout(pad=0)

        canvas_frame = ttk.Frame(self.root)
        canvas_frame.pack(fill=tk.BOTH, expand=True)
        self._canvas = FigureCanvasTkAgg(self._fig, master=canvas_frame)
        self._canvas.draw()
        self._canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # ---------- status bar ----------
        self._status_var = tk.StringVar(value='Select a directory and load frames.')
        ttk.Label(self.root, textvariable=self._status_var,
                  relief=tk.SUNKEN, anchor=tk.W,
                  padding=(4, 1)).pack(side=tk.BOTTOM, fill=tk.X)

    # -------------------------------------------------------------- #
    # Directory / selection logic                                     #
    # -------------------------------------------------------------- #

    def _browse_dir(self):
        d = filedialog.askdirectory(
            title='Select plot directory',
            initialdir=self._dir_var.get() or '.')
        if d:
            self._dir_var.set(d)
            self._on_dir_change()

    def _on_dir_change(self):
        plot_dir = self._dir_var.get().strip()
        if not os.path.isdir(plot_dir):
            self._set_status(f'Directory not found: {plot_dir}')
            return

        prefixes = _detect_prefixes(plot_dir)
        if not prefixes:
            self._set_status(f'No ANUGA PNG frames found in {plot_dir}')
            self._prefix_combo['values'] = []
            return

        self._prefix_combo['values'] = prefixes
        if self._prefix_var.get() not in prefixes:
            self._prefix_var.set(prefixes[0])

        self._refresh_qty_list()
        self._load_frames()

    def _on_selection_change(self):
        self._refresh_qty_list()
        self._load_frames()

    def _refresh_qty_list(self):
        """Update the quantity combo to show only quantities with frames."""
        plot_dir = self._dir_var.get().strip()
        prefix = self._prefix_var.get()
        if not os.path.isdir(plot_dir) or not prefix:
            return
        avail = _detect_quantities(plot_dir, prefix)
        self._qty_combo['values'] = avail if avail else list(_QUANTITIES)
        if self._qty_var.get() not in (avail or list(_QUANTITIES)):
            self._qty_var.set(avail[0] if avail else 'depth')

    # -------------------------------------------------------------- #
    # Frame loading                                                   #
    # -------------------------------------------------------------- #

    def _load_frames(self):
        self._stop_playback()
        plot_dir = self._dir_var.get().strip()
        prefix   = self._prefix_var.get()
        qty      = self._qty_var.get()

        if not os.path.isdir(plot_dir) or not prefix:
            return

        self._frames = _find_frames(plot_dir, qty, prefix)
        n = len(self._frames)

        if n == 0:
            self._set_status(f'No frames found for {prefix}_{qty}_*.png in {plot_dir}')
            self._slider.configure(to=0)
            self._frame_label.config(text='-')
            self._info_label.config(text='')
            return

        self._slider.configure(to=n - 1)
        self._slider_var.set(0)
        self._current = 0
        self._show_frame(0)
        self._info_label.config(text=f'{n} frames')
        self._set_status(f'Loaded {n} frames  |  {os.path.abspath(plot_dir)}')

    # -------------------------------------------------------------- #
    # Frame display                                                   #
    # -------------------------------------------------------------- #

    def _show_frame(self, idx):
        if not self._frames:
            return
        idx = max(0, min(idx, len(self._frames) - 1))
        self._current = idx
        self._slider_var.set(idx)
        self._frame_label.config(text=f'Frame {idx + 1} / {len(self._frames)}')

        img = mpimage.imread(self._frames[idx])
        if self._im is None:
            self._im = self._ax.imshow(img, aspect='equal')
        else:
            self._im.set_data(img)
            self._im.set_extent([0, img.shape[1], img.shape[0], 0])
        self._canvas.draw_idle()

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
        fps = max(0.1, self._fps_var.get())
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
    # Helpers                                                         #
    # -------------------------------------------------------------- #

    def _set_status(self, msg):
        self._status_var.set(msg)

    def _on_close(self):
        self._stop_playback()
        plt.close(self._fig)
        self.root.destroy()


# ------------------------------------------------------------------ #
# Entry point                                                         #
# ------------------------------------------------------------------ #

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='ANUGA PNG frame animation viewer',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--dir', default='_plot',
                        help='Plot directory containing PNG frames')
    parser.add_argument('--qty', default=None,
                        choices=['depth', 'stage', 'speed'],
                        help='Quantity to display on startup')
    args = parser.parse_args()

    root = tk.Tk()
    AnimationGUI(root, initial_dir=args.dir, initial_qty=args.qty)
    root.mainloop()


if __name__ == '__main__':
    main()
