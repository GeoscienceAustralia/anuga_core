#!/usr/bin/env python3
"""
Capture documentation screenshots for anuga_sww_gui.

Usage:
    python docs/capture_gui_screenshots.py --sww path/to/file.sww

Saves to docs/source/visualisation/img/:
    gui_main.png        — main window after frame generation
    gui_timeseries.png  — timeseries panel open (two picked points)
    gui_mesh.png        — View Mesh window (with Basemap checkbox visible)
"""

import argparse
import os
import sys

OUT_DIR = os.path.join(os.path.dirname(__file__),
                       'source', 'visualisation', 'img')


def grab_window(root):
    """Capture just the main Tk window."""
    from PIL import ImageGrab
    root.update()
    root.update_idletasks()
    x = root.winfo_rootx()
    y = root.winfo_rooty()
    w = root.winfo_width()
    h = root.winfo_height()
    return ImageGrab.grab((x, y, x + w, y + h))


def grab_toplevel(win, root):
    """Capture a Toplevel window."""
    from PIL import ImageGrab
    root.update()
    win.update()
    win.update_idletasks()
    x = win.winfo_rootx()
    y = win.winfo_rooty()
    w = win.winfo_width()
    h = win.winfo_height()
    return ImageGrab.grab((x, y, x + w, y + h))


def run(sww_path):
    import tkinter as tk
    import matplotlib
    matplotlib.use('TkAgg')

    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'scripts'))
    from anuga_sww_gui import SWWAnimationGUI

    root = tk.Tk()
    gui = SWWAnimationGUI(root, initial_sww=sww_path, initial_qty='depth')
    os.makedirs(OUT_DIR, exist_ok=True)

    step = [0]

    def next_step():
        s = step[0]
        step[0] += 1
        print(f'Step {s}')

        if s == 0:
            # Trigger generation — stride=5, DPI=80 for speed
            gui._stride_var.set('5')
            gui._dpi_var.set('80')
            gui._start_generation()
            root.after(1000, next_step)

        elif s == 1:
            # Poll until generation finishes (button re-enabled AND frames loaded)
            gen_running = gui._gen_btn['state'] == 'disabled'
            no_frames   = len(gui._frames) == 0
            if gen_running or no_frames:
                step[0] -= 1      # retry
                root.after(1000, next_step)
            else:
                # Generation done — show frame 0 and let the canvas fully render
                gui._show_frame(0)
                root.update()
                root.update_idletasks()
                root.after(1500, next_step)   # extra settle time

        elif s == 2:
            # Force a canvas redraw and wait again before grabbing
            gui._canvas.draw()
            root.update()
            root.update_idletasks()
            root.after(800, next_step)

        elif s == 3:
            # Screenshot 1: main window with animation frame visible
            img = grab_window(root)
            path = os.path.join(OUT_DIR, 'gui_main.png')
            img.save(path)
            print(f'  Saved {path}')
            root.after(300, next_step)

        elif s == 4:
            # Pick two triangles to demonstrate multi-point timeseries
            sp  = gui._splotter
            n   = len(sp.xc)
            tri1 = n // 3
            tri2 = 2 * n // 3
            gui._ts_triangles = [tri1, tri2]
            gui._compute_plot_transform(gui._last_gen_dpi)
            gui._enter_pick_mode()
            gui._update_timeseries()
            gui._update_pick_overlay()
            gui._canvas.draw()
            root.update()
            root.update_idletasks()
            # Timeseries figure needs extra time to lay out
            root.after(1500, next_step)

        elif s == 5:
            # Force another timeseries redraw then wait
            gui._ts_canvas.draw()
            root.update()
            root.update_idletasks()
            root.after(800, next_step)

        elif s == 6:
            # Screenshot 2: timeseries panel open with two picked points
            img = grab_window(root)
            path = os.path.join(OUT_DIR, 'gui_timeseries.png')
            img.save(path)
            print(f'  Saved {path}')
            gui._exit_pick_mode()
            root.after(300, next_step)

        elif s == 7:
            # Open View Mesh, then find the new Toplevel after it renders
            gui._show_mesh()
            root.update()
            root.update_idletasks()
            mesh_win = None
            for child in reversed(root.winfo_children()):
                if isinstance(child, tk.Toplevel):
                    mesh_win = child
                    break
            root.after(1500, lambda w=mesh_win: _save_mesh(w))

        else:
            root.quit()

    def _save_mesh(win):
        if win is None:
            print('  Warning: mesh Toplevel not found — skipping gui_mesh.png')
            next_step()
            return
        win.update()
        win.update_idletasks()
        img = grab_toplevel(win, root)
        path = os.path.join(OUT_DIR, 'gui_mesh.png')
        img.save(path)
        print(f'  Saved {path}')
        try:
            win.destroy()
        except Exception:
            pass
        root.after(200, next_step)

    # Start after the GUI has fully rendered
    root.after(1500, next_step)
    root.mainloop()
    print('Done.')


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sww', required=True,
                        help='SWW file to use for screenshots')
    args = parser.parse_args()
    if not os.path.isfile(args.sww):
        sys.exit(f'SWW file not found: {args.sww}')
    run(args.sww)


if __name__ == '__main__':
    main()
