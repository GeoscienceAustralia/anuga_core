from Tkinter import Button, E, W
from threading import Event
from visualiser import Visualiser

class RealtimeVisualiser(Visualiser):
    """A VTK-powered realtime visualiser which runs in its own thread.
    """
    def __init__(self, source):
        """The source parameter is assumed to be a Domain.
        """
        Visualiser.__init__(self, source)

        self.running = True

        # Synchronisation Constructs
        self.sync_idle = Event()
        self.sync_idle.clear()
        self.sync_unpaused = Event()
        self.sync_unpaused.set()
        self.sync_redrawReady = Event()
        self.sync_redrawReady.clear()

    def run(self):
        Visualiser.run(self)
        self.tk_root.after(100, self.sync_idle.set)

    def setup_gui(self):
        Visualiser.setup_gui(self)
        self.tk_pauseResume = Button(self.tk_controlFrame, text="Pause", command=self.pauseResume)
        self.tk_pauseResume.grid(row=1, column=0, sticky=E+W)

    def pauseResume(self):
        if self.sync_unpaused.isSet():
            self.sync_unpaused.clear()
            self.tk_pauseResume.config(text="Resume")
        else:
            self.sync_unpaused.set()
            self.tk_pauseResume.config(text="Pause")

    def shutdown(self):
        Visualiser.shutdown(self)
        self.running = False
        self.sync_idle.set()
        self.sync_unpaused.set()
