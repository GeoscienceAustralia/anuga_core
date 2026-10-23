
from tkinter import *
import string
import time
from os.path import join

class ToolBarButton(Label):
    def __init__(self, top, parent, tag=None, image=None, command=None,
                 statushelp='', balloonhelp='', height=21, width=21,
                 bd=1, activebackground='lightgrey', padx=0, pady=0,
                 state='normal', bg='grey75', home_dir=''):
        # HiDPI: grow the button box and icon to match the app's UI scale
        # (set by Draw.appInit); 1.0 on a normal display, so unchanged there.
        scale = getattr(top, 'ui_scale', 1.0) or 1.0
        Label.__init__(self, parent,
                       height=int(round(height * scale)),
                       width=int(round(width * scale)),
                       relief='flat', bd=bd, bg=bg)


        self.bg = bg
        self.activebackground = activebackground
        if image != None:
            if image.split('.')[1] == 'bmp':
                self.Icon = BitmapImage(file=join(home_dir,'icons/%s' % image))
            else:
                self.Icon = PhotoImage(file=join(home_dir,'icons/%s' % image))
        else:
                self.Icon = PhotoImage(file=join(home_dir,'icons/blank.gif'))
        self.Icon = self._scale_icon(self.Icon, scale)
        self.config(image=self.Icon)
        self.tag = tag
        self.icommand = command
        self.command  = self.activate
        self.bind("<Enter>",           self.buttonEnter)
        self.bind("<Leave>",           self.buttonLeave)
        self.bind("<ButtonPress-1>",   self.buttonDown)
        self.bind("<ButtonRelease-1>", self.buttonUp)
        self.pack(side='left', anchor=NW, padx=padx, pady=pady)
        if balloonhelp or statushelp:
            top.balloon().bind(self, balloonhelp, statushelp)
        self.state = state

    @staticmethod
    def _scale_icon(icon, scale):
        """Upscale a PhotoImage to the UI scale via zoom()/subsample().

        BitmapImage has no zoom(), and a scale of ~1 is a no-op, so both are
        returned unchanged.
        """
        from fractions import Fraction
        if not isinstance(icon, PhotoImage):
            return icon
        fr = Fraction(scale).limit_denominator(4)
        if fr.numerator == fr.denominator:
            return icon
        img = icon.zoom(fr.numerator)
        if fr.denominator > 1:
            img = img.subsample(fr.denominator)
        return img

    def activate(self):
        self.icommand(self.tag)

    def buttonEnter(self, event):
        if self.state != 'disabled':
            self.config(relief='raised', bg=self.bg)

    def buttonLeave(self, event):
        if self.state != 'disabled':
            self.config(relief='flat', bg=self.bg)

    def buttonDown(self, event):
        if self.state != 'disabled':
            self.config(relief='sunken', bg=self.activebackground)

    def buttonUp(self, event):
        if self.state != 'disabled':
            if self.command != None:
                self.command()
            time.sleep(0.05)
            if (self in ToolBarButton.cyclelist):
                if (ToolBarButton.sunkenbutton) and ToolBarButton.sunkenbutton != self:
                    ToolBarButton.sunkenbutton.config(relief='flat', bg=self.bg)
                ToolBarButton.sunkenbutton = self
            else:
                self.config(relief='flat', bg=self.bg)

    def enable(self):
        self.state = 'normal'

    def disable(self):
        self.state = 'disabled'

    def cycle(self, name):
        """Add button to list of buttons where one is left sunken, untill
        the next button is pressed."""
        ToolBarButton.cyclelist.append(self)
    def setInitialSunkenButton(self, name):
        if not ToolBarButton.sunkenbutton:
            ToolBarButton.sunkenbutton = self
            self.config(relief='sunken', bg=self.activebackground)

    #class variable
    cyclelist = []
    sunkenbutton = None
