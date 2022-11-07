
from tkinter import *
import string, time
from os.path import join

class ToolBarButton(Label):
    def __init__(self, top, parent, tag=None, image=None, command=None,
                 statushelp='', balloonhelp='', height=21, width=21,
                 bd=1, activebackground='lightgrey', padx=0, pady=0,
                 state='normal', bg='grey75', home_dir=''):
        Label.__init__(self, parent, height=height, width=width,
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
