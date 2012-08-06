"""This module will compile the necessary c-extensions for the mandelbrot demos
"""

import os

s = 'python compile.py mandel_ext.c'
print
os.system(s)

s = 'python compile.py mandelplot_ext.c'
print
os.system(s)
