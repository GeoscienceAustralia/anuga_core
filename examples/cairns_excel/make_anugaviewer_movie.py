#!/usr/bin/python
# -*- coding: utf-8 -*-
"""

Script to automate steps in making an anuga-viewer based movie.

It removes repeated images from anuga-viewer, and then runs 'convert' to
make a gif.

This assumes you have made a directory of jpg images with anuga-viewer.

It also assumes the 'convert' program is installed [part of ImageMajik].
However, you can optionally use the jpg output from this program with some
other movie-making software (suggest you use clean_up=False in this case).
The 'gif' files output from convert are not small.

#######################
To produce the images:

1) Start anuga-viewer:

> anuga_viewer mysww.sww

and get the view oriented how you would like for the movie.

NOTE: To cull very shallow depths, you can do things like

> anuga_viewer -hmin 0.00001 mysww.sww

You need to experiment to get a good value of hmin [it is not the physical
depth]

2) Press 1 to start recording, space to start the animation, then press 1 to
finish recording, and 3 to save the movie to 'movie.swm'

3) Close anuga-viewer

4) Make a directory with images of the movie, with:

> anuga-viewer -movie MYMOVIEDIR movie.swm

Let it run until finished.

5) Look at the images. You might want to delete the first one (seems to have a
blue background irrespective of the actual value of trace set in anuga-viewer).
Don't worry about repeats -- this script will remove them.

Finally you can run the current script, setting the input parameters as
appropriate

> python makeAnugaViewerMovie.py


Gareth Davies, Geoscience Australia 2014+
"""
import os
import glob
import numpy
import shutil

# INPUT parameters

movie_dir = 'MYMOVIEDIR'  # Directory where jpg images to make movie are stored
delay = 20  # Delay between movie frames, in (1/100)ths of a second
clean_up = False # Delete the temporary jpg files? If False, keep them so they can be used with another program

# You probably don't want to change the following parameters

file_wildcard = 'frame_*.jpg'
animation_name = 'animation.gif'  # Must end in .gif

# END INPUT parameters

#################################################################

os.chdir(movie_dir)

# File names

alljpg = glob.glob(file_wildcard)

# File name integer

m = len(alljpg[0].split('_'))
alljpg_int = [i.split('_')[m - 1] for i in alljpg]
alljpg_int = [int(i.split('.')[0]) for i in alljpg_int]

# File sizes

alljpg_size = [os.stat(i)[6] for i in alljpg]

# Sort file names and sizes

new_inds = numpy.array(alljpg_int).argsort()
alljpg = [alljpg[new_inds[i]] for i in range(len(new_inds))]
alljpg_size = [alljpg_size[new_inds[i]] for i in range(len(new_inds))]

# Make a temporary directory

tmp_dir = 'TMP/'
os.mkdir(tmp_dir)

# Copy all files with unique size to the temporary directory
# To ensure the order is interpreted correctly, we number the file names

for i in range(len(alljpg)):
    if i == 0:
        shutil.copy2(alljpg[i], tmp_dir + str(10000 + i) + '_'
                     + os.path.basename(alljpg[i]))
    else:
        if alljpg_size[i] != alljpg_size[i - 1]:
            shutil.copy2(alljpg[i], tmp_dir + str(10000 + i) + '_'
                         + os.path.basename(alljpg[i]))

# Run movie command
# NOTE: There are other programs which can make more efficient movies
movie_command = 'convert -delay ' + str(delay) + ' TMP/*' \
    + file_wildcard + ' ' + animation_name
print movie_command
os.system(movie_command)

# Clean up directory
if clean_up:
    shutil.rmtree(tmp_dir)
