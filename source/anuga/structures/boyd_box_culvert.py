#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$30/08/2010 10:56:39 AM$"

import culvert
import boyd_box_routine

class Boyd_box_culvert(culvert.Culvert):

    def __init__(self,
                 domain,
                 end_points,
                 width=None,
                 height=None,
                 verbose=False):

        culvert.Culvert.__init__(self, domain, end_points, width, height, verbose)

        self.routine = boyd_box_routine.Boyd_box_routine(self)

    def __call__(self):

        return self.routine()

