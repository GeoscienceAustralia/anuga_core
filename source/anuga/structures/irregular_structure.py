#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$30/08/2010 10:56:39 AM$"

import structure
import boyd_box_routine

class Irregular_structure(structure.Structure):

    def __init__(self,
                 domain,
                 end_points,
                 width,
                 height=None,
                 slice_coordinates,
                 verbose=False):

        structure.Structure.__init__(self, domain, end_points, width, height, verbose)
        
        self.slice_coordinates = slice_coordinates
        
        self.routine = boyd_irregular_routine.Boyd_irregular_routine(self)

    def __call__(self):

        return self.routine()
        
    def get_slice_coordinates(self):
    
        self.slice_coordinates

