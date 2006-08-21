#!/usr/bin/env python

from distutils.core import setup, Extension

setup (name = "pMesh",
       version = "1.XX",
       author="Duncan Gray",
       author_email="duncan.gray@ga.gov.au",
       ext_modules = [
                      Extension("triang",["triangle/triangle.c",
                                          "triangle/triangmodule.c"],
                                define_macros=[("TRILIBRARY",1),
                                               ("NO_TIMER",1)]
                                )           
                      
                      ]
       )
