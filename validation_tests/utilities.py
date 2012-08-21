#! /usr/bin/python

# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="stephen"
__date__ ="$20/08/2012 11:20:00 PM$"


def run_validation_script(script):
    from fabricate import run
    from validation_tests.parameters import alg,cfl
    run('python', script,  '-alg', alg, '-cfl', cfl)

    


