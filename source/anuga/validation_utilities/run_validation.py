#! /usr/bin/python


__author__="stephen"
__date__ ="$20/08/2012 11:20:00 PM$"


def run_validation_script(script):
    from anuga.validation_utilities.fabricate import run
    from anuga.validation_utilities.parameters import alg,cfl
    run('python', script,  '-alg', alg, '-cfl', cfl)

    


