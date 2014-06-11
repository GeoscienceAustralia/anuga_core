#! /usr/bin/python


__author__="stephen"
__date__ ="$20/08/2012 11:20:00 PM$"


def run_validation_script(script):
    #from anuga.validation_utilities.fabricate import run
    from anuga.validation_utilities.parameters import alg,cfl
    
    import subprocess
    cmd = 'python %s -cfl %s -alg %s' % (script, str(cfl), str(alg))
    print 50*'='
    print 'Run '+cmd
    print 50*'='
    subprocess.call([cmd], shell=True)
    #run('python', script,  '-alg', alg, '-cfl', cfl)

    


