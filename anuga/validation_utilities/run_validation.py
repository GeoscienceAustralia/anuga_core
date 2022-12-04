#! /usr/bin/python



from builtins import str
__author__="stephen"
__date__ ="$20/08/2012 11:20:00 PM$"

    

def run_validation_script(script, np=1, cfl=None, alg=None, verbose=False):
    """ For backwards compatiblity wrap new run_anuga_script as run_validation_script
    """
    
    from anuga import run_anuga_script

    run_anuga_script(script, np=np, cfl=cfl, alg=alg, verbose=verbose)

    
def run_validation_script_old(script,np=1, cfl=None, alg=None, verbose=False):
    #from anuga.validation_utilities.fabricate import run
    
    if cfl is None:
        from anuga.validation_utilities.parameters import cfl
    if alg is None:
        from anuga.validation_utilities.parameters import alg
        
    import subprocess
    
    if np>1:
        if verbose:
            cmd = 'mpiexec -np %s python %s -cfl %s -alg %s -v ' % (str(np), script, str(cfl), str(alg))
        else:
            cmd = 'mpiexec -np %s python %s -cfl %s -alg %s' % (str(np), script, str(cfl), str(alg))
        print(50*'=')
        print('Run '+cmd)
        print(50*'=')
        subprocess.call([cmd], shell=True)

    else:
        if verbose:
            cmd = 'python %s -cfl %s -alg %s -v ' % (script, str(cfl), str(alg))
        else:
            cmd = 'python %s -cfl %s -alg %s' % (script, str(cfl), str(alg))
        print(50*'=')
        print('Run '+cmd)
        print(50*'=')
        subprocess.call([cmd], shell=True)


    


