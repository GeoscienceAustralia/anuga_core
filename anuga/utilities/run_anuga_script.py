#! /usr/bin/python



from builtins import str
__author__="stephen"
__date__ ="$20/08/2012 11:20:00 PM$"

    
    
def run_script(script, args=None, np=1, alg=None, verbose=False, allow_parallel=True):
    #from anuga.validation_utilities.fabricate import run
    
    
    if args is None:    
        if alg is None:
            from anuga.validation_utilities.parameters import alg
    else:
        alg = args.alg
        np = args.np
        verbose = args.verbose
        
        
    #print args
    args_dict = vars(args)
    #print args_dict
    #print zip(args_dict.keys(), args_dict.values())
    
        
    import subprocess
    import os
    try:
        if np>1 and allow_parallel:
            if verbose:
                cmd = 'mpiexec -np %s python %s -alg %s -v ' % (str(np), script,  str(alg))
            else:
                cmd = 'mpiexec -np %s python %s -alg %s' % (str(np), script, str(alg))
                
            if verbose:
                print(50*'=')
                print('Run '+cmd)
                print(50*'=')


            #os.system(cmd)
            subprocess.call([cmd], shell=True)
    
        else:
            if verbose:
                cmd = 'python %s -alg %s -v ' % (script, str(alg))
            else:
                cmd = 'python %s -alg %s' % (script, str(alg))
            
            if verbose:
                print(50*'=')
                print('Run '+cmd)
                print(50*'=')
            
            #os.system(cmd)
            subprocess.call([cmd], shell=True)
        
        return 0 
    
    except:
        return 1
    
    


    


