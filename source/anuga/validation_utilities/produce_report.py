#! /usr/bin/python


__author__="stephen"
__date__ ="$20/08/2012 11:20:00 PM$"

    
from anuga import run_anuga_script
from anuga.validation_utilities import typeset_report

def produce_report(script, args=None):

    import anuga

    if args is None:
        args  = anuga.get_args()
        

    verbose = args.verbose
    
    
    
    #print args
    
    # Get the arguments from the calling script

    run_anuga_script(script, args=args)
    
    # We don't want to run plot_results in parallel
    args.np = 1

    run_anuga_script('plot_results.py', args=args)
    
    typeset_report(verbose=verbose)

    


