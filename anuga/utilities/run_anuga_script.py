#! /usr/bin/python
__author__="stephen"
__date__ ="$20/08/2012 11:20:00 PM$"



def run_script(script, args=None, np=1, alg=None, verbose=False, debug=False,
               cfl=None, allow_parallel=True):
    """Run a validation script in a child process.

    ``cfl`` is accepted only for backwards compatibility (run_validation_script
    forwards it). It is deliberately NOT passed on to the child: the standard
    parser no longer exposes -cfl, so a child script would reject it.
    """
    import subprocess

    if args is None:
        if alg is None:
            from anuga.validation_utilities.parameters import alg
    else:
        alg = args.alg
        np = args.np
        verbose = args.verbose
        debug = getattr(args, 'debug', False)

    # Build the child command as an argument list (no shell); forward the
    # standard flags. -cfl is not a valid child argument, so cfl is not passed.
    cmd = ['python', script, '-alg', str(alg)]
    if np > 1 and allow_parallel:
        cmd = ['mpiexec', '-np', str(np)] + cmd
    if verbose:
        cmd.append('-v')
    if debug:
        cmd.append('-d')

    if verbose:
        print(50*'=')
        print('Run ' + ' '.join(cmd))
        print(50*'=')

    try:
        return subprocess.call(cmd)
    except Exception:
        return 1







