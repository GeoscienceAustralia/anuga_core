
# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$10/07/2012 1:18:38 PM$"




def arithmetic_float(expr):
    """argparse ``type`` that evaluates a simple arithmetic expression to a float.

    Lets numeric options be given as expressions, e.g. ``-ft 3600*8`` or
    ``-ys 3600/4``. A plain number (``-ft 3600``) still works. Uses ``ast`` with
    a small operator whitelist — NOT ``eval`` — so it cannot execute arbitrary
    code. Supports ``+ - * / // % **`` and unary +/- on numeric literals.
    """
    import argparse
    import ast
    import operator

    ops = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
           ast.Div: operator.truediv, ast.FloorDiv: operator.floordiv,
           ast.Mod: operator.mod, ast.Pow: operator.pow,
           ast.USub: operator.neg, ast.UAdd: operator.pos}

    def ev(node):
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return node.value
        if isinstance(node, ast.BinOp) and type(node.op) in ops:
            return ops[type(node.op)](ev(node.left), ev(node.right))
        if isinstance(node, ast.UnaryOp) and type(node.op) in ops:
            return ops[type(node.op)](ev(node.operand))
        raise ValueError(node)

    try:
        return float(ev(ast.parse(expr, mode='eval').body))
    except Exception:
        raise argparse.ArgumentTypeError(
            f"{expr!r} is not a number or simple arithmetic expression")


def create_standard_parser():
    """ Creates a standard argument parser"""


    #from anuga.validation_utilities.parameters import cfl as default_cfl
    from anuga.validation_utilities.parameters import alg as default_alg

    import argparse
    parser = argparse.ArgumentParser(description='validation parse')

    #parser.add_argument('-cfl', type=float, default=default_cfl,
    #                   help='cfl condition')

    parser.add_argument('-ft', '--finaltime', type=arithmetic_float, default=argparse.SUPPRESS,
                       help='finaltime (accepts arithmetic, e.g. 3600*8)')

    parser.add_argument('-ys', '--yieldstep', type=arithmetic_float, default=argparse.SUPPRESS,
                       help='yieldstep (accepts arithmetic, e.g. 3600/4)')

    parser.add_argument('-os', '--outputstep', type=arithmetic_float, default=argparse.SUPPRESS,
                       help='outputstep (defaults to yieldstep if not set; accepts arithmetic)')

    parser.add_argument('-alg', type=str, default = default_alg,
                       help='flow algorithm')

    parser.add_argument('-np', type=int, default = 1,
                   help='number of processors to be used')

    parser.add_argument('-v', '--verbose', nargs='?', type=bool, const=True, default=False,
                   help='turn on verbosity')

    parser.add_argument('-d', '--debug', nargs='?', type=bool, const=True, default=False,
                   help='turn on debug-level reporting (e.g. mesh generation internals)')

    parser.add_argument('-l', '--long', nargs='?', type=bool, const=True, default=False,
                   help='include long-running validation tests (the HEC-RAS bridge/weir '
                        'behaviour cases) that run_auto_validation_tests.py skips by default')

    parser.add_argument('-nr', '--no-run', '--no_run', dest='no_run',
                   nargs='?', type=bool, const=True, default=False,
                   help='do not run the simulation: reuse the sww files already '
                        'in the case directory and only replot and typeset. Used '
                        'by produce_results.py and validations_produce_results.py')

    parser.add_argument('-cp', '--checkpointing', nargs='?', type=bool, const=True, default=False,
                   help='turn on checkpointing')

    parser.add_argument('-o', '--outname', type=str, default="domain",
                       help='sww name')

    parser.add_argument( '--partition_dir', type=str, default="PARTITIONS",
                       help='Directory for storing partitions')

    parser.add_argument( '--checkpoint_dir', type=str, default="CHECKPOINTS",
                       help='Directory for storing checkpoints')

    parser.add_argument('--checkpoint_time', type=float, default=-1.0,
                       help='checkpoint time')

    parser.add_argument('-mpm', '--multiprocessor_mode', type=int, default=argparse.SUPPRESS,
                       choices=[1, 2],
                       help='multiprocessor mode: 1=legacy CPU OpenMP, 2=unified gpu_ext '
                            'kernels. If unset, the domain keeps its default compute mode '
                            '(ANUGA_DEFAULT_COMPUTE_MODE, else legacy)')

    parser.add_argument('-nt', '--omp_num_threads', type=int, default=argparse.SUPPRESS,
                       help='number of OpenMP threads (process-wide). '
                            'Defaults to the OMP_NUM_THREADS environment variable if unset')

    # GPU offload for mode-2 (unified) domains, process-wide. Tri-state: leave
    # unset to follow the build (offload on a GPU build), or force on/off.
    parser.add_argument('-go', '--gpu_offload', dest='gpu_offload',
                       action='store_const', const=True, default=None,
                       help='enable GPU offload for mode-2 (unified) domains (process-wide)')
    parser.add_argument('-ngo', '--no-gpu_offload', '--no_gpu_offload', dest='gpu_offload',
                       action='store_const', const=False,
                       help='force CPU multicore for mode-2 (unified) domains '
                            '(process-wide); shortcut for the -go off case. '
                            'Default: follow the build (offload on a GPU build)')

    parser.add_argument('-ps', '--partition_scheme', type=str, default='metis',
                       choices=['metis', 'hilbert', 'morton', 'rcm'],
                       help='partition scheme used by distribute() '
                            '(metis, hilbert, morton, rcm)')

    parser.add_argument('-ro', '--reorder', type=str, default='none',
                       choices=['none', 'hilbert', 'morton', 'rcm', 'metis', 'metis_hilbert', 'metis_rcm'],
                       help='reorder triangles for cache locality before evolving '
                            '(none, hilbert, morton, rcm, metis, metis_hilbert, metis_rcm)')

    parser.add_argument('-rn', '--reorder_nprocs', type=int, default=argparse.SUPPRESS,
                       help='number of Metis partitions for metis/metis_hilbert/metis_rcm reordering '
                            '(defaults to OMP_NUM_THREADS if not set)')

    return parser


def parse_standard_args():
    """ Parse arguments for standard validation
    arguments. Returns values of

    alg

    """

    parser = create_standard_parser()

    args = parser.parse_args()

    #cfl = args.cfl
    alg = args.alg
    verbose= args.verbose
    np = args.np


    return alg











if __name__ == "__main__":
    print("Hello World")
