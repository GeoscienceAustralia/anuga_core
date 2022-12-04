
# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$10/07/2012 1:18:38 PM$"




def create_standard_parser():
    """ Creates a standard argument parser"""


    #from anuga.validation_utilities.parameters import cfl as default_cfl
    from anuga.validation_utilities.parameters import alg as default_alg

    import argparse
    parser = argparse.ArgumentParser(description='validation parse')

    #parser.add_argument('-cfl', type=float, default=default_cfl,
    #                   help='cfl condition')

    parser.add_argument('-ft', '--finaltime', type=float, default=argparse.SUPPRESS,
                       help='finaltime')

    parser.add_argument('-ys', '--yieldstep', type=float, default=argparse.SUPPRESS,
                       help='yieldstep')

    parser.add_argument('-alg', type=str, default = default_alg,
                       help='flow algorithm')

    parser.add_argument('-np', type=int, default = 1,
                   help='number of processors to be used')

    parser.add_argument('-v', '--verbose', nargs='?', type=bool, const=True, default=False,
                   help='turn on verbosity')

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
