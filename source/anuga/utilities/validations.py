# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="steve"
__date__ ="$10/07/2012 1:18:38 PM$"


def validation_parse():
    """ Parse arguments for standard validation
    arguments. Returns values of

    alg
    cfl

    """

    import argparse
    parser = argparse.ArgumentParser(description='validation parse')
    parser.add_argument('-cfl', type=float, default=1.0,
                       help='cfl condition')
    parser.add_argument('-alg', type=str, default = "1_5",
                       help='flow algorithm')
    args = parser.parse_args()

    cfl = args.cfl
    alg = args.alg


    return alg, cfl




def run_validation_script(script):

    alg, cfl = validation_parse()


    import os
    cmd = 'python %s -alg %s -cfl %s '% (script,alg,cfl)
    print 'Running: '+cmd
    try:
        os.system( cmd )
    except:
        print 'Failed running '+ script +' in '+os.getcwd()
        pass










if __name__ == "__main__":
    print "Hello World"
