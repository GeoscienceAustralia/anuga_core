#--------------------------------
# Setup Default values for basis
# algorithm parameters.
#--------------------------------
import argparse
parser = argparse.ArgumentParser(description='produce results')
parser.add_argument('-cfl', type=float, default=1.0,
                   help='cfl condition')
parser.add_argument('-alg', type=str, default = "1_5",
                   help='flow algorithm')
args = parser.parse_args()

cfl = args.cfl
alg = args.alg


import run_wave
import plotme
