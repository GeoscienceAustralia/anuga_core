"""
EA Test 6A limiter study: the first 40 s with the slope limiter's beta set
to 0, 0.5 and 1, for the report's figures of how the limiter moves the
hydraulic jump in front of the building.

Writes test6a_b<beta>_gauges.csv, test6a_b<beta>.sww and the rest of
run()'s outputs for each beta. produce_results.py runs it before the main
run (which writes parameters.tex last); the nightly validation does not.
"""

import os
import sys

import anuga

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import ea_dambreak  # noqa: E402

BETAS = (0.0, 0.5, 1.0)
FINAL_TIME = 40.0

if __name__ == '__main__':
    args = anuga.get_args()
    for beta in BETAS:
        print('limiter study: beta = %g' % beta)
        ea_dambreak.run('6A', args.alg, args.verbose, beta=beta, final_time=FINAL_TIME)
