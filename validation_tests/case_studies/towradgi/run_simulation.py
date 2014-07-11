
from simulation import Simulation
import anuga
from os.path import join



args = anuga.get_args()
alg = args.alg
verbose = args.verbose



towradgi = Simulation()

towradgi.run()

