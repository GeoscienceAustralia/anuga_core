
from simulation import Simulation
import anuga
from os.path import join


"""
#===============================================================================
# Setup arguments (and defaults) for towradgi simulation
#===============================================================================
def add_args(parser):
    parser.description = "ANUGA `{0}' benchmark".format(domain_name)
    parser.add_argument('--checkpointing', type=str, help='Output filename')
    parser.add_argument('--partition_dir', type=str, help='Partition directory')
    parser.add_argument('--channel_manning', type=float, help='Manning Friction in Channel')
    parser.set_defaults(
        alg = 'DE0',          
        outname = domain_name,
        partition_dir = ,
        channel_manning = 
    )
"""
# parser = anuga.create_standard_parser()
# 
# project_dict = vars(project)
# del project_dict['__builtins__']
# del project_dict['__doc__']
# del project_dict['__file__']
# del project_dict['__name__']
# del project_dict['__package__']
# del project_dict['join']
# 
# 
# parser.set_defaults(**project_dict)
# 
# args = parser.parse_args()

args = anuga.get_args()
alg = args.alg
verbose = args.verbose



towradgi = Simulation()





domain = towradgi.domain

d_dict = domain.__dict__
from pprint import pprint

d_keys = d_dict.keys()
d_keys.sort()
#pprint(d_keys)

    
# def wrap(t):
#     return [func(t), 0.0, 0.0]
# 
# func = anuga.file_function(join('Forcing','Tide','Pioneer.tms'), quantities='rainfall')
# 
# def wrapped(t):
#     return [func(t), 0.0, 0.0]
# 
# Bd = anuga.Dirichlet_boundary([0,0,0])
# Bw = anuga.Time_boundary(domain=domain, function=wrapped)
# #Bw = anuga.Time_boundary(domain=domain, function=func)
# 
# domain.set_boundary({'west': Bd, 'south': Bd, 'north': Bd, 'east': Bw})
# 
# 
# domain.fractional_step_operators



domain.checkpoint = True
domain.checkpoint_step = 4



towradgi.run(yieldstep=1.0, finaltime=20.0)
