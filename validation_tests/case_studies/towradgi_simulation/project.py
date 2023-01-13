from os.path import join


model_output_dir='MODEL_OUTPUTS'
partition_dir = 'PARTITIONS'
checkpoint_dir = 'CHECKPOINTS'
checkpoint_time = 30 # 30*60 # 30 minutes
checkpointing = True
verbose = True

finaltime = 30.0
yieldstep = 10.0

basename = join('DEM_bridges', 'towradgi')
outname = join('Towradgi_historic_flood')
meshname = join('DEM_bridges','towradgi.tsh')

channel_manning=0.03
maximum_triangle_area = 1000
base_friction = 0.04
alpha = 0.99

W=303517
N=6195670
E=308570
S=6193140
