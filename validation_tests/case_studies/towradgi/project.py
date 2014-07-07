from os.path import join


model_output_dir='MODEL_OUTPUTS'
partition_dir = 'PARTITIONS'

basename = join('DEM_bridges', 'towradgi')
outname = join('Towradgi_historic_flood')
meshname = join('DEM_bridges','towradgi.tsh')

channel_manning=0.03
maximum_triangle_area = 10000 #= 1000
base_friction = 0.04
alpha = 0.99

W=303517
N=6195670
E=308570
S=6193140




