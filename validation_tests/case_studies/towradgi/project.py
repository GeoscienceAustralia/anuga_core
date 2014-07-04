from os.path import join


model_output_dir='MODEL_OUTPUTS'
partition_dir = 'Partitions'

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

def read_polygon_list(poly_list):
    # Alternative to read_polygon_dir -- allows us to control order of polygons
    from anuga import read_polygon
    
    result = []
    for i in range(len(poly_list)):
        result.append((read_polygon(poly_list[i][0]) , poly_list[i][1]))
    return result


