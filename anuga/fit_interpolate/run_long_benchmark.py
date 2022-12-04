"""
This runs some benchmark tests, using the data here as input.

The output is a txt file with timings and memory use info.

Check pyvolution.run_profile for an example of how to use the python profile
module.

"""


from builtins import str
from .benchmark_least_squares import BenchmarkLeastSquares

ben = BenchmarkLeastSquares()

delimiter = ','

use_least_squares_list = [False]
is_fit_list = [True] #[True, False]

# a maxArea of 0.00001 gives 155195 triangles
# Simulating Cairns. 
 
#maxArea_list = [0.00001]
#num_of_points_list = [1863558]

# Simulating 1 tenth of Cairns. 
 
#maxArea_list = [0.0001]
#num_of_points_list = [186558]

# a maxArea of 0.0001 gives 15568 triangles
#maxArea_list = [0.0001]
#num_of_points_list = [46704]

# a maxArea of 0.001 gives 1539 triangles
#   3 points/tri is 4617
#  20 points/tri is 30780
# 132 points/tri is 203148
#maxArea_list = [0.001] 
#num_of_points_list = [4617,30780,203148]
#num_of_points_list = [203148]

# a maxArea of 0.005 gives 319 triangles
#   3 points/tri is 957
#  20 points/tri is 6380
# 132 points/tri is 42108
#maxArea_list = [0.005] 
num_of_points_list = [957,6380,42108]

# a maxArea of 0.01 gives 150 triangles
#   3 points/tri is 450
#  20 points/tri is 3000
# 132 points/tri is 19800
#maxArea_list = [0.01] 
#num_of_points_list = [450,3000] #,19800]

# Quick check
#maxArea_list = [0.61] 
#num_of_points_list = [4] 

# the auto-validate benchmark fit
maxArea_list = [0.0001] 
num_of_points_list = [1000]



max_points_per_cell_list = [4]
use_file_type_list = ['pts']
run_profile =  False # True #True # False #
gridded_list = [True, False]
geo_ref_list = [True, False]

if run_profile is True:
    ofile = 'profiling_lbm_results.csv'
else:
    ofile = 'lbm_results.csv'
fd = open(ofile,'a')
# write the title line

fd.write("use_file_type" + delimiter +
    "num_of_points" + delimiter +
         "maxArea" + delimiter +
         "num_of_triangles" + delimiter +
         "max_points_per_cell" + delimiter +
         "is_fit" + delimiter +
         "is_gridded" + delimiter +
         "has_geo_ref" + delimiter +
         "search_one_cell_time" + delimiter +
         "search_more_cells_time" + delimiter +
         "build_quadtree_time" + delimiter +
         "mem"  + delimiter +
         "time" + delimiter + "\n")


for is_fit in is_fit_list:
    for gridded in gridded_list:
        for geo_ref in geo_ref_list:
            for maxArea in maxArea_list:
                for use_file_type in use_file_type_list:
                    for num_of_points in num_of_points_list:
                        for max_points_per_cell in max_points_per_cell_list:
    
                            time, mem, num_tri, one_t, more_t, quad_t = ben.trial(
                                num_of_points=num_of_points
                                ,maxArea=maxArea
                                ,max_points_per_cell=max_points_per_cell
                                ,is_fit=is_fit
                                ,segments_in_mesh=False
                                ,use_file_type=use_file_type
                                ,save=True
                                ,verbose=False
                                ,run_profile=run_profile
                                ,gridded=gridded
                                ,geo_ref=geo_ref
                                )
                            print("time",time)
                            print("mem", mem)
                            print("num_tri", num_tri)
                            fd.write(str(use_file_type) + delimiter +
                                     str(num_of_points) + delimiter +
                                     str(maxArea) + delimiter +
                                     str(num_tri) + delimiter +
                                     str(max_points_per_cell) + delimiter +
                                     str(is_fit) + delimiter +
                                     str(gridded) + delimiter +
                                     str(geo_ref) + delimiter +
                                     str(one_t) + delimiter +
                                     str(more_t) + delimiter +
                                     str(quad_t) + delimiter +
                                     str(mem)  + delimiter +
                                     str(time) + delimiter + "\n")
fd.close()                         
