"""
This runs some benchmark tests, using the data here as input.

The output is a txt file with timings and memory use info.

Check pyvolution.run_profile for an example of how to use the python profile
module.

"""


from builtins import str
from .benchmark_least_squares import BenchmarkLeastSquares

ben = BenchmarkLeastSquares()

ofile = 'lbm_resultsII.csv'
delimiter = ','
run_profile = False #True
is_fit_list = [True, False]
num_of_points_list = [3, 200, 600, 2000, 6000, 10000, 20000] 
#maxArea_list = [ 0.008, 0.0016, 0.0008]
max_points_per_cell_list = [2,4,8,16,30,64]
use_file_type_list = ['pts'] #'pts'
#num_of_points_list = [10]
maxArea_list = [ 0.008, 0.0016]
max_points_per_cell_list = [4]

fd = open(ofile,'a')
# write the title line

fd.write("use_file_type" + delimiter +
    "num_of_points" + delimiter +
         "maxArea" + delimiter +
         "num_of_triangles" + delimiter +
         "max_points_per_cell" + delimiter +
         "is_fit" + delimiter +
         "is_profiling" + delimiter +
         "mem"  + delimiter +
         "time" + delimiter + "\n")


for maxArea in maxArea_list:
    for use_file_type in use_file_type_list:
        for is_fit in is_fit_list:
            for num_of_points in num_of_points_list:
                for max_points_per_cell in max_points_per_cell_list:
    
                    time, mem, num_tri = ben.trial(num_of_points=num_of_points
                                                   ,maxArea=maxArea
                                                   ,max_points_per_cell=max_points_per_cell
                                                   ,is_fit=is_fit
                                                   ,segments_in_mesh=False
                                                   ,use_file_type=use_file_type
                                                   ,save=True
                                                   ,run_profile=run_profile
                                               )
                    print("time",time)
                    print("mem", mem)
                    fd.write(str(use_file_type) + delimiter +
                             str(num_of_points) + delimiter +
                             str(maxArea) + delimiter +
                             str(num_tri) + delimiter +
                             str(max_points_per_cell) + delimiter +
                             str(is_fit) + delimiter +
                             str(run_profile) + delimiter +
                             str(mem)  + delimiter +
                             str(time) + delimiter + "\n")
fd.close()                         
