"""
This runs some benchmark tests, using the data here as input.

The output is a txt file with timings and memory use info.

Check pyvolution.run_profile for an example of how to use the python profile
module.

"""
from benchmark_least_squares import BenchmarkLeastSquares

ben = BenchmarkLeastSquares()

ofile = 'lbm_results.csv'
delimiter = ','

use_least_squares_list = [False]
is_fit_list = [True]
num_of_points_list = [1000] #, 500, 10000, 100000] #, 10000000]
maxArea_list = [0.0001] #,0.00001] #, 0.0000001] #, 0.06, 0.00001, 0.0000001]
max_points_per_cell_list = [8]
use_file_type_list = ['pts']

fd = open(ofile,'a')
# write the title line

fd.write("use_file_type" + delimiter +
    "num_of_points" + delimiter +
         "maxArea" + delimiter +
         "num_of_triangles" + delimiter +
         "max_points_per_cell" + delimiter +
         "is_fit" + delimiter +
         "mem"  + delimiter +
         "time" + delimiter + "\n")


for is_fit in is_fit_list:
    for maxArea in maxArea_list:
        for use_file_type in use_file_type_list:
            for num_of_points in num_of_points_list:
                for max_points_per_cell in max_points_per_cell_list:
    
                    time, mem, num_tri = ben.trial(
                        num_of_points=num_of_points
                        ,maxArea=maxArea
                        ,max_points_per_cell=max_points_per_cell
                        ,is_fit=is_fit
                        ,segments_in_mesh=False
                        ,use_file_type=use_file_type
                        ,save=True
                        ,run_profile=True
                        )
                    print "time",time
                    print "mem", mem
                    print "num_tri", num_tri
                    fd.write(str(use_file_type) + delimiter +
                             str(num_of_points) + delimiter +
                             str(maxArea) + delimiter +
                             str(num_tri) + delimiter +
                             str(max_points_per_cell) + delimiter +
                             str(is_fit) + delimiter +
                             str(mem)  + delimiter +
                             str(time) + delimiter + "\n")
fd.close()                         
