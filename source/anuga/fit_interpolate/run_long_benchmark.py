"""
This runs some benchmark tests, using the data here as input.

The output is a txt file with timings and memory use info.

Check pyvolution.run_profile for an example of how to use the python profile
module.

"""
from benchmark_least_squares import BenchmarkLeastSquares

ben = BenchmarkLeastSquares()

ofile = 'lbm_results.txt'
delimiter = ','

use_least_squares_list = [False]
is_fit_list = [True]
#num_of_points_list = [10]
#maxArea_list = [0.1, 0.001]
num_of_points_list = [10, 100] #, 10000, 100000] #, 10000000]
maxArea_list = [0.1, 0.001] #, 0.00001, 0.0000001]
max_points_per_cell_list = [8]
use_file_type_list = [None,'txt']

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


for use_file_type in use_file_type_list:
    for is_fit in is_fit_list:
        for num_of_points in num_of_points_list:
            for maxArea in maxArea_list:
                for max_points_per_cell in max_points_per_cell_list:
    
                    time, mem, num_tri = ben.trial(num_of_points=num_of_points
                                                   ,maxArea=maxArea
                                                   ,max_points_per_cell=max_points_per_cell
                                                   ,is_fit=is_fit
                                                   ,use_file_type=use_file_type
                                               )
                    print "time",time
                    print "mem", mem
                    fd.write(str(use_file_type) + delimiter +
                             str(is_fit) + delimiter +
                             str(num_of_points) + delimiter +
                             str(maxArea) + delimiter +
                             str(num_tri) + delimiter +
                             str(max_points_per_cell) + delimiter +
                             str(mem)  + delimiter +
                             str(time) + delimiter + "\n")
fd.close()                         
