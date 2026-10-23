

import pstats
import os
import csv
import numpy as np


#=================================
# Collect timings
#=================================

def create_benchmark_csvfile(pstat_basename, openmp_threads, verbose=True):

    """
    Create a CSV file summarizing the benchmark timings from pstat files.
    Args:
        pstat_basename (str): The base name for the pstat files.
        openmp_threads (list): List of OpenMP thread counts to analyze.
        verbose (bool): If True, prints the DataFrame to the console.

    Example usage:
        create_benchmark_csvfile('profile_example', [1, 2, 4, 8], verbose=True)
        This will create a CSV file named 'profile_example.csv' with the benchmark results.
        The verbose flag will print the DataFrame to the console.
    Note: Ensure that the pstat files are in the same directory as this script or provide the full path.
    """

    output_csv = pstat_basename+'.csv'

    table_contents = []

    # FIXME SR: Need to update script to cover cases with evolve_one_euler_step,
    # evolve_one_rk2_step, etc.
    myfuncs = ['OMP_NUM_THREADS',
    'total_time',
    'evolve',
    'evolve_one_euler_step',
    'evolve_one_rk2_step',
    'compute_fluxes',
    'distribute_to_vertices_and_edges',
    'update_conserved_quantities',
    'compute_forcing_terms',
    'update_boundary',
    'backup_conserved_quantities',
    'saxpy_conserved_quantities',
    'apply_fractional_steps'
    ]

    column_header = ['THREADS',
    'total_time',
    'evolve',
    'rk1_step',
    'rk2_step',
    'fluxes',
    'dist',
    'update',
    'force',
    'bound',
    'backup',
    'saxpy',
    'operators'
    ]

    table_contents.append(column_header)

    for threads in openmp_threads:

        pstat_file = pstat_basename+f'_omp_{threads}.pstat'

        print(f'Analysing file {pstat_file}')

        stats = pstats.Stats(pstat_file)

        benchmark_dict = {}

        for func_key, stat in stats.stats.items():
            filename, line, funcname = func_key
            for myfuncname in myfuncs:
                if funcname.endswith(myfuncname):
                    #print(f'{funcname}, {stat[3]:.3g}')
                    benchmark_dict[funcname] = stat[3]


        benchmark_dict['total_time']=stats.total_tt
        benchmark_dict['OMP_NUM_THREADS'] = threads


        # for key in myfuncs:
        #     try:
        #         print(f'{key} {benchmark_dict[key]:.3g}')
        #     except KeyError:
        #         print(f'{key} not found')

        table_line =[]

        for key in myfuncs:
            try:
                table_line.append(f'{benchmark_dict[key]:.3g}')
            except KeyError:
                table_line.append('nan')

        table_contents.append(table_line)



    # Write the cumulative times to a CSV file
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(table_contents[0])  # Write the header
        for row in table_contents[1:]:
            writer.writerow(row)  # Write the data rows


    if verbose:
        print(f'Wrote benchmark results to {output_csv}')
        import pandas as pd

        df = pd.read_csv(f'{pstat_basename}.csv')

        print('')
        print(80 * '=')
        print(f'Benchmark results for {pstat_basename}')
        print(80 * '=')
        print('')
        print(df.to_string(index=False))  # Prints the entire DataFrame without row numbers

if __name__ == '__main__':
    import sys

    if len(sys.argv) != 3:
        print("Usage: python create_benchmark_csvfile.py <pstat_basename> <openmp_threads>")
        sys.exit(1)

    pstat_basename = sys.argv[1]
    openmp_threads = list(map(int, sys.argv[2].split(',')))

    create_benchmark_csvfile(pstat_basename, openmp_threads, verbose=True)
    # Example usage:
    # python create_benchmark_csvfile.py pstat_basename 1,2,4,8
    # This will create a CSV file named 'pstat_basename.csv' with the benchmark results.

