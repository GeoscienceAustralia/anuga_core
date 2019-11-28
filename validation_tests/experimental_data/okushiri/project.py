"""Common filenames for Okushiri Island validation
Formats are given as ANUGA native netCDF where applicable.

"""

# Given boundary wave
boundary_filename = 'Benchmark_2_input.tms'

# Observed timeseries
validation_filename = 'output_ch5-7-9.txt'

# Digital Elevation Model
bathymetry_filename_stem = 'Benchmark_2_Bathymetry'
bathymetry_filename = bathymetry_filename_stem + '.pts'

# Triangular mesh
mesh_filename = 'Benchmark_2.msh'

# Model output
output_filename = 'okushiri_auto_validation.sww'

# Evolve variables
finalTime = 25.0
yieldStep = 0.05




