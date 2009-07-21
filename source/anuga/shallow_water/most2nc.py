"""This has to do with creating elevation data files for use with ferret2sww.
It reads a bathymetry ascii file and creates a NetCDF (nc) file similar to
MOSTs output.

 $Author: Peter Row
 
"""

import sys
from Scientific.IO.NetCDF import NetCDFFile
from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a
import anuga.utilities.log as log


##
# @brief Convert a MOST file to NetCDF format.
# @param input_file The input file to convert.
# @param output_file The name of the oputput NetCDF file to create,
# @param inverted_bathymetry ??
# @param verbose True if the function is to be verbose.
def most2nc(input_file, output_file, inverted_bathymetry=False, verbose=True):
    # variable names
    long_name = 'LON'
    lat_name = 'LAT'
    elev_name = 'ELEVATION'

    # set up bathymetry
    if inverted_bathymetry:
	up = -1.
    else:
	up = +1.

    # read data from the MOST file
    in_file = open(input_file,'r')

    if verbose: log.critical('reading header')

    nx_ny_str = in_file.readline()
    nx_str,ny_str = nx_ny_str.split()
    nx = int(nx_str)
    ny = int(ny_str)
    h1_list=[]
    for i in range(nx):
	h1_list.append(float(in_file.readline()))

    h2_list=[]
    for j in range(ny):
	h2_list.append(float(in_file.readline()))

    h2_list.reverse()

    if verbose: log.critical('reading depths')

    in_depth_list = in_file.readlines()
    in_file.close()

    out_depth_list = [[]]

    if verbose: log.critical('processing depths')

    k=1
    for in_line in in_depth_list:
	for string in in_line.split():
	    #j = k/nx
	    out_depth_list[(k-1)/nx].append(float(string)*up)
	    if k==nx*ny:
		break
	    if k-(k/nx)*nx ==0:
		out_depth_list.append([])
	    k+=1

    in_file.close()
    out_depth_list.reverse()
    depth_list = out_depth_list

    # write the NetCDF file
    if verbose: log.critical('writing results')

    out_file = NetCDFFile(output_file, netcdf_mode_w)

    out_file.createDimension(long_name,nx)

    out_file.createVariable(long_name,'d',(long_name,))
    out_file.variables[long_name].point_spacing='uneven'
    out_file.variables[long_name].units='degrees_east'
    out_file.variables[long_name].assignValue(h1_list)

    out_file.createDimension(lat_name,ny)
    out_file.createVariable(lat_name,'d',(lat_name,))
    out_file.variables[lat_name].point_spacing='uneven'
    out_file.variables[lat_name].units='degrees_north'
    out_file.variables[lat_name].assignValue(h2_list)

    out_file.createVariable(elev_name,'d',(lat_name,long_name))
    out_file.variables[elev_name].point_spacing='uneven'
    out_file.variables[elev_name].units='meters'
    out_file.variables[elev_name].assignValue(depth_list)

    out_file.close()
