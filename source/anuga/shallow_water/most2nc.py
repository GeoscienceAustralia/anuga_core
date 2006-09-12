"""This has to do with creating elevation data files for use with ferret2sww.
It reads a bathymetry ascii file and creates a NetCDF (nc) file similar to MOSTs output.

 $Author: Peter Row
 
"""


def most2nc(input_file=None,output_file=None,inverted_bathymetry = False,\
            verbose = True):
    #input_file = 'small.txt'
    #output_file = 'small_e.nc'

    long_name = 'LON'
    lat_name = 'LAT'
    if inverted_bathymetry:
	up = -1.
    else:
	up = +1.

    from Scientific.IO.NetCDF import NetCDFFile
    import sys

    try:
	if input_file is None:
	    input_file = sys.argv[1]
	if output_file is None:
	    output_file = sys.argv[2]
    except:
	raise 'usage is: most2nc input output'

    in_file = open(input_file,'r')
    if verbose: print 'reading header'
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

    if verbose: print 'reading depths'

    in_depth_list = in_file.readlines()
    in_file.close()

    out_depth_list = [[]]

    if verbose: print 'processing depths'

    k=1
    for in_line in in_depth_list:
	for string in in_line.split():
	    #j = k/nx
	    out_depth_list[(k-1)/nx].append(float(string)*up)
	    #print k,len(out_depth_list),(k-1)/nx,out_depth_list[(k-1)/nx][-1],len(out_depth_list[(k-1)/nx])
	    if k==nx*ny:
		break
	    if k-(k/nx)*nx ==0:
		out_depth_list.append([])
	    k+=1

    in_file.close()
    out_depth_list.reverse()
    depth_list = out_depth_list

    if verbose: print 'writing results'

    out_file = NetCDFFile(output_file,'w')

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

    out_file.createVariable('ELEVATION','d',(lat_name,long_name))
    out_file.variables['ELEVATION'].point_spacing='uneven'
    out_file.variables['ELEVATION'].units='meters'
    out_file.variables['ELEVATION'].assignValue(depth_list)

    out_file.close()
