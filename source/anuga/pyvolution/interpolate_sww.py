""" Interpolation of a sww file.

Used to interpolate height and other information from sww files.

When using as a stand alone application,
Input;
 - The sww file with stage and bed elevation information
 - An xya file specifying the points (x,y) where the height values
   (w.r.t. time) are needed.

Ouput;
An xya file with x,y and height values w.r.t. time.

The first row of the output xya file has the time value for each height
Each following row has the format x,y,[height,]

NOTE: stage = bed elevation + height

   Ole Nielsen, Stephen Roberts, Duncan Gray, Christopher Zoppou
   Geoscience Australia, 2004.   
"""

##FIXME (DSG-DSG)  no sww file? give a better error message.

from Numeric import transpose
from fit_interpolate.interpolate import Interpolate 

DEFAULT_QUANTITY = "depth"

def  interpolate_sww2xya(sww_file, quantity_name,point_file_in,
                         point_file_out, verbose = True):
    """
    This function catches exceptions.
    """
    try:
        interp = Interpolate_sww(sww_file, quantity_name)
        interp.interpolate_xya(point_file_in)
        interp.write_depth_xya(point_file_out)
    except IOError,e: #need to convert key error to ioerror
        if verbose:
            print "Could not load bad file. ", e
        import sys; sys.exit()

        # FIXME (DSG-DSG): how are bad quantities caught?
        #try:
        #    interp = Interpolate_sww(sww_file, quantity_name)
        #except KeyError:
        #    print "Error: Unknown quantity"
        #    sys.exit(1)

        interp.interpolate_xya(point_file_in)
        interp.write_depth_xya(point_file_out)



class Interpolate_sww:
    def __init__(self, file_name, quantity_name):

        #It's bad to have the quantity_name passed in here.
        # But it works with how the program is currently used.
        # Refactor when necessary. - DSG

        x, y, volumes, time, quantity = self.read_sww(file_name, quantity_name)
        vertex_coordinates = transpose([x,y])
        
        if False:
            print "****************************"
            print "x "
            print x
            print "****************************"
            print "Y "
            print y
            print "****************************"
            print "V "
            print volumes
            print "****************************"
            print "Time "
            print time
            print "****************************"
            print "quantity "
            print quantity
            print "****************************"
            print "vertex_coordinates"
            print vertex_coordinates
            print "****************************"
        self.vertex_coordinates = vertex_coordinates
        self.volumes = volumes
        self.time = time

        self.quantity_name = quantity_name
        self.quantity = quantity 
        
    def interpolate_xya(self, file_name):
        """
        Given a point file, interpolate the height w.r.t. time at the points
        specified in the point file.

        Input;
        file_name - the xya file
        """
        
        from load_mesh.loadASCII import import_points_file
        
        interp = Interpolate(self.vertex_coordinates, self.volumes)
        point_dict = import_points_file(file_name)
        self.point_coordinates = point_dict['pointlist']
        self.interpolated_quantity_raw = interp.interpolate(transpose(self.quantity),
                                                            point_coordinates = self.point_coordinates)
        #self.interpolated_quantity_raw = self.interp.interpolate(transpose(self.quantity))
        self.interpolated_quantity = {}
        for i,time_slice in enumerate(self.time):
            self.interpolated_quantity[str(time_slice)] = self.interpolated_quantity_raw[:,i]
        
            
    def read_sww(self, file_name, quantity_name):
        """
        Read in an sww file.

        Input;
        file_name - the sww file

        Output;
        x - Vector of x values
        y - Vector of y values
        z - Vector of bed elevation
        volumes - Array.  Each row has 3 values, representing
                  the vertices that define the volume
        time - Vector of the times where there is stage information
        stage - array with respect to time and vertices (x,y)
        """
        
        #FIXME Have this reader as part of data_manager?

        from Scientific.IO.NetCDF import NetCDFFile     
        import tempfile
  	import sys
        import os
            
        #Check contents
        #Get NetCDF
    
        # see if the file is there.  Throw a QUIET IO error if it isn't
        fd = open(file_name,'r')
        fd.close()

        #throws prints to screen if file not present
        
        junk = tempfile.mktemp(".txt")
        fd = open(junk,'w')
        stdout = sys.stdout
        sys.stdout = fd
        fid = NetCDFFile(file_name, 'r') 
        sys.stdout = stdout
        fd.close()
        #clean up
        os.remove(junk)    	
      
        # Get the variables
        x = fid.variables['x'][:]
        y = fid.variables['y'][:]
        volumes = fid.variables['volumes'][:] 
        time = fid.variables['time'][:]
        try:
            if quantity_name == 'depth':
                z = fid.variables['elevation'][:]
                stage = fid.variables['stage'][:,:]
                quantity = stage - z  # 2D, using broadcasting
                #print quantity
            else:
                quantity = fid.variables[quantity_name][:,:]   # 2D 
        except (KeyError, IndexError),e:
            fid.close()
            raise KeyError
            
        fid.close()
        return x, y, volumes, time, quantity

    def write_depth_xya(self,point_file, delimiter = ','):
        """
        pre condition:
        point attributes have been determined
        (interpolate_xya has been called)
        The time list is defined
        """
        from load_mesh.loadASCII import export_points_file

        xya_dict = {}
        xya_dict['pointlist'] = self.point_coordinates
        xya_dict['attributelist'] = self.interpolated_quantity
        export_points_file(point_file, xya_dict)
        
        
#-------------------------------------------------------------


if __name__ == "__main__":
    """
    Load in an sww file and an xya file and return an xya file 
    """
    import os, sys
    usage = "usage: %s pyvolution_results.sww points.xya depth.xya [depth|stage|(other quantities)]" %         os.path.basename(sys.argv[0])
    if len(sys.argv) < 4:
        print usage
    else:
        sww_file = sys.argv[1]
        point_file_in = sys.argv[2]
        point_file_out = sys.argv[3]
        if len(sys.argv) == 5:
            quantity_name = sys.argv[4]
        else:
            quantity_name = DEFAULT_QUANTITY   
        #print "quantity",quantity
        interpolate_sww2xya(sww_file, quantity_name,point_file_in,
                            point_file_out)







