"""Least squares smooting and interpolation.

   measure the speed of least squares.

   ________________________
   General comments
   
   The max_points_per_cell does effect the time spent solving a
   problem.  The best value to use is probably dependent on the number
   of triangles.  Maybe develop a simple imperical algorithm, based on
   test results.
   
   Duncan Gray
   Geoscience Australia, 2004.
"""

import os
import sys
import time
from random import seed, random
import tempfile
import profile , pstats
from math import sqrt

from anuga.fit_interpolate.interpolate import Interpolate
from anuga.fit_interpolate.fit import Fit
from anuga.pmesh.mesh import Mesh
from anuga.geospatial_data.geospatial_data import Geospatial_data
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.fit_interpolate.fit import Fit, fit_to_mesh
from anuga.fit_interpolate.interpolate import benchmark_interpolate
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.fit_interpolate.general_fit_interpolate import \
     get_build_quadtree_time


"""

Code from the web;

from ctypes import *
from ctypes.wintypes import DWORD

SIZE_T = c_ulong

class _MEMORYSTATUS(Structure):
_fields_ = [("dwLength", DWORD),
("dwMemoryLength", DWORD),
("dwTotalPhys", SIZE_T),
("dwAvailPhys", SIZE_T),
("dwTotalPageFile", SIZE_T),
("dwAvailPageFile", SIZE_T),
("dwTotalVirtual", SIZE_T),
("dwAvailVirtualPhys", SIZE_T)]
def show(self):
for field_name, field_type in self._fields_:
print field_name, getattr(self, field_name)

memstatus = _MEMORYSTATUS()
windll.kernel32.GlobalMemoryStatus(byref(memstatus ))
memstatus.show()


_______________________________

from ctypes import *
from ctypes.wintypes import *

class MEMORYSTATUS(Structure):
_fields_ = [
('dwLength', DWORD),
('dwMemoryLoad', DWORD),
('dwTotalPhys', DWORD),
('dwAvailPhys', DWORD),
('dwTotalPageFile', DWORD),
('dwAvailPageFile', DWORD),
('dwTotalVirtual', DWORD),
('dwAvailVirtual', DWORD),
]

def winmem():
x = MEMORYSTATUS()
windll.kernel32.GlobalMemoryStatus(byref(x))
return x

"""

def mem_usage():
    '''
    returns the rss.

  RSS  The total amount of physical memory used by  the  task,  in  kilo-
            bytes,  is  shown  here.  For ELF processes used library pages are
            counted here, for a.out processes not.
            
    Only works on nix systems.
    '''

    # GC stuff should be uncommented when debugging memory 'shift' problems,
    # when (if?) the amount of memory used takes a sudden leap upwards.
    #import gc
    import string

    #gc.collect()
    #print('Ran a garbage collection')
    p=os.popen('ps uwp %s'%os.getpid()) 
    lines=p.readlines()
    status=p.close() 
    if status or len(lines)!=2 or sys.platform == 'win32': 
        return None 
    return int(string.split(lines[1])[4]) 


class BenchmarkLeastSquares(object):
    r"""

    Note(DSG-DSG): If you are interested in benchmarking fitting, before
    and after blocking O:\1\dgray\before_blocking_subsandpit is before blocking

    """

    def __init__(self):
        pass

    def trial(self,
              num_of_points=20000,
              maxArea=1000,
              max_points_per_cell=13,
              is_fit=True,
              use_file_type=None,
              blocking_len=500000,
              segments_in_mesh=True,
              save=False,
              verbose=False,
              run_profile=False,
              gridded=True,
              geo_ref=True):
        '''
        num_of_points 
        '''
        if geo_ref is True:
            geo = Geo_reference(xllcorner = 2.0, yllcorner = 2.0)
        else:
            geo = None
        mesh_dict = self._build_regular_mesh_dict(maxArea=maxArea,
                                                  is_segments=segments_in_mesh,
                                                  save=save,
                                                  geo=geo)
        points_dict = self._build_points_dict(num_of_points=num_of_points,
                                              gridded=gridded,
                                              verbose=verbose)

        if is_fit is True:
            op = "Fit_"
        else:
            op = "Interp_"
        profile_file = op + "P" + str(num_of_points) + \
                       "T" + str(len(mesh_dict['triangles'])) + \
                       "PPC" + str(max_points_per_cell) + \
                       ".txt"
                    
        # Apply the geo_ref to the points, so they are relative
        # Pass in the geo_ref
        
        domain = Domain(mesh_dict['vertices'], mesh_dict['triangles'],
                        use_cache=False, verbose=verbose,
                                     geo_reference=geo)
        #Initial time and memory
        t0 = time.time()
        #m0 = None on windows
        m0 = mem_usage()
        
        # Apply the geo_ref to the points, so they are relative
        # Pass in the geo_ref
        geospatial = Geospatial_data(points_dict['points'],
                                     points_dict['point_attributes'],
                                     geo_reference=geo)
        del points_dict
        if is_fit is True:

            if use_file_type is None:
                points = geospatial
                filename = None
            else:
                #FIXME (DSG) check that the type
                fileName = tempfile.mktemp("." + use_file_type)
                geospatial.export_points_file(fileName, absolute=True)
                points = None
                filename = fileName
            if run_profile is True:
                    
                s = """domain.set_quantity('elevation',points,filename=filename,use_cache=False)"""
                pobject = profile.Profile()
                presult = pobject.runctx(s,
                                         vars(sys.modules[__name__]),
                                         vars())
                prof_file = tempfile.mktemp(".prof")
                presult.dump_stats(prof_file)
                #
                # Let process these results
                S = pstats.Stats(prof_file)
                saveout = sys.stdout 
                pfile = open(profile_file, "w")
                sys.stdout = pfile
                s = S.sort_stats('cumulative').print_stats(60)
                sys.stdout = saveout 
                pfile.close()
                os.remove(prof_file)
            else:
                domain.set_quantity('elevation',points,filename=filename,
                                    use_cache=False, verbose=verbose)
            if not use_file_type is None:
                os.remove(fileName)
                    
        else:
            # run an interploate problem.
            
            if run_profile:
                # pass in the geospatial points
                # and the mesh origin
                 
                s="""benchmark_interpolate(mesh_dict['vertices'],mesh_dict['vertex_attributes'],mesh_dict['triangles'],geospatial,max_points_per_cell=max_points_per_cell,mesh_origin=geo)"""
                pobject = profile.Profile()
                presult = pobject.runctx(s,
                                         vars(sys.modules[__name__]),
                                         vars())
                prof_file = tempfile.mktemp(".prof")
                presult.dump_stats(prof_file)
                #
                # Let process these results
                S = pstats.Stats(prof_file)
                saveout = sys.stdout 
                pfile = open(profile_file, "w")
                sys.stdout = pfile
                s = S.sort_stats('cumulative').print_stats(60)
                sys.stdout = saveout 
                pfile.close()
                os.remove(prof_file)
                    
            else:
                # pass in the geospatial points
                 benchmark_interpolate(mesh_dict['vertices'],
                                       mesh_dict['vertex_attributes'],
                                       mesh_dict['triangles'],
                                       geospatial,
                                       mesh_origin=geo,
                                       max_points_per_cell=max_points_per_cell,
                                       verbose=verbose)
        time_taken_sec = (time.time()-t0)
        m1 = mem_usage()
        if m0 is None or m1 is None:
            memory_used = None
        else:
            memory_used = (m1 - m0)
        #print 'That took %.2f seconds' %time_taken_sec

        # return the times spent in first cell searching and
        # backing up.
        
        #search_one_cell_time, search_more_cells_time = search_times()
        #reset_search_times()
        #print "bench - build_quadtree_time", get_build_quadtree_time()
        return time_taken_sec, memory_used, len(mesh_dict['triangles']), \
               get_build_quadtree_time()
    

    def _build_regular_mesh_dict(self,
                                 maxArea=1000,
                                 is_segments=True,
                                 save=False,
                                 geo=None):
      # make a normalised mesh
        # pretty regular size, with some segments thrown in.

        # don't pass in the geo ref.
        # it is applied in domain
        m = Mesh() #geo_reference=geo)
        m.addUserVertex(0,0)
        m.addUserVertex(1.0,0)
        m.addUserVertex(0,1.0)
        m.addUserVertex(1.0,1.0)
        
        m.auto_segment(alpha = 100 )

        if is_segments:
            dict = {}
            dict['points'] = [[.10,.10],[.90,.20]]
            dict['segments'] = [[0,1]] 
            dict['segment_tags'] = ['wall1']   
            m.addVertsSegs(dict)
    
            dict = {}
            dict['points'] = [[.10,.90],[.40,.20]]
            dict['segments'] = [[0,1]] 
            dict['segment_tags'] = ['wall2']   
            m.addVertsSegs(dict)
        
            dict = {}
            dict['points'] = [[.20,.90],[.60,.60]]
            dict['segments'] = [[0,1]] 
            dict['segment_tags'] = ['wall3'] 
            m.addVertsSegs(dict)
        
            dict = {}
            dict['points'] = [[.60,.20],[.90,.90]]
            dict['segments'] = [[0,1]] 
            dict['segment_tags'] = ['wall4']   
            m.addVertsSegs(dict)

        m.generateMesh(mode = "Q", maxArea = maxArea, minAngle=20.0)       
        if save is True:
            m.export_mesh_file("aaaa.tsh")
        mesh_dict =  m.Mesh2IOTriangulationDict()

        #Add vert attribute info to the mesh
        mesh_dict['vertex_attributes'] = []
        # There has to be a better way of doing this..
        for vertex in mesh_dict['vertices']:
            mesh_dict['vertex_attributes'].append([10.0])

        return mesh_dict

    def _build_points_dict(self, num_of_points=20000
                           , gridded=True, verbose=False):
        
        points_dict = {}
        points = []
        point_atts = []

        if gridded is True:
            grid = int(sqrt(num_of_points))
        
        for point in range(num_of_points):
            if gridded is True:

                # point starts at 0.0
                # the 2 and 0.25 is to make sure all points are in the
                # range 0 - 1
                points.append([float(point/grid)/float(grid*1.1)+0.0454,
                               float(point%grid)/float(grid*1.1)+0.0454])
            else:
                points.append([random(), random()])
            point_atts.append(10.0)

        points_dict['points'] = points
        points_dict['point_attributes'] = point_atts
        
        for point in points:
            assert point[0] < 1.0
            assert point[1] < 1.0
            assert point[0] > 0.0
            assert point[1] > 0.0
            
        if verbose is True:
            pass
            #print "points", points
            #import sys; sys.exit() 
            
                                    
            
        return points_dict


#-------------------------------------------------------------
if __name__ == "__main__":
        b._build_points_dict()
