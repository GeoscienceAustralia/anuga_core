"""This script tests that caching works for set_quantity using point 
data from a file.
First cache is cleared, then set_quantity is run twice checking that
fitting is evaluated only first time and that the result from cache on the
second round is correct.

This script depends on  Benchmark_2.msh and Benchmark_2_Bathymetry.pts
"""

# Module imports
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.caching import cache
from anuga.fit_interpolate.fit import _fit_to_mesh
import project
import numpy as num
import time

internal_verbose = True # Verbosity used within this function

filename=project.bathymetry_filename
alpha=0.02
from anuga.config import points_file_block_line_size as max_read_lines

#-------------------------
# Create Domain from mesh
#-------------------------
domain = cache(Domain, 
               (project.mesh_filename, 
                {'verbose': True}), 
               verbose=False)

# Clear caching of underlying function                                            
args = (filename, )
kwargs = {'vertex_coordinates': None,
          'triangles': None,
          'mesh': domain.mesh,
          'point_attributes': None,
          'alpha': alpha,
          'verbose': internal_verbose,
          'mesh_origin': None,
          'data_origin': None,
          'max_read_lines': max_read_lines,
          'attribute_name': None 
          }


cache(_fit_to_mesh,
      args, 
      kwargs,
      verbose=False,
      dependencies=[filename],
      clear=True)

# Check that cache is empty      
flag = cache(_fit_to_mesh,
             args, 
             kwargs,
             verbose=False,
             dependencies=[filename],
             test=True)
assert flag is None


#-------------------------
# Initial Conditions
#-------------------------
t0 = time.time()
domain.set_quantity('elevation',
                    filename=filename,
                    alpha=0.02,                    
                    verbose=internal_verbose,
                    use_cache=True)
compute_time = time.time()-t0
                    
ref = domain.get_quantity('elevation').get_values()                    

# Check that cache is now present (and correct)
flag = cache(_fit_to_mesh,
             args, 
             kwargs,
             verbose=False,
             dependencies=[filename],
             test=True)
             
assert flag is not None
res = domain.get_quantity('elevation').get_values()
assert num.allclose(res, ref)

# Now check this using the high level call
print('Try to read in via cache')
t0 = time.time()
domain.set_quantity('elevation',
                    filename=filename,
                    alpha=0.02,                    
                    verbose=internal_verbose,
                    use_cache=True)
cache_time = time.time()-t0                    
                    
res = domain.get_quantity('elevation').get_values() 
assert num.allclose(res, ref)

print('cache_time', cache_time)
print('compute_time', compute_time)

msg = 'Caching did not speed things up as expected'
msg += 'Compute time = %.f, Cache time = %.f' % (compute_time,
                                                 cache_time)
assert cache_time < compute_time/10, msg
