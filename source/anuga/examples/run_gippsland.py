from pyvolution.data_manager import asc_csiro2sww
import time

t0 = time.time()


asc_csiro2sww('bathymetry_expanded','elev_expanded','ucur_expanded',
              'vcur_expanded', 'test_extent_true_38.25__37.7__147__148.25.sww',zscale=1,
              mean_stage = 0,
              fail_on_NaN = False,
              elevation_NaN_filler = 0
              ,
              minlat = -38.25, maxlat = -37.7, ###-37.75,
              minlon = 147.0, maxlon = 148.25
              , verbose = True
              )

print 'That took %.2f seconds' %(time.time()-t0)
