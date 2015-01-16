"""Run a tsunami inundation scenario for Patong, as in Jakeman et al (2010)

The scenario is defined by a triangular mesh created from project.polygon, the
elevation data is compiled into a pts file through build_elevation.py and a
simulated tsunami is generated through an sts file from build_boundary.py.

Ole Nielsen and Duncan Gray, GA - 2005, Jane Sexton, Nick Bartzis, GA - 2006
Ole Nielsen, Jane Sexton and Kristy Van Putten - 2008
Gareth Davies -- 2013
Steve Roberts -- 2014
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

# Standard modules
import os
import os.path
import time
from time import localtime, strftime, gmtime

# Related major packages
from anuga.file.netcdf import NetCDFFile

import anuga




import anuga.utilities.log as log

# Parallel routines
from anuga import distribute, myid, numprocs, finalize, barrier

if myid == 0 and not os.path.isdir('topographies'):
    msg = """
################################################################################
#
# Could not the find data directories
#
# You can download these directories using the data_download.py script.
# This will download over 300 MB of data!
#
################################################################################
"""
    raise Exception(msg)    
    
    
import project  

# Application specific imports
import build_urs_boundary as bub

from anuga.file.csv_file import load_csv_as_building_polygons

import shutil

#--------------------------------------------------
# Pick up useful command line arguments (which over rides 
# values set before
#--------------------------------------------------
args = anuga.get_args()

alg = args.alg
verbose = args.verbose


if myid == 0 and verbose and numprocs == 1:
    print 80*'#'
    print '#'
    print '# Long Validation Test, takes 30 minutes on my desktop'
    print '#'
    print '# Consider running in parallel'
    print '#'
    print 80*'#'


#-------------------------------------------------------------------------------
# Copy scripts to time stamped output directory and capture screen
# output to file. Copy script must be before screen_catcher
#-------------------------------------------------------------------------------
if(myid==0):
    # Make output dir and set log filename before anything is logged
    try:
        shutil.rmtree(project.output_run)
    except:
        pass

    os.mkdir(project.output_run)


# Tell log module to store log file in output dir
log.log_filename = os.path.join(project.output_run, 'anuga_P_%g.log'%myid)
if myid == 0:
    log.console_logging_level = log.CRITICAL
else:
    log.console_logging_level = log.CRITICAL+1

if(myid==0):
    log.critical('Log filename: %s' % log.log_filename)
    
    
    output_dirname = os.path.dirname(project.__file__)
    anuga.copy_code_files(project.output_run,
                    [__file__, 
                     os.path.join(output_dirname, project.__name__+'.py')],
                    verbose=True
                    )


    #-------------------------------------------------------------------------------
    # Create the computational domain based on overall clipping polygon with
    # a tagged boundary and interior regions defined in project.py along with
    # resolutions (maximal area of per triangle) for each polygon
    #-------------------------------------------------------------------------------

    log.critical('Create computational domain')

    # Create the STS file
    # FIXME (Ole): This is deadly dangerous if buildcode changes (as was the case 24th March 2009)
    # We need to use caching instead!
    log.critical( 'project.mux_data_folder=%s' % project.mux_data_folder)
    if not os.path.exists(project.event_sts + '.sts'):
        bub.build_urs_boundary(project.mux_input_filename, project.event_sts)
else:
    print 'Hello from processor ', myid

barrier()
    

# Read in boundary from ordered sts file
event_sts = anuga.create_sts_boundary(project.event_sts)

# Reading the landward defined points, this incorporates the original clipping
# polygon minus the 100m contour
landward_boundary = anuga.read_polygon(project.landward_boundary)

# Combine sts polyline with landward points
bounding_polygon_sts = event_sts + landward_boundary

# Number of boundary segments
num_ocean_segments = len(event_sts) - 1
# Number of landward_boundary points
num_land_points = anuga.file_length(project.landward_boundary)


#======================================
# Create sequential domain
#======================================
if(myid==0):
    # Boundary tags refer to project.landward_boundary
    # 4 points equals 5 segments start at N
    boundary_tags={'back': range(num_ocean_segments+1,
                             num_ocean_segments+num_land_points),
               'side': [num_ocean_segments,
                        num_ocean_segments+num_land_points],
               'ocean': range(num_ocean_segments)}

    # Build mesh and domain
    domain = anuga.create_domain_from_regions(bounding_polygon_sts,
                                    boundary_tags=boundary_tags,
                                    maximum_triangle_area=project.bounding_maxarea,
                                    interior_regions=project.interior_regions,
                                    mesh_filename=project.meshes,
                                    use_cache=False,
                                    verbose=verbose)
    log.critical(domain.statistics())

    # FIXME(Ole): How can we make this more automatic?
    domain.geo_reference.zone = project.zone


    domain.set_name(project.scenario_name)
    domain.set_datadir(project.output_run) 

    domain.set_flow_algorithm(alg)

    #-------------------------------------------------------------------------------
    # Setup initial conditions
    #-------------------------------------------------------------------------------

    log.critical('Setup initial conditions')

    # Set the initial stage in the offcoast region only
    if project.land_initial_conditions:
        print '**********  IC ***********************************'
        IC = anuga.Polygon_function(project.land_initial_conditions,
                              default=project.tide,
                              geo_reference=domain.geo_reference)
    else:
        IC = project.tide


    domain.set_quantity('friction', project.friction) 
    
    import zipfile
    import os
    txt_filename =  os.path.basename(project.combined_elevation+'.txt')
    zipfile.ZipFile(project.combined_elevation+'.zip').extract(txt_filename,'topographies')

    if verbose: print 'Reading pts from txt'
    anuga.xya2pts(project.combined_elevation+'.txt', verbose = verbose)

    domain.set_quantity('elevation', 
                        filename=project.combined_elevation+'.pts',
                        use_cache=True,
                        verbose=True,
                        alpha=project.alpha)
    Stage = domain.quantities['stage']
    Elev  = domain.quantities['elevation']
    
    domain.set_quantity('stage', 0.8, use_cache=False, verbose=True)
    
    

    if project.use_buildings:
        # Add buildings from file
        log.critical('Reading building polygons')
        
        building_polygons, building_heights = load_csv_as_building_polygons(project.building_polygon, floor_height=3.)
        #clipping_polygons=project.building_area_polygons)

        log.critical('Creating %d building polygons' % len(building_polygons))
        def create_polygon_function(building_polygons, geo_reference=None):
            L = []
            for i, key in enumerate(building_polygons):
                if i%100==0: log.critical(i)
                poly = building_polygons[key]
                elev = building_heights[key]
                L.append((poly, elev))
                
                buildings = anuga.Polygon_function(L, default=0.0,
                                             geo_reference=geo_reference)
            return buildings

        log.critical('Creating %d building polygons' % len(building_polygons))
        #buildings = create_polygon_function(building_polygons)
        buildings = anuga.cache(create_polygon_function,
                          building_polygons,
                          {'geo_reference': domain.geo_reference},
                          verbose=True)

        log.critical('Adding buildings')
        domain.add_quantity('elevation',
                            buildings,
                            use_cache=False,
                            verbose=True)


else:
    domain=None
    #print 'Hello from Processor ', myid

barrier()


domain=distribute(domain) 
   
#-------------------------------------------------------------------------------
# Setup boundary conditions 
#-------------------------------------------------------------------------------

log.critical('Set boundary P_%g - available tags: %s' % (myid, domain.get_boundary_tags()))

Br = anuga.Reflective_boundary(domain)
Bs = anuga.Transmissive_stage_zero_momentum_boundary(domain)

if myid == 0 and verbose:
    verbose_bf = True
else:
    verbose_bf = False
    
Bf = anuga.Field_boundary(project.event_sts+'.sts',
                    domain,
                    mean_stage=project.tide,
                    time_thinning=20,
                    default_boundary=anuga.Dirichlet_boundary([0, 0, 0]),
                    boundary_polygon=bounding_polygon_sts,                    
                    use_cache=False,
                    verbose=verbose_bf)

domain.set_boundary({'back': Br,
                     'side': Bs,
                     'ocean': Bf}) 

#-------------------------------------------------------------------------------
# Evolve system through time
#-------------------------------------------------------------------------------

t0 = time.time()


## Jump up to time 400 sec to cover missing boundary data before 302 sec
for t in domain.evolve(yieldstep=400, finaltime=400):
    log.critical(domain.timestepping_statistics())
#    log.critical(domain.boundary_statistics(tags='ocean'))


import time

# Start detailed model
for t in domain.evolve(yieldstep=project.yieldstep,
                       finaltime=project.finaltime):

    log.critical(domain.timestepping_statistics())
    #log.critical(domain.boundary_statistics(tags='ocean'))


## Final print out
barrier()
if(myid==0):    
    log.critical('Simulation took %.2f seconds' %(time.time()-t0))


domain.sww_merge(delete_old=True)

finalize()


      
