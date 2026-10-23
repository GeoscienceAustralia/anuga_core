#!/usr/bin/env python3
"""
Modified run_model_3.py to use pre-partitioned mesh

CHANGES FROM ORIGINAL:
1. Removed domain creation on rank 0 (lines 159-221)
2. Removed distribute() call (line 233)
3. Added sequential_distribute_load() - all ranks load in parallel
4. Everything else stays the same

Author: Girishchandra Y. (modified for partitioned mesh)
"""
# ------------------------------------------------------------------------------
# Import necessary modules
# ------------------------------------------------------------------------------
import os
import time
import sys
import numpy as np
import linecache
import itertools
import datetime
from discharge_data import setup_discharge_3
from tide_data import setup_tide_3
import anuga
from anuga import Inlet_operator
from anuga import file_function
from anuga import Rate_operator
from anuga import Quantity
from anuga import myid, numprocs, finalize, barrier
from anuga.utilities import log as log
from anuga.abstract_2d_finite_volumes import quantity

t0 = time.time()
# --------------------------------------------------------------------------
# Setup parameters
# --------------------------------------------------------------------------
if len(sys.argv) > 1:
    input = str(sys.argv[1]).split("-")
    Current_Date = datetime.datetime(int(input[2]), int(input[1]), int(input[0])).replace(minute=0, hour=0, second=0, microsecond=0)
    Previous_Date = datetime.datetime(int(input[2]), int(input[1]), int(input[0])).replace(minute=0, hour=0, second=0, microsecond=0) - datetime.timedelta(days=1)
    Previous_hotstart_Date = datetime.datetime(int(input[2]), int(input[1]), int(input[0])).replace(minute=0, hour=0, second=0, microsecond=0) - datetime.timedelta(days=2)
else:
    #Current_Date = datetime.datetime(2021, 9, 10).replace(minute=0, hour=0, second=0, microsecond=0)
    Current_Date = datetime.datetime.today().replace(minute=0, hour=0, second=0, microsecond=0)
    Previous_Date = datetime.datetime.today().replace(minute=0, hour=0, second=0, microsecond=0) - datetime.timedelta(days=1)
    Previous_hotstart_Date = datetime.datetime.today().replace(minute=0, hour=0, second=0, microsecond=0) - datetime.timedelta(days=2)

bypass = True
TIMEFORMATCU = "%d_%m_%Y"

# CONFIGURATION: Set mesh resolution here (100, 300, or 900)
MESH_SIZE = 900 # <-- CHANGE THIS for different mesh sizes

name_stem = f'delta11372sqkm_uniform_mesh_{MESH_SIZE}sqm_Chilka{MESH_SIZE}sqm'
partition_name = f'mahanadi_delta_{MESH_SIZE}sqm'  # Must match partition_mahanadi.py output
partition_dir = './bench_partitions_check'

poly = anuga.read_polygon('polygons/Chilka_Modified_Poly.csv')
poly_bounding = anuga.read_polygon('polygons/Delta_11372_sqkm.csv')
chilka_stage = 0.15
cache = True
verbose = True
intermediate_result_extraction = True
domain_name = 'mahanadi_delta_' + Current_Date.strftime(TIMEFORMATCU)
output_dir = 'output/' + Current_Date.strftime(TIMEFORMATCU)
log.log_filename = output_dir + '/log_' + domain_name
checkpoint_dir = 'checkpoints'
#yieldstep = 50
#finaltime = 200
yieldstep = 1000
finaltime = 3000
#yieldstep = 10800
#finaltime = 43200.0 # 12 hours
#finaltime = 86400.0 # only for running 21_07_2021 simulation
#finaltime = 3600.0
min_allowed_height = 0.008
max_allowed_speed = 1.0
checkpoint_time = 60 * 60
flow_algorithm = 'DE1'
useCheckpointing = False

# File paths (some used later for hotstart)
friction_filename = 'friction_data/LULC_Apr2020_Norm.asc'
friction_location = "centroids"
#elevation_filename = f'elevation_data/Elev_{name_stem}_LIDAR1m_ALOS30m.csv'
elevation_filename = 'elevation_data/Elev_delta11372sqkm_uniform_mesh_300sqm_Chilka_300sqm_LIDAR1m_ALOS30m.csv'
elevation_location = "vertices"
previous_hotstart_filename = 'output/' + Previous_Date.strftime(TIMEFORMATCU) + '/hotstart_' + Previous_hotstart_Date.strftime(TIMEFORMATCU) + '.sww'
previous_output_filename = 'output/' + Previous_Date.strftime(TIMEFORMATCU) + '/mahanadi_delta_' + Previous_Date.strftime(TIMEFORMATCU) + '.sww'
tobe_hotstart_filename = 'output/' + Current_Date.strftime(TIMEFORMATCU) + '/hotstart_' + Previous_Date.strftime(TIMEFORMATCU) + '.sww'
sww_filename = output_dir + '/' + domain_name + '.sww'
tobe_currentday_filename = 'output/' + Current_Date.strftime(TIMEFORMATCU) + '/current_day_' + Current_Date.strftime(TIMEFORMATCU) + '_' + Current_Date.strftime(TIMEFORMATCU) + '.sww'
tide_filename = 'tide_data/output/paradip_tide_' + Current_Date.strftime(TIMEFORMATCU) + '.tms'

number_of_inlets = 2
coastal_tag_start = 289
coastal_tag_end = 303
line_naraj_barrage = [[374437.043804152, 2263860.85985414], [374665.046154188, 2264351.06490672]]
discharge_filename_naraj_barrage = 'discharge_data/output/naraj_barrage_discharge_' + Current_Date.strftime(TIMEFORMATCU) + '.tms'
line_mahanadi_barrage = [[387363.196978616, 2263420.51264919], [388385.893737998, 2265788.45575452]]
discharge_filename_mahanadi_barrage = 'discharge_data/output/mahanadi_barrage_discharge_' + Current_Date.strftime(TIMEFORMATCU) + '.tms'

imd_rainfall_factor_pt_bhub = 0.00000001157407
imd_rainfall_factor_rgdata_rain25 = 0.00000001157407
imd_rainfall_factor_wrf = 0.00000009259259
imd_rainfall_factor_gfs = 0.00000009259259
gpm_rainfall_factor = 0.00000000925925926
gfs_rainfall_factor = 0.00000009259259259
gauge_filename = 'gauage_locations_data/gauges.csv'
data_path = '/scratch/samir/15-sep/mahanadi-delta_1'
output_frequency = 432000.0
intermediate_outputs = np.arange((yieldstep * output_frequency), finaltime, (yieldstep * output_frequency))
if 86400.0 in intermediate_outputs:
    np.delete(intermediate_outputs, np.where(86400.0))
daily_time = np.arange(86400.0, finaltime, 86400.0)

t1 = time.time()
if myid == 0:
    print('Loading of data: Time', t1 - t0, flush=True)

# --------------------------------------------------------------------------
# Setup procedures
# --------------------------------------------------------------------------
def check_for_dir(release_dir):
    if os.path.isdir(release_dir):
        pass
    else:
        print(("Creating directory " + str(release_dir)))
        os.mkdir(release_dir)
    pass

# --------------------------------------------------------------------------
# LOAD PRE-PARTITIONED DOMAIN
# This replaces the old rank 0 domain creation + distribute() pattern
# --------------------------------------------------------------------------
if myid == 0:
    print('Starting domain loading from partitions', flush=True)
    print(f'Partition name: {partition_name}', flush=True)
    print(f'Partition dir: {partition_dir}', flush=True)

# Setup discharge and tide on rank 0
if myid == 0:
    tsub0 = time.time()
    if len(sys.argv) > 1:
        setup_discharge_3.setup_discharge_old(sys.argv[1])
        setup_tide_3.setup_tide_old(sys.argv[1])
    else:
        setup_discharge_3.setup_discharge()
        setup_tide_3.setup_tide()
    check_for_dir(output_dir)
    tsub1 = time.time()
    print('Setup discharge/tide: Time', tsub1 - tsub0, flush=True)

barrier()

# ALL RANKS load their partition in parallel
t_load_start = time.time()
domain = anuga.sequential_distribute_load(partition_name, partition_dir=partition_dir, verbose=verbose)
t_load = time.time() - t_load_start

if myid == 0:
    print(f'Loaded partitioned domain in {t_load:.2f}s', flush=True)

# Set domain parameters (needed after loading)
domain.set_name(domain_name)
domain.set_datadir(output_dir)

# Handle hotstart stage if previous output exists
# NOTE: This modifies stage AFTER loading partition
if myid == 0:
    tsub0 = time.time()

if os.path.exists(previous_hotstart_filename):
    if myid == 0:
        print(f'Loading hotstart from {previous_hotstart_filename}')
    pc = anuga.plot_utils.get_centroids(previous_hotstart_filename, timeSlices='last')
    domain.set_quantity('stage', numeric=pc.stage[0], location='centroids')
elif os.path.exists(previous_output_filename):
    if myid == 0:
        print(f'Loading previous output from {previous_output_filename}')
    pc = anuga.plot_utils.get_centroids(previous_output_filename)
    domain.set_quantity('stage', numeric=pc.stage[8], location='centroids')
else:
    # Stage already set in partition (elevation + Chilka polygon)
    if myid == 0:
        print('Using initial stage from partition')

if myid == 0:
    tsub1 = time.time()
    print('Stage setting: Time', tsub1 - tsub0, flush=True)

barrier()

t1 = time.time()
if myid == 0:
    print('Domain loading complete: Time', t1 - t0)
    domain.print_statistics()

barrier()

# ------------------------------------------------------------------------------
# Setup boundary conditions
# This must happen *after* domain has been loaded
# ------------------------------------------------------------------------------
Br = anuga.Reflective_boundary(domain)
Bw = anuga.Reflective_boundary(domain)
if os.path.exists(tide_filename):
    wave_function = anuga.file_function(tide_filename, quantities='stage', verbose=verbose)
    #Bw = anuga.Time_boundary(domain=domain, function=lambda t: [float(wave_function(t)), 0.0, 0.0])
    Bw = anuga.Time_boundary(domain=domain, function=lambda t: [wave_function(t).item(), 0.0, 0.0])
else:
    if myid == 0 and verbose:
        print("The paradip tide file %s does not exist !!" % tide_filename)

dict1 = {str(n): Bw for n in range(coastal_tag_start, coastal_tag_end)}
dict2 = {'exterior': Br}
dict3 = dict(itertools.chain(list(dict2.items()), list(dict1.items())))
domain.set_boundary(dict3)


domain.fixed_flux_timestep = 0.21853963
domain.set_multiprocessor_mode(1)
domain.use_c_rk2_loop = True
barrier()

# Input Discharge with Time Series Naraj Barrage
if os.path.exists(discharge_filename_naraj_barrage):
    Q0 = file_function(discharge_filename_naraj_barrage, quantities='flow')
    inlet0 = Inlet_operator(domain, line_naraj_barrage, Q0, logging=True, description='Naraj Barrage Discharge', verbose=False)
else:
    if myid == 0 and verbose:
        print("The Naraj Barrage Discharge file %s does not exist !!" % discharge_filename_naraj_barrage)

# Input Discharge with Time Series Mahanadi Barrage
if os.path.exists(discharge_filename_mahanadi_barrage):
    Q1 = file_function(discharge_filename_mahanadi_barrage, quantities='flow')
    inlet1 = Inlet_operator(domain, line_mahanadi_barrage, Q1, logging=True, description='Mahanadi Barrage Discharge', verbose=False)
else:
    if myid == 0 and verbose:
        print("The Mahanadi Barrage Discharge file %s does not exist !!" % discharge_filename_mahanadi_barrage)

# Initialize Rainfall parameter
Q = Quantity(domain, name='Rain', register=True)
rain_opertor = anuga.Rate_operator(domain, rate=Q, default_rate=0.0)

# Turn on checkpointing
if useCheckpointing:
    domain.set_checkpointing(checkpoint_time=checkpoint_time, checkpoint_dir=checkpoint_dir)

# ------------------------------------------------------------------------------
# Evolution (unchanged from original)
# ------------------------------------------------------------------------------
if myid == 0 and verbose:
    print('EVOLVE')

barrier()
#domain.gpu.flop_counters_reset()
#domain.gpu.flop_counters_start_timer()
import numpy as np
rank = domain.processor
n = domain.number_of_elements
stage = domain.quantities['stage'].centroid_values
bed = domain.quantities['elevation'].centroid_values
depth = stage - bed
n_wet = np.sum(depth > domain.minimum_allowed_height)
print(f"[Rank {rank}] elements={n}, wet={n_wet} ({100*n_wet/n:.1f}%)")
import time

t0 = time.time()
import cProfile
import pstats
profiler = cProfile.Profile()
profiler.enable()
rain_set_zero = True

for t in domain.evolve(yieldstep=yieldstep, finaltime=finaltime):
    if myid == 0: domain.write_time()
    fltStr = str(t)
    rplStr = fltStr.replace(".", "_")
    imd_daily_rain25 = "rainfall_data/imd/daily/rgdata_rain_25/" + Current_Date.strftime(
        TIMEFORMATCU) + "/imd_" + Current_Date.strftime(TIMEFORMATCU) + "_%s.csv" % rplStr
    imd_daily_daily_pt_bhub = "rainfall_data/imd/daily/pt_data_bhubaneshwar/" + Current_Date.strftime(
        TIMEFORMATCU) + "/imd_" + Current_Date.strftime(TIMEFORMATCU) + "_%s.csv" % rplStr
    imd_rain_file_wrf = "rainfall_data/imd/wrf/" + Current_Date.strftime(
        TIMEFORMATCU) + "/imd_" + Current_Date.strftime(
        TIMEFORMATCU) + "_%s.csv" % rplStr
    imd_rain_file_gfs = "rainfall_data/imd/gfs/" + Current_Date.strftime(
        TIMEFORMATCU) + "/imd_" + Current_Date.strftime(
        TIMEFORMATCU) + "_%s.csv" % rplStr
    gpm_rain_file = "rainfall_data/gpm/" + Current_Date.strftime(TIMEFORMATCU) + "/gpm_" + Current_Date.strftime(
        TIMEFORMATCU) + "_%s.csv" % rplStr
    gfs_rain_file = "rainfall_data/gfs/" + Current_Date.strftime(TIMEFORMATCU) + "/gfs_" + Current_Date.strftime(
        TIMEFORMATCU) + "_%s.csv" % rplStr
    if not rain_set_zero and len(np.where(daily_time == t)[0]) == 1:
        rain_set_zero = True;
    if rain_set_zero:
        if os.path.exists(imd_daily_daily_pt_bhub):
            if myid == 0: print(
                    "Setting up imd daily rainfall file %s !!" % imd_daily_daily_pt_bhub)
            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            for tri_index in range(len(domain)):
                vertices = domain.get_triangles(tri_index)
                triangle_index_rain.append(tri_index)
                triangle_index_elevation = []
                triangle_index_elevation.insert(0, np.double(
                    linecache.getline(imd_daily_daily_pt_bhub, vertices[0] + 1)) * imd_rainfall_factor_pt_bhub)
                triangle_index_elevation.insert(1, np.double(
                    linecache.getline(imd_daily_daily_pt_bhub, vertices[1] + 1)) * imd_rainfall_factor_pt_bhub)
                triangle_index_elevation.insert(2, np.double(
                    linecache.getline(imd_daily_daily_pt_bhub, vertices[2] + 1)) * imd_rainfall_factor_pt_bhub)
                triangle_index_elvation_main_rain.append(triangle_index_elevation)
            domain.set_quantity('Rain',
                                numeric=triangle_index_elvation_main_rain,
                                use_cache=cache,
                                verbose=True,
                                alpha=0.1, indices=triangle_index_rain,
                                location='vertices')

            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            linecache.clearcache()
            rain_opertor.set_rate(rate=Q)
            rain_set_zero = False
        elif os.path.exists(imd_daily_rain25):
            if myid == 0: print(
                    "Setting up imd daily 25 rainfall file %s !!" % imd_daily_rain25)
            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            for tri_index in range(len(domain)):
                vertices = domain.get_triangles(tri_index)
                triangle_index_rain.append(tri_index)
                triangle_index_elevation = []
                triangle_index_elevation.insert(0, np.double(
                    linecache.getline(imd_daily_rain25, vertices[0] + 3)) * imd_rainfall_factor_rgdata_rain25)
                triangle_index_elevation.insert(1, np.double(
                    linecache.getline(imd_daily_rain25, vertices[1] + 3)) * imd_rainfall_factor_rgdata_rain25)
                triangle_index_elevation.insert(2, np.double(
                    linecache.getline(imd_daily_rain25, vertices[2] + 3)) * imd_rainfall_factor_rgdata_rain25)
                triangle_index_elvation_main_rain.append(triangle_index_elevation)
            domain.set_quantity('Rain',
                                numeric=triangle_index_elvation_main_rain,
                                use_cache=cache,
                                verbose=True,
                                alpha=0.1, indices=triangle_index_rain,
                                location='vertices')

            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            linecache.clearcache()
            rain_opertor.set_rate(rate=Q)
            rain_set_zero = False
        elif os.path.exists(imd_rain_file_wrf):
            if myid == 0: print(
                    "Setting up imd wrf rainfall file %s !!" % imd_rain_file_wrf)
            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            for tri_index in range(len(domain)):
                vertices = domain.get_triangles(tri_index)
                triangle_index_rain.append(tri_index)
                triangle_index_elevation = []
                triangle_index_elevation.insert(0, np.double(
                    linecache.getline(imd_rain_file_wrf, vertices[0] + 1)) * imd_rainfall_factor_wrf)
                triangle_index_elevation.insert(1, np.double(
                    linecache.getline(imd_rain_file_wrf, vertices[1] + 1)) * imd_rainfall_factor_wrf)
                triangle_index_elevation.insert(2, np.double(
                    linecache.getline(imd_rain_file_wrf, vertices[2] + 1)) * imd_rainfall_factor_wrf)
                triangle_index_elvation_main_rain.append(triangle_index_elevation)
            domain.set_quantity('Rain',
                                numeric=triangle_index_elvation_main_rain,
                                use_cache=cache,
                                verbose=True,
                                alpha=0.1, indices=triangle_index_rain,
                                location='vertices')

            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            linecache.clearcache()
            rain_opertor.set_rate(rate=Q)
        elif os.path.exists(imd_rain_file_gfs):
            if myid == 0: print(
                    "Setting up imd gfs rainfall file %s !!" % imd_rain_file_gfs)
            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            for tri_index in range(len(domain)):
                vertices = domain.get_triangles(tri_index)
                triangle_index_rain.append(tri_index)
                triangle_index_elevation = []
                triangle_index_elevation.insert(0, np.double(
                    linecache.getline(imd_rain_file_gfs, vertices[0] + 1)) * imd_rainfall_factor_gfs)
                triangle_index_elevation.insert(1, np.double(
                    linecache.getline(imd_rain_file_gfs, vertices[1] + 1)) * imd_rainfall_factor_gfs)
                triangle_index_elevation.insert(2, np.double(
                    linecache.getline(imd_rain_file_gfs, vertices[2] + 1)) * imd_rainfall_factor_gfs)
                triangle_index_elvation_main_rain.append(triangle_index_elevation)
            domain.set_quantity('Rain',
                                numeric=triangle_index_elvation_main_rain,
                                use_cache=cache,
                                verbose=True,
                                alpha=0.1, indices=triangle_index_rain,
                                location='vertices')

            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            linecache.clearcache()
            rain_opertor.set_rate(rate=Q)
        elif os.path.exists(gpm_rain_file):
            if myid == 0: print(
                    "Setting up NOAA GPM rainfall file %s !!" % gpm_rain_file)
            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            for tri_index in range(len(domain)):
                vertices = domain.get_triangles(tri_index)
                triangle_index_rain.append(tri_index)
                triangle_index_elevation = []
                triangle_index_elevation.insert(0, np.double(
                    linecache.getline(gpm_rain_file, vertices[0] + 1)) * gpm_rainfall_factor)
                triangle_index_elevation.insert(1, np.double(
                    linecache.getline(gpm_rain_file, vertices[1] + 1)) * gpm_rainfall_factor)
                triangle_index_elevation.insert(2, np.double(
                    linecache.getline(gpm_rain_file, vertices[2] + 1)) * gpm_rainfall_factor)
                triangle_index_elvation_main_rain.append(triangle_index_elevation)
            domain.set_quantity('Rain',
                                numeric=triangle_index_elvation_main_rain,
                                use_cache=cache,
                                verbose=True,
                                alpha=0.1, indices=triangle_index_rain,
                                location='vertices')
            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            linecache.clearcache()
            rain_opertor.set_rate(rate=Q)
        elif os.path.exists(gfs_rain_file):
            if myid == 0: print(
                    "Setting up NOAA GFS rainfall file %s !!" % gfs_rain_file)
            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            for tri_index in range(len(domain)):
                vertices = domain.get_triangles(tri_index)
                triangle_index_rain.append(tri_index)
                triangle_index_elevation = []
                triangle_index_elevation.insert(0, np.double(
                    linecache.getline(gfs_rain_file, vertices[0] + 1)) * gfs_rainfall_factor)
                triangle_index_elevation.insert(1, np.double(
                    linecache.getline(gfs_rain_file, vertices[1] + 1)) * gfs_rainfall_factor)
                triangle_index_elevation.insert(2, np.double(
                    linecache.getline(gfs_rain_file, vertices[2] + 1)) * gfs_rainfall_factor)
                triangle_index_elvation_main_rain.append(triangle_index_elevation)
            domain.set_quantity('Rain',
                                numeric=triangle_index_elvation_main_rain,
                                use_cache=cache,
                                verbose=True,
                                alpha=0.1, indices=triangle_index_rain,
                                location='vertices')

            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            linecache.clearcache()
            rain_opertor.set_rate(rate=Q)
        elif rain_set_zero:
            #print("The Rainfall IMD/GPM/GFS files does not exist setting rain to 0.004 !!")
            domain.set_quantity('Rain', 0.000) #set rainfall with hardcoded value = 4mm #modifed by RK
            rain_opertor.set_rate(rate=Q)
    else:
        if myid == 0: print("Using previously set Daily Rainfall!!")
    import re
    volume = domain.compute_total_volume()
    stats = domain.timestepping_statistics()
    rainstats = rain_opertor.timestepping_statistics()
    maxInundation = Q.get_maximum_value()
    indices = domain.get_wet_elements()
    element_count = len(indices)
    file1 = open("rain_data.txt", "a+")
    file2 = open("wet_elements.txt", "a+")
    file3 = open("max_inandation.txt", "a+")
    try:
      rain_arr = re.findall(r"\d+\.\d+",rainstats)
      total_rain = float(rain_arr[0])
      file1.write(str(total_rain))
      file1.writelines("\n")
      file2.write(str(element_count))
      file2.writelines("\n")
      file3.write(str(maxInundation))
      file3.writelines("\n")
    except:
      #print("there is no Q")
      file1.write(str(0.0))
      file1.writelines("\n")
      file2.write(str(0.0))
      file2.writelines("\n")
      file3.write(str(0.0))
      file3.writelines("\n")

    if myid == 0:
        print("Total Volume =: %s Time Stepping Statistics =: %s Rain Operator Statistics=: %s" % (
            str(volume), str(stats),str(rainstats )))

    # Close file handles
    file1.close()
    file2.close()
    file3.close()

profiler.disable()
barrier()

#domain.gpu.flop_counters_stop_timer()
#domain.gpu.flop_counters_print_global()

if myid == 0:
    print("\n" + "="*80)
    print("PROFILING RESULTS - Top 40 by cumulative time")
    print("="*80)
    stats = pstats.Stats(profiler)
    stats.sort_stats('cumulative')
    stats.print_stats(40)

   # Also save to file for detailed analysis
    profiler.dump_stats('profile.prof')
    print(f"\nProfile saved to profile.prof - view with: python -m pstats profile.prof")

for p in range(numprocs):
    if myid == p:
        print('Processor %g ' % myid)
        print('That took %.2f seconds' % (time.time() - t0))
        print('Communication time %.2f seconds' % domain.communication_time)
        print('Reduction Communication time %.2f seconds' % domain.communication_reduce_time)
        print('Broadcast time %.2f seconds' % domain.communication_broadcast_time)
    else:
        pass
    barrier()

# Merge the individual sww files
#domain.sww_merge(verbose=verbose, delete_old=True)
finalize()
