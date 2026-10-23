#!/usr/bin/env python3
"""Script for running a may prediction on Mahanadi Delta regions

Source data such as elevation and boundary data is assumed to be available in
directories specified by project.py
The output sww file is stored in directory named with date after the output

Author Girishchandra Y.
"""
# ------------------------------------------------------------------------------
# Import necessary modules
# ------------------------------------------------------------------------------
# Standard modules
import os
import time
import sys
import re
import numpy as np
import linecache
import itertools
import datetime
# -----------------------------------------------------------------------------
# Import Own Data Preparation Tide and discharge Modules
# -----------------------------------------------------------------------------
from discharge_data import setup_discharge_3
from tide_data import setup_tide_3

# Related major packages
import anuga
# ------------------------
# ANUGA Modules
# ------------------------
from anuga import create_domain_from_file
from anuga import Inlet_operator  # inflow line
from anuga import file_function  # for tms file reading
from anuga import Rate_operator
from anuga import Quantity
from anuga import distribute, myid, numprocs, finalize, barrier  # for parallel
from anuga.utilities import log as log  # for logging
from anuga.abstract_2d_finite_volumes import quantity


t0 = time.time()
# --------------------------------------------------------------------------
# Setup parameters
# --------------------------------------------------------------------------
# Get key data for the simulation
if len(sys.argv) > 1:
    # Config file
    input = str(sys.argv[1]).split("-")
    Current_Date = datetime.datetime(int(input[2]), int(input[1]), int(input[0])).replace(minute=0, hour=0, second=0,
                                                                                          microsecond=0)
    Previous_Date = datetime.datetime(int(input[2]), int(input[1]), int(input[0])).replace(minute=0, hour=0, second=0,
                                                                                           microsecond=0) - datetime.timedelta(
        days=1)
    Previous_hotstart_Date = datetime.datetime(int(input[2]), int(input[1]), int(input[0])).replace(minute=0, hour=0,
                                                                                                    second=0,
                                                                                                    microsecond=0) - datetime.timedelta(
        days=2)

else:
    Current_Date = datetime.datetime.today().replace(minute=0, hour=0, second=0, microsecond=0)
    Previous_Date = datetime.datetime.today().replace(minute=0, hour=0, second=0, microsecond=0) - datetime.timedelta(
        days=1)
    Previous_hotstart_Date = datetime.datetime.today().replace(minute=0, hour=0, second=0,
                                                               microsecond=0) - datetime.timedelta(days=2)
TIMEFORMATCU = "%d_%m_%Y"
name_stem = 'delta11372sqkm_uniform_mesh_300sqm_Chilka300sqm'  # name of the ascii file
#name_stem = 'Delta_11381_sqkm_900sqm'
mesh_filename = 'mesh_file/' + name_stem + '.tsh'
poly = anuga.read_polygon('polygons/Chilka_Modified_Poly.csv')
poly_bounding = anuga.read_polygon('polygons/Delta_11372_sqkm.csv')
chilka_stage = 0.15
cache = True  # default
verbose = True  # status print on cmd prompt
intermediate_result_extraction = True
domain_name = 'mahanadi_delta_' + Current_Date.strftime(TIMEFORMATCU)
output_dir = 'output/' + Current_Date.strftime(TIMEFORMATCU)
log.log_filename = output_dir + '/log_' + domain_name
checkpoint_dir = 'checkpoints'
#yieldstep = 500
#finaltime = 1500
#finaltime = 600
# for it to load the rain data, this is the yieldstep
#yieldstep=10
#finaltime=20
yieldstep = 10800
#finaltime = 172800.0
#finaltime = 21600.0
#finaltime = 432000.0 # 432000.0 for 5 days simulation
#finaltime = 302400.0 #for 4 days simulation | 302400.0 for output at 3rd day 6 PM
#finaltime = 21600.0
finaltime = 86400.0 # only for running 21_07_2021 simulation
min_allowed_height = 0.008
max_allowed_speed = 1.0
checkpoint_time = 60 * 60
flow_algorithm = 'DE1'
useCheckpointing = False
friction_filename = 'friction_data/LULC_Apr2020_Norm.asc'
friction_location = "centroids"
elevation_filename = 'elevation_data/Elev_delta11372sqkm_uniform_mesh_300sqm_Chilka_300sqm_LIDAR1m_ALOS30m.csv'
elevation_location = "vertices"
previous_hotstart_filename = 'output/' + Previous_Date.strftime(
    TIMEFORMATCU) + '/hotstart_' + Previous_hotstart_Date.strftime(TIMEFORMATCU) + '.sww'
previous_output_filename = 'output/' + Previous_Date.strftime(
    TIMEFORMATCU) + '/mahanadi_delta_' + Previous_Date.strftime(TIMEFORMATCU) + '.sww'
tobe_hotstart_filename = 'output/' + Current_Date.strftime(TIMEFORMATCU) + '/hotstart_' + Previous_Date.strftime(
    TIMEFORMATCU) + '.sww'
sww_filename = output_dir + '/' + domain_name + '.sww'
tobe_currentday_filename = 'output/' + Current_Date.strftime(TIMEFORMATCU) + '/current_day_' + Current_Date.strftime(
    TIMEFORMATCU) + '_' + Current_Date.strftime(TIMEFORMATCU) + '.sww'
tide_filename = 'tide_data/output/paradip_tide_' + Current_Date.strftime(TIMEFORMATCU) + '.tms'
number_of_inlets = 2
coastal_tag_start = 289
coastal_tag_end = 303
line_naraj_barrage = [[374437.043804152,2263860.85985414],
                      [374665.046154188, 2264351.06490672]]  # inlet location Naraj Barrage location updated 21 Feb 2022 [[374227.523182097, 2264003.63492516], [374326.724682891, 2264397.66032303]] ORG LOC
discharge_filename_naraj_barrage = 'discharge_data/output/naraj_barrage_discharge_' + Current_Date.strftime(
    TIMEFORMATCU) + '.tms'
line_mahanadi_barrage = [[387363.196978616, 2263420.51264919],
                         [388385.893737998, 2265788.45575452]]  # inlet location Mahanadi Barrage [[387363.196978616, 2263420.51264919], [388385.893737998, 2265788.45575452]] ORG LOC
discharge_filename_mahanadi_barrage = 'discharge_data/output/mahanadi_barrage_discharge_' + Current_Date.strftime(
    TIMEFORMATCU) + '.tms'
imd_rainfall_factor_pt_bhub = 0.00000001157407
imd_rainfall_factor_rgdata_rain25 = 0.00000001157407
imd_rainfall_factor_wrf = 0.00000009259259
imd_rainfall_factor_gfs = 0.00000009259259
gpm_rainfall_factor = 0.00000000925925926
gfs_rainfall_factor = 0.00000009259259259
gauge_filename = 'gauage_locations_data/gauges.csv'
data_path = '/scratch/samir/15-sep/mahanadi-delta_1'
daily_time = np.arange(86400.0, finaltime, 86400.0)

t1 = time.time()
if myid == 0 :
    print ('Loading of data: Time',t1-t0, flush=True)

# --------------------------------------------------------------------------
# Setup procedures
# --------------------------------------------------------------------------
# a function to check for and create our release_dir if necessary
def check_for_dir(release_dir):
    if os.path.isdir(release_dir):
        pass
    else:
        print(("Creating directory " + str(release_dir)))
        os.mkdir(release_dir)
    pass


# ------------------------------------------------------------------------------
# Use a try statement to read in previous checkpoint file and if not possible
# just go ahead as normal and produce domain as usual.
#
# Though in the part of the code where you create the domain as normal,
# remember to turn on checkpointing via
# domain.set_checkpointing(checkpoint_time = 5) (see code below)
# ------------------------------------------------------------------------------

if myid == 0 :
    print ('Starting domain creation', flush=True)
if myid == 0:
    tsub0 = time.time()
    if len(sys.argv) > 1:
        setup_discharge_3.setup_discharge_old(sys.argv[1])
        setup_tide_3.setup_tide_old(sys.argv[1])
    else:
        setup_discharge_3.setup_discharge()
        setup_tide_3.setup_tide()
    check_for_dir(output_dir)
    domain = create_domain_from_file(mesh_filename)
    domain.set_name(domain_name)
    domain.set_flow_algorithm(flow_algorithm)
    domain.set_CFL(2.0)
    domain.set_datadir(output_dir)
    domain.set_minimum_allowed_height(min_allowed_height)
    domain.set_maximum_allowed_speed(max_allowed_speed)
    #print domain.statistics()
    domain.set_quantity('friction', filename=friction_filename, location=friction_location)
    triangle_index = []
    triangle_index_elvation_main = []
    tsub1 = time.time()
    print ('Domain create functioanlit: Time',tsub1-tsub0, flush=True)

    tsub0 = time.time()
    for tri_index in range(len(domain)):
        vertices = domain.get_triangles(tri_index)
        triangle_index.append(tri_index)
        triangle_index_elvation = []
        triangle_index_elvation.insert(0, np.double(linecache.getline(elevation_filename, vertices[0] + 1)))
        triangle_index_elvation.insert(1, np.double(linecache.getline(elevation_filename, vertices[1] + 1)))
        triangle_index_elvation.insert(2, np.double(linecache.getline(elevation_filename, vertices[2] + 1)))
        triangle_index_elvation_main.append(triangle_index_elvation)
    tsub1 = time.time()
    print ('For tri index loop : Time',tsub1-tsub0, flush=True)
    tsub0 = time.time()
    domain.set_quantity('elevation',
                        numeric=triangle_index_elvation_main,  # triangle_elevation
                        use_cache=cache,
                        verbose=True,
                        alpha=0.1, indices=triangle_index,
                        location=elevation_location)  # alpha refer user manual page 98 | removed ,location='centroids'
    tsub1 = time.time()
    print ('Set quantity: Time',tsub1-tsub0, flush=True)
    triangle_index = []
    triangle_index_elvation_main = []
    linecache.clearcache()
    tsub0 = time.time()
    if os.path.exists(previous_hotstart_filename):
        pc = anuga.plot_utils.get_centroids(previous_hotstart_filename, timeSlices='last')
        domain.set_quantity('stage', numeric=pc.stage[0], location='centroids')
    elif os.path.exists(previous_output_filename):
        pc = anuga.plot_utils.get_centroids(previous_output_filename)
        domain.set_quantity('stage', numeric=pc.stage[8], location='centroids') #change to 8 for regular hotstart
    else:
        domain.set_quantity('stage', expression='elevation')  # initial dry condition
        domain.add_quantity('stage', numeric=chilka_stage, polygon=poly)  # initial stage inside poly
    tsub1 = time.time()
    print ('Stage setting: Time',tsub1-tsub0, flush=True)
    domain.print_statistics()
else:
    domain = None

t1 = time.time()
if myid == 0 :
    print ('Creation of sequential domain: Time',t1-t0)

barrier()
    # --------------------------------------------------------------------------
    # Distribute sequential domain on processor 0 to other processors
    # --------------------------------------------------------------------------

if myid == 0 and verbose: print('DISTRIBUTING DOMAIN')
domain = distribute(domain, verbose=verbose)

t2 = time.time()

if myid == 0 :
    print ('Domain Distribution  Time ',t2-t1)
barrier()

#To make sure about domain distribution happening well. Samir

    # ------------------------------------------------------------------------------
    # Setup boundary conditions
    # This must currently happen *after* domain has been distributed
    # ------------------------------------------------------------------------------
Br = anuga.Reflective_boundary(domain)  # Reflective wall
Bw = anuga.Reflective_boundary(domain)
if os.path.exists(tide_filename):
    wave_function = anuga.file_function(tide_filename, quantities='stage', verbose=verbose)
    #Bw = anuga.Time_boundary(domain=domain, function=lambda t: [wave_function(t) * 1, 0.0, 0.0]) #original
    #Bw = anuga.Time_boundary(domain=domain,function=lambda t: [float(wave_function(t)), 0.0, 0.0])
    Bw = anuga.Time_boundary(domain=domain, function=lambda t: [wave_function(t).item(), 0.0, 0.0])
else:
    if myid == 0 and verbose: print("The paradip tide file %s does not exist !!" % tide_filename)


dict1 = {str(n): Bw for n in range(coastal_tag_start, coastal_tag_end)}
dict2 = {'exterior': Br}  # specify boundary condition for all other tags
dict3 = dict(itertools.chain(list(dict2.items()), list(dict1.items())))
domain.set_boundary(dict3)  # set boundary

#To make sure that boundry conditions managed well Samir
barrier()

#domain.fixed_flux_timestep = 0.46536518
# 900 sqm step size -> 0.21853963
domain.set_multiprocessor_mode(2)
domain.use_c_rk_loop = True

    ########################################################################
    # Input Discharge with Time Series Naraj Barrage
    ########################################################################
if os.path.exists(discharge_filename_naraj_barrage):
    Q0 = file_function(discharge_filename_naraj_barrage,
                      quantities='flow')  # quantities = 'flow' is the field name of discharge in TMS file so
    inlet0 = Inlet_operator(domain, line_naraj_barrage, Q0, logging=True, description='Naraj Barrage Discharge',
                                verbose=False)
else:
    if myid == 0 and verbose: print(
            "The Naraj Barrage Discharge file %s does not exist !!" % discharge_filename_naraj_barrage)
    ########################################################################
    # Input Discharge with Time Series Mahanadi Barrage
    ########################################################################
    # Q1 = 200 #constant discharge value to be specified here
if os.path.exists(discharge_filename_mahanadi_barrage):
    Q1 = file_function(discharge_filename_mahanadi_barrage,
                        quantities='flow')  # quantities = 'flow' is the field name of discharge in TMS file so disable_linear_interpolation=True
    inlet1 = Inlet_operator(domain, line_mahanadi_barrage, Q1, logging=True,
                            description='Mahanadi Barrage Discharge', verbose=False)
else:
    if myid == 0 and verbose: print(
            "The Mahanadi Barrage Discharge file %s does not exist !!" % discharge_filename_mahanadi_barrage)
    ########################################################################
    # Initialize Rainfall parameter
    ########################################################################
Q = Quantity(domain, name='Rain', register=True)
rain_operator = anuga.Rate_operator(domain, rate=Q,
                                   default_rate=0.0)  # rate factor 0.00009 for 3 hr rain & 0.1 for GPM factor & mm to m conv factor
    # -----------------------------------------------------------------------------
    # Turn on checkpointing
    # -----------------------------------------------------------------------------
if useCheckpointing:
    #print("Using Checkpointing - Samir run_model_3.py:258")
    domain.set_checkpointing(checkpoint_time=checkpoint_time, checkpoint_dir=checkpoint_dir)
# ------------------------------------------------------------------------------
# Evolution
# ------------------------------------------------------------------------------
if myid == 0 and verbose: print('EVOLVE')

barrier()

rank = domain.processor
n = domain.number_of_elements
stage = domain.quantities['stage'].centroid_values
bed = domain.quantities['elevation'].centroid_values
depth = stage - bed
n_wet = np.sum(depth > domain.minimum_allowed_height)
print(f"[Rank {rank}] elements={n}, wet={n_wet} ({100*n_wet/n:.1f}%)")

t0 = time.time()
import cProfile
import pstats
#profiler = cProfile.Profile()
#profiler.enable()
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
                global_vertices = domain.node_l2g[vertices]
                triangle_index_rain.append(tri_index)
                triangle_index_elevation = []
                triangle_index_elevation.insert(0, np.double(
                    linecache.getline(imd_daily_daily_pt_bhub, global_vertices[0] + 1)) * imd_rainfall_factor_pt_bhub)
                triangle_index_elevation.insert(1, np.double(
                    linecache.getline(imd_daily_daily_pt_bhub, global_vertices[1] + 1)) * imd_rainfall_factor_pt_bhub)
                triangle_index_elevation.insert(2, np.double(
                    linecache.getline(imd_daily_daily_pt_bhub, global_vertices[2] + 1)) * imd_rainfall_factor_pt_bhub)
                triangle_index_elvation_main_rain.append(triangle_index_elevation)  # print triangle_index_elevation
            domain.set_quantity('Rain',
                                numeric=triangle_index_elvation_main_rain,  # triangle_elevation
                                use_cache=cache,
                                verbose=True,
                                alpha=0.1, indices=triangle_index_rain,
                                location='vertices')  # alpha refer user manual page 98 | removed ,location='cen$

            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            linecache.clearcache()
            rain_operator.set_rate(rate=Q)
            rain_set_zero = False
        elif os.path.exists(imd_daily_rain25):
            if myid == 0: print(
                    "Setting up imd daily 25 rainfall file %s !!" % imd_daily_rain25)
            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            for tri_index in range(len(domain)):
                vertices = domain.get_triangles(tri_index)
                global_vertices = domain.node_l2g[vertices]
                triangle_index_rain.append(tri_index)
                triangle_index_elevation = []
                triangle_index_elevation.insert(0, np.double(
                    linecache.getline(imd_daily_rain25, global_vertices[0] + 3)) * imd_rainfall_factor_rgdata_rain25)
                triangle_index_elevation.insert(1, np.double(
                    linecache.getline(imd_daily_rain25, global_vertices[1] + 3)) * imd_rainfall_factor_rgdata_rain25)
                triangle_index_elevation.insert(2, np.double(
                    linecache.getline(imd_daily_rain25, global_vertices[2] + 3)) * imd_rainfall_factor_rgdata_rain25)
                triangle_index_elvation_main_rain.append(triangle_index_elevation)  # print triangle_index_elevation
            domain.set_quantity('Rain',
                                numeric=triangle_index_elvation_main_rain,  # triangle_elevation
                                use_cache=cache,
                                verbose=True,
                                alpha=0.1, indices=triangle_index_rain,
                                location='vertices')  # alpha refer user manual page 98 | removed ,location='cen$

            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            linecache.clearcache()
            rain_operator.set_rate(rate=Q)
            rain_set_zero = False
        elif os.path.exists(imd_rain_file_wrf):
            if myid == 0: print(
                    "Setting up imd wrf rainfall file %s !!" % imd_rain_file_wrf)
            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            for tri_index in range(len(domain)):
                vertices = domain.get_triangles(tri_index)
                global_vertices = domain.node_l2g[vertices]
                triangle_index_rain.append(tri_index)
                triangle_index_elevation = []
                triangle_index_elevation.insert(0, np.double(
                    linecache.getline(imd_rain_file_wrf, global_vertices[0] + 1)) * imd_rainfall_factor_wrf)
                triangle_index_elevation.insert(1, np.double(
                    linecache.getline(imd_rain_file_wrf, global_vertices[1] + 1)) * imd_rainfall_factor_wrf)
                triangle_index_elevation.insert(2, np.double(
                    linecache.getline(imd_rain_file_wrf, global_vertices[2] + 1)) * imd_rainfall_factor_wrf)
                triangle_index_elvation_main_rain.append(triangle_index_elevation)  # print triangle_index_elevation
            domain.set_quantity('Rain',
                                numeric=triangle_index_elvation_main_rain,  # triangle_elevation
                                use_cache=cache,
                                verbose=True,
                                alpha=0.1, indices=triangle_index_rain,
                                location='vertices')  # alpha refer user manual page 98 | removed ,location='cen$

            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            linecache.clearcache()
            rain_operator.set_rate(rate=Q)
            rain_set_zero = False
        elif os.path.exists(imd_rain_file_gfs):
            if myid == 0: print(
                    "Setting up imd gfs rainfall file %s !!" % imd_rain_file_gfs)
            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            for tri_index in range(len(domain)):
                vertices = domain.get_triangles(tri_index)
                global_vertices = domain.node_l2g[vertices]
                triangle_index_rain.append(tri_index)
                triangle_index_elevation = []
                triangle_index_elevation.insert(0, np.double(
                    linecache.getline(imd_rain_file_gfs, global_vertices[0] + 1)) * imd_rainfall_factor_gfs)
                triangle_index_elevation.insert(1, np.double(
                    linecache.getline(imd_rain_file_gfs, global_vertices[1] + 1)) * imd_rainfall_factor_gfs)
                triangle_index_elevation.insert(2, np.double(
                    linecache.getline(imd_rain_file_gfs, global_vertices[2] + 1)) * imd_rainfall_factor_gfs)
                triangle_index_elvation_main_rain.append(triangle_index_elevation)  # print triangle_index_elevation
            domain.set_quantity('Rain',
                                numeric=triangle_index_elvation_main_rain,  # triangle_elevation
                                use_cache=cache,
                                verbose=True,
                                alpha=0.1, indices=triangle_index_rain,
                                location='vertices')  # alpha refer user manual page 98 | removed ,location='cen$

            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            linecache.clearcache()
            rain_operator.set_rate(rate=Q)
            rain_set_zero = False
        elif os.path.exists(gpm_rain_file):
            if myid == 0: print(
                    "Setting up NOAA GPM rainfall file %s !!" % gpm_rain_file)
            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            for tri_index in range(len(domain)):
                vertices = domain.get_triangles(tri_index)
                global_vertices = domain.node_l2g[vertices]
                triangle_index_rain.append(tri_index)
                triangle_index_elevation = []
                triangle_index_elevation.insert(0, np.double(
                    linecache.getline(gpm_rain_file, global_vertices[0] + 1)) * gpm_rainfall_factor)
                triangle_index_elevation.insert(1, np.double(
                    linecache.getline(gpm_rain_file, global_vertices[1] + 1)) * gpm_rainfall_factor)
                triangle_index_elevation.insert(2, np.double(
                    linecache.getline(gpm_rain_file, global_vertices[2] + 1)) * gpm_rainfall_factor)
                triangle_index_elvation_main_rain.append(triangle_index_elevation)  # print triangle_index_elevation
            domain.set_quantity('Rain',
                                numeric=triangle_index_elvation_main_rain,  # triangle_elevation
                                use_cache=cache,
                                verbose=True,
                                alpha=0.1, indices=triangle_index_rain,
                                location='vertices')  # alpha refer user manual page 98 | removed ,location='cen$
            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            linecache.clearcache()
            rain_operator.set_rate(rate=Q)
            rain_set_zero = False
        elif os.path.exists(gfs_rain_file):
            if myid == 0: print(
                    "Setting up NOAA GFS rainfall file %s !!" % gfs_rain_file)
            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            for tri_index in range(len(domain)):
                vertices = domain.get_triangles(tri_index)
                global_vertices = domain.node_l2g[vertices]
                triangle_index_rain.append(tri_index)
                triangle_index_elevation = []
                triangle_index_elevation.insert(0, np.double(
                    linecache.getline(gfs_rain_file, global_vertices[0] + 1)) * gfs_rainfall_factor)
                triangle_index_elevation.insert(1, np.double(
                    linecache.getline(gfs_rain_file, global_vertices[1] + 1)) * gfs_rainfall_factor)
                triangle_index_elevation.insert(2, np.double(
                    linecache.getline(gfs_rain_file, global_vertices[2] + 1)) * gfs_rainfall_factor)
                triangle_index_elvation_main_rain.append(triangle_index_elevation)  # print triangle_index_elevation
            domain.set_quantity('Rain',
                                numeric=triangle_index_elvation_main_rain,  # triangle_elevation
                                use_cache=cache,
                                verbose=True,
                                alpha=0.1, indices=triangle_index_rain,
                                location='vertices')  # alpha refer user manual page 98 | removed ,location='cen$

            triangle_index_rain = []
            triangle_index_elvation_main_rain = []
            linecache.clearcache()
            rain_operator.set_rate(rate=Q)
            rain_set_zero = False
        elif rain_set_zero:
            print("The Rainfall IMD/GPM/GFS files does not exist setting rain to zero !!")
            domain.set_quantity('Rain', 0.00) #set rainfall with hardcoded value = 4mm #modifed by RK
            rain_operator.set_rate(rate=Q)
    else:
        if myid == 0: print("Using previously set Daily Rainfall!!")

    if myid == 0:
        volume = domain.compute_total_volume()
        stats = domain.timestepping_statistics()
        rainstats = rain_operator.timestepping_statistics()
        maxInundation = Q.get_maximum_value()
        indices = domain.get_wet_elements()
        element_count = len(indices)
        try:
            rain_arr = re.findall(r"\d+\.\d+", rainstats)
            total_rain = float(rain_arr[0])
        except (IndexError, ValueError):
            total_rain = 0.0
            element_count = 0
            maxInundation = 0.0
        with open("rain_data.txt", "a+") as f:
            f.write(str(total_rain) + "\n")
        with open("wet_elements.txt", "a+") as f:
            f.write(str(element_count) + "\n")
        with open("max_inandation.txt", "a+") as f:
            f.write(str(maxInundation) + "\n")
        print("Total Volume =: %s Time Stepping Statistics =: %s Rain Operator Statistics=: %s" % (
            str(volume), str(stats), str(rainstats)))


#profiler.disable()
barrier()


#if myid == 0:
#    print("\n" + "="*80)
#    print("PROFILING RESULTS - Top 40 by cumulative time")
#    print("="*80)
#    stats = pstats.Stats(profiler)
#    stats.sort_stats('cumulative')
#    stats.print_stats(40)

   # Also save to file for detailed analysis
#    profiler.dump_stats('profile.prof')
#    print(f"\nProfile saved to profile.prof - view with: python -m pstats profile.prof")


for p in range(numprocs):
    if myid == p:
        print('Processor %g ' % myid)
        print('That took %.2f seconds' % (time.time() - t0))
        print('Communication time %.2f seconds' % domain.communication_time)
        print('Reduction Communication time %.2f seconds' \
              % domain.communication_reduce_time)
        print('Broadcast time %.2f seconds' % domain.communication_broadcast_time)
    else:
        pass

    barrier()

# --------------------------------------------------
# Merge the individual sww files into one file
# --------------------------------------------------
domain.sww_merge(verbose=verbose, delete_old=True)  # merging all the .sww part files
'''if myid == 0:
    os.system('ssh login02')
    os.system("cd " + data_path)
    out_name = 'output_final_' + str(t0).replace(".", "_") + '_'
    cmd = 'sbatch export_result.sh ' + str(sww_filename) + ' ' + out_name + ' ' + str(finaltime / yieldstep)
    os.system(cmd)
    os.system('exit')'''
finalize()

