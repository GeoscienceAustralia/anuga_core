""" 
Towradgi Creek 17 August 1998 Storm Event Calibration
Ubuntu Linux 12.04 LTS 64 bit
By Petar Milevski, minor edits by Gareth Davies
"""
#------------------------------------------------------------------------------
# IMPORT NECESSARY MODULES
#------------------------------------------------------------------------------
print ' ABOUT to Start Simulation:- Importing Modules'

import anuga
import time
from pylab import figure, plot, axis, quiver, quiverkey, show, title, axhline
from pylab import cos, sin, pi
import numpy
import csv
import os
import zipfile
from os.path import join
import time
import sys
import os.path
import glob
from anuga.abstract_2d_finite_volumes.util import file_function
from anuga.geometry.polygon_function import Polygon_function
from anuga.geometry.polygon import read_polygon, plot_polygons, polygon_area
from anuga.abstract_2d_finite_volumes.generic_domain import Quantity
from anuga.pmesh.mesh_interface import create_mesh_from_regions
from anuga.shallow_water.shallow_water_domain import Domain
from anuga.shallow_water.forcing import Rainfall, Inflow
from anuga.structures.inlet_operator import Inlet_operator
from anuga.utilities import sww_merge

#------------------------------------------------------------------------------
# PARALLEL INTERFACE
#------------------------------------------------------------------------------

from anuga_parallel import distribute, myid, numprocs, finalize,barrier
from anuga_parallel.parallel_operator_factory import Inlet_operator, Boyd_box_operator
from anuga.operators.rate_operators import Polygonal_rate_operator
#------------------------------------------------------------------------------
# TEST FUNCTION AND DICTIONARY
#------------------------------------------------------------------------------

def read_polygon_dir(weight_dict, directory, filepattern='*.csv'):
    # get a list of matching filenames in the directory
    pattern = os.path.join(directory, filepattern)
    files = glob.glob(pattern)

    # check that the dictionary contains *all* the files
    
    errors = []
    for f in files:
        try:
            _ = weight_dict[f]
        except KeyError:
            errors.append(f)
            
    if errors:
        msg = ''
        for f in errors:
            msg = msg + ', ' + f
        raise KeyError, 'Files not defined in dictionary: %s' % msg[2:]

    # now get the result list
    result = []
    for f in files:
        result.append((read_polygon(f), weight_dict[f]))
    return result

#------------------------------------------------------------------------------
#  ADD CATCHMENT INFORMATION HERE 
#------------------------------------------------------------------------------

CatchmentDictionary = {
    join('Model', 'Bdy', 'CreekBanks.csv'): 16.0, 
	join('Model', 'Bdy', 'FineCatchment.csv'): 36.0,
	join('Model', 'Bdy', 'Catchment.csv'): 100.0
    }

ManningDictionary = {
    join('Model', 'Mannings', '1.csv'):	0.03, #park
    join('Model', 'Mannings', '2.csv'):	0.15,
    join('Model', 'Mannings', '3.csv'):	0.15,
    join('Model', 'Mannings', '4.csv'):	0.04,
    join('Model', 'Mannings', '5.csv'):	0.15,
    join('Model', 'Mannings', '6.csv'):	0.15,
    join('Model', 'Mannings', '7.csv'):	0.15,
    join('Model', 'Mannings', '8.csv'):	0.15,
    join('Model', 'Mannings', '9.csv'):	0.03, #park
    join('Model', 'Mannings', '10.csv'): 0.15,
    join('Model', 'Mannings', '11.csv'): 0.15,
    join('Model', 'Mannings', '12.csv'): 0.15,
    join('Model', 'Mannings', '13.csv'): 0.025,
    join('Model', 'Mannings', '14.csv'): 0.15,
    join('Model', 'Mannings', '15.csv'): 0.15,
    join('Model', 'Mannings', '16.csv'): 0.15,
    join('Model', 'Mannings', '17.csv'): 0.15,
    join('Model', 'Mannings', '18.csv'): 0.045,
    join('Model', 'Mannings', '18a.csv'): 0.15,
    join('Model', 'Mannings', '18b.csv'): 0.15,
    join('Model', 'Mannings', '18c.csv'): 0.15,
    join('Model', 'Mannings', '18d.csv'): 0.15,
    join('Model', 'Mannings', '18e.csv'): 0.08, #cokeworks site
    join('Model', 'Mannings', '19.csv'): 0.15,
    join('Model', 'Mannings', '20.csv'): 0.15,
    join('Model', 'Mannings', '21.csv'): 0.15,
    join('Model', 'Mannings', '22.csv'): 0.15,
    join('Model', 'Mannings', '23.csv'): 0.15,
    join('Model', 'Mannings', '24.csv'): 0.05,
    join('Model', 'Mannings', '25.csv'): 0.15,
    join('Model', 'Mannings', '26.csv'): 0.15,
    join('Model', 'Mannings', '27.csv'): 0.15,
    join('Model', 'Mannings', '28.csv'): 0.15,
    join('Model', 'Mannings', '29.csv'): 0.15,
    join('Model', 'Mannings', '30.csv'): 0.15,
    join('Model', 'Mannings', '31.csv'): 0.15,
    join('Model', 'Mannings', '32.csv'): 0.15,
    join('Model', 'Mannings', '33.csv'): 0.15,
    join('Model', 'Mannings', '34.csv'): 0.15,
    join('Model', 'Mannings', '35.csv'): 0.15,
    join('Model', 'Mannings', '36.csv'): 0.05,
    join('Model', 'Mannings', '37.csv'): 0.15,
    join('Model', 'Mannings', '38.csv'): 0.15,
    join('Model', 'Mannings', '39.csv'): 0.15,
    join('Model', 'Mannings', '40.csv'): 0.15,
    join('Model', 'Mannings', '41.csv'): 0.15,
    join('Model', 'Mannings', '42.csv'): 0.15,
    join('Model', 'Mannings', '43.csv'): 0.15,
    join('Model', 'Mannings', '44.csv'): 0.15,
    join('Model', 'Mannings', '45.csv'): 0.15,
    join('Model', 'Mannings', '46.csv'): 0.15,
    join('Model', 'Mannings', '47.csv'): 0.15,
    join('Model', 'Mannings', '48.csv'): 0.15,
    join('Model', 'Mannings', '49.csv'): 0.15,
    join('Model', 'Mannings', '50.csv'): 0.15,
    join('Model', 'Mannings', '51.csv'): 0.15,
    join('Model', 'Mannings', '52.csv'): 0.15,
    join('Model', 'Mannings', '53.csv'): 0.15,
    join('Model', 'Mannings', '54.csv'): 0.15,
    join('Model', 'Mannings', '55.csv'): 0.15,
    join('Model', 'Mannings', '56.csv'): 0.15,
    join('Model', 'Mannings', '57.csv'): 0.15,
    join('Model', 'Mannings', '58.csv'): 0.15,
    join('Model', 'Mannings', '59.csv'): 0.08,
    join('Model', 'Mannings', '60.csv'): 0.15,
    join('Model', 'Mannings', '61.csv'): 0.08,
    join('Model', 'Mannings', '62.csv'): 0.15,
    join('Model', 'Mannings', '63.csv'): 0.08,
    join('Model', 'Mannings', '64.csv'): 0.15,
    join('Model', 'Mannings', '65.csv'): 0.15,
    join('Model', 'Mannings', '66.csv'): 0.15,
    join('Model', 'Mannings', '67.csv'): 0.15,
    join('Model', 'Mannings', '68.csv'): 0.15,
    join('Model', 'Mannings', '69.csv'): 0.15,
    join('Model', 'Mannings', '70.csv'): 0.15,
    join('Model', 'Mannings', '71.csv'): 0.05,
    join('Model', 'Mannings', '72.csv'): 0.15,
    join('Model', 'Mannings', '73.csv'): 0.15,
    join('Model', 'Mannings', '74.csv'): 0.15,
    join('Model', 'Mannings', '75.csv'): 0.15,
    join('Model', 'Mannings', '76.csv'): 0.15,
    join('Model', 'Mannings', '77.csv'): 0.07,
    join('Model', 'Mannings', '78.csv'): 0.15,
    join('Model', 'Mannings', '79.csv'): 0.15,
    join('Model', 'Mannings', '80.csv'): 0.15,
    join('Model', 'Mannings', '81.csv'): 0.15,
    join('Model', 'Mannings', '82.csv'): 0.15,
    join('Model', 'Mannings', '83.csv'): 0.15,
    join('Model', 'Mannings', '84.csv'): 0.15,
    join('Model', 'Mannings', '85.csv'): 0.15,
    join('Model', 'Mannings', '86.csv'): 0.15,
    join('Model', 'Mannings', 'Escarpement.csv'): 0.15,
    join('Model', 'Mannings', 'Railway.csv'): 0.025,
    join('Model', 'Creeks', 'creeks.csv'): 0.025
	}

basename = join('DEM', 'towradgi')
outname = join('Towradgi_historic_flood')
meshname = join('DEM','towradgi.tsh')

W=303517
N=6195670
E=308570
S=6193140

maximum_triangle_area = 1000
#minimum_storable_height = 0.05
base_friction = 0.02
alpha = 0.99
verbose = True

#model_output_dir='MODEL_OUTPUTS_'+time.strftime('%Y%m%d_%H%M%S', time.localtime())  
model_output_dir='MODEL_OUTPUTS'
try:
    os.mkdir(model_output_dir)
except:
    pass
#--------------------------------------------------------------------------
# Setup Domain
#--------------------------------------------------------------------------
if myid == 0:
    #------------------------------------------------------------------------------
    # CREATING MESH
    #------------------------------------------------------------------------------
    
    bounding_polygon = [[W, S], [E, S], [E, N], [W, N]]
    interior_regions = read_polygon_dir(CatchmentDictionary, join('Model', 'Bdy'))
    
    create_mesh_from_regions(bounding_polygon, 
        boundary_tags={'south': [0], 'east': [1], 'north': [2], 'west': [3]},
        maximum_triangle_area=maximum_triangle_area,
        interior_regions=interior_regions,
        filename=meshname,
        use_cache=False,
        verbose=True)
    
    #------------------------------------------------------------------------------
    # SETUP COMPUTATIONAL DOMAIN
    #------------------------------------------------------------------------------
    
    domain = Domain(meshname, use_cache=False, verbose=True)

    from anuga.utilities.argparsing import parse_standard_args
    alg, cfl = parse_standard_args()
    domain.set_flow_algorithm(alg)
    domain.set_CFL(cfl)
    #domain.set_flow_algorithm('DE1') 
    #domain.set_store_vertices_smoothly()
    #domain.set_minimum_storable_height(minimum_storable_height) 
    domain.set_datadir(model_output_dir)
    try:
        os.mkdir()
    except:
        pass
    domain.set_name(outname)
    
    print domain.statistics()
    
    #------------------------------------------------------------------------------
    # APPLY MANNING'S ROUGHNESSES
    #------------------------------------------------------------------------------
    
    friction_list = read_polygon_dir(ManningDictionary, join('Model', 'Mannings'))
    domain.set_quantity('friction', Polygon_function(friction_list, default=base_friction, geo_reference=domain.geo_reference))
    
    # Set a Initial Water Level over the Domain
    domain.set_quantity('stage', 0)
   
    # Decompress the zip file to make a csv for reading 
    #zipfile.ZipFile('DEM/towradgi.zip').extract('towradgi.csv',path='DEM/')
    # Set the elevation with the new csv file
    domain.set_quantity('elevation', filename=basename+'.csv', use_cache=True, verbose=True, alpha=alpha)
    #os.remove('DEM/towradgi.csv') # Clean up csv file
else:
	domain = None

barrier()
if myid == 0 and verbose: print 'DISTRIBUTING DOMAIN'
domain = distribute(domain)
barrier()
 
if myid == 0 and verbose: print 'CREATING INLETS'
"""
#------------------------------------------------------------------------------
# ENTER CULVERT DATA
#------------------------------------------------------------------------------
# BrookerStCulvert branch 2
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[305772.982,6193988.557] , [305772.378,6193987.823]])
el1 = numpy.array([[305794.592,6193983.907] , [305793.988,6193983.173]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=0.9,
                                    exchange_lines=[el0, el1],
                                    height=0.9,
                                    apron=3.0,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '1H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)    

# MeadowStCulvert branch 2
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[305886.333,6193929.052] , [305883.172,6193922.986]])
el1 = numpy.array([[305906.553,6193910.461] , [305903.393,6193904.395]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=5.4,
                                    exchange_lines=[el0, el1],
                                    height=0.6,
                                    apron=3.0,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '2H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)    

# WilliamsStCulvert branch 2
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[305945.955,6193836.293] , [305945.125,6193835.387]])
el1 = numpy.array([[306040.565,6193827.573] , [306039.735,6193826.667]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=1.2,
                                    exchange_lines=[el0, el1],
                                    height=1.2,
                                    apron=3.0,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '3H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)     

# MeadowStCulvert Towradgi branch
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[305812.113,6193591.972] , [305809.390,6193588.820]])
el1 = numpy.array([[305834.913,6193588.382] , [305832.190,6193585.230]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=4.0,
                                    exchange_lines=[el0, el1],
                                    height=2.2,
                                    apron=3.0,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '4H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)  

# CollinsStCulverts tarra 5 branch
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[306330.608,6194817.116] , [306320.768,6194805.884]])
el1 = numpy.array([[306369.483,6194811.616] , [306359.643,6194800.384]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=3.6,
                                    exchange_lines=[el0, el1],
                                    height=0.93,
                                    apron=3.0,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '9H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)                                     
                                    
# Norther Distributor Culverts branch 5
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[306956.242,6194465.589] , [306950.446,6194457.411]])
el1 = numpy.array([[307003.711,6194446.089] , [306997.916,6194437.911]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=9.0,
                                    exchange_lines=[el0, el1],
                                    height=0.85,
                                    apron=3.0,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '10H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)                                      
                                    
#CokeWorksCulverts branch 5
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[307142.161,6194181.3065] , [307138.519,6194174.394]])
el1 = numpy.array([[307160.521,6194164.8165] , [307156.879,6194157.904]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=6.2,
                                    exchange_lines=[el0, el1],
                                    height=3.0,
                                    apron=3.1,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '11H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)                                      
                                    
#Northern Distributor Branch 6 Culverts
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[306950.758,6193454.717] , [306947.804,6193453.283]])
el1 = numpy.array([[306988.633,6193474.217] , [306985.679,6193472.783]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=3.6,
                                    exchange_lines=[el0, el1],
                                    height=1.2,
                                    apron=3.1,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '12H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)                                      
                                    
#Railway Branch 6 Culverts
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[307139.134,6193474.458] , [307138.492,6193473.542]])
el1 = numpy.array([[307150.884,6193469.458] , [307150.242,6193468.542]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=3.6,
                                    exchange_lines=[el0, el1],
                                    height=1.2,
                                    apron=3.1,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '13H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False) 

#Colgong Branch6 Culverts
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[307200.610,6193476.765] , [307199.140,6193475.235]])
el1 = numpy.array([[307224.610,6193475.765] , [307223.140,6193474.235]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=2.1,
                                    exchange_lines=[el0, el1],
                                    height=1.0,
                                    apron=3.1,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '14H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)   
                                    
#Basin Outlet Branch 3 Culverts
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[305629.639,6194408.883] , [305626.521,6194400.457]])
el1 = numpy.array([[305665.889,6194347.183] , [305662.771,6194338.757]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=6.0,
                                    exchange_lines=[el0, el1],
                                    height=0.86,
                                    apron=3.1,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '15H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)                                      
                                    
#BellambiRd Branch 3 Culverts
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[305777.182,6194305.377] , [305776.444,6194304.623]])
el1 = numpy.array([[305873.807,6194303.377] , [305873.069,6194302.623]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=1.05,
                                    exchange_lines=[el0, el1],
                                    height=1.0,
                                    apron=3.1,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '16H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)    

#MeadowSt Branch 3 Culverts
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[305914.649,6194322.375] , [305913.477,6194321.625]])
el1 = numpy.array([[305950.711,6194335.375] , [305949.539,6194334.625]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=1.5,
                                    exchange_lines=[el0, el1],
                                    height=1.0,
                                    apron=3.1,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '17H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)     

#13MeadowSt Branch 3 Culverts
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[305911.280,6194359.203] , [305910.260,6194358.017]])
el1 = numpy.array([[305946.090,6194353.573] , [305945.070,6194352.387]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=1.5,
                                    exchange_lines=[el0, el1],
                                    height=1.0,
                                    apron=3.1,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '18H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)     

#41AngelStBranch3Culverts
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[306196.779,6194028.193] , [306192.221,6194010.807]])
el1 = numpy.array([[306200.154,6194018.693] , [306195.596,6194001.307]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=10.0,
                                    exchange_lines=[el0, el1],
                                    height=0.35,
                                    apron=3.1,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '19H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)        

#CarrollSt Branch 7 Culverts
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[308002.045,6193820.163] , [308001.215,6193819.197]])
el1 = numpy.array([[308021.965,6193816.883] , [308021.135,6193815.917]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=1.22,
                                    exchange_lines=[el0, el1],
                                    height=0.3,
                                    apron=3.1,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '20H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)           

#ParkerRd Branch 7 Culverts
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[308105.832,6193803.622] , [308103.648,6193801.118]])
el1 = numpy.array([[308126.782,6193800.552] , [308124.598,6193798.048]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=3.2,
                                    exchange_lines=[el0, el1],
                                    height=0.3,
                                    apron=3.1,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '21H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)     

#LakePde Branch 7 Culverts
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[308251.257,6193614.658] , [308248.343,6193614.482]])
el1 = numpy.array([[308233.627,6193598.088] , [308230.713,6193597.912]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=3.0,
                                    exchange_lines=[el0, el1],
                                    height=0.75,
                                    apron=3.1,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '22H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)                                      
                                   
# PrincesHwyBridge Towradgi branch
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[306608.128,6193709.170] , [306600.772,6193696.949]])
el1 = numpy.array([[306627.238,6193696.530] , [306619.882,6193684.309]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=12.0,
                                    exchange_lines=[el0, el1],
                                    height=3.0,
                                    apron=3.0,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '5H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)  

# PioneerRdBridge Towradgi branch
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[307613.615,6193631.615] , [307597.940,6193608.357]])
el1 = numpy.array([[307635.865,6193620.848] , [307267.195,6193705.076]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=20.0,
                                    exchange_lines=[el0, el1],
                                    height=3.0,
                                    apron=3.0,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '8H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)                           

# NorthernDistributorBridge Towradgi branch
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[306944.945,6193760.601] , [306944.569,6193699.449]])
el1 = numpy.array([[306989.167,6193757.821] , [306987.147,6193699.441]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=60.0,
                                    exchange_lines=[el0, el1],
                                    height=6.0,
                                    apron=3.0,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '6H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False)    

# RailwayBridge Towradgi branch, THIS IS THE BRIDGE THAT BLOCKEDAND DAM BREAK OCCURRED
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
el0 = numpy.array([[307209.474,6193752.810] , [307197.942,6193716.493]])
el1 = numpy.array([[307275.853,6193735.391] , [307228.492,6193726.613]])
    
culvert = Boyd_box_operator(domain,
                                    losses=losses,
                                    width=3.0,
                                    exchange_lines=[el0, el1],
                                    height=1.0,
                                    apron=3.0,
                                    enquiry_gap=10.0,
                                    use_momentum_jet=False,
                                    use_velocity_head=False,
                                    manning=0.013,
                                    logging=False,
                                    label=join('Runs', '7H'), # this puts the culvert hydraulics output into the same dir as the sww's
                                    verbose=False) 
"""

#------------------------------------------------------------------------------
# APPLY RAINFALL
#------------------------------------------------------------------------------

Catchment_Rain_Polygon100 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '100.csv'))
rainfall100 = file_function(join('Forcing', 'Rainfall', 'Hort', '100.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall100, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon100, default_rate = 0.0)

Catchment_Rain_Polygon101 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '101.csv'))
rainfall101 = file_function(join('Forcing', 'Rainfall', 'Hort', '101.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall101, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon101, default_rate = 0.0)

Catchment_Rain_Polygon200 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '200.csv'))
rainfall200 = file_function(join('Forcing', 'Rainfall', 'Hort', '200.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall200, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon200, default_rate = 0.0)

Catchment_Rain_Polygon103 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '103.csv'))
rainfall103 = file_function(join('Forcing', 'Rainfall', 'Hort', '103.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall103, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon103, default_rate = 0.0)

Catchment_Rain_Polygon104 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '104.csv'))
rainfall104 = file_function(join('Forcing', 'Rainfall', 'Hort', '104.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall104, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon104, default_rate = 0.0)

Catchment_Rain_Polygon1200 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '1200.csv'))
rainfall1200 = file_function(join('Forcing', 'Rainfall', 'Hort', '1200.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall1200, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon1200, default_rate = 0.0)

Catchment_Rain_Polygon300 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '300.csv'))
rainfall300 = file_function(join('Forcing', 'Rainfall', 'Hort', '300.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall300, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon300, default_rate = 0.0)

Catchment_Rain_Polygon106 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '106.csv'))
rainfall106 = file_function(join('Forcing', 'Rainfall', 'Hort', '106.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall106, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon106, default_rate = 0.0)

Catchment_Rain_Polygon400 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '400.csv'))
rainfall400 = file_function(join('Forcing', 'Rainfall', 'Hort', '400.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall400, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon400, default_rate = 0.0)

Catchment_Rain_Polygon108 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '108.csv'))
rainfall108 = file_function(join('Forcing', 'Rainfall', 'Hort', '108.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall108, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon108, default_rate = 0.0)

Catchment_Rain_Polygon500 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '500.csv'))
rainfall500 = file_function(join('Forcing', 'Rainfall', 'Hort', '500.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall500, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon500, default_rate = 0.0)

Catchment_Rain_Polygon5005 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '5005.csv'))
rainfall5005 = file_function(join('Forcing', 'Rainfall', 'Hort', '5005.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall5005, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon5005, default_rate = 0.0)

Catchment_Rain_Polygon5004 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '5004.csv'))
rainfall5004 = file_function(join('Forcing', 'Rainfall', 'Hort', '5004.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall5004, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon5004, default_rate = 0.0)

Catchment_Rain_Polygon5006 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '5006.csv'))
rainfall5006 = file_function(join('Forcing', 'Rainfall', 'Hort', '5006.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall5006, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon5006, default_rate = 0.0)

Catchment_Rain_Polygon501 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '501.csv'))
rainfall501 = file_function(join('Forcing', 'Rainfall', 'Hort', '501.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall501, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon501, default_rate = 0.0)

Catchment_Rain_Polygon502 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '502.csv'))
rainfall502 = file_function(join('Forcing', 'Rainfall', 'Hort', '502.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall502, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon502, default_rate = 0.0)

Catchment_Rain_Polygon600 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '600.csv'))
rainfall600 = file_function(join('Forcing', 'Rainfall', 'Hort', '600.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall600, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon600, default_rate = 0.0)

Catchment_Rain_Polygon6005 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '6005.csv'))
rainfall6005 = file_function(join('Forcing', 'Rainfall', 'Hort', '6005.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall6005, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon6005, default_rate = 0.0)

Catchment_Rain_Polygon601 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '601.csv'))
rainfall601 = file_function(join('Forcing', 'Rainfall', 'Hort', '601.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall601, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon601, default_rate = 0.0)

Catchment_Rain_Polygon602 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '602.csv'))
rainfall602 = file_function(join('Forcing', 'Rainfall', 'Hort', '602.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall602, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon602, default_rate = 0.0)

Catchment_Rain_Polygon603 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '603.csv'))
rainfall603 = file_function(join('Forcing', 'Rainfall', 'Hort', '603.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall603, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon603, default_rate = 0.0)

Catchment_Rain_Polygon504 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '504.csv'))
rainfall504 = file_function(join('Forcing', 'Rainfall', 'Hort', '504.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall504, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon504, default_rate = 0.0)

Catchment_Rain_Polygon110 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '110.csv'))
rainfall110 = file_function(join('Forcing', 'Rainfall', 'Hort', '110.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall110, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon110, default_rate = 0.0)

Catchment_Rain_Polygon700 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '700.csv'))
rainfall700 = file_function(join('Forcing', 'Rainfall', 'Hort', '700.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall700, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon700, default_rate = 0.0)

Catchment_Rain_Polygon701 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '701.csv'))
rainfall701 = file_function(join('Forcing', 'Rainfall', 'Hort', '701.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall701, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon701, default_rate = 0.0)

Catchment_Rain_Polygon702 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '702.csv'))
rainfall702 = file_function(join('Forcing', 'Rainfall', 'Hort', '702.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall702, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon702, default_rate = 0.0)

Catchment_Rain_Polygon7021 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '7021.csv'))
rainfall7021 = file_function(join('Forcing', 'Rainfall', 'Hort', '7021.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall7021, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon7021, default_rate = 0.0)

Catchment_Rain_Polygon703 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '703.csv'))
rainfall703 = file_function(join('Forcing', 'Rainfall', 'Hort', '703.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall703, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon703, default_rate = 0.0)

Catchment_Rain_Polygon112 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '112.csv'))
rainfall112 = file_function(join('Forcing', 'Rainfall', 'Hort', '112.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall112, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon112, default_rate = 0.0)

Catchment_Rain_Polygon800 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '800.csv'))
rainfall800 = file_function(join('Forcing', 'Rainfall', 'Hort', '800.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall800, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon800, default_rate = 0.0)

Catchment_Rain_Polygon801 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '801.csv'))
rainfall801 = file_function(join('Forcing', 'Rainfall', 'Hort', '801.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall801, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon801, default_rate = 0.0)

Catchment_Rain_Polygon8002 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '8002.csv'))
rainfall8002 = file_function(join('Forcing', 'Rainfall', 'Hort', '8002.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall8002, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon8002, default_rate = 0.0)

Catchment_Rain_Polygon802 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '802.csv'))
rainfall802 = file_function(join('Forcing', 'Rainfall', 'Hort', '802.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall802, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon802, default_rate = 0.0)

Catchment_Rain_Polygon8021 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '8021.csv'))
rainfall8021 = file_function(join('Forcing', 'Rainfall', 'Hort', '8021.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall8021, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon8021, default_rate = 0.0)

Catchment_Rain_Polygon803 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '803.csv'))
rainfall803 = file_function(join('Forcing', 'Rainfall', 'Hort', '803.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall803, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon803, default_rate = 0.0)

Catchment_Rain_Polygon900 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '900.csv'))
rainfall900 = file_function(join('Forcing', 'Rainfall', 'Hort', '900.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall900, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon900, default_rate = 0.0)

Catchment_Rain_Polygon901 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '901.csv'))
rainfall901 = file_function(join('Forcing', 'Rainfall', 'Hort', '901.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall901, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon901, default_rate = 0.0)

Catchment_Rain_Polygon805 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '805.csv'))
rainfall805 = file_function(join('Forcing', 'Rainfall', 'Hort', '805.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall805, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon805, default_rate = 0.0)

Catchment_Rain_Polygon114 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '114.csv'))
rainfall114 = file_function(join('Forcing', 'Rainfall', 'Hort', '114.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall114, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon114, default_rate = 0.0)

Catchment_Rain_Polygon1000 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '1000.csv'))
rainfall1000 = file_function(join('Forcing', 'Rainfall', 'Hort', '1000.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall1000, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon1000, default_rate = 0.0)

Catchment_Rain_Polygon1001 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '1001.csv'))
rainfall1001 = file_function(join('Forcing', 'Rainfall', 'Hort', '1001.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall1001, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon1001, default_rate = 0.0)

Catchment_Rain_Polygon1002 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '1002.csv'))
rainfall1002 = file_function(join('Forcing', 'Rainfall', 'Hort', '1002.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall1002, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon1002, default_rate = 0.0)

Catchment_Rain_Polygon1003 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '1003.csv'))
rainfall1003 = file_function(join('Forcing', 'Rainfall', 'Hort', '1003.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall1003, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon1003, default_rate = 0.0)

Catchment_Rain_Polygon116 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '116.csv'))
rainfall116 = file_function(join('Forcing', 'Rainfall', 'Hort', '116.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall116, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon116, default_rate = 0.0)

Catchment_Rain_Polygon117 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '117.csv'))
rainfall117 = file_function(join('Forcing', 'Rainfall', 'Hort', '117.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall117, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon117, default_rate = 0.0)

Catchment_Rain_Polygon11001 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '11001.csv'))
rainfall11001 = file_function(join('Forcing', 'Rainfall', 'Hort', '11001.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall11001, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon11001, default_rate = 0.0)

Catchment_Rain_Polygon11002 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '11002.csv'))
rainfall11002 = file_function(join('Forcing', 'Rainfall', 'Hort', '11002.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall11002, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon11002, default_rate = 0.0)

Catchment_Rain_Polygon1100 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '1100.csv'))
rainfall1100 = file_function(join('Forcing', 'Rainfall', 'Hort', '1100.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall1100, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon1100, default_rate = 0.0)

Catchment_Rain_Polygon119 = read_polygon(join('Forcing', 'Rainfall', 'Gauge', '119.csv'))
rainfall119 = file_function(join('Forcing', 'Rainfall', 'Hort', '119.tms'), quantities=['rate'])
op1 = Polygonal_rate_operator(domain, rate=rainfall119, factor=1.0e-3, \
                      polygon=Catchment_Rain_Polygon119, default_rate = 0.0)

barrier()
                      
#------------------------------------------------------------------------------
# BOUNDARY CONDITIONS
#------------------------------------------------------------------------------
    
print 'Available boundary tags', domain.get_boundary_tags()

func = file_function(join('Forcing','Tide','Pioneer.tms'), quantities='rainfall')
Bd = anuga.Dirichlet_boundary([0,0,0])
Bw = anuga.Time_boundary(domain=domain, function=lambda t: [func(t)[0], 0.0, 0.0])

domain.set_boundary({'west': Bd, 'south': Bd, 'north': Bd, 'east': Bw})

#----------------------------------------------------------------
# Change conditions during the run to simulate Pioneer Rd Bridge
#----------------------------------------------------------------
#
#ManningDictionary = {join('Model', 'Creeks', 'creeks.csv'):	1.0}
#friction_list = read_polygon_dir(ManningDictionary, join('Model', 'Creeks'))
#
#ManningDictionary = {join('Model', 'Creeks', 'creeks.csv'):	0.025}
#friction_list = read_polygon_dir(ManningDictionary, join('Model', 'Creeks'))
#
#import pdb
#pdb.set_trace()
#
#for i in range(10):
#    print 'stalling'

#------------------------------------------------------------------------------
# EVOLVE SYSTEM THROUGH TIME
#------------------------------------------------------------------------------
barrier()

if myid == 0 and verbose: print 'EVOLVE'
    
t0 = time.time()
    
for t in domain.evolve(yieldstep = 300., finaltime = 83700.):
    #if t == 37800.0: #time when bridge deck starts to get submerged, increase n to act as bridge deck, handrail and blockage effects
    #    print 'Reset friction...'
    #    domain.set_quantity('friction',Polygon_function(friction_list, default=base_friction, geo_reference=domain.geo_reference))       
    #if t == 44100.0: #time when water level drops below bridge deck, bring n back down to existing conditions
    #    print 'Reset friction...'
    #    domain.set_quantity('friction', Polygon_function(friction_list, default=base_friction, geo_reference=domain.geo_reference))           
    if myid == 0:
        domain.write_time()

barrier()
if myid == 0:
    print 'Number of processors %g ' %numprocs
    print 'That took %.2f seconds' %(time.time()-t0)
    print 'Communication time %.2f seconds'%domain.communication_time
    print 'Reduction Communication time %.2f seconds'%domain.communication_reduce_time
    print 'Broadcast time %.2f seconds'%domain.communication_broadcast_time

    os.chdir(model_output_dir)
    sww_merge.sww_merge_parallel(outname, np=numprocs, verbose=True, delete_old=False)

finalize()
