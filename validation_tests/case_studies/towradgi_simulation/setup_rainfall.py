
from anuga import read_polygon
from anuga import file_function
from os.path import join

from anuga import Rate_operator

rain_polygon_dir = join('Forcing', 'Rainfall', 'Gauge')
rain_rate_dir  = join('Forcing', 'Rainfall', 'Hort')



def new_setup_rainfall(simulation):

    domain  = simulation.domain

    rate = 0.0002

    return Rate_operator(domain, rate=rate, factor=1.0e-3)


def setup_rainfall(simulation):
    """
    Setup all the rain operators 
    """
    
    # pull out info from the simulation object
    domain  = simulation.domain

    #===============================================================================
    # helper function to encode rainfall operators
    #===============================================================================
    def create_a_rain_operator(domain, base_filename):
        
        polygon = read_polygon(join(rain_polygon_dir, base_filename+'.csv'))
        rate = file_function(join(rain_rate_dir, base_filename+'.tms'), quantities=['rate'])
        
        return Rate_operator(domain, rate=rate, factor=1.0e-3, \
                             polygon = polygon)
       
        
    #=======================================================================
    # The base_filenames of the polygon and rain gauge files
    #=======================================================================
    base_filenames = ['100', '101', '200', '103', '104', '1200', '300', '106', '400', '108', 
                      '500', '5005', '5004', '5006', '501', '502', 
                      '600', '6005', '601', '602', '603', '504', '110', 
                      '700', '701', '702', '7021', '703', '112', 
                      '800', '801', '8002', '802', '8021', '803', '900', '901', '805', 
                      '114', '1000', '1001', '1002', '1003', '116', '117',
                      '11001', '11002', '1100', '119']
    
    #===========================================================================
    # Create all the operators
    #===========================================================================
    for base_filename in base_filenames:
        create_a_rain_operator(domain, base_filename)
        
    
