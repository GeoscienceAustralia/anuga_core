import sys
from . import system_tools

TEST_CON = 'test_constants'

test_constants = {'tri_a_T':1, 'tri_b_T': 1,
                  'tim_a_T':1,'fil_a_T':1.,'cons_T':1,
                  'tri_a_S':1,'cons_S':1}

# These constants come from the Major Variables script that ran on 
# tornado in serial
system_constants = {'tornado.agso.gov.au':{'tri_a_T':0.0000395, 
                                           'tri_b_T': 0.29575152,
                    'tim_a_T':0.03804736,'fil_a_T':0.005928693,
                    'cons_T':-135.0661178,
                    'tri_a_S':0.00369572,'cons_S':331.7128095},
                    TEST_CON:test_constants}

DEFAULT_HOST = 'tornado.agso.gov.au'


def estimate_time_mem(domain, yieldstep, finaltime, halt=False, 
                      log_results=True, use_test_constants=False):
    """
    Predict the time in seconds and memory in ?? that the simulation
    will need.
    
    params:
      domain: a Domain instance, used to get number of triangles
      yieldstep: the yieldstep of the simulation
      finaltime: The final time used in the simulation.
      halt: Set to True if you want ANUGA to stop after the prediction
      log_results: Add the predictions to the log file.
      use_test_constants: Use artificial test constants.
      
     Example use:
     anuga.estimate_time_mem(domain, yieldstep=yieldstep, finaltime=finaltime, 
                        halt=True)
    """
    time, memory = whole_equation(num_tri=len(domain),
                                  yieldstep=yieldstep,
                                  finaltime=finaltime, 
                                  use_test_constants=use_test_constants)
    
    if log_results: #FIXME, not loging results yet
        print("This program will run for: " + str(time) + " (s)")
        print("This program will use: " + str(memory)   + " (MB)")  
     
    if halt:
        sys.exit()
        
    return time, memory


def whole_equation(halt = False, **kwargs):
    """
    num_tri = None,
    tri_area =  None,
    time_length =  None,
    time_step =  None,
    water_depth =  None,
    velocity =  None,
    per_water_cover =  None, 
    cpus =  None,
    cpu_speed =  None,
    halt = False
    """
    if not kwargs['use_test_constants']:
        host_name = system_tools.get_host_name()
    else:
        host_name = TEST_CON
        
    constants = system_constants.get(host_name, system_constants[DEFAULT_HOST])
    kwargs['constants'] = constants
    
    time = time_equation(**kwargs) 
    memory = space_equation(**kwargs)
    
    result = (time, memory)
    
    return result
       
# Using the constants from the experiments into 
# memory and time the Time and Memory are estimated
def time_equation(**kwargs):
    
    time = kwargs['constants']['tri_a_T'] * (kwargs['num_tri']) ** 2 + \
           kwargs['constants']['tri_b_T'] * kwargs['num_tri'] + \
           kwargs['constants']['tim_a_T'] * kwargs['finaltime'] + \
           kwargs['constants']['fil_a_T'] * \
           ((kwargs['finaltime'] / kwargs['yieldstep'])) + \
           kwargs['constants']['cons_T']
 
    return time


def space_equation(**kwargs):
    memory = kwargs['constants']['tri_a_S'] * kwargs['num_tri'] + \
             kwargs['constants']['cons_S']
    return memory

    
################################################################################

if __name__ == "__main__":
    whole_equation(num_tri = 7)
