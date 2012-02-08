import sys
import system_tools

TEST_CON = 'test_constants'

test_constants = {'tri_a':1.0, 'tri_b': 10.0}
system_constants = {'tornado.agso.gov.au':{'tri_a':100.0, 'tri_b': 20.0},
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
        print "This program will run for: " + str(time) + " (s)"
        print "This program will use: " + str(memory)   + " (MB)"  
     
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
       

def time_equation(**kwargs):
    
    time = kwargs['constants']['tri_a'] * kwargs['num_tri']
    #time = (1*num_trit + 2*tri_areat + 3*time_lengtht + 
    #       4*time_stept + 5*water_deptht + 6*velocityt 
    # + 7*per_water_covert +8*cpust + 9*cpu_speedt + 10)
    return time


def space_equation(**kwargs):
    memory = 10
    return memory

    
################################################################################

if __name__ == "__main__":
    whole_equation(num_tri = 7)
