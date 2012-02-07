import sys

def get_values_from_domain(domain, stop = False):
   pass 
   wholeequation(halt=stop):


def whole_equation(num_tri = 0, tri_area = 0, time_length = 0, time_step = 0, water_depth = 0,
                  velocity = 0, per_water_cover = 0, cpus = 1, cpu_speed = 12, halt = False):
   
    time = time_equation(num_trit = num_tri, tri_areat = tri_area, time_lengtht = time_length, time_stept = time_step,       water_deptht = water_depth,velocityt = velocity, per_water_covert = per_water_cover, cpust = cpus, cpu_speedt = cpu_speed) 

    memory = space_equation(num_tris = num_tri, tri_areas = tri_area, time_lengths = time_length, time_steps = time_step,   water_depths = water_depth,velocitys = velocity, per_water_covers = per_water_cover, cpuss = cpus, cpu_speeds = cpu_speed)
    result = (time,memory)

    print "This program will run for: " + str(time) + " (s)"
    print "This program will use: " + str(memory)   + " (MB)"
 
    if halt :
       sys.exit()


def time_equation(num_trit = 0, tri_areat = 0, time_lengtht = 0, time_stept = 0, water_deptht = 0,
                  velocityt = 0, per_water_covert = 0, cpust = 1, cpu_speedt = 12):

    return (1*num_trit + 2*tri_areat + 3*time_lengtht + 4*time_stept + 5*water_deptht + 6*velocityt 
            + 7*per_water_covert +8*cpust + 9*cpu_speedt + 10)


def space_equation(num_tris = 0, tri_areas = 0, time_lengths = 0, time_steps = 0, water_depths = 0,
                  velocitys = 0, per_water_covers = 0, cpuss = 1, cpu_speeds = 12):

    return (1*num_tris + 2*tri_areas + 3*time_lengths + 4*time_steps + 5*water_depths + 6*velocitys 
            + 7*per_water_covers +8*cpuss + 9*cpu_speeds + 10)
