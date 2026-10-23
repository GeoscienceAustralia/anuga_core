def test_create_culvert_bridge_operator_boyd_pipe():
    """
    Test creation of Boyd pipe operator (based on diameter being set)
    """

    import os
    import anuga
    import numpy as num

    from anuga import Boyd_pipe_operator

    ft = 1.0
    bd = 4.0

    file_contents = \
"""
#--------------------------
# CULVERT DATA FOR: Nrth Bagot Rd
#--------------------------
#blockage = 0.0
label = 'Nrth_Bagot_6x600'
diameter = 0.6
#width=3.04
#height=2.45
apron = 0.5
enquiry_points = [ [9.,5.], [14.,5.] ]
manning=0.013
invert_elevations=[0.0,0.1]
#losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
losses = 1.0
#end_points = [[10,5], [12,5]]

exchange_lines = [ [[10, 1], [10, 9]],
                [[12, 1], [12, 9]] ]
"""

    culvert_bridge_file = 'test_boyd_pipe.txt'
    with open(culvert_bridge_file, 'w') as f:
        f.write(file_contents)

    # Create domain and add culvert/bridge operator
    # using the culvert_bridge_file
    domain1 = anuga.rectangular_cross_domain(30, 10, 30, 10)
    domain1.set_name('domain1_boyd_pipe')
    domain1.set_store(False)
    Br = anuga.Reflective_boundary(domain1)
    Bd = anuga.Dirichlet_boundary([bd, 0, 0])
    domain1.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
    anuga.Create_culvert_bridge_Operator(domain1, culvert_bridge_file)
    os.remove(culvert_bridge_file)
    domain1.evolve_to_end(ft)

    # Create domain and add culvert/bridge operator
    # using the local_vars dictionary
    domain2 = anuga.rectangular_cross_domain(30, 10, 30, 10)
    domain2.set_name('domain2_boyd_pipe')
    domain2.set_store(False)
    Br = anuga.Reflective_boundary(domain2)
    Bd = anuga.Dirichlet_boundary([bd, 0, 0])
    domain2.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
    local_vars= {}
    exec(file_contents, {}, local_vars)
    #print(local_vars)
    Boyd_pipe_operator(domain2, **local_vars)
    domain2.evolve_to_end(ft)

    # Check that the two domains give identical results
    s1 = domain1.get_quantity('stage').centroid_values
    s2 = domain2.get_quantity('stage').centroid_values

    # print('s1=', s1)
    # print('s2=', s2)

    # print(num.sum(s1 - s2))

    assert num.sum(s1 - s2) == 0.0

    # try:
    #     os.remove('domain1_boyd_pipe.sww')
    #     os.remove('domain2_boyd_pipe.sww')
    # except FileNotFoundError:
    #     pass

def test_create_culvert_bridge_operator_boyd_box():
    """
    Test creation of Boyd_box_operator (width and heightbeing set)
    """

    import os
    import anuga
    import numpy as num

    from anuga import Boyd_box_operator

    ft = 1.0
    bd = 4.0


    file_contents = \
"""
#--------------------------
# CULVERT DATA FOR: Nrth Bagot Rd
#--------------------------
#blockage = 0.0
label = 'Nrth_Bagot_6x600'
#diameter = 0.6
width=3.04
height=2.45
apron = 0.5
enquiry_points = [ [9.,5.], [14.,5.] ]
manning=0.013
invert_elevations=[0.0,0.1]
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
#end_points = [[10,5], [12,5]]

exchange_lines = [ [[10, 1], [10, 9]],
                [[12, 1], [12, 9]] ]
"""

    #print(file_contents)
    culvert_bridge_file = 'test_boyd_box.txt'
    with open(culvert_bridge_file, 'w') as f:
        f.write(file_contents)

    #Create domain and add culvert/bridge operator
    #using the culvert_bridge_file
    domain3 = anuga.rectangular_cross_domain(30, 10, 30, 10)
    domain3.set_name('domain3_boyd_box')
    domain3.set_store(False)
    Br = anuga.Reflective_boundary(domain3)
    Bd = anuga.Dirichlet_boundary([bd, 0, 0])
    domain3.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
    anuga.Create_culvert_bridge_Operator(domain3, culvert_bridge_file)
    os.remove(culvert_bridge_file)
    domain3.evolve_to_end(ft)

    #Create domain and add culvert/bridge operator
    #using the local_vars dictionary
    domain4 = anuga.rectangular_cross_domain(30, 10, 30, 10)
    domain4.set_name('domain4_boyd_box')
    domain4.set_store(False)
    Br = anuga.Reflective_boundary(domain4)
    Bd = anuga.Dirichlet_boundary([bd, 0, 0])
    domain4.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
    local_vars= {}
    exec(file_contents, {}, local_vars)
    #print(local_vars)
    Boyd_box_operator(domain4, **local_vars)
    domain4.evolve_to_end(ft)

    # Check that the two domains give identical results
    s3 = domain3.get_quantity('stage').centroid_values
    s4 = domain4.get_quantity('stage').centroid_values

    # print('s3=', s3)
    # print('s4=', s4)

    # print(num.sum(s3 - s4))

    assert num.sum(s3 - s4) == 0.0

    try:
        os.remove('domain3_boyd_box.sww')
        os.remove('domain4_boyd_box.sww')
    except FileNotFoundError:
        pass

def test_create_culvert_bridge_operator_weir_orifice_trapezoid():
    """
    Test the creation and equivalence of the Weir_orifice_trapezoid_operator via two methods:
    1. Using a culvert/bridge configuration file.
    2. Using a dictionary of local variables parsed from the same configuration.
    The test performs the following steps:
    - Writes a culvert/bridge configuration to a temporary file.
    - Creates a domain and adds the operator using the configuration file.
    - Evolves the domain for a short timestep.
    - Creates a second domain and adds the operator using the parsed local variables.
    - Evolves the second domain for the same timestep.
    - Asserts that the resulting 'stage' quantity is identical for both domains.
    - Cleans up temporary files generated during the test.

    Test creation of Weir_orifice_trapezoid_operator (based on z1 or z2 being set)
    """

    import os
    import anuga
    import numpy as num

    ft = 1.0
    bd = 4.0

    file_contents = \
"""
#--------------------------
# CULVERT DATA FOR: Nrth Bagot Rd
#--------------------------
#blockage = 0.0
label = 'Nrth_Bagot_6x600'
#diameter = 0.6
width=3.04
height=2.45
z1 = 1.0
z2 = 0.5
apron = 0.5
enquiry_points = [ [9.,5.], [14.,5.] ]
manning=0.013
invert_elevations=[0.0,0.1]
losses = {'inlet':0.5, 'outlet':1.0, 'bend':0.0, 'grate':0.0, 'pier': 0.0, 'other': 0.0}
#end_points = [[10,5], [12,5]]

exchange_lines = [ [[10, 1], [10, 9]],
                [[12, 1], [12, 9]] ]
"""

    culvert_bridge_file = 'test_weir_orifice_trapezoid_operator.txt'
    with open(culvert_bridge_file, 'w') as f:
        f.write(file_contents)

    #Create domain and add culvert/bridge operator
    #using the culvert_bridge_file
    domain1 = anuga.rectangular_cross_domain(30, 10, 30, 10)
    domain1.set_name('domain1_weir_orifice_trapezoid')
    domain1.set_store(False)
    Br = anuga.Reflective_boundary(domain1)
    Bd = anuga.Dirichlet_boundary([bd, 0, 0])
    domain1.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
    anuga.Create_culvert_bridge_Operator(domain1, culvert_bridge_file)
    os.remove(culvert_bridge_file)
    domain1.evolve_to_end(ft)

    # Create domain and add culvert/bridge operator
    # using the local_vars dictionary
    domain2 = anuga.rectangular_cross_domain(30, 10, 30, 10)
    domain2.set_name('domain2_weir_orifice_trapezoid')
    domain2.set_store(False)
    Br = anuga.Reflective_boundary(domain2)
    Bd = anuga.Dirichlet_boundary([bd, 0, 0])
    domain2.set_boundary({'left': Bd, 'right': Br, 'top': Br, 'bottom': Br})
    local_vars= {}
    exec(file_contents, {}, local_vars)
    #print(local_vars)
    anuga.Weir_orifice_trapezoid_operator(domain2, **local_vars)
    domain2.evolve_to_end(ft)

    # Check that the two domains give identical results
    s1 = domain1.get_quantity('stage').centroid_values
    s2 = domain2.get_quantity('stage').centroid_values

    # print('s1=', s1)
    # print('s2=', s2)

    # print(num.sum(s1 - s2))

    assert num.sum(s1 - s2) == 0.0

    try:
        os.remove('domain1_weir_orifice_trapezoid.sww')
        os.remove('domain2_weir_orifice_trapezoid.sww')
    except FileNotFoundError:
        pass


test_create_culvert_bridge_operator_boyd_pipe()
test_create_culvert_bridge_operator_boyd_box()
test_create_culvert_bridge_operator_weir_orifice_trapezoid()

