.. currentmodule:: anuga



Script Structure
=================


.. only:: html

.. toctree::
   :maxdepth: 1

   domain
   coordinate_reference
   initial_conditions
   boundaries
   operators
   evolve
   checkpointing
   logging
   
.. only:: html


The common way to run an ANUGA model is to use a script. This script will setup the
model. Running the script will build the model, evolve the model and concurrently 
save the results.

Setting up an ANUGA model involves five basic steps:

1. Define the computational domain

2. Set the initial conditions

3. Define the boundary conditions

4. Specify the operators

5. Evolve the model

A simple example of an ANUGA script is shown below:


>>> # Simple rain example
>>> import math
>>> import anuga
>>> 
>>> 
>>> # Setup the domain
>>> domain = anuga.rectangular_cross_domain(10,5)
>>> 
>>> # Set the initial conditions
>>> domain.set_quantity('elevation', function = lambda x,y : x/10)
>>> domain.set_quantity('stage', expression = "elevation + 0.2" )
>>> 
>>> # Define the boundary conditions
>>> Br = anuga.Reflective_boundary(domain)
>>> domain.set_boundary({'left' : Br, 'right' : Br, 'top' : Br, 'bottom' : Br})
>>> 
>>> # Specify the operators
>>> rain = anuga.Rate_operator(domain, rate=lambda t: math.exp( -t**2 ), factor=0.001)
>>> 
>>> # Evolve the model
>>> for t in domain.evolve(yieldstep=1.0, finaltime=10.0):
>>>     domain.print_timestepping_statistics()
>>>
Time = 0.0000 (sec), steps=0 (0s)
Time = 1.0000 (sec), delta t in [0.00858488, 0.01071429] (s), steps=111 (0s)
Time = 2.0000 (sec), delta t in [0.00832106, 0.00991988] (s), steps=111 (0s)
Time = 3.0000 (sec), delta t in [0.00900245, 0.00991667] (s), steps=106 (0s)
Time = 4.0000 (sec), delta t in [0.00862595, 0.00962196] (s), steps=109 (0s)
Time = 5.0000 (sec), delta t in [0.00885777, 0.00988816] (s), steps=107 (0s)
Time = 6.0000 (sec), delta t in [0.00926474, 0.00986530] (s), steps=105 (0s)
Time = 7.0000 (sec), delta t in [0.00903161, 0.00969116] (s), steps=107 (0s)
Time = 8.0000 (sec), delta t in [0.00915208, 0.00983972] (s), steps=106 (0s)
Time = 9.0000 (sec), delta t in [0.00924040, 0.00982541] (s), steps=105 (0s)
Time = 10.0000 (sec), delta t in [0.00930656, 0.00972225] (s), steps=106 (0s)

This script sets up a simple rectangular domain with a sloping bed, 
reflective boundaries and a rainfall operator. The model is then 
evolved for 10 seconds of simulation time, with the state of the model
being saved every 1 second.



