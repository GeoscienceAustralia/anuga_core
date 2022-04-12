
.. currentmodule:: anuga

Evolve
=======

Running a ANUGA model involves four basic steps:

* :doc:`Creating a domain <domain>`
* :doc:`Setting up the initial conditions <initial_conditions>`
* :doc:`Settting up the boundary condition <boundaries>`
* :doc:`Evolving the model <evolve>`

Here we describe the last step, how to run (evolve) the model for a specified amount of time. 

It would 
also be good to be able to interact with the evolving model. This 
is provided by the :doc:`evolve <anuga.Domain.evolve>` method 
of the :doc:`Domain </reference/anuga.Domain>` object. 

Suppose a Domain we have created, and the initial conditions and boundary conditions set. 
For example here is such a setup for a domain object called :code:`domain`:

>>> domain = anuga.rectangular_cross_domain(10,5)
>>> domain.set_quantity('elevation', function = lambda x,y : x/10)
>>> domain.set_quantity('stage', expression = "elevation + 0.2" )
>>> Br = anuga.Reflective_boundary(domain)
>>> domain.set_boundary({'left' : Br, 'right' : Br, 'top' : Br, 'bottom' : Br})


To evolve the model we would use the domain's evolve method, using the following 
code:

>>> for t in domain.evolve(yieldstep=1.0, finaltime=10.0):
>>>    pass

This will run the model from `time=0`` to the `finaltime 10.0`. By default the state 
of the simulation will be saved to a file (with `sww` format, by default named `domain.sww`) 
every yieldstep, in this case every 1 second of simulation time. 

As the `evolve` construct provides a `for` loop (via the python `yield`` construct) it is possible
to include extra code within the loop. A typical `evolve` loop can provide some printed feedback, i.e.

>>> for t in domain.evolve(yieldstep=1.0, finaltime=10.0):
>>>    domain.print_timestepping_statistics()
Time = 0.0000 (sec), steps=0 (33s)
Time = 1.0000 (sec), delta t in [0.00858871, 0.01071429] (s), steps=111 (0s)
Time = 2.0000 (sec), delta t in [0.00832529, 0.00994060] (s), steps=110 (0s)
Time = 3.0000 (sec), delta t in [0.00901413, 0.00993095] (s), steps=106 (0s)
Time = 4.0000 (sec), delta t in [0.00863985, 0.00963487] (s), steps=109 (0s)
Time = 5.0000 (sec), delta t in [0.00887345, 0.00990731] (s), steps=106 (0s)
Time = 6.0000 (sec), delta t in [0.00934142, 0.00988233] (s), steps=104 (0s)
Time = 7.0000 (sec), delta t in [0.00904828, 0.00970252] (s), steps=107 (0s)
Time = 8.0000 (sec), delta t in [0.00917360, 0.00985509] (s), steps=106 (0s)
Time = 9.0000 (sec), delta t in [0.00925747, 0.00984041] (s), steps=104 (0s)
Time = 10.0000 (sec), delta t in [0.00927581, 0.00973202] (s), steps=106 (0s)

While the yieldsteps are fixed, to maintain stability of the simulation, it is required 
to evolve with inner evolve timesteps which are generally much smaller than the yieldstep. 
The number of these inner evolve timesteps are reported as steps and the
range of the sizes of these evolve timesteps are reported as the delta t. 


It can also be convenient to evolve for a specific `duration`. In this case we replace the `finaltime`
argument with `duration`. I.e. let us continue the evolution for 7 seconds with yieldstep now
set to 2 seconds.

>>> for t in domain.evolve(yieldstep=2.0, duration=7.0):
>>>    domain.print_timestepping_statistics()
Time = 12.0000 (sec), delta t in [0.00932516, 0.00982159] (s), steps=209 (63s)
Time = 14.0000 (sec), delta t in [0.00941363, 0.00981322] (s), steps=210 (0s)
Time = 16.0000 (sec), delta t in [0.00944121, 0.00979934] (s), steps=208 (0s)
Time = 17.0000 (sec), delta t in [0.00945517, 0.00978655] (s), steps=105 (0s)


Sometimes it is necessary to interact with the evolution using a small 
`yieldstep`  but in this case the `sww` file stored at each yieldstep can 
become prohibitively large. 

Instead you can save the state every `outputstep` time interval, while still 
interacting every `yieldstep` interval. 

For instance. let us continue the evolution, but now with a smaller yieldstep 
of 0.5 seconds, but with output to  :code:`domain.sww` every 2 seconds.

>>> for t in domain.evolve(yieldstep=0.5, outputstep=2.0, duration=4.0):
>>>     domain.print_timestepping_statistics()
Time = 17.5000 (sec), delta t in [0.00964414, 0.00977317] (s), steps=52 (650s)
Time = 18.0000 (sec), delta t in [0.00946685, 0.00972477] (s), steps=53 (0s)
Time = 18.5000 (sec), delta t in [0.00953534, 0.00965620] (s), steps=53 (0s)
Time = 19.0000 (sec), delta t in [0.00955560, 0.00976215] (s), steps=52 (0s)
Time = 19.5000 (sec), delta t in [0.00947717, 0.00955428] (s), steps=53 (0s)
Time = 20.0000 (sec), delta t in [0.00955552, 0.00966630] (s), steps=53 (0s)
Time = 20.5000 (sec), delta t in [0.00951811, 0.00975266] (s), steps=52 (0s)
Time = 21.0000 (sec), delta t in [0.00948645, 0.00957223] (s), steps=53 (0s)


Actually typical situations could be `yieldstep = 1.0` and `outputstep=300`. 
In this case the `sww` file will be 300 times smaller than just using `yieldstep` for 
output. 

Start Time
----------

By default the evolution starts at time 0.0. TO set another start time, simply set the starttime
before the evolve loop, i.e.

>>> domain.set_starttime(-3600*24)

to set the start time one day in the past (from ANUGA's zero time). This can be used to allow 
the model to "burn in" before starting the evolution proper. 

More Subtle Start times
-----------------------

To work with dates, times and timezones we can use the python modules :code:`datetime`.
to setup a date and time (and timezone) associated with ANUGA's starttime time. 
Note the use of the :code:`datetime` argument for the 
:code:`print_timestepping_statisitics` procedure.

>>> import anuga
>>> from datetime import datetime
>>> domain = anuga.rectangular_cross_domain(10,5)
>>> domain.set_quantity('elevation', function = lambda x,y : x/10)
>>> domain.set_quantity('stage', expression = "elevation + 0.2" )
>>> Br = anuga.Reflective_boundary(domain)
>>> domain.set_boundary({'left' : Br, 'right' : Br, 'top' : Br, 'bottom' : Br})
>>>
>>> domain.set_timezone('Australia/Sydney')
>>> starttime = datetime(2021, 7, 21, 18, 45)
>>> domain.set_starttime(starttime)
>>> finaltime = datetime(2021, 7, 21, 19, 0)
>>>
>>> for t in domain.evolve(yieldstep=300, finaltime=finaltime):
>>>   domain.print_timestepping_statistics(datetime=True)
DateTime: 2021-07-21 18:45:00+1000, steps=0 (0s)
DateTime: 2021-07-21 18:50:00+1000, delta t in [0.00832571, 0.01071429] (s), steps=31233 (10s)
DateTime: 2021-07-21 18:55:00+1000, delta t in [0.00959070, 0.00964172] (s), steps=31205 (10s)
DateTime: 2021-07-21 19:00:00+1000, delta t in [0.00959070, 0.00964172] (s), steps=31205 (10s)


Essentially we use unix time as our absolute time. So time = 0 corresponds to Jan 1st 1970 UTC. 

For instance going back to an earlier example, and use the argument :code:`datetime` for the 
:code:`print_timestepping_statisitics` procedure

>>> import anuga
>>> from datetime import datetime
>>> domain = anuga.rectangular_cross_domain(10,5)
>>> domain.set_quantity('elevation', function = lambda x,y : x/10)
>>> domain.set_quantity('stage', expression = "elevation + 0.2" )
>>> Br = anuga.Reflective_boundary(domain)
>>> domain.set_boundary({'left' : Br, 'right' : Br, 'top' : Br, 'bottom' : Br})
>>>
>>> for t in domain.evolve(yieldstep=1, finaltime=5):
>>>   domain.print_timestepping_statistics(datetime=True)
DateTime: 1970-01-01 00:00:00+0000, steps=0 (10s)
DateTime: 1970-01-01 00:00:01+0000, delta t in [0.00858871, 0.01071429] (s), steps=111 (0s)
DateTime: 1970-01-01 00:00:02+0000, delta t in [0.00832529, 0.00994060] (s), steps=110 (0s)
DateTime: 1970-01-01 00:00:03+0000, delta t in [0.00901413, 0.00993095] (s), steps=106 (0s)
DateTime: 1970-01-01 00:00:04+0000, delta t in [0.00863985, 0.00963487] (s), steps=109 (0s)
DateTime: 1970-01-01 00:00:05+0000, delta t in [0.00887345, 0.00990731] (s), steps=106 (0s)


Domain.evolve
-------------

.. autosummary::
   :toctree:  
   
   Domain.evolve