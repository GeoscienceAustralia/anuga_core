
.. currentmodule:: anuga

Setting up the Evolve loop
==========================

Running a ANUGA model involves six basic steps:

* Creating a domain
* Setting up the initial conditions
* Settting up the boundary condition
* Setting up any necessary stuctures (culverts etc)
* Setting up any necessary operators (rainfall etc)
* Evolving the model

Here we describe the last step, how to run (evolve) the model for a specified amount of time. 

Evolving the Model
------------------

In addition to evolving the model, it would good to be able to interact with the evolving model. This 
is all provided by the :meth:`evolve <Domain.evolve>` method 
of the :doc:`Domain <domain.rst>` object. 

Suppose we have created and set up a Domain by completing the first 4 basic steps.
For example here is such a setup for a domain object called `domain`:

>>> domain = anuga.rectangular_cross_domain(10,5)
>>> domain.set_quantity('elevation', function = lambda x,y : x/10)
>>> domain.set_quantity('stage', expression = "elevation + 0.2" )
>>> Br = anuga.Reflective_boundary(domain)
>>> domain.set_boundary({'left' : Br, 'right' : Br, 'top' : Br, 'bottom' : Br})


To evolve the model we would use the domain's evolve method, using the following 
code:

>>> for t in domain.evolve(yieldstep=1.0, finaltime=10.0):
>>>    pass

This will run the model from `time=0` to the `finaltime = 10.0`. The method will `yield` 
to the for loop every `yieldstep = 1`. By default the state 
of the simulation will be saved to a file (by default named `domain.sww`) 
every `yieldstep`, in this case every 1 second of simulation time. 

As the `evolve` construct provides a `for` loop (via the python `yield` construct) it is possible
to include extra code within the loop. A typical `evolve` loop can provide some printed feedback
using the :meth:`print_timestepping_statistics <anuga.Domain.print_timestepping_statistics>` method, i.e.,

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

During the evolution the yieldsteps are fixed but to maintain stability of the simulation, the 
underlying computation uses inner evolve timesteps which are generally much smaller than the yieldstep. 
The number of these inner evolve timesteps are reported as steps and the
range of the sizes of these evolve timesteps are reported as the delta t. 

Duration instead of finaltime
-----------------------------

It can also be convenient to evolve for a specific `duration`. In this case we replace the `finaltime`
argument with `duration`. I.e. let us continue the evolution for 7 seconds with yieldstep now
set to 2 seconds.

>>> for t in domain.evolve(yieldstep=2.0, duration=7.0):
>>>    domain.print_timestepping_statistics()
Time = 12.0000 (sec), delta t in [0.00932516, 0.00982159] (s), steps=209 (63s)
Time = 14.0000 (sec), delta t in [0.00941363, 0.00981322] (s), steps=210 (0s)
Time = 16.0000 (sec), delta t in [0.00944121, 0.00979934] (s), steps=208 (0s)
Time = 17.0000 (sec), delta t in [0.00945517, 0.00978655] (s), steps=105 (0s)



Outputstep
----------


Sometimes it is necessary to interact with the evolution using a small 
`yieldstep` (such as controlling a hydraulic structure). In this case the `sww` file stored 
at each yieldstep can become prohibitively large. 

Instead you can save the state every `outputstep` time interval, while still 
interacting every `yieldstep` interval. 

For instance. let us continue the evolution, but now with a smaller yieldstep 
of 0.5 seconds, but with output to  `domain.sww` every 2 seconds.

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


Typical situations could be `yieldstep = 1.0` and `outputstep=300`. 
In this case the `sww` file will be 300 times smaller than just using `yieldstep` for 
output. 

Start Time
----------

By default the evolution starts at time 0.0. To set another start time, simply use 
:meth:`set_starttime <anuga.Domain.set_starttime>` 
before the evolve loop, i.e.

>>> domain.set_starttime(-3600*24)

to set the start time one day in the past (from ANUGA's zero time). This can be used to allow 
the model to "burn in" before starting the evolution proper. 

Start times with DateTime and Timezones
---------------------------------------

To work with dates, times and timezones we can use the python module `datetime`.
to setup a date and time (and timezone) associated with ANUGA's starttime time. 


Once again let's suppose we have setup a domain via:

>>> import anuga
>>> from datetime import datetime
>>> domain = anuga.rectangular_cross_domain(10,5)
>>> domain.set_quantity('elevation', function = lambda x,y : x/10)
>>> domain.set_quantity('stage', expression = "elevation + 0.2" )
>>> Br = anuga.Reflective_boundary(domain)
>>> domain.set_boundary({'left' : Br, 'right' : Br, 'top' : Br, 'bottom' : Br})

By default ANUGA uses a UTC as the default timezone for the domain. 
We can change it via :meth:`set_timezone <anuga.Domain.set_timezone>` 

>>> domain.set_timezone('Australia/Sydney')

A list of timezones names can be found on `Wikipedia <https://en.wikipedia.org/wiki/List_of_tz_database_time_zones>`_.

Suppose we want to start the model at 18:45 on the 21st July 2021. Use the `datetime` module 
to setup this date, and the set the start time, as follows:

>>> from datetime import datetime
>>> starttime = datetime(2021, 7, 21, 18, 45)
>>> domain.set_starttime(starttime)

Suppose we want to evolve until 19:00 on the 21st July 2021. Use `datetime` to setup
this `finaltime`:

>>> finaltime = datetime(2021, 7, 21, 19, 0)

And now evolve the model. Note the use of the `datetime = True` argument for the 
:meth:`print_timestepping_statisitics <Domain.print_timestepping_statistics>` procedure.

>>> for t in domain.evolve(yieldstep=300, finaltime=finaltime):
>>>   domain.print_timestepping_statistics(datetime=True)
DateTime: 2021-07-21 18:45:00+1000, steps=0 (0s)
DateTime: 2021-07-21 18:50:00+1000, delta t in [0.00832571, 0.01071429] (s), steps=31233 (10s)
DateTime: 2021-07-21 18:55:00+1000, delta t in [0.00959070, 0.00964172] (s), steps=31205 (10s)
DateTime: 2021-07-21 19:00:00+1000, delta t in [0.00959070, 0.00964172] (s), steps=31205 (10s)


Default zero time
-----------------

We use unix timestamp as our underlying absolute time. So `time = 0` actually corresponds to Jan 1st 1970 UTC. 

For instance going back to an earlier example which uses the default timezone (UTC) and 0 start time.

(Compare the output with `datetime = True` and `datetime = False` in the 
:meth:`print_timestepping_statistics <anuga.Domain.print_timestepping_statistics>` procedure.)

>>> import anuga
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


Note that the date is 1st Jan 1970, starting at time 0:00, incrementing by 1 sec and 
the UTC offset is +0000 (ie the timezone is UTC). 


Useful Domain methods
---------------------

.. autosummary::
   :toctree: generated
   
   Domain.evolve
   Domain.print_timestepping_statistics
   Domain.set_starttime
   Domain.set_timezone

