
import time


import open_channel as model

domain = model.get_domain()

finaltime = 100.0
yieldstep = 1.0

model.animate_domain(domain, yieldstep, finaltime)
print "finished"
model.plot_domain(domain)
