
import time


import culvert as model

domain = model.get_domain()

finaltime = 100.0
yieldstep = 0.1

model.animate_domain(domain, yieldstep, finaltime)
print "finished"
model.plot_domain(domain)
