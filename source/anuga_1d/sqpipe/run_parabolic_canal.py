
import time


import parabolic_canal as model

domain = model.get_domain()
finaltime = 1000.0
yieldstep = 10.0


model.animate_domain(domain, yieldstep, finaltime)
print "finished"
model.plot_domain(domain)
