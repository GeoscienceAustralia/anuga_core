"""
Script to measure how long pmesh spends doing various methods
"""
from mesh import *
from pmesh import *
import time


draw = Draw()
#draw.run()
n = 400
times = []
tinitial = time.time()
times.append(("initial",time.time()-tinitial))
#------------------------------------------
mesh = Mesh()
id = 0
for i in range(n):
    for j in range(n):
       v = mesh.addUserVertex(i,j)
       v.guiID = id
       id += 1
v1 = mesh.addUserVertex(-1,-1)
v1.guiID = id
id += 1
v2 = mesh.addUserVertex(-1,n)
v2.guiID = id
id += 1
v3 = mesh.addUserVertex(n,n)
v3.guiID = id
id += 1
v4 = mesh.addUserVertex(n,-1)
v4.guiID = id
id += 1
mesh.addUserSegment(v1,v2)
mesh.addUserSegment(v2,v3)
mesh.addUserSegment(v3,v4)
mesh.addUserSegment(v4,v1)
times.append(("user_outline_created",time.time() - tinitial ))
#------------------------------------------
#mesh.auto_segment()
#times.append(("mesh_auto_segmented",time.time() - tinitial ))
#------------------------------------------
mesh.generateMesh(mode = "Q",maxArea = 1)
times.append(("mesh_generated",time.time() - tinitial ))
#------------------------------------------
draw.mesh = mesh
draw.selectZoom(1.0)
times.append(("zoomed",time.time() - tinitial ))
#------------------------------------------


#---------------------
print "N is ", n
for time in times:
    print "%s  %0.12f" %(time[0],  time[1])


#Results - mesh.py ver   1.84 	       	1.85	
# N	400	
# initial			0		0
# user_outline_created	1.467999935	1.609999895
# mesh_generated       	21.70300007	22.3440001
#zoomed			32.81299996	33.5940001



