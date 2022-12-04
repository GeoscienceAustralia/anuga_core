"""
Script to measure how long pmesh spends doing various methods
"""


from builtins import range
from .mesh import *
from anuga.pmesh import *
import time
 
def mem_usage():
    '''
    returns the rss.

  RSS  The total amount of physical memory used by  the  task,  in  kilo-
            bytes,  is  shown  here.  For ELF processes used library pages are
            counted here, for a.out processes not.
            
    Only works on nix systems.
    '''
    import string
    p=os.popen('ps uwp %s'%os.getpid()) 
    lines=p.readlines()
    #print "lines", lines
    status=p.close() 
    if status or len(lines)!=2 or sys.platform == 'win32': 
        return None 
    return int(string.split(lines[1])[4]) 



#draw = Draw()
#draw.run()
n = 2
maxArea = 0.00005
#maxArea = 0.1
times = []
tinitial = time.time()
mem_initial = mem_usage()
print("mem_initial", mem_initial)
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
mem_now =mem_usage()
#print "mem_now", mem_now
if mem_now is not None:
    mem = mem_now - mem_initial
else:
    mem = 0.0
times.append(("user_outline_created",time.time() - tinitial, mem ))
#------------------------------------------
mesh.generateMesh(mode = "Q",maxArea = maxArea)
mem_now =mem_usage()
#print "mem_now", mem_now
if mem_now is not None:
    mem = mem_now - mem_initial
else:
    mem = 0.0
times.append(("mesh_generated",time.time()- tinitial - times[0][1], mem))
#------------------------------------------
mesh.export_mesh_file("dump.msh")
mem_now =mem_usage()
#print "mem_now", mem_now
if mem_now is not None:
    mem = mem_now - mem_initial
else:
    mem = 0.0
times.append(("export_mesh_file",time.time()- tinitial - times[1][1], mem))

#---------------------
print("Number of user verts. ", n)
print("maxArea",maxArea)
print("funtion     time   memory usage, cumulative, for nix machines")
for time in times:
    print("%s   %0.2f   %0.2f" %(time[0],  time[1], time[2]))

"""
#Results - mesh.py ver   1.84 	       	1.85	
# N	400	
# initial			0		0
# user_outline_created	1.467999935	1.609999895
# mesh_generated       	21.70300007	22.3440001
#zoomed			32.81299996	33.5940001

# results on dsg's pc dell precision 390.
# anuga version 4897
Number of user verts.  2
maxArea 0.0001
user_outline_created   0.00
mesh_generated   4.92
export_mesh_file   16.30

Number of user verts.  2
maxArea 0.0001
user_outline_created   0.00
mesh_generated   4.92
export_mesh_file   16.30

Number of user verts.  2
maxArea 5e-005
user_outline_created   0.00
mesh_generated   15.53
export_mesh_file   41.39

version 4903
mem_initial None
Number of user verts.  2
maxArea 5e-005
funtion     time   memory usage, cumulative, for nix machines
user_outline_created   0.49   0.00
mesh_generated   1.94   0.00
export_mesh_file   2.75   0.00

# anuga version 4897 
Results of tornado head node
Number of user verts.  2
maxArea 5e-05
funtion     time   memory usage, cumulative, for nix machines
user_outline_created   0.14   32.00
mesh_generated   18.81   378304.00
export_mesh_file   17.79   380620.00

version 4903
mem_initial 77816
Number of user verts.  2
maxArea 5e-05
funtion     time   memory usage, cumulative, for nix machines
user_outline_created   0.15   32.00
mesh_generated   1.26   27224.00
export_mesh_file   0.82   27224.00

"""
