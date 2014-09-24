"""
    Quick plot of the dam break outputs
"""
import anuga.utilities.plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
import csv, pprint
from numpy import zeros

#--------------------------------------------
# Recall the reference solution
# The reference solution is computed using 10^4 cells
# with FVM on varying bottom and width
# See: Roberts and Wilson, ANZIAM Journal (CTAC2010).
#--------------------------------------------
#################Opening the reference solution################
infile = open('Ver_numerical_2.000000.csv', 'r')
table = []
for row in csv.reader(infile):
    table.append(row)
infile.close()
for r in range(len(table)):
    for c in range(len(table[0])):
        table[r][c] = float(table[r][c])
Vertices = zeros(len(table))        
VerW = zeros(len(table))
VerP = zeros(len(table))
VerZ = zeros(len(table))
VerH = zeros(len(table))
VerU = zeros(len(table))
for r in range(len(table)):
    for c in range(len(table[0])):
        if c==0:
            Vertices[r] = table[r][c]             
        elif c==1:
            VerW[r] = table[r][c]
        elif c==2:
            VerP[r] = table[r][c]
        elif c==3:
            VerZ[r] = table[r][c]
        elif c==4:
            VerH[r] = table[r][c]
        elif c==5:
            VerU[r] = table[r][c]
Vertices_left = -1.0*Vertices[::-1]
VerW_left = VerW[::-1]
VerP_left = -1.0*VerP[::-1]
VerZ_left = VerZ[::-1]
VerH_left = VerH[::-1]
VerU_left = -1.0*VerU[::-1]

p_st = util.get_output('radial_dam.sww')
p2_st=util.get_centroids(p_st)

v = p2_st.y[79598]
v2=(p2_st.y==v)

#Plot stages
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.stage[-1,v2],'b.-', label='numerical stage')
pyplot.plot(Vertices, VerW,'r-', label='reference stage')
pyplot.plot(Vertices_left, VerW_left,'r-')
pyplot.title('Stage at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Radial position')
pyplot.ylabel('Stage')
pyplot.savefig('stage_plot.png')
#pyplot.show()


#Plot rmomentum
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xmom[-1,v2], 'b.-', label='numerical')
pyplot.plot(Vertices, VerP,'r-', label='reference')
pyplot.plot(Vertices_left, VerP_left,'r-')
pyplot.title('Radial momentum at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Radial position')
pyplot.ylabel('Radial momentum')
pyplot.savefig('rmom_plot.png')


#Plot rvelocities
pyplot.clf()
pyplot.plot(p2_st.x[v2], p2_st.xvel[-1,v2], 'b.-', label='numerical')
pyplot.plot(Vertices, VerU,'r-', label='reference')
pyplot.plot(Vertices_left, VerU_left,'r-')
pyplot.title('Radial velocity at an instant in time')
pyplot.legend(loc='best')
pyplot.xlabel('Radial position')
pyplot.ylabel('Radial velocity')
pyplot.savefig('rvel_plot.png')
