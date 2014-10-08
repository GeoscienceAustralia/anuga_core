"""
    Quick plot of the Yeh-Petroff dam break outputs
"""
#---------------
# Import Modules
#---------------
import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
from anuga.utilities import plot_utils as util
from numpy import sqrt, sin, cos, pi
import csv, pprint
from numpy import zeros


#################Opening the reference solution################
infile = open('for_exp.txt', 'r')
table = []
for row in csv.reader(infile, delimiter="\t"):
    table.append(row)
infile.close()

tF = zeros((len(table),2)) #t for time, F for Force
for r in range(len(table)):
    if len(table[r])==2:
        for c in range(len(table[r])):        
            tF[r][c] = float(table[r][c])
    elif len(table[r])==3:
        for c in range(len(table[r])):        
            if c==0:
                tF[r][0] = float(table[r][c])
            elif c==2:
                tF[r][1] = float(table[r][c])
            #else:
            #    print "table[r][c]=",table[r][c]
    else:
        STOP
        
t = zeros(len(tF))        
F = zeros(len(tF))
for r in range(len(tF)):
    for c in range(len(tF[0])):
        if c==0:
            t[r] = tF[r][c]
        else:
            F[r] = tF[r][c]


indicesL = [21330,
21334,
21338,
21342,
21346,
21350,
21354,
21358,
21362,
21366,
21370,
21374]
##
##indicesR = [24988,
##24992,
##24996,
##25000,
##25004,
##25008,
##25012,
##25016,
##25020,
##25024,
##25028,
##25032]


##indicesL = [
##21086,
##21090,
##21094,
##21098,
##21102,
##21106,
##21110,
##21114,
##21118,
##21122,
##21126,
##21130
##    ]

indicesR = [
25232,
25236,
25240,
25244,
25248,
25252,
25256,
25260,
25264,
25268,
25272,
25276
    ]


# Get the sww file
p = util.get_output('dam_break.sww')
p2=util.get_centroids(p)


N = len(p2.time)
Sh = zeros(N)
for i in range(N):
    for k in range(len(indicesL)):
        Sh[i] += p2.stage[i,indicesL[k]]**2
        Sh[i] -= p2.stage[i,indicesR[k]]**2

fudge_factor = 1.0
length = 0.12
g = 9.8
rho_w = 1024
average_pressure = 0.5*g*Sh/len(indicesL)*rho_w

force = fudge_factor*average_pressure*length
    

#Plot stages
pyplot.clf()
pyplot.plot(t,F,'bo', label='experimental')
pyplot.plot(p2.time+0.1,force,'k-', label='numerical')
pyplot.legend(loc='best')
pyplot.title('Force acting to the column over time')
pyplot.xlabel('Time')
pyplot.ylabel('Force')
pyplot.savefig('force.png')
#pyplot.show()


