"""View results of runup.py
"""
#---------------
# Import Modules
#---------------
import anuga
import numpy
import scipy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
#import util # Routines to read in and work with ANUGA output
#from bal_and import plot_utils as util
from anuga.utilities import plot_utils as util

p2=util.get_output('runup_riverwall.sww', minimum_allowed_height=1.0e-03)
p=util.get_centroids(p2, velocity_extrapolation=True)

#p=util.get_output('runup_v2.sww', minimum_allowed_height=1.0e-03)

#------------------
# Select line
#------------------
#py_central=p.y[scipy.argmin(abs(p.y-50.))]
#v=(p.y==p.y[py_central])
v=util.near_transect(p, [0., 50.], [100., 50.], 10.)
v=v[0]
v=v[p.x[v].argsort()]

#--------------------
# Make plot animation
#--------------------
pyplot.close() #If the plot is open, there will be problems

if False:
    pyplot.ion()
    line, = pyplot.plot( (p.x[v].min(),p.x[v].max()) ,(p.xvel[:,v].min(),p.xvel[:,v].max() ) )
    for i in range(p.xmom.shape[0]):
        line.set_xdata(p.x[v])
        line.set_ydata(p.xvel[i,v])
        pyplot.draw()
        pyplot.plot( (0,1),(0,0), 'r' )
        pyplot.title(str(i)) # : velocity does not converge to zero' )
        pyplot.xlabel('x')
        pyplot.ylabel('Velocity (m/s)')

    pyplot.savefig('runup_x_velocities.png')

## Get reference stage points
seaLev=(abs(p.x-51.)+abs(p.y-50.)).argmin() # Index near levee on sea side
landLev=(abs(p.x-49.)+abs(p.y-50.)).argmin() # Index near levee on land side
heightDiff=p.stage[:,seaLev]-p.stage[:,landLev]
## Get indices on 'landward' side of riverwall
landInds=(p.x<50.).nonzero()
# Get volume on landward side of riverwall
landVol=util.water_volume(p=p2,p2=p,subset=landInds)
l=len(landVol)
# Compute rate of influx of water
d_landVol_dt=(landVol[1:l]-landVol[0:(l-1)])/(p.time[1:l]-p.time[0:(l-1)])

# Estimated flux rate using C*L*H^(3/2)
# '1' indicates the region where the wall = -0.2 
# '2' indicates the region where the wall = 0.0
refStage1=p.stage[:,seaLev]
refStage2=p.stage[:,landLev]
w1=-0.2
w2=0.0
H1_sea =(refStage1-w1)*(refStage1> w1) 
H2_sea =(refStage1-w2)*(refStage1> w2)
H1_land=(refStage2-w1)*(refStage2> w1)
H2_land=(refStage2-w2)*(refStage2> w2)
# Compute 'uncorrected' flow rates, then villemonte correction, then add
# Note the '_land' discharge values, which are as suggested by Villemonte for
# sharp crested weirs. See the paper in my hydraulics references

ft_to_m=0.3048 # 1 foot = ft_to_m metres
#C=3.2 # Standard constant for sharp crested weir from e.g. HecRas, in units ft^(0.5)/s
#C=C*(ft_to_m)**0.5 # Convert to m^0.5/s
C=2./3.*(9.81*2./3.)**0.5 # Another standard coefficient.
L1=9.0 # Lengths of the w1/w2 regions
L2=91.
def simple_weir(H, C, L):
    return L*C*H**(1.5)

InRate1= simple_weir(H1_sea, C, L1) #L1*C*(H1_sea)**1.5
InRate1_land=simple_weir(H1_land,C,L1) #L1*C*(H1_land)**1.5
InRate2= simple_weir(H2_sea, C, L2) #L2*C*(H2_sea)**1.5
InRate2_land= simple_weir(H2_land, C, L2) #L2*C*(H2_land)**1.5

def Vil_cor(InRate_land, InRate, n=1.0):
    """
    Apply Villemonte correction to InRate
    If InRate_land>InRate, then the submergence_adjusted discharge is computed from InRate_land. In that case it is negative
    Otherwise, the submergence_adjusted discharge is computed from InRate, and is positive
     
    """
    #Fact=InRate_land/(InRate+1.0e-10)*(InRate_land<=InRate) + 1.0*(InRate_land>InRate)
    Fact=InRate_land/(InRate+1.0e-10)*(InRate_land<=InRate) + InRate/(InRate_land+1.0e-10)*(InRate<InRate_land)
    assert Fact.max()<=1.0, 'Error Fact > 1'
    out=InRate*(InRate_land<=InRate)*(1.0 - Fact**n)**0.385 -\
        InRate_land*(InRate_land>InRate)*(1.0 - Fact**n)**0.385
    return out 

#n=1.0
Villemonte_correction1=Vil_cor(InRate1_land, InRate1) #(1.0 - (InRate1_land/(InRate1+1.0e-06))**n)**0.385
Villemonte_correction2=Vil_cor(InRate2_land, InRate2) #(1.0 - (InRate2_land/(InRate2+1.0e-06))**n)**0.385


#InRate=InRate1*Villemonte_correction1+InRate2*Villemonte_correction2
InRate=Villemonte_correction1+Villemonte_correction2

pyplot.close()
pyplot.figure(figsize=(12,8))
lab2='ANUGA weir computations [yieldstep-averaged]'
pyplot.plot(0.5*(p.time[0:(l-1)]+p.time[1:l]), d_landVol_dt, color='red', label=lab2)
lab1='Simple weir with Villemonte submergence correction'
pyplot.plot(p.time[1:l], InRate[1:l],color='blue', label=lab1)
#pyplot.clf()
#pyplot.plot(InRate[1:l]/ d_landVol_dt)
#pyplot.ylim((0.,3.))

# Try the Hager approach
n=1.
c=1.
h=p.height[:, seaLev] #refStage1-p.elev[seaLev] # Depth
H=p.height[:, seaLev] #refStage1-p.elev[seaLev] # Energy head ~ depth
w=w1-p.elev[seaLev] # Weir height - bed elevation on channel side
def Hager_weir(h, H, w, n, c, L1):
    # Terms in Hager, eq 12
    y=h/(H+1.0e-30)
    W=w/(H+1.0e-30)
    # Compute first terms in inrate (which have dimension)
    InRate_H1=3./5.*n*c*(9.8*H**3)**0.5
    # Adjust with further dimensionless terms
    InRate_H1=InRate_H1*((y-W)*(y>W))**(1.5)*( (1.-W)/(3.-2.*y-W))**0.5*(1.-0.)
    # Multiply by length
    InRate_H1=InRate_H1*L1
    return(InRate_H1)

# Compute 'raw' haegar weir outflow
InRate_H1=Hager_weir(h,H,w,n,c,L1)
InRate_H1_land=Hager_weir( p.height[:, landLev], p.height[:, landLev], w1-p.elev[landLev], n, c, L1)
Vil_H1_cor=Vil_cor(InRate_H1_land,InRate_H1)
#InRate_H1=InRate_H1*Vil_H1_cor
InRate_H1=Vil_H1_cor
lab3='Hager weir with Villemonte submergence correction'
pyplot.plot(p.time[1:l], InRate_H1[1:l],color='green', label=lab3)

# Another correction, eqn 5, Fritz and Hager
yt=H1_land/(H1_sea+1.0e-30)
weir_flat_width=0.
eps=(H1_sea)/(H1_sea+weir_flat_width+1.0e-30)
yl=0.85-0.5*eps
def Hager_cor(yt, yl):
    # 0<=yt<=1
    yt = yt*(yt<=1.0)*(yt>=0.) + 1.0*(yt>1.0)
    # 
    Yt=(yt-yl)*(yt>yl)/(1.-yl)
    n=6.0
    cor=(1.-Yt)**(1./n)
    return cor

InRate_H1Raw=Hager_weir(h,H,w,n,c,L1)
HagerCor1= Hager_cor(yt,yl)
InRate_H1B = InRate_H1Raw*HagerCor1
lab4='Hager weir with Fritz and Hager submergence correction [positive flux only]'
pyplot.plot(p.time[1:l], InRate_H1B[1:l],color='pink', label=lab4)
pyplot.legend( (lab2, lab1, lab3, lab4), loc='upper right', prop={'size':10},
                title='Non-ANUGA computations assume a spatially constant \n headwater/tailwater stage, and use submergence corrections \n from the engineering literature')
pyplot.title('Flux over riverwall (m^3/s)',fontsize=20)
pyplot.xlabel('Time (s)')
pyplot.ylabel('Flux (m^3/s)')
pyplot.savefig('Fluxes.png')
pyplot.close()

pyplot.figure(figsize=(12,8))
lab1='Stage minus wall top on "sea" side of riverwall'
pyplot.plot(p.time, H1_sea, label=lab1)
lab2='Stage minus wall top on "land" side of riverwall'
pyplot.plot(p.time, H1_land, label=lab2)
#pyplot.plot(p.time[[0,l-1]], 0.*p.time[[0,l-1]]-0.2,'-')
pyplot.legend((lab1,lab2), loc='lower right')
pyplot.title('Stage above wall top at 2 points near the overflowing riverwall',fontsize=20)
pyplot.savefig('Stage.png')
pyplot.close()
# NOTE: HecRas reduces submerged weir flow OVER BRIDGES according to the
# submergence ratio = (depth of water above min weir elevation on downstream
# side)/(height of energy grade line above min weir elevation on upstream
# side). 
# The discharge reduction factor is > 0.9 for submergence ratio < 0.9, and gets
# small quickly after that. It is nearly 1.  for ratios around 0.75. Greater
# reductions are suggested by this version of Villemonte's equations. 
# On the other hand, the WA roads manual indicates that for ratios < 0.75, the
# flow becomes critical over the bridge, so there is no issue with submergence.

# HOWEVER note presentation 'Weir_submergence_issues.pdf' which discusses different relations
