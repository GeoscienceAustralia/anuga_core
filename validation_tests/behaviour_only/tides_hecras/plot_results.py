import numpy
from anuga.utilities import plot_utils as util
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot

######################################################
# Get ANUGA
p = util.get_output('channel_floodplain1.sww', 0.001)
pc=util.get_centroids(p, velocity_extrapolation=True)

# Indices in the central channel areas
v = (pc.x>10.0)*(pc.x<30.0)

######################################################
# Get hecras info
rasFile='hecras_tides_case/gauges.csv'
rasGauges=numpy.genfromtxt(rasFile,skip_header=3,delimiter=",")

# Output at 1 minute intervals
rasTime=numpy.linspace(0., 60.*(rasGauges.shape[0]-1), rasGauges.shape[0])

# Get station information for hecras
rasFileO=open(rasFile)
rasStations=rasFileO.readline().split(',')
rasStations[-1]=rasStations[-1].replace('\r\n','')
rasFileO.close()

######################################################

def get_corresponding_series(reach, station):
    """

    Convenience function to get ANUGA/hec-ras series at the same station

    """

    if(reach=='REACH1'):
        anuga_x=20.
    #elif(reach=='LEFT'):
    #    anuga_x=35.
    #elif(reach=='RIGHT'):
    #    anuga_x=5.
    else:
        raise Exception, 'reach not recognized'


    # Get station string in hecras gauges
    if(station>0. and station<1000.):
        station_str=str(int(station))+'.*'
    else:
        station_str=str(int(station))

    ras_string='RIVER1 '+reach+' '+station_str

    ras_inds=rasStations.index(ras_string)

    #if(len(ras_inds) is not 1):
    #    print ras_inds
    #    raise Exception, 'ras_inds has wrong length'

    ras_data=numpy.vstack([rasTime, rasGauges[:,ras_inds]]).transpose()

    anuga_index=((pc.x-anuga_x)**2+(pc.y-(1000.-station))**2).argmin()

    anuga_data=numpy.vstack([pc.time, pc.stage[:,anuga_index]]).transpose()

    return [ras_data, anuga_data]

def compare_reach(reach):
    stations=[0, 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000]
    stations.reverse() # This looks nicer
    colz=['red', 'blue', 'green', 'purple', 'orange', 'brown', 'skyblue', 'pink', 'yellow', 'grey']

    for i, station in enumerate(stations):
        try:
            x=get_corresponding_series(reach, station)
            pyplot.plot(x[0][:,0],x[0][:,1]+station/1000.,'-', label=str(station), color=colz[i], linewidth=2)
            # NOTE: vertical datum in ANUGA is offset from HECRAS by 1.0m
            pyplot.plot(x[1][:,0],x[1][:,1]+station/1000.+1.0, '--', color='black',linewidth=1.0)
        except:
            msg = 'Missing reach/station '+ reach + '/'+str(station)
            print msg

        pyplot.xlim((0, 15000.))
        pyplot.xlabel('Time (s)',fontsize=20)
        pyplot.ylabel('Stage with vertical offset (m)',fontsize=20)
        pyplot.title('Stage in ANUGA (dashed black) and HECRAS (solid colors). Stations are offset vertically (by station/1000.)',fontsize=20)
        pyplot.legend(loc=2)
   
    
pyplot.figure(figsize=(20,8))
compare_reach('REACH1')
pyplot.savefig('CENTRAL_CHANNEL.png') 
#pyplot.clf()
#compare_reach('RIGHT')
#pyplot.savefig('RIGHT_REACH.png') 
#pyplot.clf()
#compare_reach('LEFT')
#pyplot.savefig('LEFT_REACH.png') 

