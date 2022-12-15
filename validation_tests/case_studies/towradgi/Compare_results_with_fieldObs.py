"""

Compare the Towradgi model runs with various field observations

"""
import scipy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as pyplot
from anuga.utilities import plot_utils as util
import gdal


swwdir='MODEL_OUTPUTS/'
swwname='Towradgi_historic_flood.sww' 

p=util.get_output(swwdir+swwname)
p2=util.get_centroids(p,velocity_extrapolation=True)
floodLevels=scipy.genfromtxt('Validation/historic_1998_flood_levels_towradgi_ck.csv', delimiter=',',skip_header=1)
pioneerLevel=scipy.genfromtxt('Validation/pioneer_timeseries.txt', 
                              skip_header=1)

# Extract modelled peak at the coordinates in floodLevels
modelled_ind=floodLevels[:,0]*0
modelled_level=floodLevels[:,0]*0
for i in range(len(modelled_level)):
    # Convert the physical coordinate to ANUGA coordinate system
    pt_x=floodLevels[i,0]-p.xllcorner
    pt_y=floodLevels[i,1]-p.yllcorner

    # Find the nearest index of the physical coordinate in the centroids
    myInd=( (p2.x-pt_x)**2 + (p2.y-pt_y)**2).argmin()
    modelled_ind[i]=myInd
    modelled_level[i]=p2.stage[:,myInd].max()

    if(i==0):
        pyplot.clf()
        pyplot.plot((pioneerLevel[:,0])/60., pioneerLevel[:,2],color='black', label="Recorded Stage")
        pyplot.plot(p2.time/3600.,p2.stage[:, myInd],color='blue', label="ANUGA")
        pyplot.legend(bbox_to_anchor=(1.05, 1), borderaxespad=0.)
        pyplot.xlabel('Time (Hours)')
        pyplot.ylabel('Stage (mAHD)')
        pyplot.title('Pioneer Road Bridge Stage')
        pyplot.savefig('Pioneer_Bridge_Stage_'+swwname[:-4]+'.png')
    

pyplot.clf()
pyplot.scatter(floodLevels[:,3],modelled_level)
print (floodLevels[:,2],modelled_level)
pyplot.xlabel('Observed flood level (mAHD)')
pyplot.ylabel('Modelled flood level (mAHD)')
pyplot.title('Modelled vs Observed Flood Levels\nTowradgi Creek Historic Flood 1998')
pyplot.savefig('Modelled_vs_Observed_peakStage_'+swwname[:-4]+'.png')

pyplot.clf()
pyplot.hist(floodLevels[:,3]-modelled_level)
pyplot.xlabel('Observed minus Modelled flood level')
pyplot.title('Difference in Modelled and Observed flood levels\nTowradgi Creek Historic Flood 1998')
pyplot.savefig('Error_peakstage_'+swwname[:-4]+'.png')


# Make a bunch of GIS outputs
try:
    tif_outdir='OUTPUT_TIFS'
    CellSize=5.0
    print ('Making tifs')
    util.Make_Geotif(swwdir+swwname,
                      ['depth','velocity','depthIntegratedVelocity','elevation', 'friction'],'max',
                      CellSize=CellSize,EPSG_CODE=32756,output_dir=tif_outdir)
    print ('Made tifs')
except:
    print ('Cannot make GIS plot -- perhaps GDAL etc are not installed?')

    
# Plot depth raster with discrepency between model and data
depthFile=tif_outdir+'/Towradgi_historic_flood_depth_max.tif'
#myDepth=scipy.misc.imread(depthFile)
raster = gdal.Open(depthFile)
myDepth = scipy.array(raster.ReadAsArray())


X=scipy.arange(p.xllcorner, p.xllcorner+myDepth.shape[1]*CellSize, CellSize)
Y=scipy.arange(p.yllcorner, p.yllcorner+myDepth.shape[0]*CellSize, CellSize)
X,Y=scipy.meshgrid(X,Y)
pyplot.clf()
pyplot.figure(figsize=(12,6))
pyplot.plot([X.min(),X.max()],[Y.min(),Y.max()],' ')
pyplot.imshow(scipy.flipud(myDepth),extent=[X.min(),X.max(),Y.min(),Y.max()],origin='lower',cmap=pyplot.get_cmap('Greys'))
pyplot.gca().set_aspect('equal')
pyplot.colorbar(orientation='horizontal').set_label('Peak Depth in model (m)')
er1=floodLevels[:,3]-modelled_level
pyplot.scatter(floodLevels[:,0], floodLevels[:,1], c=er1,s=20,cmap=pyplot.get_cmap('spectral'))
pyplot.colorbar().set_label(label='Field observation - Modelled Peak Stage (m)')
pyplot.xlim([p.x.min()+p.xllcorner,p.x.max()+p.xllcorner])
pyplot.ylim([p.y.min()+p.yllcorner,p.y.max()+p.yllcorner])
pyplot.savefig('Spatial_Depth_and_Error_'+swwname[:-4]+'.png')

