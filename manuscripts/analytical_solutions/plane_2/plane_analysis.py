import sys
oldsysstdout = sys.stdout

from matplotlib import pyplot as plt
from anuga.file import netcdf
import matplotlib.tri as tri
import scipy.interpolate as sp
from glob import glob
import numpy as np

sys.stdout = oldsysstdout
from scipy.signal import savgol_filter as filter



swwFile = 'plane.sww'



fid = netcdf.NetCDFFile(swwFile)
x = fid.variables['x'][:]
y = fid.variables['y'][:]
elev = fid.variables['elevation'][:]
stage = fid.variables['stage'][:]
conc = fid.variables['concentration'][:]

depth = stage - elev

xmom = fid.variables['xmomentum'][:]
xvel = xmom / (depth + 0.000001)

fid.close()


y2 = (2<y)
y8 = y < 8
indices = y2 * y8

        
triang = tri.Triangulation(x,y)

        
# variables        
g = 9.8
rho_w = 1000.
tau_c = 0.126126
phi = 0.3
vs = 0.0138099680617

Co = 0.005       
Ke = 0.2e-6/(tau_c**0.5)

# d_star
filenames = glob('dstar/dstar*')

ds = np.load(filenames[0])
xx = np.linspace(0.5,14,45)

numbers = []
for f in filenames[1:]:
    a = f.replace('dstar/dstar','')
    b = a.replace('.0.npy','')
    numbers.append(int(b))
    
d_star = np.zeros((len(numbers),len(xx)))


fl = sp.interp1d(x[indices], elev[0][indices], kind='linear')
Dxo = filter(fl(xx), 7, 3)

for i in range(1,len(filenames)):

    ds = np.load(filenames[i])
    fl = sp.interp1d(ds[0], ds[1], kind='linear')
    d_star[i-1] = filter(fl(xx), 7, 3)
    

Cx = np.zeros((len(numbers),len(xx)))
Dx = np.zeros((len(numbers),len(xx)))


for i in numbers:
    
    Co = conc[i][x < 1].max()
    
    j = i - 2

    d_star_u = d_star[j]
    
    fl = sp.interp1d(x[indices], depth[i][indices], kind='linear')
    dxx = filter(fl(xx), 7, 3)
    
    fl = sp.interp1d(x[indices], xmom[i][indices], kind='linear')
    qxx = filter(fl(xx), 7, 3)
    
    
    S = (-np.diff(elev[i][y == 5])/np.diff(x[y ==5]))
    S_ = np.hstack((S,S[-1]))

    fl = sp.interp1d(x[y ==5], S_, kind='linear')
    slope = filter(fl(xx), 7, 3)
    

    u_star = np.sqrt(g * slope * dxx)
    tau = rho_w * u_star**2
    dot_E = Ke * (tau - tau_c)

    Exx = dot_E

    

    Cx[j] = (Co - Exx/(d_star_u * vs)) * np.exp(- d_star_u * vs * xx / qxx) + Exx/(d_star_u * vs)    
 
    Dx[j] = i / (1 - phi) * (d_star_u * vs * Co - Exx) * np.exp(- d_star_u * vs * xx / qxx) + Dxo


# check data

Cx_m = np.zeros((len(numbers),len(xx)))
Dx_m = np.zeros((len(numbers),len(xx)))

for i in numbers:
    
    j = i - 2
    
    fl = sp.interp1d(x[indices], elev[i][indices], kind='linear')
    Dx_m[j] = filter(fl(xx), 7, 3)
    
    fl = sp.interp1d(x[indices], conc[i][indices], kind='linear')
    Cx_m[j] = filter(fl(xx), 7, 3)

    
# root mean square errors

rmse_Cx = np.zeros((len(Cx_m),len(xx)))
rmse_Dx = np.zeros((len(Cx_m),len(xx)))

for i in range(len(Cx_m)):

    rmse_Cx[i] = np.sqrt((Cx_m[i] - Cx[i])**2)
    rmse_Dx[i] = np.sqrt((Dx_m[i] - Dx[i])**2)
    
    
# save outputs

np.save('rmse_concentration.npy', [xx,rmse_Cx])
np.save('rmse_bed.npy', [xx,rmse_Dx])
np.save('analytical_concentration.npy',Cx)
np.save('analytical_bed.npy',Dx)
np.save('model_concentration.npy',Cx)
np.save('model_bed.npy',Cx)

plt.figure()
plt.plot(xx, Cx_m[-1], xx, Cx[-1])
plt.legend(['model','analytical'])
plt.ylabel('Concentration')
plt.savefig('fig_concentration.png')
plt.close()

plt.figure()
plt.plot(xx, Dx_m[-1], xx, Dx[-1])
plt.legend(['model','analytical'])
plt.ylabel('Bed elevation')
plt.savefig('fig_bed.png')
plt.close()

plt.figure()
plt.plot(xx, rmse_Dx.T)
plt.ylabel('RMSE bed elevation')
plt.savefig('fig_rmse_bed.png')
plt.close()

plt.figure()
plt.plot(xx, rmse_Cx.T)
plt.ylabel('RMSE concentration')
plt.savefig('fig_rmse_concentration.png')
plt.close()
