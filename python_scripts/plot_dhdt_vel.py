import matplotlib.pyplot as plt; from mds import rdmds; from scipy.io import loadmat; import numpy as np
import sys
from IPython import embed

numyrs = float(sys.argv[1])
print(numyrs)

vavg1 = np.zeros((592,440));
uavg1 = np.zeros((592,440));
vavgn = np.zeros((592,440));
uavgn = np.zeros((592,440));
dhdt0 = np.fromfile('../dhdtcpom.bin',dtype='float64').byteswap().reshape((592,440))
dhdterr = np.fromfile('../dhdtcpomerr.bin',dtype='float64').byteswap().reshape((592,440))
rlow = rdmds('R_low_siinit')

for i in range(2,18,3):
    q = rdmds('land_ice',i+1);
    uavg1=uavg1 + 1./18.*q[0,:,:]
    vavg1=vavg1 + 1./18.*q[1,:,:]
    q = rdmds('land_ice',(numyrs-1)*18+i+1);
    uavgn=uavgn + 1./18.*q[0,:,:]
    vavgn=vavgn + 1./18.*q[1,:,:]

c = np.sqrt(uavg1**2+vavg1**2)
cn = np.sqrt(uavgn**2+vavgn**2)
h0=rdmds('H_streamiceinit')
q=rdmds('land_ice',numyrs*9.)
h1=q[2,:,:]
q=rdmds('land_ice',numyrs*18.)
h2=q[2,:,:]

haf1 = np.copy(h1);
haf1[rlow<0] = haf1[rlow<0] + 1027/917 * rlow[rlow<0]
haf2 = np.copy(h2);
haf2[rlow<0] = haf2[rlow<0] + 1027/917 * rlow[rlow<0]


plt.figure()
plt.subplot(2,2,1)
plt.imshow(c,origin='lower',vmin=0,vmax=5000); plt.colorbar();
plt.subplot(2,2,2)
dhdt = (h1-h0)/(numyrs/2.0);
dhdt[(dhdt0<-999) | (haf1<5)] = np.nan
#plt.contourf(dhdt,np.linspace(-8,8,31),cmap='bwr',extend='both'); plt.colorbar();
#plt.imshow(dhdt,origin='lower',cmap='bwr',vmin=-6,vmax=6); plt.colorbar();
plt.imshow(cn,origin='lower',vmin=0,vmax=5000); plt.colorbar();
plt.subplot(2,2,3)
dhdt = (h2-h0)/numyrs;
dhdt[(dhdt0<-999) | (haf2<5)] = np.nan
#plt.contourf(dhdt,np.linspace(-8,8,31),cmap='bwr',extend='both'); plt.colorbar();
plt.imshow(dhdt,origin='lower',cmap='bwr',vmin=-6,vmax=6); plt.colorbar();
plt.subplot(2,2,4)
dhdt0[(dhdt0<-999) | (haf2<5)] =np.nan
#plt.contourf(dhdt0,np.linspace(-8,8,31),cmap='bwr',extend='both'); plt.colorbar();
plt.imshow(dhdt0,origin='lower',cmap='bwr',vmin=-6,vmax=6); plt.colorbar();

#plt.close('all')

#plt.subplot(1,2,1)
#plt.imshow(dhdt,origin='lower',cmap='bwr',vmin=-6,vmax=6); plt.colorbar();
#plt.subplot(1,2,2)
#plt.imshow(dhdt-dhdt0,origin='lower',cmap='bwr',vmin=-6,vmax=6); plt.colorbar();

plt.savefig('dhdterr.png')
