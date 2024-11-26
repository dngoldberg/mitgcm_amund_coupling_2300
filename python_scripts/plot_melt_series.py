import matplotlib.pyplot as plt; from mds import rdmds; from scipy.io import loadmat; import numpy as np
import sys

q,x,m=rdmds('surfDiag',np.inf,returnmeta=True)

nmonths = int(x[0]/25920)

mPig = np.empty(nmonths)
mThw = np.empty(nmonths)
mDot = np.empty(nmonths)
times = np.empty(nmonths)

for i in range(25920,x[0]+1,25920):
    q=rdmds('surfDiag',i)
    ind = int(i/25920-1)
    mPig[ind] = np.sum(-1*q[1,280:,:]*1250**2*31104000e-12)
    mThw[ind] = np.sum(-1*q[1,150:280,:]*1250**2*31104000e-12)
    mDot[ind] = np.sum(-1*q[1,:150,:]*1250**2*31104000e-12)
    times[ind] = 2006 + float(ind)/12.


plt.plot(times,mPig,label='pig')
plt.plot(times,mThw,label='thw')
plt.plot(times,mDot,label='dot')
plt.legend()
plt.grid()
plt.savefig('melt_series.png')
exit()
