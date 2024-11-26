import numpy as np
import matplotlib.pyplot as plt
from mds import rdmds
from scipy.io import loadmat
from IPython import embed
import os
from scipy.interpolate import NearestNDInterpolator
from matplotlib.path import Path
import pickle as pkl

rho = 917.
rhow = 1027.


sliding = ['Snap', 'TC']
meltP = [1,3]

iceParmName = []

for i in meltP:
 iceParmName.append('iceParm' + str(i))
# Load meshcoords data, assuming it's a .mat file


q = loadmat('../../../MITgcm_forinput/amund_couple/mitgcm_amund_coupling/input/start_2009_input/meshcoords.mat')
x_mesh_mid = q['x_mesh_mid'][0]
y_mesh_mid = q['y_mesh_mid'][0]
x_mesh_mid, y_mesh_mid = np.meshgrid(x_mesh_mid, y_mesh_mid)

x_mesh_oce_mid = q['x_mesh_oce_mid'][0]
y_mesh_oce_mid = q['y_mesh_oce_mid'][0]
x_mesh_oce_mid, y_mesh_oce_mid = np.meshgrid(x_mesh_oce_mid, y_mesh_oce_mid)


q = loadmat('../../../MITgcm_forinput/amund_couple/mitgcm_amund_coupling/plotting_scripts/outlines.mat')
xpig = q['xpig'].flatten().tolist()
ypig = q['ypig'].flatten().tolist()
xthw = q['xthw'].flatten().tolist()
ythw = q['ythw'].flatten().tolist()
xsm = q['xsm'].flatten().tolist()
ysm = q['ysm'].flatten().tolist()
#q=loadmat('../../../MITgcm_forinput/thwaites_trans_cal/input_tc/bounds1.mat')
#xthw = q['x'].flatten().tolist()
#ythw = q['y'].flatten().tolist()




xgrid,ygrid = np.meshgrid(np.arange(np.min(x_mesh_mid)/1000,np.max(x_mesh_mid)/1000+.1,1),np.arange(np.min(y_mesh_mid)/1000,np.max(y_mesh_mid)/1000+.1,1))
nxgrid, nygrid = np.shape(xgrid)
xgrid, ygrid = xgrid.flatten(), ygrid.flatten()
points =np.vstack((xgrid,ygrid)).T

xgrido,ygrido = np.meshgrid(np.arange(np.min(x_mesh_oce_mid)/1000,np.max(x_mesh_oce_mid)/1000+.1,1),np.arange(np.min(y_mesh_oce_mid)/1000,np.max(y_mesh_oce_mid)/1000+.1,1))
nxgrid, nygrid = np.shape(xgrido)
xgrido, ygrido = xgrido.flatten(), ygrido.flatten()
points_oce =np.vstack((xgrido,ygrido)).T

for sh in ['pig','thw','sm']:
 exec('xs = ((np.array(x' + sh + ') - 500 + 500)/1000).tolist()') 
 exec('ys = ((np.array(y' + sh + ') - 500 + 500)/1000).tolist()') 
 poly = [(xs[i],ys[i]) for i in range(0,len(xs))]
 path = Path(poly)
 grid = path.contains_points(points).astype(float)
 NDint = NearestNDInterpolator(points,grid)
 exec('mask_' + sh + ' = NDint(x_mesh_mid/1000.,y_mesh_mid/1000.)')
 grid = path.contains_points(points_oce).astype(float)
 NDint = NearestNDInterpolator(points_oce,grid)
 exec('mask_oce_' + sh + ' = NDint(x_mesh_oce_mid/1000.,y_mesh_oce_mid/1000.)')

def getvaf(thick,area,bed,mask):

  haf = thick.copy();
  haf[mask!=1.0]=0.0
  haf[bed<0.0] = haf[bed<0.0] + bed[bed<0.0] * rhow/rho;
  haf[haf<0.0] = 0.0
  vaf = np.sum(area * haf * mask)/1e9
  return vaf

for islid in range(2):
 for PAS in [0,10,20,30]:
  for name in iceParmName:
#    PAS = 30
   slid = sliding[islid]
   if (os.path.exists(f'run_ice_2009_{slid}_{PAS}_{name}_coul')):
    print(slid)
    folder = f'run_ice_2009_{slid}_{PAS}_{name}_coul'
    ofolder = f'run_oce_2009_{slid}_{PAS}_{name}_coul'
    rA = rdmds(folder+ '/RAC')
    bed = rdmds(folder + '/R_low_siinit')
    MP = []
    MT = []
    MS = [] 
    MPdeep = []
    MTdeep = []
    MSdeep = [] 
    T = []
    VP = []
    VT = []
    VS = []
    TV = []

    q, x5, m = rdmds(folder+ '/land_ice', np.inf, rec=[5, 3, 4],returnmeta=True)
    print(type(q))
    mdiff = q[0,:, :] 
    g2 = q[2, :, :]

    if (x5[0]>33):
      num0 = x5[0]-36
    else:
      num0 = 0
    
    q2, x2, m = rdmds(folder+ '/land_ice', num0, rec=[5, 3, 4],returnmeta=True)
    mdiff2 = q2[0, :, :] 
    g22 = q2[2, :, :]

    q = rdmds(folder+ '/land_ice', 0, rec=[5, 0, 4, 3])
    mdiff0 = q[0, :, :] 
    g0 = q[2, :, :]
    msk = q[3, :, :]

    q3, x3, m = rdmds(ofolder+ '/surfDiag', np.inf, returnmeta=True)

    if (PAS>0):
      for i in range(5 * 311040, x3[0], 25920):
        i
        q = rdmds(ofolder+ '/surfDiag', i, rec=1)
        rsh = rdmds(ofolder+ '/surfDiag', i, rec=4)
        eta = rdmds(ofolder+ '/surfDiag', i, rec=0)
        draft = rsh + eta

        ml = -31104000 * q * mask_oce_pig
        MP.append(np.sum(np.sum(ml)) * (1250 ** 2))
        ml[draft>-500] = 0.
        MPdeep.append(np.sum(np.sum(ml)) * (1250 ** 2))

        ml = -31104000 * q * mask_oce_thw
        MT.append(np.sum(np.sum(ml)) * (1250 ** 2))
        ml[draft>-500] = 0.
        MTdeep.append(np.sum(np.sum(ml)) * (1250 ** 2))

        ml = -31104000 * q * mask_oce_sm
        MS.append(np.sum(np.sum(ml)) * (1250 ** 2))
        ml[draft>-500] = 0.
        MSdeep.append(np.sum(np.sum(ml)) * (1250 ** 2))

        T.append(2008 + i / 311040)
        print(i)


    import pandas as pd
    MP = pd.Series(MP).rolling(12,center=True,min_periods=1).mean().values
    MPdeep = pd.Series(MPdeep).rolling(12,center=True,min_periods=1).mean().values
    MT = pd.Series(MT).rolling(12,center=True,min_periods=1).mean().values
    MTdeep = pd.Series(MTdeep).rolling(12,center=True,min_periods=1).mean().values
    MS = pd.Series(MS).rolling(12,center=True,min_periods=1).mean().values
    MSdeep = pd.Series(MSdeep).rolling(12,center=True,min_periods=1).mean().values


    for i in range(0,x5[0]+1,3):
        print(i)
        q = rdmds(folder+ '/land_ice', i, rec=2)
        VP.append(getvaf(q,rA,bed,mask_pig))
        VT.append(getvaf(q,rA,bed,mask_thw))
        VS.append(getvaf(q,rA,bed,mask_sm))
        TV.append(i/36+2013)

    if(PAS>0):
     exec(f'MPig{slid}{PAS}{name} = MP')
     exec(f'MThw{slid}{PAS}{name} = MT')
     exec(f'MSm{slid}{PAS}{name} = MS')
     exec(f'MPigdeep{slid}{PAS}{name} = MPdeep')
     exec(f'MThwdeep{slid}{PAS}{name} = MTdeep')
     exec(f'MSmdeep{slid}{PAS}{name} = MSdeep')
    
    exec(f'T{slid}{PAS}{name} = T')
    exec(f'VPig{slid}{PAS}{name} = VP')
    exec(f'VThw{slid}{PAS}{name} = VT')
    exec(f'VSm{slid}{PAS}{name} = VS')
    exec(f'TV{slid}{PAS}{name} = TV')

    rang = [-5,5]
    d_mdiff = (mdiff - mdiff2)/(x5[0]-num0)*36
    d_mdiff[d_mdiff<rang[0]]=rang[0]
    d_mdiff[d_mdiff>rang[1]]=rang[1]
    print(np.shape(d_mdiff))
    print(np.shape(msk))
    
    d_mdiff[msk != 1] = np.nan


    # Plot the data
    fig,ax = plt.subplots()

    plt.contourf(x_mesh_mid/ 1000, y_mesh_mid / 1000, d_mdiff, np.linspace(rang[0],rang[1],31),cmap='bwr')
    cbar = plt.colorbar()
    plt.xlabel('x (km)')
    plt.ylabel('y (km)')
    cbar.set_label('dz_sdt (m/a)');
    x = x_mesh_mid
    y = y_mesh_mid
    CS = plt.contour(x/1000,y/1000,g0,[0.5],colors='k')
    CS.collections[0].set_label('2013 Grounding line')
    CS = plt.contour(x/1000,y/1000,g2,[0.5],colors='m')
    CS.collections[0].set_label(str(2013+(x5[0])/36) + ' Grounding line')
    ax.axis('equal')
    yr = 2013+(x5[0])/36;
    yr = round(yr*10)/10.
    plt.title(str(yr) + ' Grounding line')

    plt.plot(np.array(xpig)/1000,np.array(ypig)/1000,'-g',linewidth=1)
    plt.plot(np.array(xthw)/1000,np.array(ythw)/1000,'-g',linewidth=1)
    plt.plot(np.array(xsm)/1000,np.array(ysm)/1000,'-g',linewidth=1)
    plt.savefig('/home/dgoldber/www/public_html/amund_couple_figs/thin' + slid + str(PAS) + str(name) + '.png')
    os.system('chmod 755 /home/dgoldber/www/public_html/amund_couple_figs/thin' + slid + str(PAS) + str(name) + '.png')


##########

for gl in ['Pig','Thw','Sm']:
 plt.figure()
 dumpstr=''
 for slid in sliding:
  for PAS in [0,10,20,30]:
   for name in iceParmName:
    if (os.path.exists(f'run_ice_2009_{slid}_{PAS}_{name}_coul')):
     exec(f"plt.plot(T{slid}{PAS}{name},np.array(M{gl}deep{slid}{PAS}{name})/1e12,label='{slid}'+'{PAS}'+'{name}');")
     dumpstr = dumpstr+f'T{slid}{PAS}{name},M{gl}deep{slid}{PAS}{name},'
 plt.grid(which='both')
 plt.title(f'{gl} melt < 500m')
 plt.ylabel('Gt/a')
 plt.legend()

 file=open(f'{gl}_Deepmelt_arrays.pkl','wb')
 dumpstr=dumpstr[:-1]
 exec(f'pkl.dump(({dumpstr}),file)')
 file.close()

 plt.savefig(f'/home/dgoldber/www/public_html/amund_couple_figs/{gl}Deepmelt.png')
 os.system(f'chmod 755 /home/dgoldber/www/public_html/amund_couple_figs/{gl}Deepmelt.png')

###########

for gl in ['Pig','Thw','Sm']:
 plt.figure()
 dumpstr=''
 for slid in sliding:
  for PAS in [0,10,20,30]:
   for name in iceParmName:
    if (os.path.exists(f'run_ice_2009_{slid}_{PAS}_{name}_coul')):
     print(f"plt.plot(T{slid}{PAS}{name},np.array(M{gl}{slid}{PAS}{name})/1e12,label='{slid}'+'{PAS}'+'{name}');")
     exec(f"plt.plot(T{slid}{PAS}{name},np.array(M{gl}{slid}{PAS}{name})/1e12,label='{slid}'+'{PAS}'+'{name}');")
     dumpstr = dumpstr+f'T{slid}{PAS}{name},M{gl}{slid}{PAS}{name},'
 plt.grid(which='both')
 plt.title(f'{gl} melt')
 plt.ylabel('Gt/a')
 plt.legend()

 file=open(f'{gl}_Deepmelt_arrays.pkl','wb')
 dumpstr=dumpstr[:-1]
 exec(f'pkl.dump(({dumpstr}),file)')
 file.close()

 plt.savefig(f'/home/dgoldber/www/public_html/amund_couple_figs/{gl}melt.png')
 os.system(f'chmod 755 /home/dgoldber/www/public_html/amund_couple_figs/{gl}melt.png')

#########

for gl in ['Pig','Thw','Sm']:
 plt.figure()
 dumpstr=''
 for slid in sliding:
  for PAS in [0,10,20,30]:
   for name in iceParmName:
    if (os.path.exists(f'run_ice_2009_{slid}_{PAS}_{name}_coul')):
     exec(f"plt.plot(TV{slid}{PAS}{name}[1:],np.diff(np.array(V{gl}{slid}{PAS}{name}))*12,label='{slid}'+'{PAS}'+'{name}');")
 plt.grid(which='both')
 plt.title(f'{gl} dVAFdt')
 plt.ylabel('Gt/a')
 plt.legend()


 plt.savefig(f'/home/dgoldber/www/public_html/amund_couple_figs/{gl}dVafdt.png')
 os.system(f'chmod 755 /home/dgoldber/www/public_html/amund_couple_figs/{gl}dVafdt.png')

#########

for gl in ['Pig','Thw','Sm']:
 plt.figure()
 dumpstr=''
 for slid in sliding:
  for PAS in [0,10,20,30]:
   for name in iceParmName:
    if (os.path.exists(f'run_ice_2009_{slid}_{PAS}_{name}_coul')):
     #print(f"plt.plot(TV{slid}{PAS}{name}[:],(np.array(V{gl}{slid}{PAS}{name})-V{gl}{slid}{PAS}{name}[0])/400,label='{slid}'+'{PAS}'+'{name}');")
     exec(f"plt.plot(TV{slid}{PAS}{name}[:],(np.array(V{gl}{slid}{PAS}{name})-V{gl}{slid}{PAS}{name}[0])/400,label='{slid}'+'{PAS}'+'{name}');")
     dumpstr = dumpstr+f'TV{slid}{PAS}{name},V{gl}{slid}{PAS}{name},'
 plt.grid(which='both')
 plt.title(f'{gl} VAF')
 plt.ylabel('Gt/a')
 plt.legend()

 file=open(f'{gl}_VAF.pkl','wb')
 dumpstr=dumpstr[:-1]
 exec(f'pkl.dump(({dumpstr}),file)')
 file.close()

 plt.savefig(f'/home/dgoldber/www/public_html/amund_couple_figs/{gl}_VAF.png')
 os.system(f'chmod 755 /home/dgoldber/www/public_html/amund_couple_figs/{gl}_VAF.png')

#########
