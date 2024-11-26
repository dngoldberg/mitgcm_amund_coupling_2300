#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 17:08:34 2020

@author: danielgoldberg
"""

import numpy as np
import sys
import os
import shutil
sys.path.append(os.environ['ROOTDIRWORK'] + '/utils/python/MITgcmutils/MITgcmutils/');
from mds import rdmds

run_ice='run_ice_' + sys.argv[2] + '_' + sys.argv[3] + '_' + sys.argv[4] + '_' + sys.argv[5] + '_' + sys.argv[6]
run_oce='run_oce_' + sys.argv[2] + '_' + sys.argv[3] + '_' + sys.argv[4] + '_' + sys.argv[5] + '_' + sys.argv[6]

Ix = np.fromfile('../' + run_oce + '/Ix_overlap.bin').byteswap().astype(int);
Iy = np.fromfile('../' + run_oce + '/Iy_overlap.bin').byteswap().astype(int);

nzice = len(rdmds('../' + run_ice + '/RC'))

numiter = str(sys.argv[1])
data_ice,x,m = rdmds('../' + run_ice + '/pickup_streamice.ckptA',returnmeta=True)
print('vel transfer ice input iter:' + numiter)
print('vel transfer pickup iter:' + str(m['timestepnumber'][0]))

if (int(sys.argv[1]) != m['timestepnumber'][0]):
 print('ERROR -- time step of pickup does not match')
 exit()

if (np.int32(numiter)>1):
 data_ice_prev,x,m2 = rdmds('../' + run_ice + '/pickup_streamice.ckptAlast',returnmeta=True)
 if ((int(sys.argv[1])-1) != m2['timestepnumber'][0]):
  print('ERROR -- time step of pickup does not match')
  exit()
 thick_ice_prev = data_ice_prev[nzice+4,:,:]
else:
 thick_ice_prev = np.fromfile('../' + run_ice + '/BedMachineThickMod.bin').byteswap().reshape(m['dimlist'][1],m['dimlist'][0]); 

shutil.copyfile('../' + run_ice + '/pickup_streamice.ckptA.meta','../' + run_ice + '/pickup_streamice.ckptAlast.meta')
shutil.copyfile('../' + run_ice + '/pickup_streamice.ckptA.data','../' + run_ice + '/pickup_streamice.ckptAlast.data')

u_ice = np.zeros((len(Iy),len(Ix)))
v_ice = np.zeros((len(Iy),len(Ix)))
h_ice = np.zeros((len(Iy),len(Ix)))

u_ice[0:len(Iy),0:len(Ix)] = data_ice[nzice+2,(Iy[0]-1):Iy[-1],(Ix[0]-1):Ix[-1]]
v_ice[0:len(Iy),0:len(Ix)] = data_ice[nzice+3,(Iy[0]-1):Iy[-1],(Ix[0]-1):Ix[-1]]
h_ice[0:len(Iy),0:len(Ix)] = thick_ice_prev[(Iy[0]-1):Iy[-1],(Ix[0]-1):Ix[-1]]

u_ice.byteswap().tofile('../transfer_files_' + sys.argv[2] + '_' + sys.argv[3] + '_' + sys.argv[4] + '_' + sys.argv[5] + '_' + sys.argv[6] + '/Uice' + numiter + '.bin')
v_ice.byteswap().tofile('../transfer_files_' + sys.argv[2] + '_' + sys.argv[3] + '_' + sys.argv[4] + '_' + sys.argv[5] + '_' + sys.argv[6] + '/Vice' + numiter + '.bin')

u_ice.byteswap().tofile('../' + run_oce + '/Uice_oce.bin')
v_ice.byteswap().tofile('../' + run_oce + '/Vice_oce.bin')

uthickbdry = np.zeros((len(Iy),len(Ix)))
vthickbdry = np.zeros((len(Iy),len(Ix)))

ufacemask = np.fromfile('../' + run_oce + '/ufacemask_oce.bin').byteswap().reshape(len(Iy),len(Ix));
vfacemask = np.fromfile('../' + run_oce + '/vfacemask_oce.bin').byteswap().reshape(len(Iy),len(Ix));
uthickbdry[ufacemask==3] = h_ice[ufacemask==3]
vthickbdry[vfacemask==3] = h_ice[vfacemask==3]
uthickbdry.byteswap().tofile('HBCx_oce.bin')
vthickbdry.byteswap().tofile('HBCy_oce.bin')

exit()
