#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 15:47:49 2020

@author: danielgoldberg
"""

import numpy as np
import sys
import os
sys.path.append(os.environ['ROOTDIRWORK'] + '/utils/python/MITgcmutils/MITgcmutils/');
from mds import rdmds
import time
from IPython import embed

print('argument 1: ' + sys.argv[1])
print('argument 2: ' + sys.argv[2])
run_ice='run_ice_' + sys.argv[3] + '_' + sys.argv[4] + '_' + sys.argv[5] + '_' + sys.argv[6] + '_' + sys.argv[7]
run_oce='run_oce_' + sys.argv[3] + '_' + sys.argv[4] + '_' + sys.argv[5] + '_' + sys.argv[6] + '_' + sys.argv[7]

Ix = np.fromfile('../' + run_oce + '/Ix_overlap.bin').byteswap().astype(int)
Iy = np.fromfile('../' + run_oce + '/Iy_overlap.bin').byteswap().astype(int)


q = str(sys.argv[1])
niter_ice = int(q)
numiter_ice = ('0' * (10-len(q))) + q
data_ice,x,m = rdmds('../' + run_ice + '/pickup_streamice.ckptA',returnmeta=True)
print('thick transfer ice input iter:' + numiter_ice)
print('thick transfer ice pickup step:' + str(m['timestepnumber'][0]))

if (niter_ice != m['timestepnumber'][0]):
    print('ERROR -- time step of pickup does not match')
    exit()
nzice = len(rdmds('RC'))


numiter_oce = str(sys.argv[2])
niter_oce = int(sys.argv[2])
data_oce,x,m = rdmds('../' + run_oce + '/pickup_streamice.ckptA',returnmeta=True)

print('thick transfer oce input iter:' + numiter_oce)
print('thick transfer oce pickup step:' + str(m['timestepnumber'][0]))

if (niter_oce != m['timestepnumber'][0]):
    print('ERROR -- time step of pickup does not match')
    exit()
nzoce = len(rdmds('../' + run_oce + '/RC'))

thick = data_oce[nzoce+4,:len(Iy),:len(Ix)]
mask  = data_oce[nzoce+1,:len(Iy),:len(Ix)]

thick_ice_overlap = data_ice[nzice+4,(Iy[0]-1):Iy[-1],(Ix[0]-1):Ix[-1]]
thick_ice_overlap[mask==1] = thick[mask==1]

data_ice[nzice+4,(Iy[0]-1):Iy[-1],(Ix[0]-1):Ix[-1]] = thick_ice_overlap;
data_ice.byteswap().tofile('../' + run_ice + '/pickup_streamice.ckptA.data');
thick_ice = data_ice[nzice+4,:,:]
thick_ice.byteswap().tofile('../transfer_files_' + sys.argv[3] + '_' + sys.argv[4] + '_' + sys.argv[5] + '_' + sys.argv[6] + '_' + sys.argv[7] + '/thick_ice_mod' + numiter_ice + '.bin')

sys.exit()
