from mds import rdmds 
from scipy.io import loadmat
import numpy as np
import sys
import glob
from IPython import embed

YR = sys.argv[1]
PAS = sys.argv[2]
invName = sys.argv[3]
invIter = int(sys.argv[4])
sliding_law=(sys.argv[5])
tc=(sys.argv[6])

if (tc=='TC'):  


 direc='../../../mitgcm_amund_coupling_2300/ice_init_tc/run_ad_' + YR + '_' + PAS + '_' + invName + '_' + sliding_law
 print(direc)

 bglen_control=rdmds(direc + '/gradcontrol/xx_bglen',invIter)
 bglen0 = rdmds(direc + '/runoptiter000/B_glen_sqrt')
 beta_control =rdmds(direc + '/gradcontrol/xx_beta',invIter)
 beta0 = rdmds(direc + '/runoptiter000/C_basal_fric')
 
 bglen = bglen_control + bglen0
 beta = beta_control + beta0

 beta.byteswap().tofile('BetaTC' + sliding_law + PAS +'.bin')
 bglen.byteswap().tofile('BglenTC' + sliding_law + PAS + '.bin')

else:

 direc='../../../mitgcm_amund_coupling_2300/ice_init/run_ad_' + YR + '_' + sliding_law
 print(direc)

 bglen_control=rdmds(direc + '/runoptiter' + str(invIter) + '/xx_bglen',invIter)
 beta_control=rdmds(direc + '/runoptiter' + str(invIter) + '/xx_beta',invIter)
 bglen0 = rdmds(direc + '/runoptiter000/B_glen_sqrt')
 beta0 = rdmds(direc + '/runoptiter000/C_basal_fric')

 bglen = bglen_control + bglen0
 beta = beta_control + beta0

 beta.byteswap().tofile('BetaSnap' + sliding_law + '.bin')
 bglen.byteswap().tofile('BglenSnap' + sliding_law + '.bin')


