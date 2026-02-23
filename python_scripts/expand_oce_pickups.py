import numpy as np
import sys
from scipy.io import loadmat
from mds import rdmds, parsemeta
import os
from shutil import copyfile
from IPython import embed
import fileinput
import configparser



#exptfolder = sys.argv[1]
#oce_pickup = int(sys.argv[2])
#old_mesh = sys.argv[3]
#new_mesh = sys.argv[4]
#input_folder = sys.argv[5]
#newpickup_location = sys.argv[6]

def main(cfg_path):

    cp = configparser.ConfigParser()
    cp.read(cfg_path)

    cfg = cp['general']

    exptfolder = cfg.get('experiment_folder')
    oce_pickup = cfg.getint('ocean_pickup_step')
    old_mesh = cfg.get('old_mesh_matfile')
    new_mesh = cfg.get('new_mesh_matfile')
    input_folder = cfg.get('input_folder')
    newpickup_location = cfg.get('pickup_save_location')
    python_scripts_folder = cfg.get('python_scripts_folder')

    sys.path.append('../python_scripts')
    from average_timesteps import writemeta


    chkpt_count = 0

    os.makedirs(newpickup_location, exist_ok=True)
    copyfile(exptfolder + '/pickups/pickup.' + str(oce_pickup).zfill(10) + '.meta',newpickup_location + '/pickup.ckptA.meta')
    copyfile(exptfolder + '/pickups/pickup_streamice.' + str(oce_pickup).zfill(10) + '.meta',newpickup_location + '/pickup_streamice.ckptA.meta')
    copyfile(exptfolder + '/pickups/pickup_shelfice.' + str(oce_pickup).zfill(10) + '.meta',newpickup_location + '/pickup_shelfice.ckptA.meta')

    meta_main = parsemeta(exptfolder + '/pickups/pickup.' + str(oce_pickup).zfill(10) + '.meta')
    meta_shelfice = parsemeta(exptfolder + '/pickups/pickup_shelfice.' + str(oce_pickup).zfill(10) + '.meta')
    meta_streamice = parsemeta(exptfolder + '/pickups/pickup_streamice.' + str(oce_pickup).zfill(10) + '.meta')



    meshmatold = loadmat(old_mesh)
    meshmatnew = loadmat(new_mesh)
    xmeshold = meshmatold['x_mesh_oce_mid'][0]
    xmeshnew = meshmatnew['x_mesh_oce_mid'][0]
    nxnew = len(xmeshnew)
    internal_gridx_old = meshmatold['internal_grid_x'][0]-1
    internal_gridx_new = meshmatnew['internal_grid_x'][0]-1
    internal_gridy = meshmatnew['internal_grid_y'][0]-1
    zMod = meshmatnew['zMod'][0]

    dims = meta_main['dimList'].copy()
    dims[0] = nxnew
    dims[2] = nxnew

    meta_main['dimList'] = dims
    meta_shelfice['dimList'] = dims
    meta_streamice['dimList'] = dims

    ice_exptfolder = exptfolder
    ice_exptfolder = ice_exptfolder.replace('run_oce','run_ice')

    pickupfile = exptfolder + '/pickups/pickup'
    shelficepickupfile = exptfolder + '/pickups/pickup_shelfice'
    streamicepickupfile = exptfolder + '/pickups/pickup_streamice'
    large_streamicepickupfile = ice_exptfolder + '/pickups/pickup_streamice'

    qbase,x,mbase = rdmds(pickupfile,int(oce_pickup),returnmeta=True)

    nxold = np.size(qbase,2)
    ny = np.size(qbase,1)
    nstack = np.size(qbase,0)
    nz = int((nstack-3)/6)
    initTheta = np.fromfile(input_folder + '/theta.init.2010.41').byteswap().reshape(nz,ny,nxnew)
    initSalt = np.fromfile(input_folder + '/salt.init.2010.41').byteswap().reshape(nz,ny,nxnew)
    bathy = np.fromfile(input_folder + '/bathy_mod_dig_80.bin').byteswap().reshape(ny,nxnew)
    etainit = np.fromfile(input_folder + '/etainit.41.bin').byteswap().reshape(ny,nxnew)
    topoinit = np.fromfile(input_folder + '/shelftopo.41.bin').byteswap().reshape(ny,nxnew)


    fldList = mbase['fldlist']
    DataNew = np.zeros((np.size(qbase,0),np.size(qbase,1), nxnew))

    for i in range(len(fldList)):
        fldName = fldList[i]
        if fldName in ['Uvel', 'Vvel', 'GuNm1', 'GvNm1','Theta', 'Salt']:
            fldOld = qbase[i*nz:(i+1)*nz,:,:]
            fldNew = np.zeros((nz,ny,nxnew))
            fldNew[:,:,:nxold] = fldOld

            if fldName in ['Theta','Salt']:

                K = (-1-np.floor(bathy/25.)).astype('int')
                i_idx, j_idx = np.indices(K.shape)
                
                if fldName == 'Theta':

                    initField = initTheta.copy()

                if fldName == 'Salt':

                    initField = initSalt.copy()

                fldNew[K[:,(nxold-1):],i_idx[:,(nxold-1):],j_idx[:,(nxold-1):]] = \
                       initField[K[:,(nxold-1):],i_idx[:,(nxold-1):],j_idx[:,(nxold-1):]]
    #            embed()

    #            fldNew[:,:,(nxold-1):] = initField[:,:,(nxold-1):]

            DataNew[i*nz:(i+1)*nz,:,:] = fldNew

        if fldName in ['EtaN', 'dEtaHdt', 'EtaH']:

            i_e = i-9
            fldOld = qbase[i_e,:,:]
            fldNew = np.zeros((ny,nxnew))
            fldNew[:,:nxold] = fldOld[:,:nxold]
            
            if fldName in ['EtaN','EtaH']:

                 fldNew[:,nxold:] = etainit[:,nxold:]
            
            DataNew[i_e,:,:] = fldNew

    DataNew.byteswap().tofile(newpickup_location + '/pickup.ckptA.data')

    ###

    qSI,x,mSI = rdmds(streamicepickupfile,oce_pickup,returnmeta=True)
    qLSI,x,mLSI = rdmds(large_streamicepickupfile,int((oce_pickup-1555200)/8640),returnmeta=True)
    nzLSI = np.size(qLSI,0)-9
    icethickLSI = qLSI[nzLSI+4,:,:]

    fldList = mSI['fldlist']
    DataNewSI = np.zeros((np.size(qSI,0),np.size(qSI,1), nxnew))


    for i in range(len(fldList)):
        fldName = fldList[i]
        if (i >= 0) & (i < 5):

            fldNew = np.zeros((ny,nxnew))
            fldLSI = qLSI[nzLSI+i,:,:]
            fldNew[:,:] = fldLSI[internal_gridy[:,None],internal_gridx_new]
            fldNew[:,:(nxold-1)] = qSI[i+nz,:,:(nxold-1)]
            DataNewSI[nz+i,:,:] = fldNew

    DataNewSI.byteswap().tofile(newpickup_location + '/pickup_streamice.ckptA.data')

    ###

    qShI,x,mShI = rdmds(shelficepickupfile,oce_pickup,returnmeta=True)
    DataNewShI = np.zeros((np.size(qShI,0),np.size(qShI,1), nxnew))

    fldNew = np.zeros((ny,nxnew))
    fldNew[:,(nxold-1):] = 917 * icethickLSI[internal_gridy[:,None],internal_gridx_new[(nxold-1):]]
    fldNew[:,:nxold-1] = qShI[0,:,:-1]
    DataNewShI[0,:,:] = fldNew

    fldNew = np.zeros((ny,nxnew))
    fldNew[:,:] = topoinit.copy()
    fldNew[:,:(nxold-1)] = qShI[1,:,:(nxold-1)]
    DataNewShI[1,:,:] = fldNew

    DataNewShI.byteswap().tofile(newpickup_location + '/pickup_shelfice.ckptA.data')

    ###

    for nm in ['pickup_shelfice.ckptA.meta', 'pickup.ckptA.meta', 'pickup_streamice.ckptA.meta']:

     with open(newpickup_location+'/'+nm,'r') as file:
      filedata = file.read()
     
     filedata = filedata.replace(str(nxold),str(nxnew))
     
     with open(newpickup_location+'/'+nm, 'w') as file:
      file.write(filedata)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        cfg_file = 'params.ini'
    else:
        cfg_file = sys.argv[1]
    main(cfg_file)



#writemeta(meta_main, newpickup_location + '/pickup.ckptA.meta')
#writemeta(meta_shelfice, newpickup_location + '/pickup_shelfice.ckptA.meta')
#writemeta(meta_streamice, newpickup_location + '/pickup_streamice.ckptA.meta')

#  ['SHI_mass', 'R_Shelfi']








