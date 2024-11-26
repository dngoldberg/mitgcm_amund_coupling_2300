import numpy as np
from scipy.io import savemat, loadmat
from netCDF4 import Dataset
from IPython import embed
import sys
import time
from scipy.interpolate import RegularGridInterpolator

def readOceDataRegional(PAS, input_dir, prefix, ystart, yend, sub_dir, calc_bounds, calc_init, coords_file, owner):

    PASint = int(PAS)
    ystart = int(ystart)
    yend = int(yend)
    if ( calc_bounds == 'True' ):
        calc_bounds = True
    else:
        calc_bounds = False
    if ( calc_init == 'True' ):
        calc_init = True
    else:
        calc_init = False
    data = loadmat(f'{sub_dir}/zmod.mat',variable_names=['zMod'])
    zMod = data['zMod']


    start_year = int(ystart);
    end_year = int(yend);
    if (PASint>0):
     start_month = (start_year - 1955) * 12
     end_month = (end_year - 1955) * 12 + 11
    else:
     start_month = (start_year - 1991) * 12
     end_month = (end_year - 1991 + 1) * 12 + 1
   
    if PASint > 100:
        nc = f'/exports/geos.ed.ac.uk/iceocean/dgoldber/pahol_output/PAS_{PAS}/run/'
        lens = False
    elif PASint > 0 and PASint <= 100:
        ncbase = f'{input_dir}/{prefix}_{PAS}_'
        lens = True
    else:
        nc = '/home/dgoldber/network_links/ice_data/toshi_data/ASE/'
        lens = False

    if PASint > 100:
        lon = np.double(Dataset(nc + 'stateUvel.nc')['LONGITUDE'][:])
        lat = np.double(Dataset(nc + 'stateUvel.nc')['LATITUDE'][:])
        z = np.double(Dataset(nc + 'stateUvel.nc')['DEPTH'][:])
    else:
        lon = np.double(Dataset(ncbase + '2015.nc')['XG'][:])
        lat = np.double(Dataset(ncbase + '2015.nc')['YG'][:])
        z = np.double(Dataset(ncbase + '2015.nc')['Z'][:])

    # Load other variables like x_mesh_mid, y_mesh_mid, etc.

    data = loadmat(f'{sub_dir}/{coords_file}') 
#               variable_names=['x_mesh_mid', 'y_mesh_mid', 'x_mesh_oce_mid', 'y_mesh_oce_mid', 
#                              'internal_grid_x', 'internal_grid_y', 'diffx', 'diffy'])

#    x_mesh_mid = data['x_mesh_mid'].astype('float64')
#    y_mesh_mid = data['y_mesh_mid'].astype('float64')
    x_mesh_oce_mid = data['x_mesh_oce_mid'].astype('float64')
    y_mesh_oce_mid = data['y_mesh_oce_mid'].astype('float64')
#    internal_grid_x = data['internal_grid_x'].astype('float64')
#    internal_grid_y = data['internal_grid_y'].astype('float64')
#    diffx = data['diffx'].astype('float64')
#    diffy = data['diffy'].astype('float64')
    
    x_mesh_mid, y_mesh_mid = np.meshgrid(x_mesh_oce_mid, y_mesh_oce_mid)

    delx = x_mesh_mid[0, 1] - x_mesh_mid[0, 0]
    dely = y_mesh_mid[1, 0] - y_mesh_mid[0, 0]
    x_mesh_mid = x_mesh_mid[0,:]
    y_mesh_mid = y_mesh_mid[:,0]


    x_botbdry = x_mesh_mid
    y_botbdry = np.ones(len(x_mesh_mid)) * y_mesh_mid[0]
    x_lbdry = np.ones(len(y_mesh_mid)) * x_mesh_mid[0]
    y_lbdry = y_mesh_mid
    x_topbdry = x_mesh_mid
    y_topbdry = np.ones(len(x_mesh_mid)) * y_mesh_mid[-1]

    start_year = ystart
    end_year = yend
    

    if calc_bounds:
 
        lon_ps, lat_ps = np.meshgrid(lon,lat)

        Tleftbdry = np.zeros((end_month-start_month+1,len(x_lbdry),len(z)));
        Sleftbdry = np.zeros((end_month-start_month+1,len(x_lbdry),len(z)));
        Uleftbdry = np.zeros((end_month-start_month+1,len(x_lbdry),len(z)));
        Vleftbdry = np.zeros((end_month-start_month+1,len(x_lbdry),len(z)));
        phileft,lamleft = polarstereo_inv(x_lbdry,y_lbdry)
        lamleft += 360
        ntot = end_month-start_month+1

        for n in range(ntot):

          print ('month ' + str(n));

          if (lens):
            year = start_year + int(np.floor(n/12));
            mon = np.mod(n,12);
            nc = ncbase + str(year) + '.nc'
            with Dataset(nc, "r") as dataset:                
                Ttemp = dataset.variables["THETA"][mon, :, :, :]
                Stemp = dataset.variables["SALT"][mon, :, :, :]
                Utemp = dataset.variables["UVEL"][mon, :, :, :]
                Vtemp = dataset.variables["VVEL"][mon, :, :, :]
          else:
            nc_theta = f"{nc}stateTheta.nc"
            nc_salt = f"{nc}stateSalt.nc"
            nc_uvel = f"{nc}stateUvel.nc"
            nc_vvel = f"{nc}stateVvel.nc"
            with Dataset(nc_theta, "r") as dataset_theta, Dataset(nc_salt, "r") as dataset_salt:
                    Ttemp = dataset_theta.variables["THETA"][start_month+n, :, :, :]
                    Stemp = dataset_salt.variables["SALT"][start_month+n, :, :, :]
            with Dataset(nc_uvel, "r") as dataset_uvel, Dataset(nc_vvel, "r") as dataset_vvel:
                    Utemp = dataset_uvel.variables["UVEL"][start_month+n, :, :, :]
                    Vtemp = dataset_vvel.variables["VVEL"][start_month+n, :, :, :]

          Ttemp[Stemp==0]=np.nan;
          Vtemp[Stemp==0]=np.nan;
          Utemp[Stemp==0]=np.nan;
          Stemp[Stemp==0]=np.nan;

          for k in range(len(z)):
              
              f_T = RegularGridInterpolator((lat, lon), Ttemp[k, :, :],bounds_error=False)
              f_S = RegularGridInterpolator((lat, lon), Stemp[k, :, :],bounds_error=False)
              
              Tleftbdry[n,:,k] = f_T(np.array([phileft, lamleft]).T)
              Sleftbdry[n,:,k] = f_S(np.array([phileft, lamleft]).T)

              xdummy,ydummy,Ups,Vps = vec_ll2ps(lon_ps,lat_ps,Utemp[k,:,:],Vtemp[k,:,:])

              f_U = RegularGridInterpolator((lat, lon), Ups, bounds_error=False)
              f_V = RegularGridInterpolator((lat, lon), Vps, bounds_error=False)

              Uleftbdry[n,:,k] = f_U(np.array([phileft, lamleft]).T)
              Vleftbdry[n,:,k] = f_V(np.array([phileft, lamleft]).T)

              print ( 'lft ' + str(k) + ' ' + str(n) + ' ' + str(Tleftbdry[n,10,k]) ) 

        ########################################              

        Tbotbdry = np.zeros((end_month-start_month+1,len(x_botbdry),len(z)));
        Sbotbdry = np.zeros((end_month-start_month+1,len(x_botbdry),len(z)));
        Ubotbdry = np.zeros((end_month-start_month+1,len(x_botbdry),len(z)));
        Vbotbdry = np.zeros((end_month-start_month+1,len(x_botbdry),len(z)));
        phibot,lambot = polarstereo_inv(x_botbdry,y_botbdry)
        lambot += 360
        ntot = end_month-start_month+1

        for n in range(ntot):

          print ('month ' + str(n));

          if (lens):
            year = start_year + int(np.floor(n/12));
            mon = np.mod(n,12);
            nc = ncbase + str(year) + '.nc'
            with Dataset(nc, "r") as dataset:
                Ttemp = dataset.variables["THETA"][mon, :, :, :]
                Stemp = dataset.variables["SALT"][mon, :, :, :]
                Utemp = dataset.variables["UVEL"][mon, :, :, :]
                Vtemp = dataset.variables["VVEL"][mon, :, :, :]
          else:
            nc_theta = f"{nc}stateTheta.nc"
            nc_salt = f"{nc}stateSalt.nc"
            nc_uvel = f"{nc}stateUvel.nc"
            nc_vvel = f"{nc}stateVvel.nc"
            with Dataset(nc_theta, "r") as dataset_theta, Dataset(nc_salt, "r") as dataset_salt:
                    Ttemp = dataset_theta.variables["THETA"][start_month+n, :, :, :]
                    Stemp = dataset_salt.variables["SALT"][start_month+n, :, :, :]
            with Dataset(nc_uvel, "r") as dataset_uvel, Dataset(nc_vvel, "r") as dataset_vvel:
                    Utemp = dataset_uvel.variables["UVEL"][start_month+n, :, :, :]
                    Vtemp = dataset_vvel.variables["VVEL"][start_month+n, :, :, :]

          Ttemp[Stemp==0]=np.nan;
          Vtemp[Stemp==0]=np.nan;
          Utemp[Stemp==0]=np.nan;
          Stemp[Stemp==0]=np.nan;

          for k in range(len(z)):

              f_T = RegularGridInterpolator((lat, lon), Ttemp[k, :, :],bounds_error=False)
              f_S = RegularGridInterpolator((lat, lon), Stemp[k, :, :],bounds_error=False)

              Tbotbdry[n,:,k] = f_T(np.array([phibot, lambot]).T)
              Sbotbdry[n,:,k] = f_S(np.array([phibot, lambot]).T)

              tic = time.time()
              xdummy,ydummy,Ups,Vps = vec_ll2ps(lon_ps,lat_ps,Utemp[k,:,:],Vtemp[k,:,:])
              toc = time.time()

              f_U = RegularGridInterpolator((lat, lon), Ups, bounds_error=False)
              f_V = RegularGridInterpolator((lat, lon), Vps, bounds_error=False)

              Ubotbdry[n,:,k] = f_U(np.array([phibot, lambot]).T)
              Vbotbdry[n,:,k] = f_V(np.array([phibot, lambot]).T)

              print ( 'bot ' + str(k) + ' ' + str(n) + ' ' + str(Tbotbdry[n,10,k]) )


        bdry_data = {
             'Tleftbdry': np.transpose(Tleftbdry,(2,1,0)),
             'Sleftbdry': np.transpose(Sleftbdry,(2,1,0)),
             'Uleftbdry': np.transpose(Uleftbdry,(2,1,0)),
             'Vleftbdry': np.transpose(Vleftbdry,(2,1,0)),
             'Tbotbdry': np.transpose(Tbotbdry,(2,1,0)),
             'Sbotbdry': np.transpose(Sbotbdry,(2,1,0)),
             'Ubotbdry': np.transpose(Ubotbdry,(2,1,0)),
             'Vbotbdry': np.transpose(Vbotbdry,(2,1,0)),
             'ystart': ystart,
             'yend': yend,
             'z': z
             }
        savemat(f'{sub_dir}/bdryDatapaholPy{prefix}{PAS}{owner}.mat', bdry_data)








    
    
    # Save data to .mat file
    if calc_init:
    
        if (lens):
            nc = ncbase + str(start_year) + '.nc'
#            with Dataset(nc_theta, "r") as dataset_theta, Dataset(nc_salt, "r") as dataset_salt:
#                    Ttemp = dataset_theta.variables["THETA"][start_month+n, :, :, :]
#                    Stemp = dataset_salt.variables["SALT"][start_month+n, :, :, :]

        xgridp, ygridp, zgridp = np.meshgrid(x_mesh_mid,y_mesh_mid,zMod);   
        phi,lam = polarstereo_inv(xgridp,ygridp)
        lam = lam + 360
        ntot=12
        mSz, nSz, o = np.shape(phi);
        
        Tinit = np.zeros((len(z),nSz,mSz))
        Sinit = np.zeros((len(z),nSz,mSz))
        latInterp, lonInterp = np.meshgrid(lat,lon)
        
        for n in range(ntot):
                    
        
            if lens:
                nc = f"{ncbase}{start_year}.nc"
                with Dataset(nc, "r") as dataset:
                    Ttemp = dataset.variables["THETA"][n, :, :, :]
                    Stemp = dataset.variables["SALT"][n, :, :, :]
            else:
                nc_theta = f"{nc}stateTheta.nc"
                nc_salt = f"{nc}stateSalt.nc"
                with Dataset(nc_theta, "r") as dataset_theta, Dataset(nc_salt, "r") as dataset_salt:
                    Ttemp = dataset_theta.variables["THETA"][start_month+n, :, :, :]
                    Stemp = dataset_salt.variables["SALT"][start_month+n, :, :, :]

            Ttemp[Stemp == 0] = np.nan
            Stemp[Stemp == 0] = np.nan

            for k in range(len(z)):
                f_T = RegularGridInterpolator((lat, lon), Ttemp[k, :, :],bounds_error=False)
                f_S = RegularGridInterpolator((lat, lon), Stemp[k, :, :],bounds_error=False)
                Tinit[k, :, :] = Tinit[k, :, :] + 1 / ntot * f_T(np.array([phi[:, :, 0].ravel(), lam[:, :, 1].ravel()]).T).reshape((nSz,mSz))
                Sinit[k, :, :] = Sinit[k, :, :] + 1 / ntot * f_S(np.array([phi[:, :, 0].ravel(), lam[:, :, 1].ravel()]).T).reshape((nSz,mSz))

                tlayer = Tinit[k, :, :]
                I = ~np.isnan(tlayer)
                print(f"month {n}; level {k} {np.nanmean(tlayer[I])}")
    
        init_data = {
            'Sinit': np.transpose(Sinit,(2,1,0)),
            'Tinit': np.transpose(Tinit,(2,1,0)),
            'z': z
        }
        savemat(f'{sub_dir}/initpaholPy{prefix}{PAS}{owner}.mat', init_data)

        
def polarstereo_inv(x, y, a=6378137.0, e=0.08181919, phi_c=-71, lambda_0=0):
        # Convert to radians
        phi_c = np.deg2rad(phi_c)
        lambda_0 = np.deg2rad(lambda_0)

        # If the standard parallel is in the Southern Hemisphere, switch signs
        if phi_c < 0:
            pm = -1
            phi_c = -phi_c
            lambda_0 = -lambda_0
            x = -x
            y = -y
        else:
            pm = 1

        t_c = np.tan(np.pi / 4 - phi_c / 2) / ((1 - e * np.sin(phi_c)) / (1 + e * np.sin(phi_c)))**(e / 2)
        m_c = np.cos(phi_c) / np.sqrt(1 - e**2 * (np.sin(phi_c))**2)
        rho = np.sqrt(x**2 + y**2)
        t = rho * t_c / (a * m_c)

        chi = np.pi / 2 - 2 * np.arctan(t)
        phi = chi + (e**2 / 2 + 5 * e**4 / 24 + e**6 / 12 + 13 * e**8 / 360) * np.sin(2 * chi) \
            + (7 * e**4 / 48 + 29 * e**6 / 240 + 811 * e**8 / 11520) * np.sin(4 * chi) \
            + (7 * e**6 / 120 + 81 * e**8 / 1120) * np.sin(6 * chi) \
            + (4279 * e**8 / 161280) * np.sin(8 * chi)

        lambda_val = lambda_0 + np.arctan2(x, -y)

        # Correct the signs and phasing
        phi = pm * phi
        lambda_val = pm * lambda_val
        lambda_val = (lambda_val + np.pi) % (2 * np.pi) - np.pi  # Want longitude in the range -pi to pi

        # Convert back to degrees
        phi = np.rad2deg(phi)
        lambda_val = np.rad2deg(lambda_val)

        return phi, lambda_val

def vec_ll2ps(lon, lat, U, V, Sp=-71, Cm=0):
        Ups = np.zeros(U.shape)
        Vps = np.zeros(V.shape)
        delt = 1
        
        def polarstereo_fwd(phi, lon, a=6378137.0, e=0.08181919, phi_c=-71, lambda_0=0):
        # Convert to radians
            phi = np.deg2rad(phi)
            phi_c = np.deg2rad(phi_c)
            lon = np.deg2rad(lon)
            lambda_0 = np.deg2rad(lambda_0)

            # If the standard parallel is in the Southern Hemisphere, switch signs
            if (phi_c < 0):
                pm = -1
                phi = -phi
                phi_c = -phi_c
                lon = -lon
                lambda_0 = -lambda_0
            else:
                pm = 1

            t = np.tan(np.pi / 4 - phi / 2) / ((1 - e * np.sin(phi)) / (1 + e * np.sin(phi)))**(e / 2)
            t_c = np.tan(np.pi / 4 - phi_c / 2) / ((1 - e * np.sin(phi_c)) / (1 + e * np.sin(phi_c)))**(e / 2)
            m_c = np.cos(phi_c) / np.sqrt(1 - e**2 * (np.sin(phi_c))**2)
            rho = a * m_c * t / t_c

            m = np.cos(phi) / np.sqrt(1 - e**2 * (np.sin(phi))**2)
            x = pm * rho * np.sin(lon - lambda_0)
            y = -pm * rho * np.cos(lon - lambda_0)
            k = rho / (a * m)

            return x, y

#        for t in range(U.shape[2]):
        endlon = lon + U / (2 * np.pi * np.cos(np.deg2rad(lat)) * 6378137 / 360) * delt
        endlat = lat + V / (2 * np.pi * 6378137 / 360) * delt

            # convert end points to lat-lon

        tic = time.time()
        x,y = polarstereo_fwd(lat, lon, phi_c=Sp, lambda_0=Cm)
        toc = time.time()
        tic = time.time()
        endx, endy = polarstereo_fwd(endlat, endlon, phi_c=Sp, lambda_0=Cm)
        toc = time.time()

            # calculate lat-lon displacements and convert to metres and then velocities
        Ups = (endx - x) / delt
        Vps = (endy - y) / delt

        return x, y, Ups, Vps
        
if __name__ == "__main__":
    readOceDataRegional(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9],sys.argv[10])
