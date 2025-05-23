function gen_mesh(use_2004,PAS,weert,dig_depth,gamma_depth,gammaT,cd,bc_setup)

% gen_mesh(use_2004,PAS,weert,dig_depth,gamma_depth,gammaT,cd,bc_setup)
% gen_mesh(false,80,false,80,200,0.00014,.006,false)

ice_setup = true;
oce_setup = true;





global sub_dir

dx0 = 1250;
dy0 = 1250;

if (use_2004)
 sub_dir = 'start_2004_input'
else
 sub_dir = 'start_2009_input'
end

if(length(dir(sub_dir))==0);
     eval(['!mkdir ' sub_dir]);

     
end
     if(length(dir([sub_dir '/oce']))==0);
         eval(['!mkdir ' sub_dir '/oce']);
     end
     if(length(dir([sub_dir '/ice']))==0);
         eval(['!mkdir ' sub_dir '/ice']);
     end

addpath('/home/dgoldber/network_links/datastore/ice_data/ThwaitesDataProphet/CODE');
addpath('/home/dgoldber/network_links/datastore/ice_data/bamber/');

unif_bounds = [-1750e3 -1100e3 -750e3 50e3];

Xuniform = unif_bounds(1):1e3:unif_bounds(2);
Yuniform = unif_bounds(3):1e3:unif_bounds(4);

[Xuniform Yuniform] = meshgrid(Xuniform,Yuniform);


npx_oce = 16;
npy_oce = 16;

npx_ice = 8;
npy_ice = 16;

oce_extent_y_init = [-7.12e5 -2e5];
ny_init = diff(oce_extent_y_init) / npy_oce / dy0;
nyoce = ceil(ny_init) * npy_oce;

oce_extent_x_init = [-1.7e6 -1.4e6];
nx_init = diff(oce_extent_x_init) / npx_oce / dx0;
nxoce = ceil(nx_init) * npx_oce;

x_mesh_oce = oce_extent_x_init(1) + (0:dx0:((nxoce)*dx0));
y_mesh_oce = oce_extent_y_init(1) + (0:dy0:((nyoce)*dy0));

y_mesh_oce_mid = .5 * (y_mesh_oce(1:end-1)+y_mesh_oce(2:end));

res_slope = 7.5;

dy_var = dy0 + res_slope * (1:400);
dx_var = dx0 + res_slope * (1:400);
y_var = cumsum(dy_var);
x_var = cumsum(dx_var);

x_mesh = [(x_mesh_oce(1)-fliplr(x_var)) x_mesh_oce (x_mesh_oce(end)+x_var)];
y_mesh = [(y_mesh_oce(1)-fliplr(y_var)) y_mesh_oce (y_mesh_oce(end)+y_var)];

k = max(find(y_mesh<unif_bounds(3)));
k2 = min(find(y_mesh>unif_bounds(4)));
k2 = k + npy_ice * ceil((k2-k)/npy_ice);
y_mesh = y_mesh(k:k2);

l = max(find(x_mesh<unif_bounds(1)));
l2 = min(find(x_mesh>unif_bounds(2)));
l2 = l + npx_ice * ceil((l2-l)/npx_ice);
x_mesh = x_mesh(l:l2);

% now do new ocean grid

oce_extent_y_init = [-7.12e5 -1.35e5];
oce_extent_x_init = [-1.7e6 -1.3e6];

% revert to old
%oce_extent_y_init = [-7.12e5 0e5];
%oce_extent_x_init = [-1.7e6 -1.2e6];

npx_oce = 16;
npy_oce = 24;

k2_oce = min(find(y_mesh>oce_extent_y_init(2)))
l2_oce = min(find(x_mesh>oce_extent_x_init(2)));

ind2x = find(x_mesh==x_mesh_oce(end));
ind2y = find(y_mesh==y_mesh_oce(end))
y_mesh(k2_oce)


x_mesh_oce = [x_mesh_oce x_mesh((ind2x+1):l2_oce)];
y_mesh_oce = [y_mesh_oce y_mesh((ind2y+1):k2_oce)];

y_mesh_oce(end)

nxoce = length(x_mesh_oce)-1;
nyoce = length(y_mesh_oce)-1

nxoce = floor(nxoce/npx_oce)*npx_oce;
nyoce = floor(nyoce/npy_oce)*npy_oce

x_mesh_oce = x_mesh_oce(1:(nxoce+1));
y_mesh_oce = y_mesh_oce(1:(nyoce+1));


x_mesh_mid = .5 * (x_mesh(1:end-1)+x_mesh(2:end));
y_mesh_mid = .5 * (y_mesh(1:end-1)+y_mesh(2:end));
x_mesh_oce_mid = .5 * (x_mesh_oce(1:end-1)+x_mesh_oce(2:end));
y_mesh_oce_mid = .5 * (y_mesh_oce(1:end-1)+y_mesh_oce(2:end));



nxice = length(x_mesh_mid);
nyice = length(y_mesh_mid);

diffx = diff(x_mesh);
diffx_oce = diff(x_mesh_oce);
diffy = diff(y_mesh);
diffy_oce = diff(y_mesh_oce);

write_coupled_input('delX.bin',diffx,'oce');
write_coupled_input('delY.bin',diffy,'oce');
write_coupled_input('delX_oce.bin',diffx_oce,'oce');
write_coupled_input('delY_oce.bin',diffy_oce,'oce');


ind1x = find(x_mesh_mid==x_mesh_oce_mid(1));
ind2x = find(x_mesh_mid==x_mesh_oce_mid(end));

ind1y = find(y_mesh_mid==y_mesh_oce_mid(1));
ind2y = find(y_mesh_mid==y_mesh_oce_mid(end));




nxoce = length(x_mesh_oce_mid);
nyoce = length(y_mesh_oce_mid);

nx_grid = ceil(nxoce/npx_oce)
ny_grid = ceil(nyoce/npy_oce)

gx = npx_oce*nx_grid-nxoce;
gy = npy_oce*ny_grid-nyoce;

length(x_mesh_mid)/npx_ice
length(y_mesh_mid)/npy_ice

[xoce yoce] = meshgrid(x_mesh_oce_mid,y_mesh_oce_mid);

internal_grid_x = ind1x:ind2x;
internal_grid_y = ind1y:ind2y;

zMod = -12.5:-25:-2425;
length(zMod)

save([sub_dir '/meshcoords.mat'],'x_mesh_mid', 'y_mesh_mid', 'x_mesh_oce_mid', 'y_mesh_oce_mid', ...
    'internal_grid_x','internal_grid_y','diffx','diffy','zMod');


write_coupled_input('Ix_overlap.bin',internal_grid_x,'oce');
write_coupled_input('Iy_overlap.bin',internal_grid_y,'oce');


calc_init = true;
if (length(dir([sub_dir '/initpaholPy' num2str(PAS) '.mat']))==1);
    calc_init = false;
end
calc_bounds = true;
if (length(dir([sub_dir '/bdryDatapaholPy' num2str(PAS) '.mat']))==1);
    calc_bounds = false;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




if(use_2004)
    ystart = 1999;
else
    ystart = 2006;
end

if (PAS>100);
    yend = 2006;
elseif (PAS>10);
    ystart = 2010;
    yend = 2100;
elseif (PAS>1);
    ystart = 0;
    yend = 0;
else
    yend = 2019;
end

if (ice_setup)
setup_experiment(nxice,nyice,0,0,use_2004,PAS,weert,dig_depth,gamma_depth,gammaT,cd);
end

%readOceDataRegional (PAS, zMod, ystart, yend, calc_bounds, calc_init);

% pas=20: paris2 mean
% pas=30: rcp85 mean
% pas=300: ssp585_2300 corrected



if (PAS==20);
     exptstr = 'Paris2Mean';
elseif (PAS==30);
     exptstr = 'rcp85Mean';
elseif (PAS==10)
     exptstr = 'BaseClimMean'
elseif (PAS==80)
     exptstr = 'CESMRCP852300'
else 
     error('bad PAS code')
end

eval(['!ln -s /home/dgoldber/network_links/geosIceOcean/dgoldber/pahol_output/naughten_jas/bdryDatapahol' 'rcp85Mean' '.mat ' sub_dir '/bdryDatapahol' num2str(30) '.mat']);
eval(['!ln -s /exports/geos.ed.ac.uk/iceocean/dgoldber/MITgcm_forinput/amund_couple/mitgcm_amund_coupling_2300/prepare_2300_data/naughtenClim.mat ' sub_dir '/naughtenClim.mat']);
eval(['!ln -s /exports/geos.ed.ac.uk/iceocean/dgoldber/MITgcm_forinput/amund_couple/mitgcm_amund_coupling_2300/prepare_2300_data/slices_anom_forcing.mat ' sub_dir '/slices_anom_forcing.mat']);


if (bc_setup)
rdmds_init(PAS,nxoce,nyoce,zMod,gx,gy,2,ystart,yend);
end
if (oce_setup)
gendata(PAS,nxoce,nyoce,zMod,gx,gy,2,1030,ystart,dig_depth,gamma_depth,gammaT,cd);
end


 
return
