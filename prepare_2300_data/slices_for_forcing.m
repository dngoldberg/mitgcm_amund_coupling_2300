load /home/dgoldber/network_links/geosIceOcean/dgoldber/MITgcm_forinput/amund_couple/mitgcm_amund_coupling_2300/input/start_2009_input/meshcoords.mat
addpath('/home/dgoldber/network_links/ice_data/ThwaitesDataProphet/CODE/');

salt = ncread('/home/dgoldber/scratch/globus/climatology_from_obs_1995-2017/obs_salinity_1995-2017_8km_x_60m_clipped.nc','salinity');
temp = ncread('/home/dgoldber/scratch/globus/climatology_from_obs_1995-2017/obs_temperature_1995-2017_8km_x_60m_clipped.nc','temperature');
X = ncread('/home/dgoldber/scratch/globus/ccsm4_RCP85/1995-2300/CCSM4_RCP85_salinity_8km_x_60m_clipped_1995.nc','x');
Y = ncread('/home/dgoldber/scratch/globus/ccsm4_RCP85/1995-2300/CCSM4_RCP85_salinity_8km_x_60m_clipped_1995.nc','y');
Z = ncread('/home/dgoldber/scratch/globus/ccsm4_RCP85/1995-2300/CCSM4_RCP85_salinity_8km_x_60m_clipped_1995.nc','z');

x_mesh_oce_mid=x_mesh_oce_mid(1:240);
y_mesh_oce_mid=y_mesh_oce_mid(1:416);
[x y]=meshgrid(x_mesh_oce_mid,y_mesh_oce_mid);
[X Y] = meshgrid(X,Y);

bed = interpBedmachineAntarctica(X,Y,'bed');
thick = interpBedmachineAntarctica(X,Y,'thickness');
surf = interpBedmachineAntarctica(X,Y,'surface');
base = surf-thick;
zmask = interpBedmachineAntarctica(X,Y,'mask');


bed_naught = interpBedmachineAntarctica(x,y,'bed');

Zmat = zeros(1,1,length(Z));
Zmat(1,1,:)=Z;
Zmat = repmat(Zmat,[size(X) 1]);
mask = logical(zeros(size(Zmat)));
for i=1:length(Z);
    mask(:,:,i) = Zmat(:,:,i)<base & Zmat(:,:,i)>bed;
end;




%pcolor(X,Y,bed); shading flat; colorbar;
%hold on;
%plot(x(1,:),y(1,:),'k')
%plot(x(:,1),y(:,1),'k')
%contour(X,Y,zmask); 

%xlim([-2e6 -1.2e6]);
%ylim([-8e5 -1e5]);
%caxis([-1800 0]);


%plot(x_dt,[y_pitw(1) y_pitw(1)],'r');
%plot(x_dt,[y_pitw(2) y_pitw(2)],'r');
%plot([x_dt(1) x_pitw(1)],y_pitw,'r');
%plot([x_dt(2) x_pitw(2)],y_pitw,'r');

%plot(x_dt,[y_dt(1) y_dt(1)],'r');
%plot(x_dt,[y_dt(2) y_dt(2)],'r');
%plot([x_dt(1) x_dt(1)],y_dt,'r');
%plot([x_dt(2) x_dt(2)],y_dt,'r');
%hold off



years = 2005:2299;

SleftAnomslice = zeros(length(y_mesh_oce_mid),length(znaught),length(years));
TleftAnomslice = zeros(length(y_mesh_oce_mid),length(znaught),length(years));
SbotAnomslice = zeros(length(x_mesh_oce_mid),length(znaught),length(years));
TbotAnomslice = zeros(length(x_mesh_oce_mid),length(znaught),length(years));



X = ncread('/home/dgoldber/scratch/globus/ccsm4_RCP85/1995-2300/CCSM4_RCP85_salinity_8km_x_60m_clipped_1995.nc','x');
Y = ncread('/home/dgoldber/scratch/globus/ccsm4_RCP85/1995-2300/CCSM4_RCP85_salinity_8km_x_60m_clipped_1995.nc','y');
Z = ncread('/home/dgoldber/scratch/globus/ccsm4_RCP85/1995-2300/CCSM4_RCP85_salinity_8km_x_60m_clipped_1995.nc','z');
[X Y Z] = meshgrid(X,Y,Z);
[xleft yleft zleft] = meshgrid(x_mesh_oce_mid(1),y_mesh_oce_mid,znaught);
[xbot ybot zbot] = meshgrid(x_mesh_oce_mid,y_mesh_oce_mid(1),znaught);


folder='cesm2_ssp585/1995-2299';
for i=1:length(years);
%  [scol, tcol] = get_profile(folder, years(i), X, Y, Z, x, y, I_pitw, mask);
%  Scol_pitw(:,i) = scol;
%  Tcol_pitw(:,i) = tcol;
%  [scol, tcol] = get_profile(folder, years(i), X, Y, Z, x, y, I_dt, mask);
%  Scol_dt(:,i) = scol;
%  Tcol_dt(:,i) = tcol;
%  i
    
    disp(num2str(years(i)))
    salt = ncread(['/home/dgoldber/scratch/globus/cesm2_ssp585/1995-2299/CESM2-WACCM_SSP585_salinity_8km_x_60m_clipped_' num2str(years(i)) '.nc'],'salinity');
    theta = ncread(['/home/dgoldber/scratch/globus/cesm2_ssp585/1995-2299/CESM2-WACCM_SSP585_temperature_8km_x_60m_clipped_' num2str(years(i)) '.nc'],'temperature');  
    salt(~mask)=nan;
    theta(~mask)=nan;

    SleftAnomslice(:,:,i) = squeeze(griddata(X(mask),Y(mask),Z(mask),salt(mask),xleft,yleft,zleft,'nearest')) - (i>1)*SleftAnomslice(:,:,1);
    TleftAnomslice(:,:,i) = squeeze(griddata(X(mask),Y(mask),Z(mask),theta(mask),xleft,yleft,zleft,'nearest')) - (i>1)*TleftAnomslice(:,:,1);
    SbotAnomslice(:,:,i) = squeeze(griddata(X(mask),Y(mask),Z(mask),salt(mask),xbot,ybot,zbot,'nearest')) - (i>1)*SbotAnomslice(:,:,1);
    TbotAnomslice(:,:,i) = squeeze(griddata(X(mask),Y(mask),Z(mask),theta(mask),xbot,ybot,zbot,'nearest')) - (i>1)*TbotAnomslice(:,:,1);

    %salt_slice_left = interp3(X,Y,Z,salt,xleft,yleft,zleft);
    %theta_slice_left = interp3(X,Y,Z,theta,xleft,yleft,zleft);
    %salt_slice_bot = interp3(X,Y,Z,salt,xbot,ybot,zbot);
    %theta_slice_left = interp3(X,Y,Z,theta,xbot,ybot,zbot);
    
end


SleftAnomslice(:,:,1) = 0;
TleftAnomslice(:,:,1) = 0; 
SbotAnomslice(:,:,1) = 0;
TbotAnomslice(:,:,1) = 0;


save slices_anom_forcing.mat years SleftAnomslice TleftAnomslice  SbotAnomslice TbotAnomslice
