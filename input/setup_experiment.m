function setup_experiment(nx,ny,gx,gy,use_2004,PAS,weert,dig_depth,gamma_depth,gammaT,cd);

global sub_dir

deltaT_tc = 31104000/12;
i_deltaT_tc = 31104000/deltaT_tc;

addpath('/exports/csce/datastore/geos/users/dgoldber/ice_data/ThwaitesDataProphet/CODE');
addpath('/exports/csce/datastore/geos/users/dgoldber/ice_data/bamber');

cpom_file = '/home/dgoldber/network_links/ice_data/cryosat_data/dhdt_cpom/antarctic_dhdt_5km_grid_1992_2017.nc'

pwd_dir=[pwd '/'];
streamicedatafile = [pwd_dir sub_dir '/ice/data.streamice'];
streamiceocedatafile = [pwd_dir sub_dir '/oce/data.streamice'];
datafile = [pwd_dir sub_dir '/ice/data'];
ocedatafile = [pwd_dir sub_dir '/oce/data'];

if (weert)
	chars='w'
	chars2='weert'
else
	chars='c'
	chars2='coul'
end

create_inputs_snap = false;
create_inputs_couplesnap = false;
create_inputs_coupletc = false;

if (use_2004)
 sub_dir = 'start_2004_input'
 oce_outdir = ['/home/dgoldber/network_links/geosIceOcean/dgoldber/archer_output/AMUND_COUPLE/OCE_SPINUP/run_2004_' num2str(PAS) '_' ...
     num2str(gamma_depth) '_' num2str(cd) '_' num2str(gammaT) '_' num2str(dig_depth)]
 ice_optimdir = ['/home/dgoldber/network_links/geosIceOcean/dgoldber/archer_output/AMUND_COUPLE/ICE_INIT/run_ad_2004']
 optimnumber = 200;
 s = dir([ice_optimdir '/runoptiter*']);
 if (length(s)>0);
  create_inputs_couplesnap = true;
  sname = s(end).name;
  optimnumbercouplesnap=str2num(sname(end-2:end));
 end
 ice_optimdirtc = ['/home/dgoldber/network_links/geosIceOcean/dgoldber/archer_output/AMUND_COUPLE/ICE_INIT_TC/run_ad_2004_' num2str(PAS) '_' chars 'NoCpomErr4yr']
 s = dir([ice_optimdirtc '/runoptiter*'])
 if (length(s)>0);
  create_inputs_coupletc = true;
  sname = s(end).name;
  optimnumbercoupletc=str2num(sname(end-2:end));
 end
else
 sub_dir = 'start_2009_input'
 oce_outdir = ['/home/dgoldber/network_links/geosIceOcean/dgoldber/archer_output/AMUND_COUPLE/OCE_SPINUP/run_2009_' num2str(PAS) '_' ...
     num2str(gamma_depth) '_' strrep(num2str(cd),'0.','.') '_' strrep(num2str(gammaT),'0.','.') '_' num2str(dig_depth)]
 if (weert)
     ice_optimdir = ['/home/dgoldber/network_links/geosIceOcean/dgoldber/archer_output/AMUND_COUPLE/ICE_INIT/run_ad_2009_weert']
 else
     ice_optimdir = ['/home/dgoldber/network_links/geosIceOcean/dgoldber/archer_output/AMUND_COUPLE/ICE_INIT/run_ad_2009_coul']
 end
 optimnumber = 200;
 s = dir([ice_optimdir '/runoptiter*'])
 if (length(s)>0);
  create_inputs_couplesnap = true;
  sname = s(end).name;
  optimnumbercouplesnap=str2num(sname(end-2:end));
 end
 ice_optimdirtc = ['/home/dgoldber/network_links/geosIceOcean/dgoldber/archer_output/AMUND_COUPLE/ICE_INIT_TC/run_ad_2009_' num2str(PAS) '_myiter_' chars2]
 s = dir([ice_optimdirtc '/runoptiter*'])
 if (length(s)>0);
  create_inputs_coupletc = true;
  sname = s(end).name;
  optimnumbercoupletc=str2num(sname(end-2:end));
 end
end
disp(['optimnumbercoupletc:'  num2str(optimnumbercoupletc) ' ' ice_optimdirtc]);

disp ('gx :')
disp(gx)
disp('gy :')
disp(gy)

useBM=false

load([sub_dir '/meshcoords.mat'],'x_mesh_mid', 'y_mesh_mid', 'x_mesh_oce_mid', 'y_mesh_oce_mid', ...
    'internal_grid_x','internal_grid_y','diffx','diffy');

density_ice = 917;
density_oce = 1027;


replace_param(streamicedatafile, 'streamice_density_ocean_avg', density_oce);


air_depth = 0;

%nx=360;
%ny=480;

x_mesh_d = 1000;
y_mesh_d = 1000;

[X Y] = meshgrid(x_mesh_mid,y_mesh_mid);

clear I J x y x_mesh x_mesh_mid y_mesh y_mesh_mid spd

if (use_2004)

    [vx vy]= interpMouginotAmund(X,Y,8,1);
    [v_err]= interpMouginotAmund(X,Y,8,2);
    
else
   
    %[vx vy]= interpMouginotAnt2017(X,Y,1);
    %[v_err]= interpMouginotAnt2017(X,Y,2);
    load ../../../amund/exps_thwaites_remove_shelf/input/PROPHET_Inversion_VelocityMeasurements.mat
    vx = VX(X,Y);
    vy = VY(X,Y);
    sdx = SDX(X,Y);
    sdy = SDY(X,Y);
    v_err = sqrt(sdx.^2+sdy.^2);
    
end
    
if (use_2004)
    dhdt1=ncread(cpom_file,'dhdt_2002_2006')';
    err1=ncread(cpom_file,'uncert_2002_2006')';
%    err = .5 *(err1+err2)';
%    dhdt = .5 *(dhdt1+dhdt2)';
    err = err1;
    dhdt = dhdt1;
    gaussFilter = fct_GaussianFilter([4 4], 1, 0);
    [err,im_conv,count,NaNcount] = fct_convNaN(err, gaussFilter, 'same', .2);
    gaussFilter = fct_GaussianFilter([2 2], 1, 0);
    [dhdt,im_conv,count,NaNcount] = fct_convNaN(dhdt, gaussFilter, 'same', .5);
    err(:,:)=1;
    
else
    dhdt1=ncread(cpom_file,'dhdt_2007_2011');
    dhdt2=ncread(cpom_file,'dhdt_2012_2016');
    err1=ncread(cpom_file,'uncert_2007_2011');
    err2=ncread(cpom_file,'uncert_2012_2016');
    err = .5 *(err2)';
    dhdt = (dhdt2)';
    err(:,:) = 1;
end
xcpom = ncread(cpom_file,'x');
ycpom = ncread(cpom_file,'y');
size(dhdt)
size(xcpom)
size(ycpom)
dhdtmap=interp2(xcpom,ycpom,dhdt,X,Y);
dhdterrmap=interp2(xcpom,ycpom,err,X,Y);

dhdtmap(isnan(dhdterrmap))=nan;
if(use_2004)
dhdtmap(dhdtmap>-.3)=nan;
else
dhdtmap(dhdtmap>-.5 | dhdtmap<-10)=nan;
end
%pcolor(dhdtmap); shading flat; colorbar; caxis([-10 10]);
%dhdtmap(isnan(dhdtmap))=-9999;
%dhdterrmap(isnan(dhdterrmap))=-9999;

bed = interpBedmachineAntarctica(X,Y,'bed');
thick = interpBedmachineAntarctica(X,Y,'thickness');
mask0 = interpBedmachineAntarctica(X,Y,'mask');
firn = interpBedmachineAntarctica(X,Y,'firn');
geoid = interpBedmachineAntarctica(X,Y,'geoid');

if (use_2004)
    open_bamber
    surf = interp2(bamber_x,bamber_y,bamber_dem,X,Y,'linear');
    mask = zeros(size(surf));
    mask(surf>0)=2;
    mask(mask0==1)=1;
    mask(surf<0)=0;
    surf(surf==0)=nan;
    surf(surf<-8)=nan;
else
    surf = interpBedmachineAntarctica(X,Y,'surface');
    mask = mask0;
end
   
if(~use_2004)
    surf2 = surf;
    surf2(thick<5 & surf2>0)=.5;
    surf2(surf2==0 | mask==1 | mask==0)=nan;
    gaussFilter = fct_GaussianFilter([1 1], 1, 0);
    [surf2,im_conv,count,NaNcount] = fct_convNaN(surf2, gaussFilter, 'same', .5);
    surf2(mask==0 | mask==1)=nan;
    surf2(isnan(surf2))=surf(isnan(surf2));
    surf = surf2;
%    surf(isnan(surf))=0;
end


    



if (use_2004)
hasfirn = (firn>0);
nofirn = (firn==0 & ~isnan(surf));
firn(nofirn) = InvDistWeighting(X(hasfirn),Y(hasfirn),firn(hasfirn),...
    X(nofirn),Y(nofirn),2,5e4);
surfadj = surf - firn - geoid;
surfadj(isnan(surf)) = 0;
else
surfadj = surf;
end

thick_floatation = (density_ice * air_depth - density_oce*surfadj) / (density_ice - density_oce);
thick_floatation(surfadj==0)=0;
base_floatation = surfadj - thick_floatation;
base = max(bed,base_floatation);
thick_mod = surfadj-base;
thick_mod(surfadj==0)=0;
mask(thick_mod<0 & mask~=1) = 1;
save surfadj.mat surfadj

%%%%%%% DIG HERE %%%%%%%%%%%%%


bed2 = bed;

if(dig_depth>0);
for i=1:ny;
 for j=1:nx;
  if (X(i,j)>-1.65e6 & X(i,j)<1.548e6 & Y(i,j)<-2e5 & Y(i,j)>-3.5e5);
   if (base(i,j)>bed(i,j)+10);
    rlow=bed(i,j);
    rshelf=base(i,j);
    bedmin = min(rlow,dig_depth);
    for k=[-1 1];
     if (base(i+k,j)>bed(i+k,j)+10);
      bedmin = min(bedmin,base(i+k,j)-dig_depth);
     end
     if (base(i,j+k)>bed(i,j+k)+10);
      bedmin = min(bedmin,base(i,j+k)-dig_depth);
     end
    end
    bed2(i,j) = bedmin;
   end  
  end
 end
end
end
       
save temp.mat bed bed2
bed = bed2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x_m=[];
y_m=[];
load /exports/csce/datastore/geos/users/dgoldber/ice_data/ThwaitesDataProphet/RACMO2.3/basin_Antarctica_1km_RACMO2
smb = InterpFromGrid(.5*(x_m(1:end-1)+x_m(2:end)),.5*(y_m(1:end-1)+y_m(2:end)), ...
    accumulation,double(X),double(Y),'cubic');



Xhere = X;
Yhere = Y;
load /exports/geos.ed.ac.uk/iceocean/dgoldber/pattyn/Temp_2013
load /exports/geos.ed.ac.uk/iceocean/dgoldber/pattyn/Zeta
Zeta = ones(1,1,length(zeta));
Zeta(1,1,:) = zeta;
Zeta = repmat(Zeta,[size(temp507,1) size(temp507,2) 1]);
Xp = X;
Yp = Y;

X = Xhere;
Y = Yhere;

Aglen = apaterson(temp507-273.15) * 31104000;
Bglen = Aglen.^(-1/3);

Bbar = .5*sum((Bglen(:,:,1:end-1)+Bglen(:,:,2:end)).*diff(Zeta,1,3),3);

Bbar = sqrt(InterpFromGrid(Xp(1,:)'*1000,Yp(:,1)*1000,Bbar,double(X),double(Y),'linear'));
Bbar2 = Bbar;
Bbar(isnan(Bbar))=0;
Bbar(Bbar==0)=700;
Bbar2(isnan(Bbar2))=0;


hmask = ones(ny,nx);


oce_mask = (surf==0);
hmask(mask==1)=-1;
hmask(mask==0)=0;
hmask([1 end],:) = -1;
hmask(:,[1 end]) = -1;
hmask(thick_mod<10 & (hmask==1 | hmask==2))= -1;
hmask(surf>2000)=-1;

% hmask = ones(ny,nx);
% hmask(thick==0 & bed<0) = 0;
% hmask(thick<10 & bed>=0) = -1;
% hmask(surf>2200)=-1;
% hmask([1 end],:) = -1;
% hmask(:,[1 end]) = -1;
% hmask(mask==1)=-1;



%hmask(192,74)=0;
%nytemp = max(find(Y(:,1)-Y(1,1)<1.48e5));
%s1 = surf(1:nytemp,:); 
%h1 = hmask(1:nytemp,:); 
%h1(s1>2000)=-1; 
%hmask(1:nytemp,:)=h1;
%s1 = surf(27:55,320:340); h1 = hmask(27:55,320:340); h1(s1>2000)=-1; hmask(27:55,320:340)=h1;


% hmask(I==1) = -1;



load([sub_dir '/meshcoords.mat'],'x_mesh_mid', 'y_mesh_mid', 'x_mesh_oce_mid', 'y_mesh_oce_mid', ...
    'internal_grid_x','internal_grid_y','diffx','diffy');
%x_mesh_d = 1000;
%y_mesh_d = 1000;

faketopog = -1000*ones(ny,nx);
faketopog([1 end],:) = 0;
faketopog(:,[1 end]) = 0;

con_mask = thick_mod~=0 & hmask==1;  
CC = bwconncomp(con_mask,4);
pp = CC.PixelIdxList;
for i=1:length(pp);
    length_pp(i) = length(pp{i});
end
[xx isort] = sort(length_pp,'descend');
if (length(pp)>1);
 for i=2:length(pp);
  for k=1:length(pp{isort(i)});
   hmask(pp{isort(i)}(k)) = -1;
  end;
 end;
end;

faketopog(hmask~=1) = 10;

if(use_2004)
snip = hmask(180:230,100:150); snip(snip==-1)=0;
hmask(180:230,100:150)=snip;
end

%%% melt rate from oce model  %%%%%%%

avgmelt = zeros(ny,nx);

s = dir([oce_outdir '/surfDiag*']);
if (length(s)>0);

 for i=((311040*2+25920):25920:(11*311040));
    i
    ml=rdmds([oce_outdir '/surfDiag'],i,'rec',2)';
    avgmelt(internal_grid_y,internal_grid_x(1:size(ml,2))) = ...
        avgmelt(internal_grid_y,internal_grid_x(1:size(ml,2))) + ...
        -31104000/917 * ml / (9 * 12);
 end
end

% todo

%%% snapshot parameters

s = dir([ice_optimdir '/runoptiter' appNum(optimnumber,3)]);

if (length(s)>0); 
create_inputs_snap = true;
beta0 = rdmds([ice_optimdir '/runoptiter' appNum(optimnumber,3) '/C_basal_fric']);
betap = rdmds([ice_optimdir '/runoptiter' appNum(optimnumber,3) '/xx_genarr2d2']);
Bglen0 = rdmds([ice_optimdir '/runoptiter' appNum(optimnumber,3) '/B_glen_sqrt']);
Bglenp = rdmds([ice_optimdir '/runoptiter' appNum(optimnumber,3) '/xx_genarr2d1']);
betaS = beta0 + betap;
BglenS = Bglen0 + Bglenp;
q = rdmds([ice_optimdir '/runoptiter' appNum(optimnumber,3) '/land_ice']);
grd = q(:,:,5);
end
%betaS(grd==0)=80;

if (create_inputs_couplesnap);
beta0 = rdmds([ice_optimdir '/runoptiter' appNum(optimnumbercouplesnap,3) '/C_basal_fric']);
betap = rdmds([ice_optimdir '/runoptiter' appNum(optimnumbercouplesnap,3) '/xx_genarr2d2']);
Bglen0 = rdmds([ice_optimdir '/runoptiter' appNum(optimnumbercouplesnap,3) '/B_glen_sqrt']);
Bglenp = rdmds([ice_optimdir '/runoptiter' appNum(optimnumbercouplesnap,3) '/xx_genarr2d1']);
betaSnap = beta0 + betap;
BglenSnap = Bglen0 + Bglenp;
betaSnap(grd==0) = 1.e-6;
end

if (create_inputs_coupletc);
beta0 = rdmds([ice_optimdirtc '/runoptiter' appNum(optimnumbercoupletc,3) '/C_basal_fric']);
betap = rdmds([ice_optimdirtc '/runoptiter' appNum(optimnumbercoupletc,3) '/xx_beta'],optimnumbercoupletc);
Bglen0 = rdmds([ice_optimdirtc '/runoptiter' appNum(optimnumbercoupletc,3) '/B_glen_sqrt']);
Bglenp = rdmds([ice_optimdirtc '/runoptiter' appNum(optimnumbercoupletc,3) '/xx_bglen'],optimnumbercoupletc);
betaTC = beta0 + betap;
BglenTC = Bglen0 + Bglenp;
betaTC(grd==0) = 1.e-6;
end


%%%%

H = [[smb zeros(ny,gx)];zeros(gy,nx+gx)];
write_coupled_input('SMB_racmo.bin',H','ice');
replace_param(streamicedatafile, 'STREAMICEAdotFile', 'SMB_racmo.bin');

H = [[bed zeros(ny,gx)];zeros(gy,nx+gx)];
if (dig_depth>0);
 write_coupled_input(['topog_dig_' num2str(dig_depth) '.bin'],H','ice');
else
 write_coupled_input('topog.bin',H','ice');
end
replace_param(streamicedatafile, 'STREAMICEtopogFile', 'topog.bin');

H = [[faketopog zeros(ny,gx)];zeros(gy,nx+gx)];
write_coupled_input('faketopog.bin',H','ice');
replace_param(datafile, 'bathyFile', 'faketopog.bin');

H = [[thick_mod zeros(ny,gx)];zeros(gy,nx+gx)];
write_coupled_input('BedMachineThickMod.bin',H','ice');
replace_param(streamicedatafile, 'streamicethickFile', 'BedMachineThickMod.bin');

vx(isnan(vx))=-9999;
vy(isnan(vy))=-9999;
v_err(isnan(v_err))=-9999;

H = [[vx zeros(ny,gx)];zeros(gy,nx+gx)];
write_coupled_input('velobsu.bin',H','ice');

H = [[vy zeros(ny,gx)];zeros(gy,nx+gx)];
write_coupled_input('velobsv.bin',H','ice');

H = [[v_err zeros(ny,gx)];zeros(gy,nx+gx)];
write_coupled_input('velobserr.bin',H','ice');

%%%%

H = [[vx zeros(ny,gx)];zeros(gy,nx+gx)];
write_coupled_input(['velobs' appNum(i_deltaT_tc,10) 'u.bin'],H','ice');

H = [[vy zeros(ny,gx)];zeros(gy,nx+gx)];
write_coupled_input(['velobs' appNum(i_deltaT_tc,10) 'v.bin'],H','ice');

H = [[v_err zeros(ny,gx)];zeros(gy,nx+gx)];
write_coupled_input(['velobs' appNum(i_deltaT_tc,10) 'err.bin'],H','ice');

%%%%

Bbar = [[Bbar zeros(ny,gx)];zeros(gy,nx+gx)];
write_coupled_input('BglenPattyn.bin',(Bbar)','ice');
replace_param(streamicedatafile, 'STREAMICEglenconstfile', 'BglenPattyn.bin');

H = [[Bbar2 zeros(ny,gx)];zeros(gy,nx+gx)];
write_coupled_input('BglenPattynMask.bin',(H)','ice');

write_coupled_input('ufacemask.bin',-1*ones(size(H))','ice');
write_coupled_input('vfacemask.bin',-1*ones(size(H))','ice');
replace_param(streamicedatafile, 'STREAMICEuFaceBdryFile','ufacemask.bin');
replace_param(streamicedatafile, 'STREAMICEvFaceBdryFile','vfacemask.bin');

write_coupled_input('delX_ice.bin',[diffx ones(1,gx)],'ice')
write_coupled_input('delY_ice.bin',[diffy ones(1,gy)],'ice')
replace_param(datafile, 'delxfile','delX.bin');
replace_param(datafile, 'delyfile','delY.bin');

H = [[hmask -1*ones(ny,gx)];-1*ones(gy,nx+gx)];
write_coupled_input('hmask.bin',H','ice');
replace_param(streamicedatafile, 'streamicehmaskfile','hmask.bin');

% v = [[dhdtmap zeros(ny,gx)];zeros(gy,nx+gx)];
% write_coupled_input('dhdtcpom.bin',(v)','ice');
% v = [[dhdterrmap zeros(ny,gx)];zeros(gy,nx+gx)];
% write_coupled_input('dhdtcpomerr.bin',(v)','ice');

for i=i_deltaT_tc:i_deltaT_tc:(8*i_deltaT_tc);
    surftemp = surfadj + dhdtmap*(i/i_deltaT_tc);
    surftemp(isnan(surftemp))=-999999;
    err(isnan(surftemp))=-999999;

    v = [[surftemp zeros(ny,gx)];zeros(gy,nx+gx)];
    write_coupled_input(['dhdtcpom' appNum(i,10) '.bin'],(v)','ice');
    v = [[err zeros(ny,gx)];zeros(gy,nx+gx)];
    write_coupled_input(['dhdtcpom' appNum(i,10) 'err.bin'],(v)','ice');

end

v = [[avgmelt zeros(ny,gx)];zeros(gy,nx+gx)];
write_coupled_input(['avgmelt_spinup_' num2str(PAS) '.bin'],(v)','ice');


if (create_inputs_snap)
write_coupled_input(['BetaSnap_' num2str(optimnumber) '.bin'],betaS,'ice')
write_coupled_input(['BglenSnap_' num2str(optimnumber) '.bin'],BglenS,'ice')
end

if(chars=='w')
    sname = 'weert'
else
    sname = 'coul'
end

if (create_inputs_couplesnap)
write_coupled_input(['BetaSnap' sname '.bin'],betaSnap,'ice')
write_coupled_input(['BglenSnap' sname '.bin'],BglenSnap,'ice')
end

if (create_inputs_coupletc)
write_coupled_input(['BetaTC' sname num2str(PAS) '.bin'],betaTC,'ice')
write_coupled_input(['BglenTC' sname num2str(PAS) '.bin'],BglenTC,'ice')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set up here for streamice stub

bed = bed(internal_grid_y,internal_grid_x);
if (dig_depth>0);
	write_coupled_input(['bathy_init_oce_dig_' num2str(dig_depth) '.bin'],bed','oce');
else
	write_coupled_input('bathy_init_oce.bin',bed','oce');
end
replace_param(streamiceocedatafile, 'streamicetopogfile','bathy_init_oce.bin');

base = base(internal_grid_y,internal_grid_x);
write_coupled_input('icetopo_init_oce.bin',base','oce');

thick = thick_mod(internal_grid_y,internal_grid_x);
write_coupled_input('icethick_oce.bin',thick','oce');
replace_param(streamiceocedatafile, 'streamicethickfile','icethick_oce.bin');

Uice = 0*vx(internal_grid_y,internal_grid_x);
write_coupled_input('Uice_oce.bin',Uice','oce');
replace_param(streamiceocedatafile, 'streamice_uvel_ext_file','Uice_oce.bin');

Vice = 0*vy(internal_grid_y,internal_grid_x);
write_coupled_input('Vice_oce.bin',Vice','oce');
replace_param(streamiceocedatafile, 'streamice_vvel_ext_file','Vice_oce.bin');

Hbcx = 0*vx(internal_grid_y,internal_grid_x);
write_coupled_input('HBCx_oce.bin',Hbcx','oce');
replace_param(streamiceocedatafile, 'STREAMICEHBCxFile','HBCx_oce.bin');

Hbcy = 0*vy(internal_grid_y,internal_grid_x);
write_coupled_input('HBCy_oce.bin',Hbcy','oce');
replace_param(streamiceocedatafile, 'STREAMICEHBCyFile','HBCy_oce.bin');

hmask = hmask(internal_grid_y,internal_grid_x);
hmask(end,:)=-1;
hmask(:,end)=-1;
write_coupled_input('hmask_oce.bin',hmask','oce');
replace_param(streamiceocedatafile, 'streamicehmaskfile','hmask_oce.bin');

ufacemask = -1*ones(size(hmask));
vfacemask = -1*ones(size(hmask));
ufacemask(:,end) = 3;
vfacemask(end,:) = 3;
write_coupled_input('ufacemask_oce.bin',ufacemask','oce');
write_coupled_input('vfacemask_oce.bin',vfacemask','oce');
replace_param(streamiceocedatafile, 'STREAMICEuFaceBdryFile','ufacemask_oce.bin');
replace_param(streamiceocedatafile, 'STREAMICEvFaceBdryFile','vfacemask_oce.bin');

smb = smb(internal_grid_y,internal_grid_x);
write_coupled_input('SMB_racmo_oce.bin',smb','oce');
replace_param(streamiceocedatafile, 'STREAMICEAdotFile', 'SMB_racmo_oce.bin');


