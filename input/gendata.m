function gendata(PAS, nx_in,ny_in,zmod,gx,gy,mwct,rho0,ystart,dig_depth,gamma_depth,gammaT,cd);

global sub_dir
% Dimensions of grid
nx=nx_in
ny=ny_in
nz=length(zmod);
delz = zmod(1)-zmod(2);
time_step=100;

dirpwd=[pwd '/'];

ocedatafile = [dirpwd sub_dir '/oce/data'];
shelficedatafile = [dirpwd sub_dir '/oce/data.shelfice'];
obcsdatafile = [dirpwd sub_dir '/oce/data.obcs'];

if (ystart==2001)
  year_spinup = 2001;
  year_pickup = 2001;
  oce_outdir = ['/home/dgoldber/network_links/geosIceOcean/dgoldber/archer_output/AMUND_COUPLE/OCE_SPINUP/run_2004_1_' num2str(PAS)]
else
  year_spinup = 2008;
  year_pickup = 2013;
  oce_outdir = ['/home/dgoldber/network_links/geosIceOcean/dgoldber/archer_output/AMUND_COUPLE/OCE_SPINUP/run_2009_' num2str(PAS) '_' ...
     num2str(gamma_depth) '_' strrep(num2str(cd),'0.','.') '_' strrep(num2str(gammaT),'0.','.') '_' num2str(dig_depth)]
end

hfacMin = round(mwct/4/delz,4);

replace_param(ocedatafile, 'cg2dmincolumneps', mwct*2.5);
replace_param(ocedatafile, 'rhoconst', rho0);
replace_param(ocedatafile, 'hfacmin', hfacMin);
replace_param(ocedatafile, 'hfacinf', hfacMin);
replace_param(ocedatafile, 'delr', [num2str(nz) '*' num2str(delz)],false);

rho_ice = 917;



eos = 'jmd95z';
acc = 'real*8';
prec = 8;

if (dig_depth>0);
	bathyname = ['bathy_init_oce_dig_' num2str(dig_depth) '.bin']
	bathymodname = ['bathy_mod_dig_' num2str(dig_depth) '.bin']
else
	bathyname = 'bathy_init_oce.bin';
	bathymodname = 'bathy_mod.bin';
end

bathy0 = read_coupled_input(bathyname,'oce',nx+gx,ny+gx);

bathy = bathy0;

zgp1 = 0:-delz:(-nz*delz);
zc = .5 * (zgp1(1:end-1)+zgp1(2:end));
zg = zgp1(1:end-1);

T0 = read_coupled_input(['theta.init.' num2str(ystart)  '.' num2str(PAS)],'oce',nx+gx,ny+gy,nz);
S0 = read_coupled_input(['salt.init.' num2str(ystart)  '.' num2str(PAS)],'oce',nx+gx,ny+gy,nz);
tref = squeeze(T0(1,1,:));
sref = squeeze(S0(1,1,:));

% Gravity
gravity=9.81;
rhoConst = rho0;

% compute potential field underneath ice shelf
talpha = 2e-4;
sbeta  = 7.4e-4;
t    = tref;
s    = sref;
gravity = 9.81;
k=1;
dzm = abs([zg(1)-zc(1) .5*diff(zc)]);
dzp = abs([.5*diff(zc) zc(end)-zg(end)]);
p = abs(zc)*gravity*rhoConst*1e-4;
dp = p;
kp = 0;



while rms(dp) > 1e-13
  phiHydF(k) = 0;
  p0 = p;
  kp = kp+1
  for k = 1:nz
    switch eos
     case 'linear'
      drho = rhoConst*(1-talpha*(t(k)-tref(k))+sbeta*(s(k)-sref(k)))-rhoConst;
     case 'jmd95z'
      drho = densjmd95(s(k),t(k),p(k))-rhoConst;
      rho0(k) = drho;
     case 'mdjwf'
      drho = densmdjwf(s(k),t(k),p(k))-rhoConst;
     otherwise
      error(sprintf('unknown EOS: %s',eos))
    end
    phiHydC(k)   = phiHydF(k) + dzm(k)*gravity*drho/rhoConst;
    phiHydF(k+1) = phiHydC(k) + dzp(k)*gravity*drho/rhoConst;
  end
  switch eos
   case 'mdjwf'
    p = (gravity*rhoConst*abs(zc) + phiHydC*rhoConst)/gravity/rhoConst;
  end
  dp = p-p0;
end

p = (gravity*rhoConst*abs(zc) + phiHydC*rhoConst)/gravity;

thickTD = read_coupled_input('icethick_oce.bin','oce',nx+gx,ny+gy);
thick1=thickTD;
shelficemass=thick1*rho_ice;

masstest = zeros(nx+gx,ny+gy);
topo = zeros(nx+gx,ny+gy);
bathy(bathy>-delz)=0;
bathy2 = bathy;
for ix=1:nx
  for iy=1:ny
          cellFace = mod(abs(bathy2(ix,iy)),delz);
          if (cellFace<hfacMin*delz & cellFace>0);
              bathy2(ix,iy) = -floor(abs(bathy2(ix,iy))/delz)*delz;
          end
  end
end
bathy = bathy2;

for ix=1:(nx+gx)
  for iy=1:(ny+gy)

     mass = shelficemass (ix,iy);
     massFuncC = rhoConst * (phiHydC/gravity - zc);
    
     massFuncF = rhoConst * (phiHydF/gravity - zgp1);

     k = max (find ( massFuncF < mass ));
     if (k>nz)
         topo(ix,iy) = bathy(ix,iy)+mwct;
     else
        if (isempty(k))
         k=0;
        end
         
     
        if (k>0)
         if (mass < massFuncC(k))
          ztopo = zg(k) - (mass-massFuncF(k)) * delz/2 / (massFuncC(k)-massFuncF(k));
          topo(ix,iy) = max(ztopo,bathy(ix,iy)+mwct);
         else
          ztopo = zc(k) - (mass-massFuncC(k)) * delz/2 / (massFuncF(k+1)-massFuncC(k));
          topo(ix,iy) = max(ztopo,bathy(ix,iy)+mwct);
         end
        end
     end
     
  end
end


save([sub_dir '/masses.mat'], 'massFuncC', 'massFuncF', 'phiHydC'); 


etainit = zeros(size(topo));

% new topography: icetopo rounded to the nearest k * deltaZ
%                 eta_init set to make difference

icetopo2 = topo;
save([sub_dir '/topo.mat'],'topo')

for ix=1:nx
  for iy=1:ny
    k=max(find(abs(zg)<abs(icetopo2(ix,iy))));
    if isempty(k)
      k=0;
    else
%      k, zg(k),abs(icetopo2(ix,iy)) 
      dr = 1-(zg(k) - icetopo2(ix,iy))/delz;
      if (dr > .2 | ((zg(k+1)-bathy(ix,iy))<hfacMin))
          % bring Ro_surf *up* to closest grid face & make etainit negative
          % to compensate
          icetopo2(ix,iy) = zg(k);
          etainit(ix,iy) = (dr-1)*delz;
      else
          % bring Ro_surf *down* to closest grid face & make etainit pos
          % to compensate
          icetopo2(ix,iy) = zg(k+1);
          etainit(ix,iy) = (dr)*delz;
      end

       
    end
  end
end

topo = icetopo2;
save([sub_dir '/topo2.mat'], 'icetopo2', 'etainit');

%if (min_thickness>0);
% con_mask = bathy<0 & (topo-bathy)>0;
% CC = bwconncomp(con_mask,4);
% pp = CC.PixelIdxList;
% if (length(pp)>1);
%  for i=2:length(pp);
%   disp(['length: ' num2str(length(pp{i}))]);
%   for k=1:length(pp{i});
%    I(pp{i}(k)) = 0;
%   end;
%  end;
% end;
%end





bathy(bathy0>-20) = 0;
bathy(1,294:end)=0;
bathy(80:end,1)=0;

replace_param(obcsdatafile,'ob_jsouth',[num2str(79) '*1 ' num2str(nx-79) '*0,'],false)
replace_param(obcsdatafile,'ob_jwest',[num2str(293) '*1 ' num2str(ny-293) '*0,'],false)

bathy(end,:) = 0;
bathy(:,end) = 0;
% bathy = bathy';
% bathytemp = bathy(141:end,1:18);
% etatemp = etainit(141:end,1:18);
% topotemp = icetopo2(141:end,1:18);
% 
% I = (etatemp+topotemp-bathytemp > 1);
% bathytemp(I)=0;
% topotemp(I)=0;
% etatemp(I)=0;
% bathy(141:end,1:18)=bathytemp;
% etainit(141:end,1:18)=etatemp;
% icetopo2(141:end,1:18)=topotemp;
 

% bathy(2,141:end)=0;
%bathy(1:12,135:end)=0;


%bathy(bathy>-20)=0;
%bathy(107:end,1)=0;
%bathy(177:end,220)=0;
%bathy(:,end-gy:end)=0;
%bathy(end-gx:end,:)=0;
%bathy(104:end,1)=0;
%bathy(1,338:end)=0;


write_coupled_input(['etainit.' num2str(PAS) '.bin'] ,etainit,'oce');
write_coupled_input(['shelftopo.' num2str(PAS) '.bin'],topo,'oce');
write_coupled_input('shelficemassinit.bin',shelficemass,'oce');
write_coupled_input(bathymodname,bathy,'oce');

replace_param(ocedatafile, 'psurfinitfile', ['etainit.bin.' num2str(PAS)]);
replace_param(ocedatafile, bathymodname, 'bathy_mod.bin');
replace_param(shelficedatafile, 'shelficetopofile', ['shelftopo.bin.' num2str(PAS)]);
replace_param(shelficedatafile, 'shelficemassfile', 'shelficemassinit.bin');
r=(etainit+topo-bathy)'; shelfmask = double(topo'<0 & r>2); 
shelfmask(shelfmask==1)=2;

temp = shelfmask(280:end,:); 
temp(temp>0) = 3;
shelfmask(280:end,:)=temp;

temp = shelfmask(1:156,:); 
temp(temp>0) = 1;
shelfmask(1:156,:)=temp;

write_coupled_input('shelf_mask.bin',shelfmask','oce');

s = dir([oce_outdir '/pickup*']);

if (length(s)>1);

pickup_count = appNum((year_pickup-year_spinup)*31104000/time_step,10);

evalstr = ['!cp ' oce_outdir '/pickup.' pickup_count '* ' sub_dir '/oce'];
eval(evalstr)
evalstr = ['!cp ' oce_outdir '/pickup_streamice.' pickup_count '* ' sub_dir '/oce'];
eval(evalstr)
evalstr = ['!cp ' oce_outdir '/pickup_shelfice.' pickup_count '* ' sub_dir '/oce'];
eval(evalstr)


% q = binread([sub_dir '/oce/pickup_streamice.' pickup_count '.data'],8,nx,ny,92); 
% mask=q(:,:,85); 
% mask(:,end)=-1; mask(end,:)=-1; q(:,:,85)=mask; 
% binwrite([sub_dir '/oce/pickup_streamice.' pickup_count '.data'],q);

%q = binread([sub_dir '/oce/pickup_shelfice.0002488320.data'],8,nx,ny);
%q2 = zeros(nx,ny,2);
%q2(:,:,2) = q;
%q2(:,:,1) = shelficemass;
%binwrite([sub_dir '/oce/pickup_shelfice.0002488320.data'],q2);

end
