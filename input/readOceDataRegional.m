function readOceDataRegional (PAS, zMod, ystart, yend, calc_bounds, calc_init)

global sub_dir 

if (PAS>100)
     nc = ['/exports/geos.ed.ac.uk/iceocean/dgoldber/pahol_output/PAS_' num2str(PAS) '/run/'];
elseif (PAS>0 & PAS<=100)
     ncbase = ['/exports/geos.ed.ac.uk/iceocean/dgoldber/pahol_output/LENS_001/'];
     lens=true;
else
     nc = ['/home/dgoldber/network_links/ice_data/toshi_data/ASE/']
end


if (isempty(calc_init))
    calc_init = true;
end
if (isempty(calc_bounds))
    calc_bounds = true;
end

if (PAS>100);
lon = double(ncread([nc 'stateUvel.nc'],'LONGITUDE'));
lat = double(ncread([nc 'stateUvel.nc'],'LATITUDE'));
z = double(ncread([nc 'stateUvel.nc'],'DEPTH'));
else
lon = double(ncread([ncbase '2015.nc'],'XG'));
lat = double(ncread([ncbase '2015.nc'],'YG'));
z = double(ncread([ncbase '2015.nc'],'Z'));
end

load([sub_dir '/meshcoords.mat'],'x_mesh_mid', 'y_mesh_mid', 'x_mesh_oce_mid', 'y_mesh_oce_mid', ...
    'internal_grid_x','internal_grid_y','diffx','diffy');

[x_mesh_mid y_mesh_mid] = meshgrid(x_mesh_oce_mid,y_mesh_oce_mid);

delx = x_mesh_mid(2)-x_mesh_mid(1);
dely = y_mesh_mid(2)-y_mesh_mid(1);
x_mesh_mid = x_mesh_mid(1,:)';
y_mesh_mid = y_mesh_mid(:,1);

x_botbdry = x_mesh_mid;
y_botbdry = (y_mesh_mid(1)) * ones(length(x_mesh_mid),1);
x_lbdry = (x_mesh_mid(1)) * ones(length(y_mesh_mid),1);
y_lbdry = y_mesh_mid;
x_topbdry = x_mesh_mid;
y_topbdry = (y_mesh_mid(end)) * ones(length(x_mesh_mid),1);


start_year = ystart;
end_year = yend;
if (PAS>100)
 start_month = (start_year - 1955) * 12 
 end_month = (end_year - 1955 + 1) * 12 + 1
elseif (PAS>0 & PAS<=100)
 start_month = 1;
 end_month = (end_year-start_year+1)*12;
else
 start_month = (start_year - 1991) * 12
 end_month = (end_year - 1991 + 1) * 12 + 1
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (calc_init)

if (lens)
	nc = [ncbase num2str(start_year) '.nc'];
end

[xgridp ygridp zgridp] = meshgrid(x_mesh_mid,y_mesh_mid,zMod);

[phi,lambda]=polarstereo_inv(xgridp,ygridp,[],[],-71,0);
lambda = lambda + 360;

gradLat = zeros([size(phi(:,:,1)) 2]);
gradLon = zeros([size(phi(:,:,1)) 2]);
[mSz nSz o] = size(phi);

[latTempA lonTemp] = polarstereo_inv(xgridp(:,:,1)+100,ygridp(:,:,1),[],[],-71,0);
[latTempB lonTemp] = polarstereo_inv(xgridp(:,:,1)-100,ygridp(:,:,1),[],[],-71,0);
[latTempC lonTemp] = polarstereo_inv(xgridp(:,:,1),ygridp(:,:,1)+100,[],[],-71,0);
[latTempD lonTemp] = polarstereo_inv(xgridp(:,:,1),ygridp(:,:,1)-100,[],[],-71,0);
gradLat(:,:,1) = (latTempA-latTempB)/200;
gradLat(:,:,2) = (latTempC-latTempD)/200;
latNorm = sqrt(sum(gradLat.^2,3));
gradLat(:,:,1) = gradLat(:,:,1) ./ latNorm;
gradLat(:,:,2) = gradLat(:,:,2) ./ latNorm;

gradLon(:,:,1) = gradLat(:,:,2) ;
gradLon(:,:,2) = -gradLat(:,:,1);

ntot = 12;

Tinit = zeros(mSz,nSz,length(z));
Sinit = zeros(mSz,nSz,length(z));
Uinit = zeros(mSz,nSz,length(z));
Vinit = zeros(mSz,nSz,length(z));

[latInterp lonInterp] = meshgrid(lat,lon);

for n=1:ntot;
    
  start_month-1+n

  if (lens)
   nc = [ncbase num2str(start_year) '.nc'];
   Ttemp = double(ncread(nc,'THETA',[1 1 1 n],[600 384 50 1]));
   Stemp = double(ncread(nc,'SALT',[1 1 1 n],[600 384 50 1]));
  else
   Ttemp = double(ncread([nc 'stateTheta.nc'],'THETA',[1 1 1 start_month-1+n],[600 384 50 1]));
   Stemp = double(ncread([nc 'stateSalt.nc'],'SALT',[1 1 1 start_month-1+n],[600 384 50 1]));
  end
  disp('got here')

  
  Ttemp(Stemp==0)=nan;
  Vtemp(Stemp==0)=nan;
  Utemp(Stemp==0)=nan;
  Stemp(Stemp==0)=nan;
  
  for k=1:length(z)
  Tinit(:,:,k) = Tinit(:,:,k) + 1/ntot * interp2(lat,lon,Ttemp(:,:,k),phi(:,:,1),lambda(:,:,1));
  Sinit(:,:,k) = Sinit(:,:,k) + 1/ntot * interp2(lat,lon,Stemp(:,:,k),phi(:,:,1),lambda(:,:,1));
  %Uinit(:,:,k) = Uinit(:,:,k) + 1/ntot * interp2(lat,lon,Utemp(:,:,k),phi(:,:,1),lambda(:,:,1));
  %Vinit(:,:,k) = Vinit(:,:,k) + 1/ntot * interp2(lat,lon,Vtemp(:,:,k),phi(:,:,1),lambda(:,:,1));
  %Uinit(:,:,k) = inpaint_nans(Uinit(:,:,k),1);
  %Vinit(:,:,k) = inpaint_nans(Vinit(:,:,k),1);
  tlayer = Tinit(:,:,k);
  I = ~isnan(tlayer);
  disp(['month ' num2str(n) '; level ' num2str(k) ' ' num2str(mean(tlayer(I)))]);
  end
  
  
  
  
 end
% for k=1:length(z)
%   Tinit(:,:,k) = inpaint_nans(Tinit(:,:,k),1);
%   Sinit(:,:,k) = inpaint_nans(Sinit(:,:,k),1);
%   k
% end

% for k=1:length(z);
%     
%     uPolar = Uinit(:,:,n) .* gradLon(:,:,1) + Vinit(:,:,n) .* gradLat(:,:,1);
%     vPolar = Uinit(:,:,n) .* gradLon(:,:,2) + Vinit(:,:,n) .* gradLat(:,:,2);
%     Uinit(:,:,n) = uPolar;
%     Vinit(:,:,n) = vPolar;
%     
% end
    

save([sub_dir '/initpahol' num2str(PAS) '.mat'], 'Sinit', 'Tinit', 'z', 'xgridp', 'ygridp')


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[latps lonps] = meshgrid(lat,lon);
if (calc_bounds)

lbdryX = x_botbdry; lbdryY = y_botbdry;
Tbotbdry = zeros(length(z),length(lbdryX),end_month-start_month+1);
Sbotbdry = zeros(length(z),length(lbdryX),end_month-start_month+1);
Ubotbdry = zeros(length(z),length(lbdryX),end_month-start_month+1);
Vbotbdry = zeros(length(z),length(lbdryX),end_month-start_month+1);
[phiLeft,lambdaLeft]=polarstereo_inv(lbdryX,lbdryY,[],[],-71,0);
lambdaLeft = lambdaLeft + 360;

ntot = end_month-start_month+1

for n=1:ntot;
  if (PAS==0 & n>ntot-2)
   n2 = ntot-2;
  else
   n2 = n;
  end
  
  disp(['month ' num2str(n2)]);

  if (lens)
   year = start_year + floor((n-1)/12);
   mon = mod((n-1),12)+1;
   nc = [ncbase num2str(year) '.nc'];
   Ttemp = double(ncread(nc,'THETA',[1 1 1 mon],[600 384 50 1]));
   Stemp = double(ncread(nc,'SALT',[1 1 1 mon],[600 384 50 1]));
   Utemp = double(ncread(nc,'UVEL',[1 1 1 mon],[600 384 50 1]));
   Vtemp = double(ncread(nc,'VVEL',[1 1 1 mon],[600 384 50 1]));
  else
   Ttemp = double(ncread([nc 'stateTheta.nc'],'THETA',[1 1 1 start_month-1+n2],[600 384 50 1]));
   Stemp = double(ncread([nc 'stateSalt.nc'],'SALT',[1 1 1 start_month-1+n2],[600 384 50 1]));
   Utemp = double(ncread([nc 'stateUvel.nc'],'UVEL',[1 1 1 start_month-1+n2],[600 384 50 1]));
   Vtemp = double(ncread([nc 'stateVvel.nc'],'VVEL',[1 1 1 start_month-1+n2],[600 384 50 1]));
  end
  
  Ttemp(Stemp==0)=nan;
  Vtemp(Stemp==0)=nan;
  Utemp(Stemp==0)=nan;
  Stemp(Stemp==0)=nan;
      
    
  for k=1:length(z);
      
      Tbotbdry(k,:,n) = interp2(lat,lon,Ttemp(:,:,k),phiLeft,lambdaLeft);
      Sbotbdry(k,:,n) = interp2(lat,lon,Stemp(:,:,k),phiLeft,lambdaLeft);

      [xdummy ydummy Ups Vps] = vec_ll2ps(lonps,latps,Utemp(:,:,k),Vtemp(:,:,k),[],[]);
      Ubotbdry(k,:,n) = interp2(lat,lon,Ups,phiLeft,lambdaLeft);
      Vbotbdry(k,:,n) = interp2(lat,lon,Vps,phiLeft,lambdaLeft);

      
      disp(['bot ' num2str([k n Tbotbdry(k,1,n)])])
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%lbdryX = X(:,1)-500; lbdryY = Y(:,1);
lbdryX = x_lbdry; lbdryY = y_lbdry;
Tleftbdry = zeros(length(z),length(lbdryX),end_month-start_month+1);
Sleftbdry = zeros(length(z),length(lbdryX),end_month-start_month+1);
Uleftbdry = zeros(length(z),length(lbdryX),end_month-start_month+1);
Vleftbdry = zeros(length(z),length(lbdryX),end_month-start_month+1);
[phiLeft,lambdaLeft]=polarstereo_inv(lbdryX,lbdryY,[],[],-71,0);
lambdaLeft = lambdaLeft + 360;


ntot = end_month-start_month+1;

for n=1:ntot;
  if (PAS==0 & n>ntot-2)
   n2 = ntot-2;
  else
   n2 = n;
  end
  
  disp(['month ' num2str(n2)]);

  if (lens)
   year = start_year + floor((n-1)/12);
   mon = mod((n-1),12)+1;
   nc = [ncbase num2str(year) '.nc'];
   Ttemp = double(ncread(nc,'THETA',[1 1 1 mon],[600 384 50 1]));
   Stemp = double(ncread(nc,'SALT',[1 1 1 mon],[600 384 50 1]));
   Utemp = double(ncread(nc,'UVEL',[1 1 1 mon],[600 384 50 1]));
   Vtemp = double(ncread(nc,'VVEL',[1 1 1 mon],[600 384 50 1]));
  else

  
  Ttemp = double(ncread([nc 'stateTheta.nc'],'THETA',[1 1 1 start_month-1+n2],[600 384 50 1]));
  Stemp = double(ncread([nc 'stateSalt.nc'],'SALT',[1 1 1 start_month-1+n2],[600 384 50 1]));
  Utemp = double(ncread([nc 'stateUvel.nc'],'UVEL',[1 1 1 start_month-1+n2],[600 384 50 1]));
  Vtemp = double(ncread([nc 'stateVvel.nc'],'VVEL',[1 1 1 start_month-1+n2],[600 384 50 1]));
  end
  
  Ttemp(Stemp==0)=nan;
  Vtemp(Stemp==0)=nan;
  Utemp(Stemp==0)=nan;
  Stemp(Stemp==0)=nan;
      
  for k=1:length(z);
%     for i=1:length(lbdryX);
      
      Tleftbdry(k,:,n) = interp2(lat,lon,Ttemp(:,:,k),phiLeft,lambdaLeft);
      Sleftbdry(k,:,n) = interp2(lat,lon,Stemp(:,:,k),phiLeft,lambdaLeft);
      [xdummy ydummy Ups Vps] = vec_ll2ps(lonps,latps,Utemp(:,:,k),Vtemp(:,:,k),[],[]);
      Uleftbdry(k,:,n) = interp2(lat,lon,Ups,phiLeft,lambdaLeft);
      Vleftbdry(k,:,n) = interp2(lat,lon,Vps,phiLeft,lambdaLeft);
      
      disp(['left ' num2str([k n Tleftbdry(k,1,n)])])
%     end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% lbdryX = x_topbdry; lbdryY = y_topbdry;
% Ttopbdry = zeros(length(z),length(lbdryX),end_month-start_month+3);
% Stopbdry = zeros(length(z),length(lbdryX),end_month-start_month+3);
% Utopbdry = zeros(length(z),length(lbdryX),end_month-start_month+3);
% Vtopbdry = zeros(length(z),length(lbdryX),end_month-start_month+3);
% [phiLeft,lambdaLeft]=polarstereo_inv(lbdryX,lbdryY,[],[],-71,0);
% lambdaLeft = lambdaLeft + 360;
% 
% ntot = end_month-start_month+1;
% 
% for n=1:ntot;
%   if (PAS==0 & n>ntot-2)
%    n2 = ntot-2;
%   else
%    n2 = n;
%   end  
%   Ttemp = double(ncread([nc 'stateTheta.nc'],'THETA',[1 1 1 start_month-1+n2],[600 384 50 1]));
%   Stemp = double(ncread([nc 'stateSalt.nc'],'SALT',[1 1 1 start_month-1+n2],[600 384 50 1]));
%   Utemp = double(ncread([nc 'stateUvel.nc'],'UVEL',[1 1 1 start_month-1+n2],[600 384 50 1]));
%   Vtemp = double(ncread([nc 'stateVvel.nc'],'VVEL',[1 1 1 start_month-1+n2],[600 384 50 1]));
%   
%   Ttemp(Stemp==0)=nan;
%   Vtemp(Stemp==0)=nan;
%   Utemp(Stemp==0)=nan;
%   Stemp(Stemp==0)=nan;
%       
%   for k=1:length(z);
% %     for i=1:length(lbdryX);
%       
%       Ttopbdry(k,:,n) = interp2(lat,lon,Ttemp(:,:,k),phiLeft,lambdaLeft);
%       Stopbdry(k,:,n) = interp2(lat,lon,Stemp(:,:,k),phiLeft,lambdaLeft);
%       [xdummy ydummy Ups Vps] = vec_ll2ps(lonps,latps,Utemp(:,:,k),Vtemp(:,:,k),[],[]);
%       Utopbdry(k,:,n) = interp2(lat,lon,Ups,phiLeft,lambdaLeft);
%       Vtopbdry(k,:,n) = interp2(lat,lon,Vps,phiLeft,lambdaLeft);
%       
%       disp(['left ' num2str([k n Ttopbdry(k,1,n)])])
% %     end
%   end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([sub_dir '/bdryDatapahol' num2str(PAS) '.mat'],'Tleftbdry','Sleftbdry', 'Uleftbdry', ...
     'Vleftbdry', 'Tbotbdry', 'Sbotbdry', 'Ubotbdry', 'Vbotbdry', 'ystart', 'yend','z')

end
