%This program creates open boundary prescription files for the PIG
%experiment based on the 10 yrs spin-up run using OBCS with 
%U,V = 0, Tref, Sref

%following variables 
%-----3D fields-----
% T Temperature (C)
% S Salinity (psu)
% U u-velocity (m/s)
% PH ocean pressure (or atm geopotential)

function rdmds_init(PAS,nx_in,ny_in,zmod,gx,gy,smoothnum,ystart_in,yend_in)

global sub_dir

ocedatafile = ['/exports/geos.ed.ac.uk/iceocean/dgoldber/MITgcm_forinput/amund_couple/mitgcm_amund_coupling_2300/input/' sub_dir '/oce/data'];
exfdatafile = ['/exports/geos.ed.ac.uk/iceocean/dgoldber/MITgcm_forinput/amund_couple/mitgcm_amund_coupling_2300/input/' sub_dir '/oce/data.exf'];
caldatafile = ['/exports/geos.ed.ac.uk/iceocean/dgoldber/MITgcm_forinput/amund_couple/mitgcm_amund_coupling_2300/input/' sub_dir '/oce/data.cal'];
obcsdatafile = ['/exports/geos.ed.ac.uk/iceocean/dgoldber/MITgcm_forinput/amund_couple/mitgcm_amund_coupling_2300/input/' sub_dir '/oce/data.obcs'];


close all

%Set the grid size;
nx = nx_in;    delx = 1;   X = nx*delx;
ny = ny_in;    dely = 1;   Y = ny*dely;
nz = length(zmod);     delz = zmod(1)-zmod(2);  Z = nz*delz;

load([sub_dir '/meshcoords.mat'])

%x axis 
x = x_mesh_mid(1,:);
%y axis
y = y_mesh_mid(:,1);
%z axis [m]
z = 0:delz:Z;
zmid = zmod;

%dzmod = [   10.0000   10.0000   10.0000   10.0000   12.0000   14.5000   17.5000   21.0000   25.0000   30.0000   35.0000  40.0000  45.0000  50 ...
%   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000 ...
%   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000 ...
%   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000   50.0000 50 50 50 50 50];
dzmod = delz*ones(1,nz);
zface = [0 cumsum(dzmod)];
%zmid = .5 * (zface(1:end-1)+zface(2:end));

%Note: 
%Ensure that the volume flux is zero at the open boundary



%% ________________________________________________________________________
%  First try: nearest velocity profile. 
%  From U = -0.025 ms-1 (outwards) to  0.025 ms-1  
%  Constant in y direction (meridionally)
%  Do not take exactly Umin and Umax from the OBC of the original
%  experiment: Take modified numbers to that Umin = -Umax and volume is conserved

% Initialize variable
V_OBC(1:nx,1:nz) = NaN;


v_sfc = -0.025;                         %U at sfc
v_bottom = 0.025;                       %U at bottom
del_v = (v_bottom - v_sfc)/1000;   %delu/delz
    
for iz = 1:nz;
    V_OBC(1:nx,iz) = v_sfc + del_v*((iz-1)*delz);    
end



%  PW Linear T profile

% Initialize variable
T_OBC(1:nx,1:nz) = NaN;

T_sfc = -1.9;


%for iz = 1:nz;
%    
%    if (zmid(iz) < 500)
%%%%        
%        T_OBC(1:nx,iz) = T_sfc + zmid(iz) * .0048;
%        
%    elseif (zmid>=500 & zmid<700)
%        
%        T_OBC(1:nx,iz) = T_sfc + 500 * .0048 + (zmid(iz)-500) * .0015;
%        
%    else
%        
%        T_OBC(1:nx,iz) = T_sfc + 500 * .0048 + (700-500) * .0015 + (zmid(iz)-700) * .000333;
%        
%    end
%end


%  PW Linear S profile

% Initialize variable
S_OBC(1:nx,1:nz) = NaN;

S_sfc = 34.2050;




%% Initial T & S conditions

% Take western open boundary conditions for T & S 
% and assume no change in x-direction


% for iy = 1:ny;
%     T_init(:,iy,:) = repmat(squeeze(T_OBC(1,:)),[nx 1]);
%     S_init(:,iy,:) = repmat(squeeze(S_OBC(1,:)),[nx 1]);
% end



%% Print OBCS files for T, S

% fid=fopen('uvel.obw','w','b');  fprintf(fid,'%10.4f',U_OBC);fclose(fid);
% fid=fopen('vvel.obs','w','b');  fwrite(fid,V_OBC,'real*8'); fclose(fid);
% fid=fopen('theta.obs','w','b');  fwrite(fid,T_OBC,'real*8'); fclose(fid);
% fid=fopen('salt.obs','w','b');  fwrite(fid,S_OBC,'real*8'); fclose(fid);
% 
% fid=fopen('uvel.obw','w','b');  fwrite(fid,V_OBC(1:135,:),'real*8'); fclose(fid);
% fid=fopen('theta.obw','w','b');  fwrite(fid,T_OBC(1:135,:),'real*8'); fclose(fid);
% fid=fopen('salt.obw','w','b');  fwrite(fid,S_OBC(1:135,:),'real*8'); fclose(fid);


% fid=fopen('theta.obw','w');  fprintf(fid,'%10.4f',T_OBC);fclose(fid);
% fid=fopen('salt.obw','w');  fprintf(fid,'%10.4f',S_OBC);fclose(fid);



%% Print init files for T, S




z2 = 0:delz:Z;
z2 = zface;

%load([sub_dir '/bdryDatapahol' num2str(PAS) '.mat']);
%load([sub_dir '/initpahol' num2str(PAS) '.mat']);
load ([sub_dir '/naughtenClim.mat']);
load ([sub_dir '/slices_anom_forcing.mat']);
load ([sub_dir '/bdryDatapahol30.mat']);

if (PAS<100)
z = -z;
end

if((ystart_in < year_array(1)) & (ystart_in~=0))
    error('regional data has wrong start year')
end

T_init = NaN(ny,nx,nz);
S_init = NaN(ny,nx,nz);

Tinit = T2005;
Sinit = S2005;

nx_forcing = size(Tinit,2);
ny_forcing = size(Tinit,1);

for i=1:ny_forcing; 
    for j=1:nx_forcing; 
        T_init(i,j,:) = interp1(z,squeeze(Tinit(i,j,:)),-zmid);
        S_init(i,j,:) = interp1(z,squeeze(Sinit(i,j,:)),-zmid);
    end
end

T_init_prof = zeros(1,1,nz);
S_init_prof = zeros(1,1,nz);
for k=1:nz;
    Tlayer = T_init(1:ny_forcing,1:nx_forcing,k);
    Slayer = S_init(1:ny_forcing,1:nx_forcing,k);
    T_init_prof(k) = mean(Tlayer(~isnan(Tlayer)));
    S_init_prof(k) = mean(Slayer(~isnan(Slayer)));
end

size(T_init_prof)
size(isnan(T_init_prof))

[m ii]=max(squeeze(S_init_prof))

T_init_prof((ii+1):end)=T_init_prof(1,1,ii);
S_init_prof((ii+1):end)=S_init_prof(1,1,ii);

T_init  = repmat(T_init_prof,[nx+gx ny+gy 1]);
S_init  = repmat(S_init_prof,[nx+gx ny+gy 1]);

%fid_T=fopen('theta_section.init','w','b');fwrite(fid_T,permute(T_init(14:17,401:403,:),[1 2 3]),'real*8');fclose(fid_T);
%fid_T=fopen('salt_section.init','w','b');fwrite(fid_T,permute(S_init(14:17,401:403,:),[1 2 3]),'real*8');fclose(fid_T);

write_coupled_input(['theta.init.' num2str(ystart_in) '.' num2str(PAS)],permute(T_init,[1 2 3]),'oce');
write_coupled_input(['salt.init.' num2str(ystart_in) '.' num2str(PAS)],permute(S_init,[1 2 3]),'oce');


    replace_param(ocedatafile, 'hydrogthetafile', ['theta.init.' num2str(ystart_in) '.' num2str(PAS)]);
    replace_param(ocedatafile, 'hydrogsaltfile', ['salt.init.' num2str(ystart_in) '.' num2str(PAS)]);



% for iz = 1:nz;
%     fprintf(fid_T,'%10.4f',T_init(:,:,iz));
%     fprintf(fid_S,'%10.4f',S_init(:,:,iz)); 
% end
% 
% fclose(fid_T);
% fclose(fid_S);

% S_col = squeeze(S_OBC(1,:));
% T_col = squeeze(T_OBC(1,:));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dstart = round(time2num(ystart_in-year_array(1)));

n_months = 12*(length(year_array) - (dstart));

z=double(-z);
z2=double(-z2);

V_south = Vbotbdry; 
U_south = Ubotbdry; 

V_southPAS = permute(V_south,[2 1 3]);
U_southPAS = permute(U_south,[2 1 3]);


% need USouth and VSoutn
% append 10-year cycle to end of series
PASyears = size(U_south,3)/12;
U_south = zeros(nx_forcing,length(z),n_months);
V_south = zeros(nx_forcing,length(z),n_months);
U_south(:,:,1:(12*PASyears-12)) = U_southPAS(:,:,1:(end-12));
V_south(:,:,1:(12*PASyears-12)) = V_southPAS(:,:,1:(end-12));
U_south(:,:,(12*PASyears-11):end) = repmat(U_southPAS(:,:,(end-119):end),[1 1 (1+length(year_array)-dstart-PASyears)/10]);
V_south(:,:,(12*PASyears-11):end) = repmat(V_southPAS(:,:,(end-119):end),[1 1 (1+length(year_array)-dstart-PASyears)/10]);


% now, to create S_south and T_south!
S_south = zeros(nx_forcing,length(z),n_months);
T_south = zeros(nx_forcing,length(z),n_months);

SbotClim = permute(SbotClim,[2 1 3]);
TbotClim = permute(TbotClim,[2 1 3]);

for i = 1:(length(year_array) - (ystart_in-year_array(1)));

    TbotAnomslice_i = TbotAnomslice(:,:,i+dstart);
    TbotAnomslice_im1 = TbotAnomslice(:,:,i+dstart-1);
    SbotAnomslice_i = SbotAnomslice(:,:,i+dstart);
    SbotAnomslice_im1 = SbotAnomslice(:,:,i+dstart-1);

    if (i==(length(year_array) - (ystart_in-year_array(1))));
        TbotAnomslice_ip1 = TbotAnomslice(:,:,i+dstart);
        SbotAnomslice_ip1 = SbotAnomslice(:,:,i+dstart);
    else
        TbotAnomslice_ip1 = TbotAnomslice(:,:,i+dstart+1);
        SbotAnomslice_ip1 = SbotAnomslice(:,:,i+dstart+1);
    end

    for j=1:12;

        if (j<=6);
           TbotAnom = (12-(6.5-j))/12*TbotAnomslice_i + ...
                      (6.5-j)/12*TbotAnomslice_im1;
           SbotAnom = (12-(6.5-j))/12*SbotAnomslice_i + ...
                      (6.5-j)/12*SbotAnomslice_im1;
        else
           TbotAnom = (12-(j-.5))/12*TbotAnomslice_i + ...
                      ((j-.5))/12*TbotAnomslice_ip1;
           SbotAnom = (12-(j-.5))/12*SbotAnomslice_i + ...
                      ((j-.5))/12*SbotAnomslice_ip1;
        end

        S_south(:,:,(i-1)*12+j) = ...
            SbotClim(:,:,j) + SbotAnom;
        T_south(:,:,(i-1)*12+j) = ...
            TbotClim(:,:,j) + TbotAnom;

    end
end

%S_south(120:end,:,:)=nan;
%T_south(120:end,:,:)=nan;
%V_south(120:end,:,:)=nan;
%U_south(120:end,:,:)=nan;

%should use nearest neighbour not S_col and T_col... will cause spurious
%t/s and convection

%Do this instead;
S_Int = zeros(nx+gx,nz,n_months);
T_Int = zeros(nx+gx,nz,n_months);
U_Int = zeros(nx+gx,nz,n_months);
V_Int = zeros(nx+gx,nz,n_months);
[X,Y] = meshgrid(linspace(1,nx,nx),z); Y=Y';X=X';
zmid = .5 * (z2(1:end-1)+z2(2:end));
[X2,Y2] = meshgrid(linspace(1,nx,nx),zmid); Y2=Y2';X2=X2';

for t=1:n_months
    
    S_time = zeros(nx,length(z)); S_time((size(S_south,1)+1):end,:)=nan;
    T_time = zeros(nx,length(z)); T_time((size(S_south,1)+1):end,:)=nan;
    U_time = zeros(nx,length(z)); U_time((size(S_south,1)+1):end,:)=nan;
    V_time = zeros(nx,length(z)); V_time((size(S_south,1)+1):end,:)=nan;

    S_time(1:size(S_south,1),:) = smooth_obc(S_south(:,:,t),smoothnum);
    T_time(1:size(S_south,1),:) = smooth_obc(T_south(:,:,t),smoothnum);
    U_time(1:size(S_south,1),:) = smooth_obc(U_south(:,:,t),smoothnum);
    V_time(1:size(S_south,1),:) = smooth_obc(V_south(:,:,t),smoothnum);
    ibad = isnan(S_time);
    
    S_interp = griddata(X(~ibad), Y(~ibad), S_time(~ibad), X2, Y2, 'nearest'); S_interp(isnan(S_interp))=-9999;
    S_Int(1:nx,:,t) = S_interp;
    T_interp = griddata(X(~ibad), Y(~ibad), T_time(~ibad), X2, Y2, 'nearest'); T_interp(isnan(T_interp))=-9999;
    T_Int(1:nx,:,t) = T_interp;
    U_interp = griddata(X(~ibad), Y(~ibad), U_time(~ibad), X2, Y2, 'nearest'); U_interp(isnan(U_interp))=-9999;
    U_Int(1:nx,:,t) = U_interp;
    
    V_interp = griddata(X(~ibad), Y(~ibad), V_time(~ibad), X2, Y2, 'nearest'); V_interp(isnan(V_interp))=-9999;
    V_Int(1:nx,:,t) = V_interp;
    t
    
    if(mod(t,12)==0);
      if (PAS>10)
        cyear = ystart_in + floor((t-1)/12)
        write_coupled_input(['vvel.obs.' num2str(PAS) '_' num2str(cyear)],V_Int(:,:,(t-11:t)),'oce');
        write_coupled_input(['uvel.obs.' num2str(PAS) '_' num2str(cyear)],U_Int(:,:,(t-11:t)),'oce');
        write_coupled_input(['salt.obs.' num2str(PAS) '_' num2str(cyear)],S_Int(:,:,(t-11:t)),'oce');
        write_coupled_input(['temp.obs.' num2str(PAS) '_' num2str(cyear)],T_Int(:,:,(t-11:t)),'oce');
      elseif (PAS==10)
        write_coupled_input(['vvel.obs.' num2str(PAS)],V_Int(:,:,(t-11:t)),'oce');
        write_coupled_input(['uvel.obs.' num2str(PAS)],U_Int(:,:,(t-11:t)),'oce');
        write_coupled_input(['salt.obs.' num2str(PAS)],S_Int(:,:,(t-11:t)),'oce');
        write_coupled_input(['temp.obs.' num2str(PAS)],T_Int(:,:,(t-11:t)),'oce');
      end
    end
    
      

    



    
end

% figure()
% pcolor(T_Int(:,:,12)'); hold on;set(gca,'YDir','reverse'); shading flat; colorbar;  %caxis([33.9,34.6]);
% plot(linspace(x_mesh_mid(1:nx)),-topo(:,1)/20, 'Color','k')

if (PAS==0 | PAS>100)
    write_coupled_input(['vvel_' num2str(ystart) '.obs'  '.' num2str(PAS)],V_Int,'oce');
    write_coupled_input(['uvel_' num2str(ystart) '.obs'  '.' num2str(PAS)],U_Int,'oce');
    write_coupled_input(['salt_' num2str(ystart) '.obs'  '.' num2str(PAS)],S_Int,'oce');
    write_coupled_input(['temp_' num2str(ystart) '.obs'  '.' num2str(PAS)],T_Int,'oce');
    replace_param(obcsdatafile, 'obstfile', ['temp_' num2str(ystart) '.obs'  '.' num2str(PAS)]);
    replace_param(obcsdatafile, 'obssfile', ['salt_' num2str(ystart) '.obs'  '.' num2str(PAS)]);
    replace_param(obcsdatafile, 'obsufile', ['uvel_' num2str(ystart) '.obs'  '.' num2str(PAS)]);
    replace_param(obcsdatafile, 'obsvfile', ['vvel_' num2str(ystart) '.obs'  '.' num2str(PAS)]);
else
    replace_param(obcsdatafile, 'obstfile', ['temp.obs.' num2str(PAS)]);
    replace_param(obcsdatafile, 'obssfile', ['salt.obs.' num2str(PAS)]);
    replace_param(obcsdatafile, 'obsufile', ['uvel.obs.' num2str(PAS)]);
    replace_param(obcsdatafile, 'obsvfile', ['vvel.obs.' num2str(PAS)]);
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


U_west = Uleftbdry; 
V_west = Vleftbdry; 

V_westPAS = permute(V_west,[2 1 3]);
U_westPAS = permute(U_west,[2 1 3]);


% need U and V
% append 10-year cycle to end of series
PASyears = size(U_west,3)/12;
U_west = zeros(ny_forcing,length(z),n_months);
V_west = zeros(ny_forcing,length(z),n_months);
U_west(:,:,1:(12*PASyears-12)) = U_westPAS(:,:,1:(end-12));
V_west(:,:,1:(12*PASyears-12)) = V_westPAS(:,:,1:(end-12));
U_west(:,:,(12*PASyears-11):end) = repmat(U_westPAS(:,:,(end-119):end),[1 1 (1+length(year_array)-dstart-PASyears)/10]);
V_west(:,:,(12*PASyears-11):end) = repmat(V_westPAS(:,:,(end-119):end),[1 1 (1+length(year_array)-dstart-PASyears)/10]);

% now, to create S_south and T_south!
S_west = zeros(ny_forcing,length(z),n_months);
T_west = zeros(ny_forcing,length(z),n_months);

SleftClim = permute(SleftClim,[2 1 3]);
TleftClim = permute(TleftClim,[2 1 3]);

for i = 1:(length(year_array) - (ystart_in-year_array(1)));

    TleftAnomslice_i = TleftAnomslice(:,:,i+dstart);
    TleftAnomslice_im1 = TleftAnomslice(:,:,i+dstart-1);
    SleftAnomslice_i = SleftAnomslice(:,:,i+dstart);
    SleftAnomslice_im1 = SleftAnomslice(:,:,i+dstart-1);

    if (i==(length(year_array) - (ystart_in-year_array(1))));
        TleftAnomslice_ip1 = TleftAnomslice(:,:,i+dstart);
        SleftAnomslice_ip1 = SleftAnomslice(:,:,i+dstart);
    else
        TleftAnomslice_ip1 = TleftAnomslice(:,:,i+dstart+1);
        SleftAnomslice_ip1 = SleftAnomslice(:,:,i+dstart+1);
    end

    for j=1:12;

        if (j<=6);
           TleftAnom = (12-(6.5-j))/12*TleftAnomslice_i + ...
                      (6.5-j)/12*TleftAnomslice_im1;
           SleftAnom = (12-(6.5-j))/12*SleftAnomslice_i + ...
                      (6.5-j)/12*SleftAnomslice_im1;
        else
           TleftAnom = (12-(j-.5))/12*TleftAnomslice_i + ...
                      ((j-.5))/12*TleftAnomslice_ip1;
           SleftAnom = (12-(j-.5))/12*SleftAnomslice_i + ...
                      ((j-.5))/12*SleftAnomslice_ip1;
        end

        S_west(:,:,(i-1)*12+j) = ...
            SleftClim(:,:,j) + SleftAnom;
        T_west(:,:,(i-1)*12+j) = ...
            TleftClim(:,:,j) + TleftAnom;

    end
end

%%%%%%


[X,Y] = meshgrid(linspace(1,ny,ny),z); Y=Y';X=X';
[X2,Y2] = meshgrid(linspace(1,ny,ny),zmid); Y2=Y2';X2=X2';
S_Int = zeros(ny+gy,nz,n_months);
T_Int = zeros(ny+gy,nz,n_months);
U_Int = zeros(ny+gy,nz,n_months);
V_Int = zeros(ny+gy,nz,n_months);
for t=1:n_months
    S_time = zeros(ny,length(z)); S_time((size(S_west,1)+1):end,:)=nan;
    T_time = zeros(ny,length(z)); T_time((size(S_west,1)+1):end,:)=nan;
    U_time = zeros(ny,length(z)); U_time((size(S_west,1)+1):end,:)=nan;
    V_time = zeros(ny,length(z)); V_time((size(S_west,1)+1):end,:)=nan;
    
    S_time(1:size(S_west,1),:) = smooth_obc(S_west(:,:,t),smoothnum);
    T_time(1:size(S_west,1),:) = smooth_obc(T_west(:,:,t),smoothnum);
    U_time(1:size(S_west,1),:) = smooth_obc(U_west(:,:,t),smoothnum);
    V_time(1:size(S_west,1),:) = smooth_obc(V_west(:,:,t),smoothnum);
    ibad = isnan(S_time);
    S_interp = griddata(X(~ibad), Y(~ibad), S_time(~ibad), X2, Y2, 'nearest'); S_interp(isnan(S_interp))=-9999;
    S_Int(1:ny,:,t) = S_interp;
    T_interp = griddata(X(~ibad), Y(~ibad), T_time(~ibad), X2, Y2, 'nearest'); T_interp(isnan(T_interp))=-9999;
    T_Int(1:ny,:,t) = T_interp;
    U_interp = griddata(X(~ibad), Y(~ibad), U_time(~ibad), X2, Y2, 'nearest'); U_interp(isnan(U_interp))=-9999;
    U_Int(1:ny,:,t) = U_interp; 
    V_interp = griddata(X(~ibad), Y(~ibad), V_time(~ibad), X2, Y2, 'nearest'); V_interp(isnan(V_interp))=-9999;
    V_Int(1:ny,:,t) = V_interp;
    t

    if(mod(t,12)==0);
      if (PAS>10)
        cyear = ystart_in + floor((t-1)/12)
        write_coupled_input(['vvel.obw.' num2str(PAS) '_' num2str(cyear)],V_Int(:,:,(t-11:t)),'oce');
        write_coupled_input(['uvel.obw.' num2str(PAS) '_' num2str(cyear)],U_Int(:,:,(t-11:t)),'oce');
        write_coupled_input(['salt.obw.' num2str(PAS) '_' num2str(cyear)],S_Int(:,:,(t-11:t)),'oce');
        write_coupled_input(['temp.obw.' num2str(PAS) '_' num2str(cyear)],T_Int(:,:,(t-11:t)),'oce');
      elseif (PAS==10)
        write_coupled_input(['vvel.obw.' num2str(PAS)],V_Int(:,:,(t-11:t)),'oce');
        write_coupled_input(['uvel.obw.' num2str(PAS)],U_Int(:,:,(t-11:t)),'oce');
        write_coupled_input(['salt.obw.' num2str(PAS)],S_Int(:,:,(t-11:t)),'oce');
        write_coupled_input(['temp.obw.' num2str(PAS)],T_Int(:,:,(t-11:t)),'oce');
      end
    end
    

end

% figure()
% pcolor(T_Int(:,:,12)'); hold on;set(gca,'YDir','reverse'); shading flat; colorbar
% plot(linspace(1,ny,ny),-topo(1,:)/20, 'Color','k')

if (PAS==0 | PAS>100)
    write_coupled_input(['vvel_' num2str(ystart) '.obw' '.' num2str(PAS)],V_Int,'oce');
    write_coupled_input(['uvel_' num2str(ystart) '.obw' '.' num2str(PAS)],U_Int,'oce');
    write_coupled_input(['salt_' num2str(ystart) '.obw' '.' num2str(PAS)],S_Int,'oce');
    write_coupled_input(['temp_' num2str(ystart) '.obw' '.' num2str(PAS)],T_Int,'oce');
    replace_param(obcsdatafile, 'obwtfile', ['temp_' num2str(ystart) '.obw'  '.' num2str(PAS)]);
    replace_param(obcsdatafile, 'obwsfile', ['salt_' num2str(ystart) '.obw'  '.' num2str(PAS)]);
    replace_param(obcsdatafile, 'obwufile', ['uvel_' num2str(ystart) '.obw'  '.' num2str(PAS)]);
    replace_param(obcsdatafile, 'obwvfile', ['vvel_' num2str(ystart) '.obw'  '.' num2str(PAS)]);
else
    replace_param(obcsdatafile, 'obwtfile', ['temp.obw.' num2str(PAS)]);
    replace_param(obcsdatafile, 'obwsfile', ['salt.obw.' num2str(PAS)]);
    replace_param(obcsdatafile, 'obwufile', ['uvel.obw.' num2str(PAS)]);
    replace_param(obcsdatafile, 'obwvfile', ['vvel.obw.' num2str(PAS)]);
end


if(ystart_in==2001);
datestr1 = '20001216';
datestr2 = '20010101';
else
datestr1 = '20051216';
datestr2 = '20060101';
end

    replace_param(exfdatafile, 'obcsWstartdate1',datestr1,false);
    replace_param(exfdatafile, 'obcssstartdate1',datestr1,false);
    replace_param(caldatafile, 'startdate_1',datestr2,false);
    
    
