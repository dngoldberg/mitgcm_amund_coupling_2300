%load ../input/start_2009_input/bdryDatapahol30.mat

% we make a 2001-2010 climatology at boundaries

SbotClim = zeros(50,240,12);
SleftClim = zeros(50,416,12);

TbotClim = zeros(50,240,12);
TleftClim = zeros(50,416,12);

S2005 = zeros(416,240,50);
T2005 = zeros(416,240,50);


UbotClim = zeros(50,240,12);
UleftClim = zeros(50,416,12);

VbotClim = zeros(50,240,12);
VleftClim = zeros(50,416,12);



for i=1:10;
    load (['bdryDatapaholPyRCP850' appNum(i,2) 'Dan01_10.mat']);
    for j=1:12;
        SbotClim(:,:,j) = SbotClim(:,:,j) + ...
            1/10 * mean(Sbotbdry(:,1:240,j:12:end),3);
        TbotClim(:,:,j) = TbotClim(:,:,j) + ...
            1/10 * mean(Tbotbdry(:,1:240,j:12:end),3);
        SleftClim(:,:,j) = SleftClim(:,:,j) + ...
            1/10 * mean(Sleftbdry(:,:,j:12:end),3);
        TleftClim(:,:,j) = TleftClim(:,:,j) + ...
            1/10 * mean(Tleftbdry(:,:,j:12:end),3);


        UbotClim(:,:,j) = UbotClim(:,:,j) + ...
            1/10 * mean(Ubotbdry(:,1:240,j:12:end),3);
        VbotClim(:,:,j) = VbotClim(:,:,j) + ...
            1/10 * mean(Vbotbdry(:,1:240,j:12:end),3);
        UleftClim(:,:,j) = UleftClim(:,:,j) + ...
            1/10 * mean(Uleftbdry(:,:,j:12:end),3);
        VleftClim(:,:,j) = VleftClim(:,:,j) + ...
            1/10 * mean(Vleftbdry(:,:,j:12:end),3);
    end
    load (['initpaholPyRCP850' appNum(i,2) 'Dan01_10.mat']);
    S2005 = S2005 + 1/10 * Sinit(:,1:240,:);
    T2005 = T2005 + 1/10 * Tinit(:,1:240,:);
    i
end

Sinit = S2005;
Tinit = T2005;
vars = {'U','V','T','S'};
for i=1:4;
    k=vars{i};
    eval([k 'botbdry = ' k 'botClim;']);
    eval([k 'leftbdry = ' k 'leftClim;']);
end

ystart=2001;


save naughtenClim.mat z SbotClim   SleftClim  TbotClim   TleftClim S2005 T2005
save naughtenBdryClim2010.mat Sbotbdry   Sleftbdry  Tbotbdry   Tleftbdry  Ubotbdry   Uleftbdry  Vbotbdry   Vleftbdry  ystart
save naughtenInit2010.mat Sinit Tinit z


znaught = z;

addpath('/home/dgoldber/scratch/globus/');

%slices_for_forcing

disp('finished')


