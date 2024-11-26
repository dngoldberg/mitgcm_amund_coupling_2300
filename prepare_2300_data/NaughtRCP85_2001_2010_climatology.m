%load ../input/start_2009_input/bdryDatapahol30.mat

% we make a 2001-2010 climatology at boundaries

SbotClim = zeros(50,240,12);
SleftClim = zeros(50,416,12);

TbotClim = zeros(50,240,12);
TleftClim = zeros(50,416,12);

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
    end
    i
end

save naughtenClim.mat z SbotClim   SleftClim  TbotClim   TleftClim

znaught = z;

addpath('/home/dgoldber/scratch/globus/');

slices_for_forcing

disp('finished')


