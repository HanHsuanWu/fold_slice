% This script prepares experimental electron ptycho. data for PtychoShelves
clear;
voltage=100; %kev
ADU=1; %intensity count to electron conversion
rbf=149.3/2; % radius of center disk in pixels
alpha=24.3 ; %convergence angle in mrad
scan_number = 2; %Ptychoshelves needs
crop_nyx=[1,40,1,40]; % start from smaller data [lower y, higher y, lower x, higer x]

data_dir = 'D:\Wuhanhsuan\20240806BTO Trial2 bulk\'; %change this
data_dir = strrep(data_dir, '\', '/'); % needed this to avoid having \U
cbed_strut=load(strcat(data_dir,'Spectrum Image EELS Image.mat'));
dp= cat(4,cbed_strut.data1,cbed_strut.data2);
clear cbed_strut;

% Step 3: go back to .../fold_slice/ptycho and pre-process data
addpath(strcat(pwd,'/utils_electron/'))
Np_p = [360, 360]; % size of diffraction patterns used during reconstruction. can also pad to 256
% pad cbed
[ndpy,ndpx,npy,npx]=size(dp);
if ndpy < Np_p(1) % pad zeros
    dp=padarray(dp,[(Np_p(1)-ndpy)/2,(Np_p(2)-ndpx)/2,0,0],0,'both');
else
    dp=crop_pad(dp,Np_p);
end

%crop in ny,nx

dp=dp(:,:,crop_nyx(1):crop_nyx(2),crop_nyx(3):crop_nyx(4));

dp = dp / ADU; % convert to electron count
dp=reshape(dp,Np_p(1),Np_p(2),[]);
Itot=mean(squeeze(sum(sum(dp,1),2))); %need this for normalizting initial probe

% calculate pxiel size (1/A) in diffraction plane
[~,lambda]=electronwavelength(voltage);


dk=alpha/1e3/rbf/lambda; %%% PtychoShelves script needs this %%%

% Step 4: save CBED in a .hdf5 file (needed by Ptychoshelves)
save_dir = strcat(data_dir,num2str(scan_number),'/');
mkdir(save_dir)
roi_label = '0_Ndp360';
saveName = strcat('data_roi',roi_label,'_dp.hdf5');
h5name=strcat(save_dir,saveName);
%%
h5create(h5name, '/dp', size(dp),'ChunkSize',[size(dp,1), size(dp,2), 1],'Deflate',4);
%%
h5write(h5name, '/dp', dp);
