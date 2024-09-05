
% This script prepares experimental electron ptycho. data for PtychoShelves
clear
%% Step 1: download the sample data from the PARADIM website:
% https://data.paradim.org/doi/ssmm-2j11/
% Note: it's a good practice to store data (and reconstructions) in a
% different folder from fold_slice 
scriptfolder = 'C:\Users\PanGroupWorkstation\Desktop\hanhsuanwu\fold_slice\ptycho';
addpath(strcat(pwd,'/utils_electron/'))

% Step 2: load data
data_dir = 'D:\Wuhanhsuan\20240806BTO Trial1 bulk\'; %change this
data_dir = strrep(data_dir,'\','/')
filename = 'Spectrum Image EELS Image.npy';
h5_suf = 'mask_bin2_crop50x50';
scan_number = 1; %Ptychoshelves needs
bin = 2;
cutoff = 180;
crop_idx=[1,50,1,50]; % start from smaller data [lower y, higher y, lower x, higer x]
% load(strcat(data_dir,'s01_R3_0_0.mat'))

% input parameters if they are not included in the mat file
exp_p = {};
exp_p.ADU =     1.0;
exp_p.voltage = 100; % kV
exp_p.alpha =   25.0; % mrad
[~,lambda] =    electronwavelength(exp_p.voltage);
%dk =            0.0367;
exp_p.defocus = -50.0; % Angst.
exp_p.scan_step_size = 0.4; % Angst.
exp_p.nv = 180; %final dp pattern size
% exp_p.rot_ang = 20.0;

% calculate pxiel size (1/A) in diffraction plane
% [~,lambda]=electronwavelength(exp_p.voltage);

exp_p.rbf=149.3/2; % radius of center disk in pixels
dk=exp_p.alpha/1e3/exp_p.rbf/lambda;
%exp_p.rbf = exp_p.alpha/1e3/dk/lambda; % radius of center disk in pixels
exp_p.scan_number = scan_number;
roi_label = strcat('0_Ndp', num2str(exp_p.nv), h5_suf);
exp_p.roi_label = roi_label;
folder = data_dir;
save_dir = strcat(folder,num2str(scan_number),'/'); % where to save
mkdir(save_dir)
% save experimental parameters
save(strcat(save_dir,'/exp_para', h5_suf, '.mat'),'exp_p');
% copyfile(strcat(scriptfolder, mfilename, '.m'), strcat(save_dir, mfilename, '.m'));

% Step 3: load raw data
% nv, nv, ny, nx
f = strcat(data_dir,filename);
if strcmp(filename(end-2:end), 'npy')
    dp = readNPY(f);
    dp = permute(dp, [3 4 1 2]);
elseif strcmp(filename(end-2:end), 'mat')
    load(f)
    dp = m;
elseif strcmp(filename(end-1:end), 'h5')
    dp = io_TEAM(f, 0, 0, 0, 0);
end
% dp = dp(2:end,2:end,:,:);
pacbed = mean(dp, [3 4]);
%figure(); imagesc(pacbed); colorbar; axis image;

dp = applyCircularCutoff(dp, cutoff);
pacbed = mean(dp, [3 4]);
%figure(); imagesc(pacbed); colorbar; axis image;
% bin / pad
%%% bin 4d
if bin > 1
    dp = bin4d(dp, bin, bin);
    warning('CBEDs are binned.')
end
Np_p = [exp_p.nv,exp_p.nv]; % size of diffraction patterns used during reconstruction. can also pad to 256
%%% crop data
dp=dp(:,:,crop_idx(1):crop_idx(2),crop_idx(3):crop_idx(4));

%%% pad cbed
[ndpy,ndpx,npy,npx]=size(dp);
if ndpy < Np_p(1) % pad zeros
    dp=padarray(dp,[(Np_p(1)-ndpy)/2,(Np_p(2)-ndpx)/2,0,0],0,'both');
else
    dp=crop_pad(dp,Np_p);
end
% calc_rotation(dp)

dp = dp / exp_p.ADU; % convert to electron count, contained in the data file
dp=reshape(dp,Np_p(1),Np_p(2),[]);
pacbed = mean(dp, 3);
figure(); 
imagesc(pacbed.^0.5); colorbar; axis image;
Itot=mean(squeeze(sum(sum(dp,1),2))); %need this for normalizting initial probe

% Get the size of the image
[rows, cols] = size(pacbed);

% Define the center and radius of the circle
centerX = cols / 2;
centerY = rows / 2;
radius = exp_p.rbf/bin;  % Adjust the radius as needed

% Draw the circle
rectangle('Position', [centerX - radius, centerY - radius, 2*radius, 2*radius], ...
          'Curvature', [1, 1], 'EdgeColor', 'r', 'LineWidth', 2);
hold off;

% Step 4: save CBED in a .hdf5 file (needed by Ptychoshelves)
% scan_number = 1; %Ptychoshelves needs
% save_dir = strcat(data_dir,num2str(scan_number),'/');
saveName = strcat('data_roi',roi_label,'_dp.hdf5');
h5create(strcat(save_dir,saveName), '/dp', size(dp),'ChunkSize',[size(dp,1), size(dp,2), 1],'Deflate',4)  % Deflate 为压缩级别
h5write(strcat(save_dir,saveName), '/dp', dp)

%% Step 5: prepare initial probe
dx=1/Np_p(1)/dk; %% pixel size in real space (angstrom)

par_probe = {};
par_probe.df = exp_p.defocus;
par_probe.voltage = exp_p.voltage;
par_probe.alpha_max = exp_p.alpha;
par_probe.plotting = true;
probe = make_tem_probe(dx, Np_p(1), par_probe);

probe=probe/sqrt(sum(sum(abs(probe.^2))))*sqrt(Itot)/sqrt(Np_p(1)*Np_p(2));
probe=single(probe);
% add parameters for PtychoShelves
p = {};
p.binning = false;
p.detector.binning = false;

%% Step 6: save initial probe
save(strcat(save_dir,'/init_probe.mat'),'probe','p')