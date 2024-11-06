
% This script prepares experimental electron ptycho. data for PtychoShelves
clear
%% Step 1: download the sample data from the PARADIM website:
% https://data.paradim.org/doi/ssmm-2j11/
% Note: it's a good practice to store data (and reconstructions) in a
% different folder from fold_slice 
scriptfolder = 'C:\Users\hanhsuan\Documents\GitHub\fold_slice\ptycho';
addpath(strcat(scriptfolder,'\utils_electron_from_huozhi'));

% Step 2: load data
data_dir = '\\PanGroupOffice4\PGO4_v1\Han-Hsuan\Ptychography_test\20241010_AlGaAs-30s_arm\trial6bin8\'; %change this
data_dir = strrep(data_dir,'\','/');
filename = 'ALGAAS30-6_bin8.npy';
h5_suf = '';
scan_number = 1; %Ptychoshelves needs
bin = 1;
%cutoff = 600;
%crop_idx = [1,100,1,100]; % start from smaller data [lower y, higher y, lower x, higer x]
% Positive is up and left.
shift_dp = [0,0]; % [shift ky, shift kx] shift the center of dp by croping kx ky pixels then pad with 0. Has to be even number
dp_size = 100; % Initial size of diffraction pattern


% load(strcat(data_dir,'s01_R3_0_0.mat'))

% input parameters if they are not included in the mat file
exp_p = {};
exp_p.ADU =     1.0;
exp_p.voltage = 300; % kV
exp_p.alpha =   25.0; % mrad
[~,lambda] =    electronwavelength(exp_p.voltage);
%dk =            0.0367;
exp_p.defocus = 0.0; % Angst.
exp_p.scan_step_size = 0.42; % Angst.
exp_p.nv = dp_size; %final dp pattern size
% exp_p.rot_ang = 20.0;

% calculate pxiel size (1/A) in diffraction plane
% [~,lambda]=electronwavelength(exp_p.voltage);

exp_p.rbf=50.0/2/bin; % radius of center disk in pixels
dk=exp_p.alpha/1e3/exp_p.rbf/lambda;
%exp_p.rbf = exp_p.alpha/1e3/dk/lambda; % radius of center disk in pixels
exp_p.scan_number = scan_number;
roi_label = strcat('0_Ndp', num2str(exp_p.nv/bin), h5_suf);
exp_p.roi_label = roi_label;
folder = data_dir;
save_dir = strcat(folder,num2str(scan_number),'/'); % where to save
mkdir(save_dir)
% save experimental parameters
save(strcat(save_dir,'/exp_para', h5_suf, '.mat'),'exp_p');
% copyfile(strcat(scriptfolder, mfilename, '.m'), strcat(save_dir, mfilename, '.m'));

%% Step 3: load raw data
% nv, nv, ny, nx
f = strcat(data_dir,filename);
if strcmp(filename(end-2:end), 'npy') 
    dp = readNPY(f); %GrandArm do not need any permute
    %dp = permute(dp, [4 3 1 2]); %EMPAD
    disp(size(dp));
elseif strcmp(filename(end-2:end), 'mat')
    dp_struct = load(f);
    fields = fieldnames(dp_struct);
    dp = dp_struct.(fields{1});

    for i = 2:length(fields)
        dp = cat(4,dp,dp_struct.(fields{i}));
    end

    dp = permute(dp, [3 4 1 2]); %after this permute [ky kx y x]
    clear dp_struct;
elseif strcmp(filename(end-1:end), 'h5')
    dp = io_TEAM(f, 0, 0, 0, 0);
end
    dp = double(dp);
% dp needs to be in [ky kx y x]
%% Check CBED center of 1 dp
dp1=dp(:,:,1,1);
%crop to center dp
if shift_dp(1) >= 0
    row_start = shift_dp(1) + 1;
    row_end = size(dp1, 1);
else
    row_start = 1;
    row_end = size(dp1, 1) + shift_dp(1);
end

if shift_dp(2) >= 0
    col_start = shift_dp(2) + 1;
    col_end = size(dp1, 2);
else
    col_start = 1;
    col_end = size(dp1, 2) + shift_dp(2);
end

Np_p = [exp_p.nv,exp_p.nv]; % size of diffraction patterns used during reconstruction. can also pad to 256
%%% crop data

dp1 = dp1(row_start:row_end, col_start:col_end, :, :);
%%% pad cbed
[ndpy,ndpx,~,~]=size(dp1);
if ndpy < Np_p(1) % pad zeros
    dp1=padarray(dp1,[(Np_p(1)-ndpy)/2,round((Np_p(2)-ndpx)/2),0,0],0,'both');
else
    dp1=crop_pad(dp1,Np_p);
end

%cutoff cbed with circular mask
%dp1 = applyCircularCutoff(dp1, cutoff);
pacbed2 = mean(dp1, [3 4]);
figure(); imagesc(pacbed2); colorbar; axis image;

% Get the size of the image
[rows, cols] = size(pacbed2);

% Define the center and radius of the circle
centerX = cols / 2;
centerY = rows / 2;
radius = exp_p.rbf*bin;  % Adjust the radius as needed

% Draw the circle and it's center
rectangle('Position', [centerX - radius, centerY - radius, 2*radius, 2*radius], ...
          'Curvature', [1, 1], 'EdgeColor', 'r', 'LineWidth', 2);

circle_radius = 2;
rectangle('Position', [centerX - circle_radius/2, centerY - circle_radius/2, circle_radius, circle_radius], ...
          'Curvature', [1, 1], 'EdgeColor', 'r', 'LineWidth', 2);

hold off;

%% Center all dp

%crop to center dp
if shift_dp(1) >= 0
    row_start = shift_dp(1) + 1;
    row_end = size(dp, 1);
else
    row_start = 1;
    row_end = size(dp, 1) + shift_dp(1);
end

if shift_dp(2) >= 0
    col_start = shift_dp(2) + 1;
    col_end = size(dp, 2);
else
    col_start = 1;
    col_end = size(dp, 2) + shift_dp(2);
end
dp = dp(row_start:row_end, col_start:col_end, :, :);

Np_p = [exp_p.nv,exp_p.nv]; % size of diffraction patterns used during reconstruction. can also pad to 256
%%% crop data
%dp=dp(:,:,crop_idx(1):crop_idx(2),crop_idx(3):crop_idx(4));

%%% pad cbed
[ndpy,ndpx,npy,npx]=size(dp);
if ndpy < Np_p(1) % pad zeros
    dp=padarray(dp,[(Np_p(1)-ndpy)/2,round((Np_p(2)-ndpx)/2),0,0],0,'both');
else
    dp=crop_pad(dp,Np_p);
end

%cutoff cbed with circular mask
%dp = applyCircularCutoff(dp, cutoff);

% bin / pad
%%% bin 4d
if bin > 1
    dp = bin4d(dp, bin, bin);
    warning('CBEDs are binned.')
end
%% Check Average CBED
pacbed = mean(dp, [3 4]);
figure(); imagesc(pacbed.^0.5); colorbar; axis image;
%% Check virtual BF image
cutoffbf = exp_p.rbf ;
bfdata = applyCircularCutoff(dp, cutoffbf);
pacbed = mean(bfdata, [3 4]);
figure(); imagesc(pacbed.^0.5); colorbar; axis image;
title('crop BF cbed', 'FontSize', 14);
saveas(gcf,strcat(data_dir,'/crop_BF_cbed.tiff'));
close();
%%
bf_image=squeeze(sum(sum(bfdata,1),2)).^0.5;
bf = imagesc(transpose(bf_image)); colorbar ; axis image;
colormap(flipud(parula));
title('reverse bf image', 'FontSize', 14);
saveas(gcf,strcat(data_dir,'/reverse_bf_image.tiff'));
close();
%% check rotation and flip 
calc_rotation(dp)
saveas(gcf,strcat(data_dir,'/curlCOM.tiff'));
close()
%%
dp = dp / exp_p.ADU; % convert to electron count, contained in the data file
dp=reshape(dp,Np_p(1)/bin,Np_p(2)/bin,[]);
pacbed = mean(dp, 3);
figure('Position', [600 200 400 400]);
%imagesc(pacbed); colorbar; axis image;
imagesc(pacbed.^0.5); colorbar; axis image;
Itot=mean(squeeze(sum(sum(dp,1),2))); %need this for normalizting initial probe

% Get the size of the image
[rows, cols] = size(pacbed);

% Define the center and radius of the circle
centerX = cols / 2;
centerY = rows / 2;
radius = exp_p.rbf;  % Adjust the radius as needed

% Draw the circle and it's center
rectangle('Position', [centerX - radius, centerY - radius, 2*radius, 2*radius], ...
          'Curvature', [1, 1], 'EdgeColor', 'r', 'LineWidth', 0.5);

circle_radius = 0.5;
rectangle('Position', [centerX - circle_radius/2, centerY - circle_radius/2, circle_radius, circle_radius], ...
          'Curvature', [1, 1], 'EdgeColor', 'r', 'LineWidth', 0.5);

hold off;

saveas(gcf,strcat(data_dir,'/cbed.tiff'));
close();
%% save file
copyfile([mfilename('fullpath'), '.m'], strcat(data_dir,num2str(scan_number),'/prepare_data.m'));
%% save dp hdf5
%caxis([0 7.5]);


% Step 4: save CBED in a .hdf5 file (needed by Ptychoshelves)
% scan_number = 1; %Ptychoshelves needs
% save_dir = strcat(data_dir,num2str(scan_number),'/');
saveName = strcat('data_roi',roi_label,'_dp.hdf5');
h5create(strcat(save_dir,saveName), '/dp', size(dp),'ChunkSize',[size(dp,1), size(dp,2), 1],'Deflate',4)  % Deflate 为压缩级别
h5write(strcat(save_dir,saveName), '/dp', dp)

%{

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
%}
