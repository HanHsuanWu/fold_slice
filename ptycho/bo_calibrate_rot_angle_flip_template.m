directory_content = dir; % contains everything of the current directory
exe_path = directory_content(1).folder; % returns the path that is currently open

scriptfolder = '/mnt/pgo4/pgo4_v1/Han-Hsuan/fold_slice/ptycho';
scriptfolder = strrep(scriptfolder,'\','/');
cd(scriptfolder);
addpath(strcat(pwd,'/utils/'))
addpath(core.find_base_package)

%% Step 0: Run the prepare_data script to generate data for ptycho reconstruction
% Step 1: Prepare data and reconstruction parameters
par = {};
par.verbose_level = 3;
par.scan_number = 1;
par.beam_source = 'electron';

base_path = '/mnt/pgo4/pgo4_v1\Han-Hsuan\Ptychography_test\20240917_AlGaAs_nion\data2\';
par.base_path = strrep(base_path,'\','/');
clear base_path;
par.scan_format = '%01d';
par.scan_string_format = '%01d';
par.roi_label = '0_Ndp360mask'; %don't forget to change this to match file name
par.Ndp = 360;
bin = 1;

num_ang = 5;
ang_range = [115.2, 235.2];
scan_custom_fliplr =    ones(1,num_ang);
scan_custom_flipud =    ones(1,num_ang);
scan_custom_transpose = ones(1,num_ang);
rot_ang =               ang_range(1):(ang_range(2)-ang_range(1))/(num_ang-1):ang_range(2);
%%
%{
%[1,1,1] means no flip
scan_custom_fliplr =    [1,1,1,1,1,1,1,1,1,1,1];
scan_custom_flipud =    [1,1,1,1,1,1,1,1,1,1,1];
scan_custom_transpose = [1,1,1,1,1,1,1,1,1,1,1];
rot_ang =               [57.9,67.9,77.9,82.9,87.9,92.9,97.9,107.9,-92.1,-102.1,-112.1];
%}

%{
scan_custom_fliplr =    [0,1,0,0,1,1,1,0];
scan_custom_flipud =    [0,0,1,0,1,1,0,1];
scan_custom_transpose = [0,0,0,1,1,0,1,1];
rot_ang =               [0,0,0,0,0,0,0,0];  1.8561  1.8848
%}

Niter=100;

par.defocus = -100; %overfocus is negative
par.energy = 60;
par.rbf = 185.8/2/bin;

par.cen_dp_y = floor(par.Ndp/2)+1;
par.cen_dp_x = floor(par.Ndp/2)+1;

par.scan_nx = 200;
par.scan_ny = 200;

par.scan_step_size_x = 0.191;
par.scan_step_size_y = 0.191;

par.detector_name = 'empad';
par.data_preparator = 'matlab_aps';
par.src_positions =  'matlab_pos';
par.scan_type = 'raster';

par.use_model_probe = true;
par.normalize_init_probe = true;

par.output_dir_base = par.base_path;
par.Niter = Niter;
par.Niter_save_results_every = Niter;
par.save.save_reconstructions = true;

par.eng_name = 'GPU_MS';
par.method = 'MLc';
par.momentum = 0.5;

par.Nprobe = 2;
par.grouping = 50;
par.apply_multimodal_update = true;

par.Nlayers = 1;
par.regularize_layers = 1;
par.variable_probe_modes = 1;
par.Ndp_presolve = 360;

% Step 1.5 (optional): Run a single reconstruction to check parameters
par.GPU_list = [1,2,3,4];

par.alpha_max = 38.0;

par.thickness = 160; %in angstroms
thickness = 160;

data_error_list = zeros(1, length(scan_custom_fliplr));

%% make it iterative
for ieng=1:length(scan_custom_fliplr)
    %These uses custom_flip the newer version
    par.scan_custom_fliplr = scan_custom_fliplr(ieng);
    par.scan_custom_flipud = scan_custom_flipud(ieng);
    par.scan_custom_transpose = scan_custom_transpose(ieng);
    par.rot_ang = rot_ang(ieng);

    output_dir_suffix_base = strcat('_flip', num2str(par.scan_custom_fliplr), num2str(par.scan_custom_flipud), num2str(par.scan_custom_transpose), '_RotAng', num2str(par.rot_ang));
    par.output_dir_suffix_base = strrep(output_dir_suffix_base,'\','/');

    data_error_list(ieng) = ptycho_recon_exp_data(par, 'thickness', thickness);
end
%{
%% plot result
scatter(rot_ang,data_error_list);
% Add labels and title
xlabel('Rotation Angle');
ylabel('Fourier Error');
title('');
saveas(gcf,strcat(par.base_path,num2str(par.scan_number),'/rot_angError.tiff'));
%}
%{
%% Step 2: Use Bayesian optimization with Gaussian processes to find experimental parameters that minimize data error
% Note: Parallel BO is generally recommended for multislice reconstructions
close all

par.GPU_list = [1];

par.scan_custom_fliplr = 1;
par.scan_custom_flipud = 1;
par.scan_custom_transpose = 1;

rot_ang = optimizableVariable('rot_ang', [0, 360]); %angstroms

N_workers = length(par.GPU_list);
if N_workers>1
    delete(gcp('nocreate'))
    c = parcluster('local');
    c.NumWorkers = N_workers;
    p = parpool(c);
end

fun = @(x)ptycho_recon_exp_data(par, 'rot_ang', x.rot_ang);
results = bayesopt(fun, [rot_ang],...
    'Verbose', 4,...
    'AcquisitionFunctionName', 'expected-improvement-plus',...
    'IsObjectiveDeterministic', false,...
    'MaxObjectiveEvaluations', 30,...
    'NumSeedPoints', N_workers,...
    'PlotFcn', {@plotObjectiveModel, @plotMinObjective}, 'UseParallel', N_workers>1);

delete(gcp('nocreate'))
%}