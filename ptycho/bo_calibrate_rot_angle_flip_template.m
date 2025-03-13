directory_content = dir; % contains everything of the current directory
exe_path = directory_content(1).folder; % returns the path that is currently open

scriptfolder = '/mnt/pgo4/pgo4_v1/Han-Hsuan/fold_slice/ptycho';
scriptfolder = strrep(scriptfolder,'\','/');
cd(scriptfolder);
%%
addpath(strcat(pwd,'/utils/'))
addpath(core.find_base_package)

% Step 0: Run the prepare_data script to generate data for ptycho reconstruction
% Step 1: Prepare data and reconstruction parameters
par = {};
par.verbose_level = 3;
par.scan_number = 2;
par.beam_source = 'electron';

base_path = '\\PanGroupOffice4\PGO4_v1\Han-Hsuan\Ptychography_test\20241010_AlGaAs-30s_arm\trial4\';
%base_path = '/mnt/pgo4/pgo4_v1/Han-Hsuan\Ptychography_test\20241009_AlGaAs-90s_arm\trial6\';
par.base_path = strrep(base_path,'\','/');
par.roi_label = '0_Ndp400mask';
par.scan_format = '%01d';
par.Ndp = 400;  % size of cbed
par.alpha0 = 25.0; % semi-convergen1e angle (mrad)
bin = 2;

Niter=100;

par.defocus = -200; %overfocus is negative
par.energy = 300;
par.rbf = 400./2/bin;

par.cen_dp_y = floor(par.Ndp/2)+1;
par.cen_dp_x = floor(par.Ndp/2)+1;

par.scan_nx = 120;
par.scan_ny = 122;

par.scan_step_size_x = 0.417;
par.scan_step_size_y = 0.417;

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
par.method = 'MLs';
par.momentum = 0;

par.Nprobe = 5;
par.grouping = 32;
par.apply_multimodal_update = false;

par.Nlayers = 2;
par.regularize_layers = 1;
par.variable_probe_modes = 1;
par.Ndp_presolve = 400;
par.alpha_max = 25.0;
par.thickness = 80;
par.beta_probe = 0.3;

%{
% Step 1.5 (optional): Run a single reconstruction to check parameters
par.GPU_list = [1];



par.thickness = 300; %in angstroms

%data_error_list = zeros(1, length(scan_custom_fliplr));

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
%}
%{
%% plot result
scatter(rot_ang,data_error_list);
% Add labels and title
xlabel('Rotation Angle');
ylabel('Fourier Error');
title('');
saveas(gcf,strcat(par.base_path,num2str(par.scan_number),'/rot_angError.tiff'));
%}

%% Step 2: Use Bayesian optimization with Gaussian processes to find experimental parameters that minimize data error
% Note: Parallel BO is generally recommended for multislice reconstructions
close all

par.rot_ang = 52.2;
par.GPU_list = [1];

par.scan_custom_fliplr = 1;
par.scan_custom_flipud = 1;
par.scan_custom_transpose = 1;

rot_ang = optimizableVariable('rot_ang', [50.5, 54.2]);
defocus = optimizableVariable('defocus', [-220, -180]); %angstroms

output_dir_suffix_base = strcat('_flip', num2str(par.scan_custom_fliplr), num2str(par.scan_custom_flipud), num2str(par.scan_custom_transpose));
par.output_dir_suffix_base = strrep(output_dir_suffix_base,'\','/');

N_workers = length(par.GPU_list);
if N_workers>1
    delete(gcp('nocreate'))
    c = parcluster('local');
    c.NumWorkers = N_workers;
    p = parpool(c);
end

fun = @(x)ptycho_recon_exp_data(par, 'rot_ang', x.rot_ang, 'defocus', x.defocus);
results = bayesopt(fun, [rot_ang, defocus],...
    'Verbose', 4,...
    'AcquisitionFunctionName', 'expected-improvement-plus',...
    'IsObjectiveDeterministic', false,...
    'MaxObjectiveEvaluations', 50,...
    'NumSeedPoints', N_workers,...
    'PlotFcn', {@plotObjectiveModel, @plotMinObjective}, ...
    'UseParallel', N_workers>1);

delete(gcp('nocreate'))