%clear variables
addpath(strcat(pwd,'/utils/'))
addpath(core.find_base_package)
utils.ccc
% Step 0: Run the prepare_data script to generate data for ptycho reconstruction
% Step 1: Prepare data and reconstruction parameters
par = {};
par.verbose_level = 3;
par.scan_number = 1;
par.beam_source = 'electron';

par.Ndp = 360;
base_path = 'D:\Wuhanhsuan\20240806BTO Trial8 defect\';
par.base_path = strrep(base_path,'\','/');

clear base_path;

par.scan_format = '%01d';
par.scan_string_format = '%01d';
par.roi_label = '0_Ndp360'; %don't forget to change this to match file name

par.cen_dp_y = floor(par.Ndp/2)+1;
par.cen_dp_x = floor(par.Ndp/2)+1;

par.energy = 100;
par.rbf = 149.3/2;

par.scan_nx = 100;
par.scan_ny = 100;

par.scan_step_size_x = 0.4;
par.scan_step_size_y = 0.4;

par.scan_custom_fliplr = 0;
par.scan_custom_flipud = 0;
par.scan_custom_transpose = 0;

par.detector_name = 'ELA';
par.data_preparator = 'matlab_aps';
par.src_positions =  'matlab_pos';
par.scan_type = 'raster';

par.use_model_probe = true;
par.normalize_init_probe = true;

par.output_dir_base = par.base_path;
par.Niter = 100;
par.Niter_save_results_every = 100;
par.save.save_reconstructions = true;

par.eng_name = 'GPU_MS';
par.method = 'MLc';
par.momentum = 0.5;

par.Nprobe = 2;
par.grouping = inf;
par.apply_multimodal_update = false;

par.Nlayers = 1;
par.regularize_layers = 1;
par.variable_probe_modes = 1;
par.Ndp_presolve = 360;

%% Step 1.5 (optional): Run a single reconstruction to check parameters
par.GPU_list = 1;
par.rot_ang = 91.5;
par.alpha_max = 24.35;

defocus = 100; %overfocus is negative
par.thickness = 164; %in angstroms

output_dir_suffix_base = '';
par.output_dir_suffix_base = strrep(output_dir_suffix_base,'\','/');

%data_error = ptycho_recon_exp_data(par, 'defocus', defocus, 'thickness', thickness);

%% Step 2: Use Bayesian optimization with Gaussian processes to find experimental parameters that minimize data error
% Note: Parallel BO is generally recommended for multislice reconstructions
close all
par.alpha_max = 24.35;
par.GPU_list = [1];

defocus = optimizableVariable('defocus', [90, 110]); %angstroms
rot_ang = optimizableVariable('rot_ang', [-180, 180]); %angstroms

N_workers = length(par.GPU_list);
if N_workers>1
    delete(gcp('nocreate'))
    c = parcluster('local');
    c.NumWorkers = N_workers;
    p = parpool(c);
end

fun = @(x)ptycho_recon_exp_data(par, 'defocus', x.defocus, 'rot_ang', x.rot_ang);
results = bayesopt(fun, [defocus, rot_ang],...
    'Verbose', 4,...
    'AcquisitionFunctionName', 'expected-improvement-plus',...
    'IsObjectiveDeterministic', false,...
    'MaxObjectiveEvaluations', 30,...
    'NumSeedPoints', N_workers,...
    'PlotFcn', {@plotObjectiveModel, @plotMinObjective}, 'UseParallel', N_workers>1);

delete(gcp('nocreate'))
