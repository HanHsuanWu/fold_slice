clear variables
directory_content = dir; % contains everything of the current directory
exe_path = directory_content(1).folder; % returns the path that is currently open

scriptfolder = '/mnt/pgo4/pgo4_v1/Han-Hsuan/fold_slice/ptycho';
scriptfolder = strrep(scriptfolder,'\','/');
cd(scriptfolder);

addpath(strcat(pwd,'/utils/'))
addpath(core.find_base_package)

% Step 0: Run the prepare_data script to generate data for ptycho reconstruction
% Step 1: Prepare data and reconstruction parameters
par = {};
par.verbose_level = 3;
par.scan_number = 1;
par.beam_source = 'electron';

base_path = '\\PanGroupOffice4\PGO4_v1\Han-Hsuan\Ptychography_test\20250304_200kV_AlGaAs_bo\P18_AlGaAs_df=-10nm_rot0_step0.4\';
par.base_path = strrep(base_path,'\','/');
par.base_path = strrep(par.base_path,'//PanGroupOffice4/PGO4_v1','/mnt/pgo4/pgo4_v1');

par.roi_label = '';
par.scan_format = '%01d';
par.Ndp = 180;  % size of cbed

Niter=150;

%par.defocus = -100; %overfocus is negative
par.energy = 200;
par.rbf = 35.0;

par.cen_dp_y = floor(par.Ndp/2)+1;
par.cen_dp_x = floor(par.Ndp/2)+1;

par.scan_nx = 100;
par.scan_ny = 100;

par.scan_step_size_x = 0.4;
par.scan_step_size_y = 0.4;

par.detector_name = 'empad';
par.data_preparator = 'matlab_aps';
par.src_positions =  'matlab_pos';
par.scan_type = 'raster';

par.use_model_probe = true;
par.normalize_init_probe = true;

par.output_dir_base = par.base_path;
par.Niter = Niter;
par.Niter_save_results_every = 5;
par.save.save_reconstructions = true;

par.eng_name = 'GPU_MS';
par.method = 'MLs';

par.Nprobe = 6;
par.grouping = 8;
par.apply_multimodal_update = false;


par.regularize_layers = 0.3;
par.variable_probe_modes = 1;
par.Ndp_presolve = par.Ndp;
par.alpha_max = 25.0;

par.beta_probe = 0.3;
par.beta_object = 0.8;
par.beta_LSQ = 0.5;
par.delta_p = 0.1;
par.rot_ang = 0;
par.probe_change_start = 5;
par.Nlayers = 20;
par.thickness = 200;
par.delta_z = 15.0;
par.preshift_ML_probe = false;
par.diff_pattern_blur = 1.3; %ELA at 200kev
par.detector_upsampling = false;

par.GPU_list = [1,2,3,4];

defocus = optimizableVariable('defocus', [-100, 0]); %angstroms
Nlayers = optimizableVariable('Nlayers', [5, 20],'Type','integer'); %angstroms

par.output_dir_suffix_base = '';

N_workers = length(par.GPU_list);
if N_workers>1
    delete(gcp('nocreate'))
    c = parcluster('local');
    c.NumWorkers = N_workers;
    p = parpool(c);
end

fun = @(x)ptycho_recon_exp_data(par, 'defocus', x.defocus,'Nlayers',x.Nlayers);
results = bayesopt(fun, [defocus,Nlayers],...
    'Verbose', 4,...
    'AcquisitionFunctionName', 'expected-improvement-plus',...
    'IsObjectiveDeterministic', false,...
    'MaxObjectiveEvaluations', 20,...
    'NumSeedPoints', N_workers,...
    'PlotFcn', {@plotObjectiveModel, @plotMinObjective}, ...
    'UseParallel', N_workers>1);

save_path = sprintf([par.base_path, '/summary/bo/']);
save_name = 'bo.mat';

% Find the next available numbered filename
counter = 1;
while exist(fullfile(save_path, save_name), 'file')
    save_name = sprintf('bo_%d.mat', counter);
    counter = counter + 1;
end

if ~exist(save_path, 'dir'); mkdir(save_path); end
save(fullfile(save_path, save_name), 'results');

delete(gcp('nocreate'))