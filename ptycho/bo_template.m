directory_content = dir; % contains everything of the current directory
exe_path = directory_content(1).folder; % returns the path that is currently open

scriptfolder = '/mnt/pgo4/pgo4_v1/Han-Hsuan/fold_slice/ptycho';
scriptfolder = strrep(scriptfolder,'\','/');
cd(scriptfolder);

addpath(strcat(pwd,'/utils/'))
addpath(core.find_base_package)
%%
par = {};
par.verbose_level = 3;
par.scan_number = 3;
par.beam_source = 'electron';

base_path = '\\PanGroupOffice4\PGO4_v1\Han-Hsuan\Ptychography_test\PT3_-20df_285CL_41Mx_n64\';
par.base_path = strrep(base_path,'\','/');
par.base_path = strrep(par.base_path,'//PanGroupOffice4/PGO4_v1','/mnt/pgo4/pgo4_v1');

par.roi_label = '0_Ndp124';
par.scan_format = '%01d';
par.Ndp = 124;  % size of cbed
par.alpha0 = 25.0; % semi-convergen1e angle (mrad)

Niter=100;

par.defocus = -150; %overfocus is negative
par.energy = 300;
par.rbf = 22.0;

par.cen_dp_y = floor(par.Ndp/2)+1;
par.cen_dp_x = floor(par.Ndp/2)+1;

par.scan_nx = 64;
par.scan_ny = 64;

par.scan_step_size_x = 0.367;
par.scan_step_size_y = 0.367;

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

par.Nprobe = 5;
par.grouping = 16;
par.apply_multimodal_update = false;

par.Nlayers = 20;
par.regularize_layers = 0.3;
par.variable_probe_modes = 1;
par.Ndp_presolve = par.Ndp;
par.alpha_max = 25.0;
par.thickness = 200;
par.beta_probe = 0.3;
par.beta_object = 0.8;
par.delta_p = 0.1;
par.probe_change_start=5;

par.rot_ang = 0;
par.GPU_list = [1,2,3,4];


defocus = optimizableVariable('defocus', [-200, 300]); %angstroms
par.output_dir_suffix_base = '';
%%
N_workers = length(par.GPU_list);
if N_workers>1
    delete(gcp('nocreate'))
    c = parcluster('local');
    c.NumWorkers = N_workers;
    p = parpool(c);
end

fun = @(x)ptycho_recon_exp_data(par, 'defocus', x.defocus);
results = bayesopt(fun, [defocus],...
    'Verbose', 4,...
    'AcquisitionFunctionName', 'expected-improvement-plus',...
    'IsObjectiveDeterministic', false,...
    'MaxObjectiveEvaluations', 32,...
    'NumSeedPoints', N_workers,...
    'PlotFcn', {@plotObjectiveModel, @plotMinObjective}, ...
    'UseParallel', N_workers>1);

delete(gcp('nocreate'))