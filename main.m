clear; clearvars -global; clc;

%% Paths
current_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder, 'utils'), '-begin');
setup_runtime_paths(current_folder);

%% Verify data files
required_files = {'calibration_data.mat', 'TFP_process.mat'};
missing_files = {};
for i = 1:numel(required_files)
    if ~exist(required_files{i}, 'file')
        missing_files{end+1} = required_files{i}; %#ok<SAGROW>
    end
end
if ~isempty(missing_files)
    error('main_IRs:MissingDataFiles', 'Missing: %s', strjoin(missing_files, ', '));
end

%% Configuration
config = runtime_config();

%% Parameters
params = params_config();
params_config_snapshot = params;
config_sources = struct( ...
    'runtime', fullfile(current_folder, 'runtime_config.m'), ...
    'params', fullfile(current_folder, 'params_config.m'));

%% Derived settings
N_SECTORS = 37;
sector_indices = config.sector_indices;
validate_sector_indices(sector_indices, N_SECTORS, 'main_IRs');
save_label = strcat(config.date, config.exp_label);

%% Experiment folder
exp_paths = setup_experiment_folder(save_label);
fprintf('\nExperiment: %s | Sectors: [%s]\n', save_label, num2str(sector_indices));

%% Load calibration
[calib_data, params] = load_calibration_data( ...
    params, sector_indices, config.model_type, config.shock_scaling, ...
    config.smooth, config.wds, config.covariance_scale);
labels = calib_data.labels;
fprintf('Calibration loaded (model: %s, smooth: %s, wds: %s, covariance_scale: %.3f)\n', ...
    config.model_type, mat2str(calib_data.smooth), mat2str(calib_data.wds), ...
    calib_data.covariance_scale);
print_empirical_targets(calib_data.empirical_targets);

%% Steady state
fprintf('\n--- Steady State ---\n');
ss_file = fullfile(exp_paths.experiment, 'steady_state.mat');
ss_cached = exist(ss_file, 'file') && ~config.force_recalibrate;

if ss_cached
    [ModData, params, cache_info] = load_or_build_cached_steady_state( ...
        params, config, exp_paths, struct('verbose', false));
    fprintf('Loaded cached SS (%s): %s\n', cache_info.source, cache_info.path);
else
    fprintf('Targets: sig_c=%.2f sig_m=%.2f sig_q=%.2f sig_y=%.2f sig_I=%.2f sig_l=%.2f\n', ...
        params.sigma_c, params.sigma_m, params.sigma_q, params.sigma_y, params.sigma_I, params.sigma_l);

    calib_opts = struct();
    calib_opts.gridpoints = config.gridpoints;
    calib_opts.verbose = true;
    calib_opts.sol_guess_file = config.sol_guess_file;
    calib_opts.fsolve_options = config.fsolve_options;

    tic;
    [ModData, params] = calibrate_steady_state(params, calib_opts);
    elapsed_calib = toc;

    if config.save_results
        save(ss_file, 'ModData', 'params');
        fprintf('SS cached: %s (%.1fs)\n', ss_file, elapsed_calib);
    else
        fprintf('SS not cached because save_results=false (%.1fs)\n', elapsed_calib);
    end
end

%% Base simulation
fprintf('\n--- Dynare Analysis ---\n');
n_shocks = numel(config.shock_values);
run_any_irs = config.run_firstorder_irs || config.run_secondorder_irs || config.run_pf_irs;
run_any_simul = config.run_firstorder_simul || config.run_secondorder_simul || config.run_pf_simul || config.run_mit_shocks_simul;

AllShockResults = struct('DynareResults', {cell(n_shocks,1)}, ...
    'ShockArtifacts', {cell(n_shocks,1)});
BaseRuntimeResults = struct();

if run_any_simul
    dynare_opts_base = build_dynare_opts(config, sector_indices, 'base');
    if config.run_mit_shocks_simul && strcmpi(config.mit_solver_mode, 'rolling')
        dynare_opts_base.mit_checkpoint_file = fullfile(exp_paths.temp, 'mit_rolling_checkpoint.mat');
        dynare_opts_base.mit_checkpoint_stride = max(1, min(50, ceil(config.simul_T / 20)));
    end
    params.IRshock = config.shock_values(1).value;
    print_simulation_config(config);

    tic;
    BaseRuntimeResults = run_dynare_analysis(ModData, params, dynare_opts_base);
    fprintf('Base simulation: %.1fs\n', toc);
end

%% IRF analysis
if run_any_irs
    [AllShockResults, BaseRuntimeResults] = run_irf_loop(config, sector_indices, ...
        ModData, params, BaseRuntimeResults, AllShockResults, exp_paths, labels);
else
    fprintf('\nIRF analysis skipped\n');
end

%% Package results
ModelData = build_ModelData(config, save_label, sector_indices, n_shocks, ...
    calib_data, labels, params, ModData, exp_paths);

simulation_artifacts = struct();
[ModelData_simulation, flags] = build_ModelData_simulation( ...
    BaseRuntimeResults, params, save_label, exp_paths, simulation_artifacts);
ModelData = attach_simulation_statistics(ModelData, ModelData_simulation);

%% Build IRF data
has_irfs = run_any_irs && isfield(AllShockResults, 'ShockArtifacts') && ...
    ~isempty(AllShockResults.ShockArtifacts) && ~isempty(AllShockResults.ShockArtifacts{1});
ModelData_IRs = build_ModelData_IRs(AllShockResults, config, save_label, sector_indices, n_shocks);

%% Diagnostics
if flags.has_1storder || flags.has_2ndorder || flags.has_pf || flags.has_mit || has_irfs
    Diagnostics = build_nonlinearity_diagnostics(ModelData_simulation, AllShockResults, params, config, ModData);
else
    Diagnostics = [];
end

[ModelData, ModelData_simulation, ModelData_IRs] = finalize_model_outputs( ...
    ModelData, ModelData_simulation, ModelData_IRs, flags, has_irfs, Diagnostics);

validate_pipeline_outputs(ModelData, ModelData_simulation, ModelData_IRs, ...
    flags, has_irfs, 'main');

%% Save
save_experiment_results(config, params_config_snapshot, exp_paths, ...
    ModelData, ModelData_simulation, ModelData_IRs, flags, has_irfs, config_sources);

%% Summary
print_summary_table(ModelData);
