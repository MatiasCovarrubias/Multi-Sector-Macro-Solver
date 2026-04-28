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
fprintf('Targets: sig_c=%.2f sig_m=%.2f sig_q=%.2f sig_y=%.2f sig_I=%.2f sig_l=%.2f\n', ...
    params.sigma_c, params.sigma_m, params.sigma_q, params.sigma_y, params.sigma_I, params.sigma_l);

[ModData, params, cache_info] = load_or_build_cached_steady_state( ...
    params, config, exp_paths, struct( ...
    'verbose', true, ...
    'save_local_cache', config.save_results));

if strcmp(cache_info.source, 'experiment_cache') || strcmp(cache_info.source, 'project_cache')
    fprintf('Loaded cached SS (%s): %s\n', cache_info.source, cache_info.path);
elseif cache_info.saved
    fprintf('SS cached: %s (%.1fs)\n', cache_info.path, cache_info.elapsed_seconds);
else
    fprintf('SS not cached because save_results=false (%.1fs)\n', cache_info.elapsed_seconds);
end

%% Base simulation
fprintf('\n--- Dynare Analysis ---\n');
n_shocks = numel(config.shock_values);
run_any_irs = config.run_firstorder_irs || config.run_secondorder_irs || config.run_pf_irs;
run_any_simul = config.run_firstorder_simul || config.run_secondorder_simul || config.run_pf_simul || config.run_mit_shocks_simul;
needs_base_dynare_run = run_any_simul || ~run_any_irs || n_shocks == 0;

AllShockResults = struct('DynareResults', {cell(n_shocks,1)}, ...
    'ShockArtifacts', {cell(n_shocks,1)});
BaseRuntimeResults = struct();

if needs_base_dynare_run
    dynare_opts_base = build_dynare_opts(config, sector_indices, 'base');
    if config.run_mit_shocks_simul && strcmpi(config.mit_solver_mode, 'rolling')
        dynare_opts_base.mit_checkpoint_file = fullfile(exp_paths.temp, 'mit_rolling_checkpoint.mat');
        dynare_opts_base.mit_checkpoint_stride = max(1, min(50, ceil(config.simul_T / 20)));
    end
    if n_shocks >= 1
        params.IRshock = config.shock_values(1).value;
    end
    if run_any_simul
        print_simulation_config(config);
    else
        fprintf('No simulation requested; running Dynare once to compute/save TheoStats.\n');
    end

    tic;
    BaseRuntimeResults = run_dynare_analysis(ModData, params, dynare_opts_base);
    if run_any_simul
        fprintf('Base simulation: %.1fs\n', toc);
    else
        fprintf('TheoStats-only Dynare pass: %.1fs\n', toc);
    end
end

%% IRF analysis
if run_any_irs
    [AllShockResults, BaseRuntimeResults] = run_irf_loop(config, sector_indices, ...
        ModData, params, BaseRuntimeResults, AllShockResults, exp_paths, labels);
else
    fprintf('\nIRF analysis skipped\n');
end

if ~isfield(BaseRuntimeResults, 'TheoStats') || ...
        ~isstruct(BaseRuntimeResults.TheoStats) || ...
        isempty(fieldnames(BaseRuntimeResults.TheoStats))
    error('main:MissingTheoStats', ...
        ['TheoStats were not produced by the Dynare workflow. ' ...
         'Expected BaseRuntimeResults.TheoStats to be a non-empty struct.']);
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

%% IRF summary diagnostics
if has_irfs
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
