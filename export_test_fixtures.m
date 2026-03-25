function ExportInfo = export_test_fixtures(cfg)
% EXPORT_TEST_FIXTURES Save reusable fixtures from a MATLAB/Dynare run

current_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder, 'testing'), '-begin');

if nargin < 1 || isempty(cfg)
    cfg = build_test_defaults();
end

old_dir = pwd;
cleanup_dir = onCleanup(@() cd(old_dir));
cd(cfg.project_root);

setup_test_environment(struct('require_dynare', true));

assert(~isempty(which('dynare')), ...
    'export_test_fixtures:DynareMissing', ...
    'Dynare must be on the MATLAB path before exporting fixtures.');

assert_required_local_files({'calibration_data.mat', 'TFP_process.mat', 'beadat_37sec.mat'});

if exist(cfg.fixture_dir, 'dir') ~= 7
    mkdir(cfg.fixture_dir);
end

runtime_config = build_test_runtime_config(struct( ...
    'model_type', cfg.model_type, ...
    'covariance_scale', cfg.covariance_scale, ...
    'sector_indices', cfg.sector_indices, ...
    'force_recalibrate', cfg.force_recalibrate, ...
    'date', cfg.fixture.date, ...
    'exp_label', cfg.fixture.exp_label, ...
    'gridpoints', cfg.fixture.gridpoints, ...
    'ir_horizon', cfg.fixture.ir_horizon, ...
    'ir_plot_length', cfg.fixture.ir_plot_length, ...
    'simul_T', cfg.fixture.simul_T, ...
    'simul_burn_in', cfg.fixture.simul_burn_in, ...
    'simul_burn_out', cfg.fixture.simul_burn_out, ...
    'shock_sizes_pct', cfg.fixture.shock_sizes_pct, ...
    'run_firstorder_simul', cfg.fixture.run_firstorder_simul, ...
    'run_secondorder_simul', false, ...
    'run_pf_simul', false, ...
    'run_mit_shocks_simul', false, ...
    'run_firstorder_irs', cfg.fixture.run_firstorder_irs, ...
    'run_secondorder_irs', false, ...
    'run_pf_irs', cfg.fixture.run_pf_irs));
save_label = strcat(runtime_config.date, runtime_config.exp_label);
fixture_exp_paths = setup_experiment_folder([save_label '_fixture_run'], cfg.fixture_dir);

params = build_default_params();
[calib_data, params] = load_calibration_data( ...
    params, cfg.sector_indices, cfg.model_type, runtime_config.shock_scaling, ...
    runtime_config.smooth, runtime_config.wds, runtime_config.covariance_scale);
labels = calib_data.labels;

[ModData, params] = load_or_build_cached_steady_state(params, runtime_config, fixture_exp_paths, struct( ...
    'project_root', cfg.project_root, ...
    'allow_project_cache', cfg.allow_project_cache, ...
    'verbose', false));

base_opts = build_dynare_opts(runtime_config, cfg.sector_indices, 'base');
base_opts.verbose = false;
params_base = params;
params_base.IRshock = runtime_config.shock_values(1).value;
BaseResults = run_dynare_analysis(ModData, params_base, base_opts);

irf_opts = build_dynare_opts(runtime_config, cfg.sector_indices, 'irf');
irf_opts.verbose = false;
params_irf = params;
params_irf.IRshock = runtime_config.shock_values(1).value;
DynareResults = run_dynare_analysis(ModData, params_irf, irf_opts);

process_opts = struct();
process_opts.plot_graphs = false;
process_opts.save_graphs = false;
process_opts.save_intermediate = false;
process_opts.exp_paths = fixture_exp_paths;
process_opts.save_label = runtime_config.shock_values(1).label;
process_opts.ir_horizon = runtime_config.ir_horizon;
process_opts.ir_plot_length = runtime_config.ir_plot_length;
process_opts.shock_description = runtime_config.shock_values(1).description;
process_opts.shock_config = runtime_config.shock_values(1);
IRFResults = process_sector_irs(DynareResults, params_irf, ModData, labels, process_opts);

AllShockResults = struct();
AllShockResults.DynareResults = {DynareResults};
AllShockResults.ShockArtifacts = {IRFResults};

flags = struct('has_1storder', true, 'has_2ndorder', false, 'has_pf', false, 'has_mit', false);
has_irfs = true;

fixture_info = struct();
fixture_info.created_at = char(datetime('now', 'Format', "yyyyMMdd'T'HHmmss"));
fixture_info.save_label = save_label;
fixture_info.sector_indices = cfg.sector_indices;
fixture_info.model_type = cfg.model_type;
fixture_info.base_results_fields = fieldnames(BaseResults);
fixture_info.irf_fields = fieldnames(IRFResults);
fixture_info.flags = flags;
fixture_info.has_irfs = has_irfs;

save(cfg.stage_inputs_fixture, 'fixture_info', 'runtime_config', 'save_label', 'fixture_exp_paths', ...
    'calib_data', 'labels');
save(cfg.steady_state_fixture, 'fixture_info', 'runtime_config', 'save_label', 'fixture_exp_paths', ...
    'calib_data', 'labels', 'params', 'ModData');
save(cfg.fixture_bundle, 'fixture_info', 'runtime_config', 'save_label', 'fixture_exp_paths', ...
    'calib_data', 'labels', 'params', 'ModData', 'BaseResults', 'AllShockResults');

ExportInfo = struct();
ExportInfo.stage_inputs_fixture = cfg.stage_inputs_fixture;
ExportInfo.steady_state_fixture = cfg.steady_state_fixture;
ExportInfo.fixture_bundle = cfg.fixture_bundle;
ExportInfo.save_label = save_label;
ExportInfo.exp_paths = fixture_exp_paths;

fprintf('Stage inputs fixture: %s\n', ExportInfo.stage_inputs_fixture);
fprintf('Steady-state fixture: %s\n', ExportInfo.steady_state_fixture);
fprintf('Packaging fixture: %s\n', ExportInfo.fixture_bundle);
end

function assert_required_local_files(required_files)
missing_files = {};
for i = 1:numel(required_files)
    if exist(required_files{i}, 'file') ~= 2
        missing_files{end+1} = required_files{i}; %#ok<AGROW>
    end
end

assert(isempty(missing_files), ...
    'export_test_fixtures:MissingLocalData', ...
    'Missing required local data files: %s', strjoin(missing_files, ', '));
end
