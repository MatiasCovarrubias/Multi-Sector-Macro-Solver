function Results = run_fixture_tests(cfg)
% RUN_FIXTURE_TESTS Replay packaging and validation steps from fixtures

current_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder, 'testing'), '-begin');

if nargin < 1 || isempty(cfg)
    cfg = build_test_defaults();
elseif ischar(cfg) || isstring(cfg)
    override_cfg = build_test_defaults();
    override_cfg.fixture_bundle = char(cfg);
    cfg = override_cfg;
end

old_dir = pwd;
cleanup_dir = onCleanup(@() cd(old_dir));
cd(cfg.project_root);

setup_test_environment();

assert(exist(cfg.fixture_bundle, 'file') == 2, ...
    'run_fixture_tests:MissingFixtureBundle', ...
    'Fixture bundle not found: %s', cfg.fixture_bundle);

log_file = start_log(cfg.log_dir, 'fixture');
cleanup_log = onCleanup(@() stop_log());

fprintf('\nRunning fixture tests from %s\n', cfg.fixture_bundle);

loaded = load(cfg.fixture_bundle);
assert_required_fixture_fields(loaded, ...
    {'runtime_config', 'save_label', 'fixture_exp_paths', 'calib_data', 'labels', ...
     'params', 'ModData', 'BaseResults', 'AllShockResults'});

Results = struct();
Results.ok = true;
Results.log_file = log_file;
Results.fixture_bundle = cfg.fixture_bundle;

replay_exp_paths = setup_experiment_folder([loaded.save_label '_fixture_replay'], cfg.replay_root);
runtime_config = loaded.runtime_config;
runtime_config.save_results = true;

ModelData = build_ModelData(runtime_config, loaded.save_label, runtime_config.sector_indices, ...
    numel(runtime_config.shock_values), loaded.calib_data, loaded.labels, loaded.params, ...
    loaded.ModData, replay_exp_paths);

[ModelData_simulation, flags] = build_ModelData_simulation( ...
    loaded.BaseResults, loaded.params, loaded.save_label, replay_exp_paths);
ModelData = attach_simulation_statistics(ModelData, ModelData_simulation);

has_irfs = isfield(loaded.AllShockResults, 'ShockArtifacts') && ...
    ~isempty(loaded.AllShockResults.ShockArtifacts) && ...
    ~isempty(loaded.AllShockResults.ShockArtifacts{1});

ModelData_IRs = build_ModelData_IRs(loaded.AllShockResults, runtime_config, ...
    loaded.save_label, runtime_config.sector_indices, numel(runtime_config.shock_values));

Diagnostics = build_nonlinearity_diagnostics( ...
    ModelData_simulation, loaded.AllShockResults, loaded.params, runtime_config, loaded.ModData);

[ModelData, ModelData_simulation, ModelData_IRs] = finalize_model_outputs( ...
    ModelData, ModelData_simulation, ModelData_IRs, flags, has_irfs, Diagnostics);

validate_pipeline_outputs(ModelData, ModelData_simulation, ModelData_IRs, ...
    flags, has_irfs, 'run_fixture_tests');

assert_validator_fails(@() fail_model_validation(ModelData), 'validate_ModelData');
assert_validator_fails(@() fail_simulation_validation(ModelData_simulation, flags), 'validate_ModelData_simulation');
if has_irfs
    assert_validator_fails(@() fail_irf_validation(ModelData_IRs), 'validate_ModelData_IRs');
end

save_experiment_results(runtime_config, params_config_defaults(), replay_exp_paths, ...
    ModelData, ModelData_simulation, ModelData_IRs, flags, has_irfs);

assert(exist(fullfile(replay_exp_paths.experiment, 'runtime_config.mat'), 'file') == 2, ...
    'run_fixture_tests:MissingRuntimeConfigSave', ...
    'Replay save pipeline did not write runtime_config.mat');
assert(exist(fullfile(replay_exp_paths.experiment, 'params_config.mat'), 'file') == 2, ...
    'run_fixture_tests:MissingParamsConfigSave', ...
    'Replay save pipeline did not write params_config.mat');

assert(exist(fullfile(replay_exp_paths.experiment, 'ModelData.mat'), 'file') == 2, ...
    'run_fixture_tests:MissingModelDataSave', ...
    'Replay save pipeline did not write ModelData.mat');
assert(exist(fullfile(replay_exp_paths.experiment, 'ModelData_simulation.mat'), 'file') == 2, ...
    'run_fixture_tests:MissingSimulationSave', ...
    'Replay save pipeline did not write ModelData_simulation.mat');
if has_irfs
    assert(exist(fullfile(replay_exp_paths.experiment, 'ModelData_IRs.mat'), 'file') == 2, ...
        'run_fixture_tests:MissingIRSave', ...
        'Replay save pipeline did not write ModelData_IRs.mat');
end

Results.replay_exp_paths = replay_exp_paths;
Results.flags = flags;
Results.has_irfs = has_irfs;
fprintf('Fixture replay completed. Log file: %s\n', Results.log_file);
end

function fail_model_validation(ModelData)
invalid_model = ModelData;
invalid_model.metadata.n_shocks = -1;
validate_ModelData(invalid_model, 'run_fixture_tests.invalid_model');
end

function fail_simulation_validation(ModelData_simulation, flags)
invalid_sim = ModelData_simulation;
invalid_sim.Shocks.T = invalid_sim.Shocks.T + 1;
validate_ModelData_simulation(invalid_sim, flags, 'run_fixture_tests.invalid_sim');
end

function fail_irf_validation(ModelData_IRs)
invalid_irf = ModelData_IRs;
invalid_irf.shocks(1).entries(1).sector_idx = -1;
validate_ModelData_IRs(invalid_irf, 'run_fixture_tests.invalid_irf');
end

function assert_validator_fails(fn, validator_name)
did_fail = false;
try
    fn();
catch
    did_fail = true;
end

assert(did_fail, ...
    'run_fixture_tests:ValidatorDidNotFail', ...
    '%s should reject the mutated artifact.', validator_name);
end

function assert_required_fixture_fields(loaded, required_fields)
for i = 1:numel(required_fields)
    assert(isfield(loaded, required_fields{i}), ...
        'run_fixture_tests:InvalidFixtureBundle', ...
        'Fixture bundle is missing required field: %s', required_fields{i});
end
end

function log_file = start_log(log_dir, prefix)
if exist(log_dir, 'dir') ~= 7
    mkdir(log_dir);
end

timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
log_file = fullfile(log_dir, sprintf('%s_%s.log', prefix, timestamp));
diary(log_file);
end

function stop_log()
try
    diary off;
catch
end
end
