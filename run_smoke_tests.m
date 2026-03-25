function run_smoke_tests()
% RUN_SMOKE_TESTS Lightweight validation for the MATLAB/Dynare pipeline

current_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder, 'utils'), '-begin');
setup_runtime_paths(current_folder);

fprintf('\nRunning MATLAB smoke tests...\n');

config = build_smoke_config();
params = build_smoke_params();
sector_indices = config.sector_indices;
save_label = strcat(config.date, config.exp_label);

[calib_data, params] = load_calibration_data( ...
    params, sector_indices, config.model_type, config.shock_scaling, ...
    config.smooth, config.wds, config.covariance_scale);
labels = calib_data.labels;
exp_paths = setup_experiment_folder(save_label);
[ModData, params] = load_or_build_smoke_steady_state(config, params, exp_paths);

run_firstorder_smoke(ModData, params, config);
run_secondorder_smoke(ModData, params, config);
run_runtime_validation_smoke(ModData, params, config);
run_pf_irf_smoke(ModData, params, config);
run_irf_pipeline_smoke(ModData, params, config, labels, exp_paths);
run_pf_simulation_smoke(ModData, params, config);
run_mit_smoke(ModData, params, config);
run_save_pipeline_smoke(ModData, params, config, calib_data, labels, exp_paths, save_label);

fprintf('Smoke tests passed.\n');
end

function config = build_smoke_config()
config = build_test_runtime_config(struct( ...
    'date', "_Mar_2026", ...
    'exp_label', "_smoke_tests", ...
    'force_recalibrate', false, ...
    'gridpoints', 8, ...
    'ir_horizon', 40, ...
    'ir_plot_length', 20, ...
    'simul_T', 20, ...
    'simul_burn_in', 5, ...
    'simul_burn_out', 5));
end

function params = build_smoke_params()
params = build_default_params();
end

function [ModData, params] = load_or_build_smoke_steady_state(config, params, exp_paths)
[ModData, params] = load_or_build_cached_steady_state(params, config, exp_paths, struct( ...
    'project_root', current_project_root(), ...
    'allow_project_cache', true, ...
    'verbose', false));
end

function project_root = current_project_root()
project_root = fileparts(mfilename('fullpath'));
end

function run_firstorder_smoke(ModData, params, config)
opts = build_dynare_opts(config, config.sector_indices, 'base');
opts.run_secondorder_simul = false;
opts.run_firstorder_irs = false;
opts.run_secondorder_irs = false;
opts.run_pf_irs = false;
opts.run_pf_simul = false;
opts.run_mit_shocks_simul = false;
params.IRshock = config.shock_values(1).value;

Results = run_dynare_analysis(ModData, params, opts);
assert(isfield(Results, 'SimulFirstOrder'), 'Smoke test failed: missing first-order simulation.');
assert(isfield(Results.SimulFirstOrder, 'summary_stats'), 'Smoke test failed: missing first-order summary_stats.');
assert(isfield(Results.SimulFirstOrder.summary_stats, 'ModelStats'), 'Smoke test failed: missing first-order model stats.');
assert_simulation_block_shape(Results.SimulFirstOrder, config.simul_burn_in, config.simul_T, config.simul_burn_out, false);
assert(strcmp(Results.SimulFirstOrder.summary_stats.ModelStats.sample_window, 'shocks_simul'), ...
    'Smoke test failed: first-order moments must use shocks_simul.');
end

function run_secondorder_smoke(ModData, params, config)
opts = build_dynare_opts(config, config.sector_indices, 'base');
opts.run_firstorder_simul = false;
opts.run_secondorder_simul = true;
opts.run_firstorder_irs = false;
opts.run_secondorder_irs = false;
opts.run_pf_irs = false;
opts.run_pf_simul = false;
opts.run_mit_shocks_simul = false;
params.IRshock = config.shock_values(1).value;

Results = run_dynare_analysis(ModData, params, opts);
assert(isfield(Results, 'SimulSecondOrder'), 'Smoke test failed: missing second-order simulation.');
assert(isfield(Results.SimulSecondOrder, 'summary_stats'), 'Smoke test failed: missing second-order summary_stats.');
assert(isfield(Results.SimulSecondOrder.summary_stats, 'ModelStats'), 'Smoke test failed: missing second-order model stats.');
assert_simulation_block_shape(Results.SimulSecondOrder, config.simul_burn_in, config.simul_T, config.simul_burn_out, false);
assert(strcmp(Results.SimulSecondOrder.summary_stats.ModelStats.sample_window, 'shocks_simul'), ...
    'Smoke test failed: second-order moments must use shocks_simul.');
end

function run_pf_irf_smoke(ModData, params, config)
opts = build_dynare_opts(config, config.sector_indices, 'irf');
opts.run_firstorder_simul = false;
opts.run_secondorder_simul = false;
opts.run_pf_simul = false;
opts.run_mit_shocks_simul = false;
opts.run_firstorder_irs = false;
opts.run_secondorder_irs = false;
opts.run_pf_irs = true;
params.IRshock = config.shock_values(1).value;

Results = run_dynare_analysis(ModData, params, opts);
assert(isfield(Results, 'IRSPerfectForesight_raw'), 'Smoke test failed: missing PF IRFs.');
assert(~isempty(Results.IRSPerfectForesight_raw), 'Smoke test failed: empty PF IRFs.');
end

function run_runtime_validation_smoke(ModData, params, config)
opts = build_dynare_opts(config, config.sector_indices, 'base');
params.IRshock = config.shock_values(1).value;

invalid_opts = opts;
invalid_opts.model_type = 'BAD_MODEL';
assert_validator_fails(@() run_dynare_analysis(ModData, params, invalid_opts), ...
    'validate_dynare_runtime_opts');

invalid_opts = opts;
invalid_opts.simul_T = 0;
assert_validator_fails(@() run_dynare_analysis(ModData, params, invalid_opts), ...
    'validate_dynare_runtime_opts');

assert_validator_fails(@() build_shock_values([100]), 'build_shock_values');
end

function run_irf_pipeline_smoke(ModData, params, config, labels, exp_paths)
opts = build_dynare_opts(config, config.sector_indices, 'irf');
opts.run_firstorder_simul = false;
opts.run_secondorder_simul = false;
opts.run_pf_simul = false;
opts.run_mit_shocks_simul = false;
params.IRshock = config.shock_values(1).value;

DynareResults = run_dynare_analysis(ModData, params, opts);
ir_opts = struct();
ir_opts.plot_graphs = false;
ir_opts.save_graphs = false;
ir_opts.save_intermediate = false;
ir_opts.exp_paths = exp_paths;
ir_opts.save_label = config.shock_values(1).label;
ir_opts.ir_horizon = config.ir_horizon;
ir_opts.ir_plot_length = config.ir_plot_length;
ir_opts.shock_description = config.shock_values(1).description;
ir_opts.shock_config = config.shock_values(1);

IRFResults = process_sector_irs(DynareResults, params, ModData, labels, ir_opts);
assert(isfield(IRFResults, 'entries'), 'Smoke test failed: missing canonical IR entries.');
assert(isfield(IRFResults, 'summary_stats'), 'Smoke test failed: missing canonical IR summary_stats.');
AllShockResults = struct('DynareResults', {{DynareResults}}, 'ShockArtifacts', {{IRFResults}});
ModelData_IRs = build_ModelData_IRs(AllShockResults, config, strcat(config.date, config.exp_label), ...
    config.sector_indices, numel(config.shock_values));

flags = struct('has_1storder', false, 'has_2ndorder', false, 'has_pf', false, 'has_mit', false);
ModelData_IRs.metadata.run_flags = struct( ...
    'has_1storder', logical(config.run_firstorder_irs), ...
    'has_2ndorder', logical(config.run_secondorder_irs), ...
    'has_pf', logical(config.run_pf_irs), ...
    'has_mit', false);
ModelData_IRs.metadata.has_irfs = true;
validate_ModelData_IRs(ModelData_IRs, 'run_irf_pipeline_smoke');

invalid_irf = ModelData_IRs;
invalid_irf.shocks(1).entries(1).sector_idx = -1;
assert_validator_fails(@() validate_ModelData_IRs(invalid_irf, 'run_irf_pipeline_smoke.invalid_irf'), ...
    'validate_ModelData_IRs');
end

function run_pf_simulation_smoke(ModData, params, config)
opts = build_dynare_opts(config, config.sector_indices, 'base');
opts.run_firstorder_simul = false;
opts.run_secondorder_simul = false;
opts.run_firstorder_irs = false;
opts.run_secondorder_irs = false;
opts.run_pf_irs = false;
opts.run_pf_simul = true;
opts.run_mit_shocks_simul = false;
params.IRshock = config.shock_values(1).value;

Results = run_dynare_analysis(ModData, params, opts);
assert(isfield(Results, 'SimulPerfectForesight'), 'Smoke test failed: missing PF simulation.');
assert(~isempty(Results.SimulPerfectForesight), 'Smoke test failed: empty PF simulation.');
assert_simulation_block_shape(Results.SimulPerfectForesight, config.simul_burn_in, config.simul_T, config.simul_burn_out, false);
end

function run_mit_smoke(ModData, params, config)
opts = build_dynare_opts(config, config.sector_indices, 'base');
opts.run_firstorder_simul = false;
opts.run_secondorder_simul = false;
opts.run_firstorder_irs = false;
opts.run_secondorder_irs = false;
opts.run_pf_irs = false;
opts.run_pf_simul = false;
opts.run_mit_shocks_simul = true;
params.IRshock = config.shock_values(1).value;

Results = run_dynare_analysis(ModData, params, opts);
assert(isfield(Results, 'SimulMITShocks'), 'Smoke test failed: missing MIT simulation.');
assert(~isempty(Results.SimulMITShocks), 'Smoke test failed: empty MIT simulation.');
assert_simulation_block_shape(Results.SimulMITShocks, 0, config.simul_T, config.simul_burn_out, true);
end

function run_save_pipeline_smoke(ModData, params, config, calib_data, labels, exp_paths, save_label)
opts = build_dynare_opts(config, config.sector_indices, 'base');
opts.run_secondorder_simul = false;
opts.run_firstorder_irs = false;
opts.run_secondorder_irs = false;
opts.run_pf_irs = false;
opts.run_pf_simul = false;
opts.run_mit_shocks_simul = false;
params.IRshock = config.shock_values(1).value;

BaseResults = run_dynare_analysis(ModData, params, opts);
ModelData = build_ModelData(config, save_label, config.sector_indices, numel(config.shock_values), ...
    calib_data, labels, params, ModData, exp_paths);
[ModelData_simulation, flags] = build_ModelData_simulation( ...
    BaseResults, params, save_label, exp_paths);
ModelData = attach_simulation_statistics(ModelData, ModelData_simulation);
Diagnostics = build_nonlinearity_diagnostics(ModelData_simulation, [], params, config, ModData);
ModelData_IRs = struct();
has_irfs = false;

[ModelData, ModelData_simulation, ModelData_IRs] = finalize_model_outputs( ...
    ModelData, ModelData_simulation, ModelData_IRs, flags, has_irfs, Diagnostics);

validate_pipeline_outputs(ModelData, ModelData_simulation, ModelData_IRs, ...
    flags, has_irfs, 'run_save_pipeline_smoke');

invalid_model = ModelData;
invalid_model.metadata.n_shocks = -1;
assert_validator_fails(@() validate_ModelData(invalid_model, 'run_save_pipeline_smoke.invalid_model'), ...
    'validate_ModelData');

invalid_sim = ModelData_simulation;
invalid_sim.FirstOrder.T_total = invalid_sim.FirstOrder.T_total + 1;
assert_validator_fails(@() validate_ModelData_simulation(invalid_sim, flags, 'run_save_pipeline_smoke.invalid_sim'), ...
    'validate_ModelData_simulation');

save_config = config;
save_config.save_results = true;
save_experiment_results(save_config, build_smoke_params(), exp_paths, ...
    ModelData, ModelData_simulation, ModelData_IRs, flags, has_irfs);
assert(exist(fullfile(exp_paths.experiment, 'runtime_config.mat'), 'file') == 2, ...
    'Smoke test failed: runtime_config.mat was not saved.');
assert(exist(fullfile(exp_paths.experiment, 'params_config.mat'), 'file') == 2, ...
    'Smoke test failed: params_config.mat was not saved.');
assert(exist(fullfile(exp_paths.experiment, 'ModelData.mat'), 'file') == 2, ...
    'Smoke test failed: ModelData.mat was not saved.');
assert(exist(fullfile(exp_paths.experiment, 'ModelData_simulation.mat'), 'file') == 2, ...
    'Smoke test failed: ModelData_simulation.mat was not saved.');
end

function assert_validator_fails(fn, validator_name)
did_fail = false;
try
    fn();
catch
    did_fail = true;
end

assert(did_fail, 'Smoke test failed: %s should reject invalid artifacts.', validator_name);
end

function assert_simulation_block_shape(simul_block, expected_burn_in, expected_active, expected_burn_out, is_mit)
assert(isstruct(simul_block), 'Smoke test failed: simulation block must be a struct.');
assert(size(simul_block.burnin_simul, 2) == expected_burn_in, ...
    'Smoke test failed: burn-in window has wrong length.');
assert(size(simul_block.shocks_simul, 2) == expected_active, ...
    'Smoke test failed: active shock window has wrong length.');
assert(size(simul_block.burnout_simul, 2) == expected_burn_out, ...
    'Smoke test failed: burn-out window has wrong length.');
assert(simul_block.T_active == expected_active, ...
    'Smoke test failed: T_active mismatch.');
assert(simul_block.T_total == expected_burn_in + expected_active + expected_burn_out, ...
    'Smoke test failed: T_total mismatch.');
if is_mit
    assert(simul_block.burn_in == 0, 'Smoke test failed: MIT burn_in must equal 0.');
else
    assert(simul_block.burn_in == expected_burn_in, ...
        'Smoke test failed: burn_in mismatch.');
end
assert(simul_block.burn_out == expected_burn_out, ...
    'Smoke test failed: burn_out mismatch.');
end
