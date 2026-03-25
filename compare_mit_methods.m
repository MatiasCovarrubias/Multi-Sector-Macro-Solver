function comparison = compare_mit_methods(overrides)
% compare_mit_methods  Compare legacy and current rolling MIT solvers.

if nargin < 1 || isempty(overrides)
    overrides = struct();
end

current_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder, 'utils'), '-begin');
setup_runtime_paths(current_folder);

config = runtime_config_defaults();
config.date = '_Mar_2026';
config.exp_label = '_matlab_smoke';
config.force_recalibrate = false;
config.save_results = false;
config.run_firstorder_simul = false;
config.run_secondorder_simul = false;
config.run_firstorder_irs = false;
config.run_secondorder_irs = false;
config.run_pf_irs = false;
config.run_pf_simul = false;
config.run_mit_shocks_simul = true;
config.simul_T = 4;
config.simul_burn_in = 0;
config.simul_burn_out = 2;
config.verbose = false;
config.rng_seed = 0;
config = apply_overrides(config, overrides);

if ~isfield(overrides, 'mit_lookahead_horizon') || isempty(overrides.mit_lookahead_horizon)
    config.mit_lookahead_horizon = config.simul_T + config.simul_burn_out;
end

save_label = strcat(config.date, config.exp_label);
ss_file = fullfile(current_folder, 'experiments', save_label(2:end), 'steady_state.mat');
if ~exist(ss_file, 'file')
    error('compare_mit_methods:MissingSteadyState', ...
        'Cached steady state not found: %s', ss_file);
end

params = params_config();
[~, params] = load_calibration_data( ...
    params, config.sector_indices, config.model_type, config.shock_scaling, ...
    config.smooth, config.wds, config.covariance_scale);
exp_paths = struct('experiment', fullfile(current_folder, 'experiments', save_label(2:end)));
[ModData, params] = load_or_build_cached_steady_state( ...
    params, config, exp_paths, struct('verbose', false));
params.IRshock = config.shock_values(1).value;

n_sectors = params.n_sectors;
if ~isempty(config.rng_seed)
    rng(config.rng_seed);
end
shared_shocks = mvnrnd(zeros([n_sectors,1]), params.Sigma_A, config.simul_T);

opts = build_dynare_opts(config, config.sector_indices, 'base');
opts.verbose = config.verbose;
opts.continue_on_failure = false;
opts.external_shocks = shared_shocks;

opts.mit_solver_mode = 'legacy';
legacy_results = run_dynare_analysis(ModData, params, opts);

opts.mit_solver_mode = 'rolling';
legacy_path = [legacy_results.SimulMITShocks.shocks_simul legacy_results.SimulMITShocks.burnout_simul];
idx = legacy_results.SimulMITShocks.variable_indices;

comparison = struct();
comparison.config = config;
comparison.shared_shocks = shared_shocks;
comparison.legacy = legacy_results;
comparison.rolling = run_rolling_case(ModData, params, opts, legacy_path, idx);

fprintf('\nMIT method comparison\n');
fprintf('  solver horizon: %d\n', config.mit_lookahead_horizon);
fprintf('  active periods: %d\n', config.simul_T);
fprintf('  burn-out:       %d\n', config.simul_burn_out);
print_case_summary('rolling', comparison.rolling);

end

function s = apply_overrides(s, overrides)
if ~isstruct(overrides)
    error('compare_mit_methods:InvalidOverrides', ...
        'overrides must be a struct.');
end

override_fields = fieldnames(overrides);
for i = 1:numel(override_fields)
    field_name = override_fields{i};
    s.(field_name) = overrides.(field_name);
end
end

function case_result = run_rolling_case(ModData, params, base_opts, legacy_path, idx)
case_opts = base_opts;
case_opts.mit_solver_mode = 'rolling';

rolling_results = run_dynare_analysis(ModData, params, case_opts);
rolling_path = [rolling_results.SimulMITShocks.shocks_simul rolling_results.SimulMITShocks.burnout_simul];
path_diff = rolling_path - legacy_path;

case_result = struct();
case_result.results = rolling_results;
case_result.path = rolling_path;
case_result.path_diff = path_diff;
case_result.max_abs_diff = max(abs(path_diff(:)));
case_result.mean_abs_diff = mean(abs(path_diff(:)));
case_result.max_abs_diff_by_period = max(abs(path_diff), [], 1);
case_result.max_abs_diff_by_variable = max(abs(path_diff), [], 2);
case_result.max_abs_diff_capital = max(abs(path_diff(idx.k(1):idx.k(2), :)), [], 'all');
case_result.max_abs_diff_tfp = max(abs(path_diff(idx.a(1):idx.a(2), :)), [], 'all');
case_result.max_abs_diff_cagg = max(abs(path_diff(idx.c_agg, :)), [], 'all');
case_result.max_abs_diff_yagg = max(abs(path_diff(idx.gdp_agg, :)), [], 'all');
case_result.max_abs_diff_iagg = max(abs(path_diff(idx.i_agg, :)), [], 'all');
if size(rolling_path, 2) >= 2
    shifted_diff = rolling_path(:, 2:end) - legacy_path(:, 1:end-1);
    case_result.max_abs_diff_shifted = max(abs(shifted_diff(:)));
else
    case_result.max_abs_diff_shifted = NaN;
end
end

function print_case_summary(label, case_result)
fprintf('\n  %s\n', label);
fprintf('    max abs diff:   %.6e\n', case_result.max_abs_diff);
fprintf('    mean abs diff:  %.6e\n', case_result.mean_abs_diff);
fprintf('    max |diff| k:   %.6e\n', case_result.max_abs_diff_capital);
fprintf('    max |diff| a:   %.6e\n', case_result.max_abs_diff_tfp);
fprintf('    max |diff| C:   %.6e\n', case_result.max_abs_diff_cagg);
fprintf('    max |diff| Y:   %.6e\n', case_result.max_abs_diff_yagg);
fprintf('    max |diff| I:   %.6e\n', case_result.max_abs_diff_iagg);
fprintf('    shifted by 1:   %.6e\n', case_result.max_abs_diff_shifted);
end
