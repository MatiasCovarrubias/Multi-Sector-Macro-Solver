function config = runtime_config()
% RUNTIME_CONFIG User-facing runtime and Dynare configuration for main.m

% Set to false to use the values below.
use_defaults = false;

config = struct();

config.run_firstorder_simul = true;
config.run_secondorder_simul = false;
config.run_firstorder_irs = true;
config.run_secondorder_irs = false;
config.run_pf_irs = true;
config.run_pf_simul = false;
config.run_mit_shocks_simul = false;
config.sector_indices = 1:37;
% config.sector_indices = [1,3,4,5,6,8,20,24,28];
config.verbose = true;
config.continue_on_failure = false;
config.ir_horizon = 200;
% Set to an integer for reproducible simulated shocks, or [] to use the current RNG state.
config.rng_seed = 1;
config.model_type = 'GO_noVA';
config.simul_T = 5000;
config.simul_burn_in = 100;
config.simul_burn_out = 150;
config.mit_lookahead_horizon = 150;
config.mit_solver_mode = 'rolling';

% Set mit_solver_mode = 'legacy' to use the previous expectation-errors solver.

config.smooth = false;
config.wds = true;
config.covariance_scale = 0.0;
config.save_results = true;
config.force_recalibrate = true;

config.date = "_April_2026";
config.exp_label = "_benchmark_allsectors";
config.force_recalibrate = false;

config.gridpoints = 4;
config.sol_guess_file = '';
config.fsolve_options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10, ...
    'MaxFunEvals',10000000,'MaxIter',10000);

config.ir_plot_length = 60;
config.plot_irs = false;

config.shock_sizes_pct = [5,15,25];
% config.shock_sizes_pct = [25];
% config.shock_scaling = struct('sectors', [1], 'factor', 2.1);
% config.shock_scaling = struct('sectors', [1], 'factor', 1.6);
config.shock_scaling = struct('sectors', [], 'factor', 1.0);

if use_defaults
    config = runtime_config_defaults();
end

config.shock_values = build_shock_values(config.shock_sizes_pct);
end
