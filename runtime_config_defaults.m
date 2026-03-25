function config = runtime_config_defaults()
% RUNTIME_CONFIG_DEFAULTS Canonical runtime and Dynare defaults for main.m

config = struct();

config.run_firstorder_simul = true;
config.run_secondorder_simul = true;
config.run_firstorder_irs = true;
config.run_secondorder_irs = false;
config.run_pf_irs = true;
config.run_pf_simul = false;
config.run_mit_shocks_simul = false;
config.sector_indices = [1];
config.verbose = true;
config.continue_on_failure = false;
config.ir_horizon = 200;
config.rng_seed = [];
config.model_type = 'VA';
config.simul_T = 1000;
config.simul_burn_in = 100;
config.simul_burn_out = 100;
config.mit_lookahead_horizon = 200;
config.mit_solver_mode = 'rolling';

config.smooth = false;
config.wds = false;
config.covariance_scale = 1.0;
config.save_results = true;
config.force_recalibrate = false;

config.date = "_March_2026";
config.exp_label = "_baseline";

config.gridpoints = 16;
config.sol_guess_file = 'SS_CDsolution_norm_permanent.mat';
config.fsolve_options = optimset('Display','iter','TolX',1e-10,'TolFun',1e-10, ...
    'MaxFunEvals',10000000,'MaxIter',10000);

config.ir_plot_length = 60;
config.plot_irs = false;

config.shock_sizes_pct = [10,20,30];
config.shock_scaling = struct('sectors', [], 'factor', 1.0);
config.shock_values = build_shock_values(config.shock_sizes_pct);
end
