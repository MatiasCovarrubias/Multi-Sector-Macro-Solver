function defaults = default_dynare_runtime_opts()
% DEFAULT_DYNARE_RUNTIME_OPTS Canonical defaults for Dynare runtime options

defaults = struct();
defaults.run_firstorder_simul = true;
defaults.run_secondorder_simul = true;
defaults.run_firstorder_irs = true;
defaults.run_secondorder_irs = false;
defaults.run_pf_irs = true;
defaults.run_pf_simul = false;
defaults.run_mit_shocks_simul = false;
defaults.sector_indices = [1];
defaults.shock_scaling = struct('sectors', [], 'factor', 1.0);
defaults.verbose = true;
defaults.continue_on_failure = false;
defaults.ir_horizon = 200;
defaults.rng_seed = [];
defaults.model_type = 'VA';
defaults.simul_T = 1000;
defaults.simul_burn_in = 100;
defaults.simul_burn_out = 100;
defaults.mit_lookahead_horizon = 200;
defaults.mit_solver_mode = 'rolling';
end
