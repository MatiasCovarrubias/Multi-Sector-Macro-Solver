function opts = build_dynare_opts(config, sector_indices, mode)
    opts = struct();
    opts.sector_indices = sector_indices;
    opts.shock_scaling = config.shock_scaling;
    opts.verbose = true;
    opts.ir_horizon = config.ir_horizon;
    opts.rng_seed = config.rng_seed;
    opts.model_type = config.model_type;
    opts.simul_T = config.simul_T;
    opts.simul_burn_in = config.simul_burn_in;
    opts.simul_burn_out = config.simul_burn_out;
    opts.mit_lookahead_horizon = config.mit_lookahead_horizon;
    opts.mit_solver_mode = config.mit_solver_mode;

    switch mode
        case 'base'
            opts.run_firstorder_simul = config.run_firstorder_simul;
            opts.run_secondorder_simul = config.run_secondorder_simul;
            opts.run_firstorder_irs = false;
            opts.run_secondorder_irs = false;
            opts.run_pf_irs = false;
            opts.run_pf_simul = config.run_pf_simul;
            opts.run_mit_shocks_simul = config.run_mit_shocks_simul;

        case 'irf'
            opts.run_firstorder_simul = false;
            opts.run_secondorder_simul = false;
            opts.run_firstorder_irs = config.run_firstorder_irs;
            opts.run_secondorder_irs = config.run_secondorder_irs;
            opts.run_pf_irs = config.run_pf_irs;
            opts.run_pf_simul = false;
            opts.run_mit_shocks_simul = false;

        otherwise
            error('build_dynare_opts:InvalidMode', 'Unknown mode: %s (expected ''base'' or ''irf'')', mode);
    end
end
