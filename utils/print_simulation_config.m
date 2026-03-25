function print_simulation_config(config)
% PRINT_SIMULATION_CONFIG Display simulation method configuration box
%
% INPUTS:
%   config - Configuration structure with simulation settings

fprintf('\n  ┌─ Simulation Configuration ────────────────────────────────────┐\n');
fprintf('  │  Approximation methods:                                        │\n');
fprintf('  │  Common windows: burn_in=%d, active=%d, burn_out=%d         │\n', ...
    config.simul_burn_in, config.simul_T, config.simul_burn_out);
if config.run_firstorder_simul
    fprintf('  │    • First-order (linear):    ON  (3-window schema)          │\n');
else
    fprintf('  │    • First-order (linear):    OFF                            │\n');
end
if config.run_secondorder_simul
    fprintf('  │    • Second-order (quadratic):ON  (3-window schema)          │\n');
else
    fprintf('  │    • Second-order (quadratic):OFF                            │\n');
end
if config.run_pf_simul
    pf_total = config.simul_burn_in + config.simul_T + config.simul_burn_out;
    fprintf('  │    • Perfect foresight (NL):  ON  (T_total=%d)              │\n', pf_total);
    fprintf('  │      Uses shared active shocks with zero-shock boundaries    │\n');
else
    fprintf('  │    • Perfect foresight (NL):  OFF                            │\n');
end
if config.run_mit_shocks_simul
    mit_total = config.simul_T + config.simul_burn_out;
    fprintf('  │    • MIT shocks (surprise):   ON  (T_total=%d, burn_in=0)   │\n', mit_total);
    fprintf('  │      Active shocks start immediately, then burn-out          │\n');
    fprintf('  │      MIT solver mode: %s                                     │\n', ...
        upper(char(config.mit_solver_mode)));
    if strcmpi(char(config.mit_solver_mode), 'rolling')
        fprintf('  │      Rolling horizon input: %d periods                      │\n', ...
            config.mit_lookahead_horizon);
        fprintf('  │      Internal rolling details use fixed horizon + tail      │\n');
    end
else
    fprintf('  │    • MIT shocks (surprise):   OFF                            │\n');
end
fprintf('  │  Moment sample: shocks_simul only                            │\n');
fprintf('  └────────────────────────────────────────────────────────────────┘\n\n');

end
