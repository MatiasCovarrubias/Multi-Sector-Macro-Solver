function Results = run_dynare_analysis(ModData, params, opts)

%% Defaults
if nargin < 3, opts = struct(); end

opts = normalize_dynare_opts(opts);
runtime_ctx = resolve_runtime_context(opts);

validate_params(params, {'n_sectors', 'IRshock', 'Sigma_A'}, 'run_dynare_analysis');
validate_dynare_runtime_opts(opts, params.n_sectors, 'run_dynare_analysis');

n_sectors = params.n_sectors;
announce_runtime_step(runtime_ctx, 'Initialize Dynare runtime');
session = initialize_dynare_session(ModData, params, opts);
v = session.verbose;

%% Generate shocks (only when at least one simulation is requested)
run_any_simul = opts.run_firstorder_simul || opts.run_secondorder_simul || ...
                opts.run_pf_simul || opts.run_mit_shocks_simul;

if run_any_simul
    T_active = opts.simul_T;

    if ~isempty(opts.rng_seed), rng(opts.rng_seed); end
    rng_state = rng;

    if isfield(opts, 'external_shocks') && ~isempty(opts.external_shocks)
        shocks_active = opts.external_shocks;
        validate_shock_matrix(shocks_active, T_active, n_sectors, 'run_dynare_analysis');
    else
        shocks_active = mvnrnd(zeros([n_sectors,1]), params.Sigma_A, T_active);
    end

    Shocks = struct('data', shocks_active, 'T', T_active, ...
        'rng_state', rng_state, 'Sigma_A', params.Sigma_A);

    if v
        if isfield(opts, 'external_shocks') && ~isempty(opts.external_shocks)
            fprintf('Shared active shocks: %d x %d (external)\n', ...
                size(shocks_active, 1), size(shocks_active, 2));
        else
            fprintf('Shared active shocks: %d x %d (seed=%d)\n', ...
                size(shocks_active, 1), size(shocks_active, 2), rng_state.Seed);
        end
    end
else
    rng_state = rng;
    shocks_active = [];
    Shocks = struct('data', [], 'T', 0, 'rng_state', rng_state, 'Sigma_A', params.Sigma_A);
    if v, fprintf('IRF-only mode, skipping shock generation\n'); end
end

idx = session.idx;
policies_ss = session.policies_ss;
endostates_ss = session.endostates_ss;

Results = struct( ...
    'params', params, ...
    'C_ss', session.steady_state_aggregates.C_ss, ...
    'L_ss', session.steady_state_aggregates.L_ss, ...
    'GDP_ss', session.steady_state_aggregates.GDP_ss, ...
    'I_ss', session.steady_state_aggregates.I_ss, ...
    'K_ss', session.steady_state_aggregates.K_ss, ...
    'utility_intratemp_ss', session.steady_state_aggregates.utility_intratemp_ss, ...
    'steady_state_aggregates', session.steady_state_aggregates);

run_any_simul = opts.run_firstorder_simul || opts.run_secondorder_simul || ...
    opts.run_pf_simul || opts.run_mit_shocks_simul;
needs_theoretical_stats = true;
needs_1st_order_solution = opts.run_firstorder_simul || opts.run_firstorder_irs || needs_theoretical_stats;

have_1st_order_solution = false;
have_2nd_order_solution = false;
oo_1st = [];
M_1st = [];
options_1st = [];
oo_2nd = [];
M_2nd = [];
options_2nd = [];

if needs_1st_order_solution
    [oo_1st, M_1st, options_1st] = load_or_reuse_dynare_solution( ...
        1, session.dynare_folder, v, opts, runtime_ctx);
    have_1st_order_solution = true;
    Results.oo_1st = oo_1st;
    Results.M_1st = M_1st;
    Results.options_1st = options_1st;
    Results.steady_state = oo_1st.steady_state;
end

if needs_theoretical_stats
    TheoStats = compute_theoretical_statistics(oo_1st, M_1st, policies_ss, n_sectors, endostates_ss);
    Results.TheoStats = TheoStats;

    if v && ~isempty(fieldnames(TheoStats))
        fprintf('TheoStats: sig_GDP=%.4f sig_C=%.4f sig_I=%.4f sig_L=%.4f sig_K=%.4f\n', ...
            TheoStats.sigma_VA_agg, TheoStats.sigma_C_agg, TheoStats.sigma_I_agg, TheoStats.sigma_L_agg, TheoStats.sigma_K_agg);
    end
end

%% 2. First-order simulation
if opts.run_firstorder_simul
    simul_cfg = get_simulation_window_config(opts, false);
    announce_runtime_step(runtime_ctx, sprintf( ...
        'Run first-order simulation (burn_in=%d, T_active=%d, burn_out=%d)', ...
        simul_cfg.burn_in, simul_cfg.T_active, simul_cfg.burn_out));
    if v
        fprintf('\n1st-order simulation (burn_in=%d, T_active=%d, burn_out=%d, total=%d)...\n', ...
            simul_cfg.burn_in, simul_cfg.T_active, simul_cfg.burn_out, simul_cfg.T_total);
    end
    tic;

    shock_path_1st = build_perturbation_shock_path(shocks_active, simul_cfg, session.n_sectors);
    dynare_simul_1st_raw = simult_(M_1st, options_1st, oo_1st.steady_state, oo_1st.dr, shock_path_1st, 1);
    simul_block_1st = build_simulation_block_from_raw( ...
        dynare_simul_1st_raw, simul_cfg, idx, 'First-order simulation');
    SolData = extract_state_space(oo_1st, M_1st, n_sectors, policies_ss);
    Shocks.usage.FirstOrder = struct('start', 1, 'end', simul_cfg.T_active);

    variables_var = var(simul_block_1st.shocks_simul, 0, 2);
    SolData.shocks_sd = sqrt(var(shocks_active, 0, 1)).';
    SolData.states_sd = sqrt(variables_var(1:idx.n_states));
    SolData.policies_sd = sqrt(variables_var(idx.n_states+1:idx.n_dynare));

    simul_block_1st.summary_stats = build_simulation_summary( ...
        simul_block_1st, idx, policies_ss, n_sectors, endostates_ss, 'First-order');

    Results.SolData = SolData;
    Results.SimulFirstOrder = simul_block_1st;

    if v, fprintf('1st-order done (%.1fs)\n', toc); end
end

%% 3. Second-order simulation
if opts.run_secondorder_simul
    if ~have_2nd_order_solution
        [oo_2nd, M_2nd, options_2nd] = load_or_reuse_dynare_solution( ...
            2, session.dynare_folder, v, opts, runtime_ctx);
        have_2nd_order_solution = true;
        Results.oo_2nd = oo_2nd;
        Results.M_2nd = M_2nd;
        Results.options_2nd = options_2nd;
    end

    simul_cfg = get_simulation_window_config(opts, false);
    announce_runtime_step(runtime_ctx, sprintf( ...
        'Run second-order simulation (burn_in=%d, T_active=%d, burn_out=%d)', ...
        simul_cfg.burn_in, simul_cfg.T_active, simul_cfg.burn_out));
    if v
        fprintf('\n2nd-order simulation (burn_in=%d, T_active=%d, burn_out=%d, total=%d)...\n', ...
            simul_cfg.burn_in, simul_cfg.T_active, simul_cfg.burn_out, simul_cfg.T_total);
    end
    tic;

    shock_path_2nd = build_perturbation_shock_path(shocks_active, simul_cfg, session.n_sectors);
    dynare_simul_2nd_raw = simult_(M_2nd, options_2nd, oo_2nd.steady_state, oo_2nd.dr, shock_path_2nd, 2);
    simul_block_2nd = build_simulation_block_from_raw( ...
        dynare_simul_2nd_raw, simul_cfg, idx, 'Second-order simulation');
    Shocks.usage.SecondOrder = struct('start', 1, 'end', simul_cfg.T_active);

    simul_block_2nd.summary_stats = build_simulation_summary( ...
        simul_block_2nd, idx, policies_ss, n_sectors, endostates_ss, 'Second-order');

    Results.SimulSecondOrder = simul_block_2nd;

    if v, fprintf('2nd-order done (%.1fs)\n', toc); end
end

%% 4. First-order IRFs
if opts.run_firstorder_irs
    announce_runtime_step(runtime_ctx, sprintf('Run first-order IRs (H=%d)', opts.ir_horizon));
    if v, fprintf('\n1st-order IRFs (H=%d)...\n', opts.ir_horizon); end
    tic;
    Results.IRSFirstOrder_raw = compute_perturbation_irfs( ...
        oo_1st, M_1st, options_1st, params, opts, 1);
    Results.ir_sector_indices = opts.sector_indices;
    if v, fprintf('1st-order IRFs done (%.1fs)\n', toc); end
end

%% 5. Second-order IRFs
if opts.run_secondorder_irs
    if ~have_2nd_order_solution
        [oo_2nd, M_2nd, options_2nd] = load_or_reuse_dynare_solution( ...
            2, session.dynare_folder, v, opts, runtime_ctx);
        have_2nd_order_solution = true;
        Results.oo_2nd = oo_2nd;
        Results.M_2nd = M_2nd;
        Results.options_2nd = options_2nd;
    end
    announce_runtime_step(runtime_ctx, sprintf('Run second-order IRs (H=%d)', opts.ir_horizon));
    if v, fprintf('\n2nd-order IRFs (H=%d)...\n', opts.ir_horizon); end
    tic;
    Results.IRSSecondOrder_raw = compute_perturbation_irfs( ...
        oo_2nd, M_2nd, options_2nd, params, opts, 2);
    if v, fprintf('2nd-order IRFs done (%.1fs)\n', toc); end
end

%% 6. Perfect foresight IRFs
if opts.run_pf_irs
    announce_runtime_step(runtime_ctx, sprintf('Run perfect-foresight IRs (H=%d)', opts.ir_horizon));
    if v, fprintf('\nPF IRFs (H=%d)...\n', opts.ir_horizon); end
    tic;
    Results.IRSPerfectForesight_raw = run_pf_irfs(session, params, opts);
    if v, fprintf('PF IRFs done (%.1fs)\n', toc); end
end

%% 7. Perfect foresight simulation
if opts.run_pf_simul
    simul_cfg = get_simulation_window_config(opts, false);
    announce_runtime_step(runtime_ctx, sprintf( ...
        'Run perfect-foresight simulation (T_active=%d, burn_in=%d, burn_out=%d)', ...
        simul_cfg.T_active, simul_cfg.burn_in, simul_cfg.burn_out));
    if v
        fprintf('\nPF simulation (T_active=%d, burn_in=%d, burn_out=%d, total=%d)...\n', ...
            simul_cfg.T_active, simul_cfg.burn_in, simul_cfg.burn_out, simul_cfg.T_total);
    end
    tic;
    [Results, Shocks] = run_pf_simulation(Results, Shocks, shocks_active, session, idx, ...
        policies_ss, endostates_ss, opts);
    if isfield(Results, 'SimulPerfectForesight') && ~isempty(Results.SimulPerfectForesight) && v
        fprintf('PF simulation done (%.1fs)\n', toc);
    end
end

%% 8. MIT shocks simulation
if opts.run_mit_shocks_simul
    simul_cfg = get_simulation_window_config(opts, true);
    announce_runtime_step(runtime_ctx, sprintf( ...
        'Run MIT shocks simulation (T_active=%d, burn_out=%d, H=%d)', ...
        simul_cfg.T_active, simul_cfg.burn_out, opts.mit_lookahead_horizon));
    if v
        fprintf('\nMIT shocks (T_active=%d, burn_out=%d, total=%d, H=%d)...\n', ...
            simul_cfg.T_active, simul_cfg.burn_out, simul_cfg.T_total, ...
            opts.mit_lookahead_horizon);
    end
    tic;
    [Results, Shocks] = run_mit_shocks_simulation(Results, Shocks, shocks_active, session, idx, ...
        policies_ss, endostates_ss, opts);
    if isfield(Results, 'SimulMITShocks') && ~isempty(Results.SimulMITShocks) && v
        fprintf('MIT shocks done (%.1fs)\n', toc);
    end
end

%% Store shocks
Results.Shocks = Shocks;
Results.rng_state = rng_state;

if v, fprintf('\nDynare analysis complete. Fields: {%s}\n', strjoin(fieldnames(Results)', ', ')); end

end


%% ==================== Local Helpers ====================

function opts = normalize_dynare_opts(opts)
defaults = default_dynare_runtime_opts();
default_fields = fieldnames(defaults);
for i = 1:numel(default_fields)
    field_name = default_fields{i};
    opts = set_default(opts, field_name, defaults.(field_name));
end
opts = promote_legacy_simulation_opts(opts);
end

function [oo, M, options] = load_or_reuse_dynare_solution(order, dynare_folder, verbose, opts, runtime_ctx)
cached_solution = get_cached_solution(opts, order);

if ~isempty(cached_solution)
    announce_runtime_step(runtime_ctx, sprintf('Reuse cached order=%d perturbation solution', order));
    oo = cached_solution.oo;
    M = cached_solution.M;
    options = cached_solution.options;
    return;
end

announce_runtime_step(runtime_ctx, sprintf('Solve order=%d perturbation solution', order));
[oo, M, options] = ensure_dynare_solution(order, dynare_folder, verbose);
end

function cached_solution = get_cached_solution(opts, order)
cached_solution = [];

if ~isfield(opts, 'solution_cache') || isempty(opts.solution_cache)
    return;
end

if order == 1
    field_name = 'first_order';
else
    field_name = 'second_order';
end

if ~isfield(opts.solution_cache, field_name)
    return;
end

candidate = opts.solution_cache.(field_name);
required_fields = {'oo', 'M', 'options'};
if all(isfield(candidate, required_fields))
    cached_solution = candidate;
end
end

function session = initialize_dynare_session(ModData, params, opts)
[dynare_folder, ~, ~] = fileparts(mfilename('fullpath'));
n_sectors = params.n_sectors;
idx = get_variable_indices(n_sectors);
policies_ss = ModData.policies_ss;
endostates_ss = ModData.endostates_ss;

session = struct();
session.dynare_folder = dynare_folder;
session.verbose = opts.verbose;
session.n_sectors = n_sectors;
session.idx = idx;
session.policies_ss = policies_ss;
session.endostates_ss = endostates_ss;
session.model_parameters = ModData.parameters;
session.steady_state_aggregates = collect_steady_state_aggregates(ModData);

setup_dynare_runtime(session, ModData, params, opts);
end

function IRSDeterm_all = run_pf_irfs(session, params, opts)
ir_simul_periods = opts.ir_horizon - 2;
generate_determ_irs_mod(session.dynare_folder, ir_simul_periods);
IRSDeterm_all = cell(numel(opts.sector_indices), 1);

for ii = 1:numel(opts.sector_indices)
    sector_idx = opts.sector_indices(ii);
    shocksim_0 = zeros([session.n_sectors, 1]);
    shocksim_0(sector_idx) = -params.IRshock;
    workspace_vars = struct();
    workspace_vars.shocksim_0 = shocksim_0;
    workspace_vars.shockssim_ir = zeros([opts.ir_horizon, session.n_sectors]);

    dynare_irf_paths = run_generated_dynare_job(session, 'determ_irs_generated', workspace_vars);
    IRSDeterm_all{ii} = normalize_simulation_matrix(dynare_irf_paths, session.idx.n_dynare, 'PF IRF simulation');
    if opts.verbose, fprintf('  Sector %d done\n', sector_idx); end
end
end

function [Results, Shocks] = run_pf_simulation(Results, Shocks, shocks_active, session, idx, ...
        policies_ss, endostates_ss, opts)
simul_cfg = get_simulation_window_config(opts, false);

Shocks.usage.PerfectForesight = struct('start', 1, 'end', simul_cfg.T_active);

shockssim_pf = build_pf_shock_path(shocks_active, simul_cfg, session.n_sectors);
simul_periods = size(shockssim_pf, 1) - 2;
generate_determ_simul_mod(session.dynare_folder, simul_periods);

try
    workspace_vars = struct('shockssim_determ', shockssim_pf);
    dynare_simul_pf_raw = normalize_simulation_matrix( ...
        run_generated_dynare_job(session, 'determ_simul_generated', workspace_vars), ...
        session.idx.n_dynare, 'PF simulation');
    simul_block_pf = build_simulation_block_from_raw( ...
        dynare_simul_pf_raw, simul_cfg, idx, 'PF simulation');
    simul_block_pf.summary_stats = build_simulation_summary( ...
        simul_block_pf, session.idx, policies_ss, session.n_sectors, endostates_ss, ...
        'Perfect foresight');

    Results.SimulPerfectForesight = simul_block_pf;
catch ME
    if opts.continue_on_failure
        if opts.verbose, fprintf('PF simulation failed: %s\n', ME.message); end
        Results.SimulPerfectForesight = [];
        Results.pf_simul_error = ME;
    else
        error('run_dynare_analysis:PFSimulationFailed', ...
            'PF simulation failed: %s', ME.message);
    end
end
end

function [Results, Shocks] = run_mit_shocks_simulation(Results, Shocks, shocks_active, session, idx, ...
        policies_ss, endostates_ss, opts)
simul_cfg = get_simulation_window_config(opts, true);
mit_solver_mode = char(opts.mit_solver_mode);

Shocks.usage.MITShocks = struct('start', 1, 'end', simul_cfg.T_active);

try
    switch mit_solver_mode
        case 'rolling'
            dynare_simul_mit_raw = run_mit_shocks_simulation_rolling( ...
                shocks_active, session, idx, endostates_ss, opts);
        case 'legacy'
            dynare_simul_mit_raw = run_mit_shocks_simulation_legacy( ...
                shocks_active, session, opts);
        otherwise
            error('run_dynare_analysis:UnsupportedMITSolverMode', ...
                'Unsupported MIT solver mode: %s', mit_solver_mode);
    end

    simul_block_mit = build_simulation_block_from_raw( ...
        dynare_simul_mit_raw, simul_cfg, idx, 'MIT simulation');
    simul_block_mit.summary_stats = build_simulation_summary( ...
        simul_block_mit, session.idx, policies_ss, session.n_sectors, endostates_ss, ...
        'MIT shocks');

    Results.SimulMITShocks = simul_block_mit;
catch ME
    if opts.continue_on_failure
        if opts.verbose, fprintf('MIT shocks failed: %s\n', ME.message); end
        Results.SimulMITShocks = [];
        Results.mit_simul_error = ME;
    else
        error('run_dynare_analysis:MITSimulationFailed', ...
            'MIT shocks simulation failed: %s', ME.message);
    end
end
end

function dynare_simul_mit_raw = run_mit_shocks_simulation_legacy(shocks_active, session, opts)
simul_cfg = get_simulation_window_config(opts, true);
generate_mit_shocks_mod(session.dynare_folder, simul_cfg.T_active, simul_cfg.burn_out, ...
    session.n_sectors, shocks_active);
dynare_simul_mit_raw = normalize_simulation_matrix( ...
    run_generated_dynare_job(session, 'mit_shocks_simul_generated', struct()), ...
    session.idx.n_dynare, 'MIT simulation');
end

function dynare_simul_mit_raw = run_mit_shocks_simulation_rolling(shocks_active, session, idx, ...
        endostates_ss, opts)
total_steps = opts.simul_T + opts.simul_burn_out;
rolling_steps = opts.simul_T;
stitched_path = zeros(session.idx.n_dynare, total_steps);
rolling_cache = initialize_mit_rolling_pf_cache();
a_curr = zeros(session.n_sectors, 1);
k_curr = endostates_ss(:);
previous_endo_simul = [];
start_step = 1;
rolling_run_tic = tic;

checkpoint_cfg = get_mit_rolling_checkpoint_config(opts);
checkpoint_signature = build_mit_rolling_checkpoint_signature( ...
    shocks_active, session, total_steps, rolling_steps, opts);
[checkpoint_state, checkpoint_loaded] = load_mit_rolling_checkpoint( ...
    checkpoint_cfg, checkpoint_signature, opts);

if checkpoint_loaded
    stitched_cols = size(checkpoint_state.stitched_path, 2);
    stitched_path(:, 1:stitched_cols) = checkpoint_state.stitched_path;

    if strcmp(checkpoint_state.status, 'complete')
        dynare_simul_mit_raw = stitched_path;
        return;
    end

    k_curr = checkpoint_state.k_curr;
    a_curr = checkpoint_state.a_curr;
    previous_endo_simul = checkpoint_state.previous_endo_simul;
    start_step = checkpoint_state.completed_step + 1;

    if opts.verbose
        fprintf('  Resuming MIT rolling checkpoint at step %d/%d\n', start_step, rolling_steps);
    end
end

for t = start_step:rolling_steps
    step_tic = tic;
    checkpoint_saved = false;

    if t <= opts.simul_T
        shock_curr = shocks_active(t, :).';
    else
        shock_curr = zeros(session.n_sectors, 1);
    end

    local_horizon = resolve_mit_rolling_horizon(opts);
    [rolling_cache, rolling_session] = get_mit_rolling_pf_session( ...
        rolling_cache, session, local_horizon);
    [oo_local, local_path] = solve_rolling_mit_step( ...
        rolling_session, idx, k_curr, a_curr, shock_curr, previous_endo_simul);
    stitched_path(:, t) = local_path(:, 2);

    if t == opts.simul_T
        stitched_path = append_mit_burnout_tail(stitched_path, local_path, total_steps, t);
        save_mit_rolling_checkpoint(checkpoint_cfg, create_mit_rolling_checkpoint( ...
            'complete', t, stitched_path, [], [], [], checkpoint_signature));
        checkpoint_saved = checkpoint_cfg.enabled;

        if opts.verbose
            print_mit_rolling_progress(t, rolling_steps, start_step, toc(step_tic), ...
                toc(rolling_run_tic), checkpoint_saved, total_steps - rolling_steps);
        end
        break;
    end

    k_curr = local_path(idx.k(1):idx.k(2), 2);
    a_curr = local_path(idx.a(1):idx.a(2), 2);
    previous_endo_simul = oo_local.endo_simul;

    if should_save_mit_rolling_checkpoint(checkpoint_cfg, t, rolling_steps)
        save_mit_rolling_checkpoint(checkpoint_cfg, create_mit_rolling_checkpoint( ...
            'in_progress', t, stitched_path(:, 1:t), k_curr, a_curr, ...
            previous_endo_simul, checkpoint_signature));
        checkpoint_saved = true;
    end

    if opts.verbose
        print_mit_rolling_progress(t, rolling_steps, start_step, toc(step_tic), ...
            toc(rolling_run_tic), checkpoint_saved, total_steps - rolling_steps);
    end
end

dynare_simul_mit_raw = stitched_path;
end

function stitched_path = append_mit_burnout_tail(stitched_path, local_path, total_steps, active_step)
tail_periods = total_steps - active_step;
if tail_periods <= 0
    return;
end

available_tail_cols = size(local_path, 2) - 2;
if available_tail_cols < tail_periods
    error('run_dynare_analysis:MITRollingTailTooShort', ...
        ['Active-only rolling solve produced only %d future columns, but %d ' ...
         'burn-out periods are required. Increase the local horizon or adjust ' ...
         'the rolling horizon mode.'], available_tail_cols, tail_periods);
end

stitched_path(:, active_step + 1:total_steps) = local_path(:, 3:2 + tail_periods);
end

function rolling_cache = initialize_mit_rolling_pf_cache()
rolling_cache = struct();
rolling_cache.horizons = [];
rolling_cache.sessions = {};
end

function local_horizon = resolve_mit_rolling_horizon(opts)
local_horizon = opts.mit_lookahead_horizon + 3;
end

function [rolling_cache, rolling_session] = get_mit_rolling_pf_session( ...
        rolling_cache, session, lookahead_horizon)
match_idx = find(rolling_cache.horizons == lookahead_horizon, 1);
if ~isempty(match_idx)
    rolling_session = rolling_cache.sessions{match_idx};
    return;
end

rolling_session = initialize_mit_rolling_pf_session(session, lookahead_horizon);
rolling_cache.horizons(end + 1) = lookahead_horizon;
rolling_cache.sessions{end + 1} = rolling_session;
end

function rolling_session = initialize_mit_rolling_pf_session(session, lookahead_horizon)
generate_mit_rolling_pf_mod(session.dynare_folder, lookahead_horizon, 1e-3);

workspace_vars = struct();
workspace_vars.a_init = zeros(session.n_sectors, 1);
workspace_vars.k_init = session.endostates_ss(:);
workspace_vars.shockssim_mit_pf = zeros(lookahead_horizon, session.n_sectors);
assign_base_workspace_vars(workspace_vars);
evalc('run_dynare_mod(session.dynare_folder, ''mit_rolling_pf_generated'');');

rolling_session = struct();
rolling_session.lookahead_horizon = lookahead_horizon;
rolling_session.zero_exo = workspace_vars.shockssim_mit_pf;
rolling_session.oo_template = read_base_workspace_var('oo_');
rolling_session.M = read_base_workspace_var('M_');
rolling_session.options = read_base_workspace_var('options_');
end

function [oo_local, local_path] = solve_rolling_mit_step(rolling_session, idx, k_curr, a_curr, ...
        shock_curr, previous_endo_simul)
oo_local = rolling_session.oo_template;
oo_local.steady_state(idx.k(1):idx.k(2)) = k_curr;
oo_local.steady_state(idx.a(1):idx.a(2)) = a_curr;

options_local = rolling_session.options;
oo_local = perfect_foresight_setup(rolling_session.M, options_local, oo_local);
oo_local.exo_simul = build_mit_local_exo_path(rolling_session.zero_exo, shock_curr);

if ~isempty(previous_endo_simul)
    oo_local.endo_simul = apply_rolling_warm_start(previous_endo_simul, oo_local.endo_simul);
end

options_local.dynatol.f = 1e-3;
options_local.endogenous_terminal_period = true;
[~, oo_local, ~] = evalc('perfect_foresight_solver(rolling_session.M, options_local, oo_local);');
local_path = extract_transition_path(oo_local.endo_simul, ...
    rolling_session.lookahead_horizon, 'MIT rolling local solve');
end

function local_exo = build_mit_local_exo_path(zero_exo, shock_curr)
local_exo = zero_exo;
if size(local_exo, 1) < 2
    error('run_dynare_analysis:InvalidMITLocalExoPath', ...
        'MIT local exogenous path must have at least 2 rows.');
end
local_exo(2, :) = shock_curr(:).';
end

function warm_start = apply_rolling_warm_start(previous_endo_simul, template_endo_simul)
warm_start = template_endo_simul;

if ~isequal(size(previous_endo_simul), size(template_endo_simul))
    return;
end

n_cols = size(template_endo_simul, 2);
if n_cols <= 2
    return;
end

warm_start(:, 2:end-1) = previous_endo_simul(:, 3:end);
warm_start(:, 1) = template_endo_simul(:, 1);
warm_start(:, end) = template_endo_simul(:, end);
end

function print_mit_rolling_progress(step_idx, total_steps, start_step, step_elapsed_seconds, ...
        total_elapsed_seconds, checkpoint_saved, burnout_periods)
completed_steps_this_run = max(1, step_idx - start_step + 1);
avg_step_seconds = total_elapsed_seconds / completed_steps_this_run;
remaining_steps = max(0, total_steps - step_idx);
eta_seconds = avg_step_seconds * remaining_steps;
progress_pct = 100 * step_idx / total_steps;

line = sprintf(['  MIT rolling step %d/%d (%.1f%%) | step %s | avg %s | ' ...
    'elapsed %s | eta %s'], ...
    step_idx, total_steps, progress_pct, ...
    format_duration_seconds(step_elapsed_seconds), ...
    format_duration_seconds(avg_step_seconds), ...
    format_duration_seconds(total_elapsed_seconds), ...
    format_duration_seconds(eta_seconds));

if checkpoint_saved
    line = sprintf('%s | checkpoint', line);
end

if step_idx == total_steps && burnout_periods > 0
    line = sprintf('%s | burn-out appended (%d periods)', line, burnout_periods);
end

fprintf('%s\n', line);
end

function text = format_duration_seconds(seconds)
seconds = max(0, seconds);

if seconds < 60
    text = sprintf('%.1fs', seconds);
    return;
end

hours = floor(seconds / 3600);
minutes = floor(mod(seconds, 3600) / 60);
secs = mod(seconds, 60);

if hours > 0
    text = sprintf('%dh%02dm%04.1fs', hours, minutes, secs);
else
    text = sprintf('%dm%04.1fs', minutes, secs);
end
end

function checkpoint_cfg = get_mit_rolling_checkpoint_config(opts)
checkpoint_cfg = struct('enabled', false, 'file', '', 'stride', 1);

if ~isfield(opts, 'mit_checkpoint_file') || isempty(opts.mit_checkpoint_file)
    return;
end

checkpoint_cfg.enabled = true;
checkpoint_cfg.file = char(opts.mit_checkpoint_file);

if isfield(opts, 'mit_checkpoint_stride') && ~isempty(opts.mit_checkpoint_stride)
    checkpoint_cfg.stride = max(1, round(opts.mit_checkpoint_stride));
end
end

function signature = build_mit_rolling_checkpoint_signature(shocks_active, session, total_steps, ...
        rolling_steps, opts)
signature = struct();
signature.version = 1;
signature.total_steps = total_steps;
signature.rolling_steps = rolling_steps;
signature.simul_T = opts.simul_T;
signature.simul_burn_out = opts.simul_burn_out;
signature.mit_lookahead_horizon = opts.mit_lookahead_horizon;
signature.n_dynare = session.idx.n_dynare;
signature.n_sectors = session.n_sectors;
signature.shocks_active = shocks_active;
signature.endostates_ss = session.endostates_ss(:);
signature.policies_ss = session.policies_ss(:);
source_file = which('run_dynare_analysis');
signature.source_file = source_file;
signature.source_datenum = get_mit_checkpoint_file_datenum(source_file);

generator_file = which('generate_mit_rolling_pf_mod');
signature.generator_file = generator_file;
signature.generator_datenum = get_mit_checkpoint_file_datenum(generator_file);

base_mod_file = resolve_base_mod_file(session.dynare_folder, opts.model_type);
signature.base_mod_file = base_mod_file;
signature.base_mod_datenum = get_mit_checkpoint_file_datenum(base_mod_file);
end

function base_mod_file = resolve_base_mod_file(dynare_folder, model_type)
switch char(model_type)
    case 'VA'
        base_mod_name = 'ProdNetRbc_base.mod';
    case {'GO', 'GO_noVA'}
        base_mod_name = 'ProdNetRbc_base_GO.mod';
    otherwise
        error('run_dynare_analysis:UnsupportedModelType', ...
            'Unsupported model_type: %s', char(model_type));
end

base_mod_file = fullfile(dynare_folder, base_mod_name);
end

function [checkpoint_state, checkpoint_loaded] = load_mit_rolling_checkpoint( ...
        checkpoint_cfg, checkpoint_signature, opts)
checkpoint_state = struct('status', '', 'completed_step', 0, 'stitched_path', [], ...
    'k_curr', [], 'a_curr', [], 'previous_endo_simul', []);
checkpoint_loaded = false;

if ~checkpoint_cfg.enabled || exist(checkpoint_cfg.file, 'file') ~= 2
    return;
end

try
    loaded = load(checkpoint_cfg.file, 'checkpoint');
catch ME
    if opts.verbose
        fprintf('  Ignoring unreadable MIT checkpoint: %s\n', ME.message);
    end
    return;
end

if ~isfield(loaded, 'checkpoint') || ~isstruct(loaded.checkpoint)
    return;
end

candidate = loaded.checkpoint;
if ~is_valid_mit_rolling_checkpoint(candidate, checkpoint_signature)
    if opts.verbose
        fprintf('  Ignoring stale MIT checkpoint: %s\n', checkpoint_cfg.file);
    end
    return;
end

checkpoint_state = candidate;
checkpoint_loaded = true;
end

function tf = is_valid_mit_rolling_checkpoint(checkpoint, checkpoint_signature)
tf = false;
required_fields = {'version', 'status', 'completed_step', 'stitched_path', ...
    'k_curr', 'a_curr', 'previous_endo_simul', 'signature'};

for i = 1:numel(required_fields)
    if ~isfield(checkpoint, required_fields{i})
        return;
    end
end

if ~isequal(checkpoint.signature, checkpoint_signature)
    return;
end

status = char(checkpoint.status);
if ~any(strcmp(status, {'in_progress', 'complete'}))
    return;
end

if ~isnumeric(checkpoint.completed_step) || ~isscalar(checkpoint.completed_step) || ...
        checkpoint.completed_step ~= floor(checkpoint.completed_step)
    return;
end

if size(checkpoint.stitched_path, 1) ~= checkpoint_signature.n_dynare
    return;
end

if strcmp(status, 'complete')
    if checkpoint.completed_step ~= checkpoint_signature.rolling_steps || ...
            size(checkpoint.stitched_path, 2) ~= checkpoint_signature.total_steps
        return;
    end
    tf = true;
    return;
end

if checkpoint.completed_step < 1 || checkpoint.completed_step >= checkpoint_signature.rolling_steps
    return;
end

if size(checkpoint.stitched_path, 2) ~= checkpoint.completed_step || ...
        numel(checkpoint.k_curr) ~= checkpoint_signature.n_sectors || ...
        numel(checkpoint.a_curr) ~= checkpoint_signature.n_sectors
    return;
end

if ~isempty(checkpoint.previous_endo_simul) && ...
        size(checkpoint.previous_endo_simul, 1) ~= checkpoint_signature.n_dynare
    return;
end

tf = true;
end

function tf = should_save_mit_rolling_checkpoint(checkpoint_cfg, step_idx, total_steps)
if ~checkpoint_cfg.enabled
    tf = false;
    return;
end

tf = step_idx == 1 || step_idx == total_steps || mod(step_idx, checkpoint_cfg.stride) == 0;
end

function file_datenum = get_mit_checkpoint_file_datenum(file_path)
file_datenum = [];
if isempty(file_path) || exist(file_path, 'file') ~= 2
    return;
end

file_info = dir(file_path);
if ~isempty(file_info)
    file_datenum = file_info.datenum;
end
end

function checkpoint = create_mit_rolling_checkpoint(status, completed_step, stitched_path, ...
        k_curr, a_curr, previous_endo_simul, signature)
checkpoint = struct();
checkpoint.version = 1;
checkpoint.status = status;
checkpoint.completed_step = completed_step;
checkpoint.stitched_path = stitched_path;
checkpoint.k_curr = k_curr;
checkpoint.a_curr = a_curr;
checkpoint.previous_endo_simul = previous_endo_simul;
checkpoint.signature = signature;
checkpoint.saved_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
end

function save_mit_rolling_checkpoint(checkpoint_cfg, checkpoint)
if ~checkpoint_cfg.enabled
    return;
end

checkpoint_dir = fileparts(checkpoint_cfg.file);
if ~isempty(checkpoint_dir) && exist(checkpoint_dir, 'dir') ~= 7
    mkdir(checkpoint_dir);
end

temp_file = [checkpoint_cfg.file, '.tmp'];
save(temp_file, 'checkpoint');

if exist(checkpoint_cfg.file, 'file') == 2
    delete(checkpoint_cfg.file);
end

movefile(temp_file, checkpoint_cfg.file);
end

function validate_shock_matrix(shocks_active, expected_rows, expected_cols, context)
assert(isnumeric(shocks_active) && ismatrix(shocks_active), ...
    'run_dynare_analysis:InvalidExternalShocks', ...
    '[%s] external_shocks must be a numeric matrix.', context);
assert(size(shocks_active, 1) == expected_rows && size(shocks_active, 2) == expected_cols, ...
    'run_dynare_analysis:ExternalShockSizeMismatch', ...
    '[%s] external_shocks must have size %d x %d.', ...
    context, expected_rows, expected_cols);
end

function IRS_all = compute_perturbation_irfs(oo_, M_, options_, params, opts, order)
    n_sectors = params.n_sectors;
    IRS_all = cell(numel(opts.sector_indices), 1);

    for ii = 1:numel(opts.sector_indices)
        sector_idx = opts.sector_indices(ii);

        ss_shocked = oo_.steady_state;
        ss_shocked(n_sectors + sector_idx) = -params.IRshock;

        shockssim_ir = zeros([opts.ir_horizon, n_sectors]);
        IRS_all{ii} = simult_(M_, options_, ss_shocked, oo_.dr, shockssim_ir, order);

        if opts.verbose, fprintf('  Sector %d done\n', sector_idx); end
    end
end

function simul = normalize_simulation_matrix(simul, expected_rows, context_label)
if ~(isnumeric(simul) && ismatrix(simul))
    error('run_dynare_analysis:InvalidSimulationOutput', ...
        '%s must be a numeric matrix.', context_label);
end

if size(simul, 1) == expected_rows
    return;
end

if size(simul, 2) == expected_rows
    simul = simul.';
    return;
end

error('run_dynare_analysis:SimulationRowMismatch', ...
    '%s has size %d x %d; expected one dimension to equal %d endogenous variables.', ...
    context_label, size(simul, 1), size(simul, 2), expected_rows);
end

function opts = promote_legacy_simulation_opts(opts)
if ~isfield(opts, 'simul_T')
    opts.simul_T = first_present(opts, ...
        {'simul_T_firstorder', 'simul_T_secondorder', 'simul_T_pf', 'simul_T_mit'}, 1000);
end
if ~isfield(opts, 'simul_burn_in')
    opts.simul_burn_in = first_present(opts, {'pf_burn_in'}, 100);
end
if ~isfield(opts, 'simul_burn_out')
    opts.simul_burn_out = first_present(opts, {'pf_burn_out', 'mit_burn_out'}, 100);
end
end

function value = first_present(s, field_names, fallback)
value = fallback;
for i = 1:numel(field_names)
    field_name = field_names{i};
    if isfield(s, field_name) && ~isempty(s.(field_name))
        value = s.(field_name);
        return;
    end
end
end

function aggregates = collect_steady_state_aggregates(ModData)
required_fields = {'C_ss', 'L_ss', 'GDP_ss', 'I_ss', 'K_ss', 'utility_intratemp_ss'};
aggregates = struct();

for i = 1:numel(required_fields)
    field_name = required_fields{i};
    assert(isfield(ModData, field_name), ...
        'run_dynare_analysis:MissingSteadyStateAggregate', ...
        'ModData is missing required steady-state field: %s', field_name);
    aggregates.(field_name) = ModData.(field_name);
end
end

function runtime_ctx = resolve_runtime_context(opts)
runtime_ctx = struct();
if isfield(opts, 'runtime_context') && isstruct(opts.runtime_context)
    runtime_ctx = opts.runtime_context;
end
end

function announce_runtime_step(runtime_ctx, message)
fprintf('%s %s\n', runtime_prefix(runtime_ctx), message);
end

function prefix = runtime_prefix(runtime_ctx)
prefix = '[Dynare]';

if isempty(runtime_ctx) || ~isstruct(runtime_ctx)
    return;
end

if isfield(runtime_ctx, 'stage') && strcmpi(runtime_ctx.stage, 'IRF')
    prefix = '[Dynare|IRF]';
    if isfield(runtime_ctx, 'shock_idx') && isfield(runtime_ctx, 'n_shocks') ...
            && isfield(runtime_ctx, 'shock_label')
        prefix = sprintf('[Dynare|IRF %d/%d|%s]', ...
            runtime_ctx.shock_idx, runtime_ctx.n_shocks, runtime_ctx.shock_label);
    end
elseif isfield(runtime_ctx, 'stage') && ~isempty(runtime_ctx.stage)
    prefix = sprintf('[Dynare|%s]', char(runtime_ctx.stage));
end
end

function simul_cfg = get_simulation_window_config(opts, is_mit)
simul_cfg = struct();
simul_cfg.T_active = opts.simul_T;
simul_cfg.burn_in = opts.simul_burn_in;
simul_cfg.burn_out = opts.simul_burn_out;
if is_mit
    simul_cfg.burn_in = 0;
end
simul_cfg.T_total = simul_cfg.burn_in + simul_cfg.T_active + simul_cfg.burn_out;
end

function shock_path = build_perturbation_shock_path(shocks_active, simul_cfg, n_sectors)
shock_path = zeros(simul_cfg.T_total, n_sectors);
shock_path(simul_cfg.burn_in + 1:simul_cfg.burn_in + simul_cfg.T_active, :) = shocks_active;
end

function shock_path = build_pf_shock_path(shocks_active, simul_cfg, n_sectors)
required_rows = simul_cfg.burn_in + simul_cfg.T_active + 1;
shock_path = zeros(max(simul_cfg.T_total, required_rows), n_sectors);
active_start = simul_cfg.burn_in + 2;
active_end = active_start + simul_cfg.T_active - 1;
shock_path(active_start:active_end, :) = shocks_active;
end

function simul_block = build_simulation_block_from_raw(simul_raw, simul_cfg, idx, context_label)
simul_path = extract_transition_path(simul_raw, simul_cfg.T_total, context_label);
simul_block = build_simulation_block(simul_path, simul_cfg, idx);
end

function simul_path = extract_transition_path(simul_raw, expected_total, context_label)
actual_cols = size(simul_raw, 2);
if actual_cols < expected_total
    error('run_dynare_analysis:SimulationColumnMismatch', ...
        '%s has %d columns, but expected at least %d periods.', ...
        context_label, actual_cols, expected_total);
end

if actual_cols == expected_total
    simul_path = simul_raw;
else
    simul_path = simul_raw(:, 2:expected_total+1);
end
end

function simul_block = build_simulation_block(simul_path, simul_cfg, idx)
burn_in = simul_cfg.burn_in;
T_active = simul_cfg.T_active;
burn_out = simul_cfg.burn_out;
expected_total = burn_in + T_active + burn_out;

assert(size(simul_path, 2) == expected_total, ...
    'run_dynare_analysis:SimulationWindowMismatch', ...
    'Simulation path has %d columns, but expected %d.', size(simul_path, 2), expected_total);

simul_block = struct();
simul_block.burnin_simul = slice_simulation_window(simul_path, 1, burn_in);
simul_block.shocks_simul = slice_simulation_window(simul_path, burn_in + 1, burn_in + T_active);
simul_block.burnout_simul = slice_simulation_window(simul_path, burn_in + T_active + 1, expected_total);
simul_block.variable_indices = idx;
simul_block.burn_in = burn_in;
simul_block.burn_out = burn_out;
simul_block.T_active = size(simul_block.shocks_simul, 2);
simul_block.T_total = simul_block.burn_in + simul_block.T_active + size(simul_block.burnout_simul, 2);
end

function window = slice_simulation_window(simul_path, start_col, end_col)
if end_col < start_col
    window = simul_path(:, 1:0);
else
    window = simul_path(:, start_col:end_col);
end
end

function summary_stats = build_simulation_summary(simul_block, idx, policies_ss, ...
        n_sectors, endostates_ss, label)
states_simul = simul_block.shocks_simul(1:idx.n_states, :);
policies_simul = simul_block.shocks_simul(idx.n_states+1:end, :);

summary_stats = struct();
summary_stats.sample_window = 'shocks_simul';
summary_stats.aggregate_definition = 'exact_logdev_to_deterministic_ss';
summary_stats.burn_in = simul_block.burn_in;
summary_stats.burn_out = simul_block.burn_out;
summary_stats.T_active = simul_block.T_active;
summary_stats.T_total = simul_block.T_total;
summary_stats.states_mean = mean(states_simul, 2);
summary_stats.states_std = std(states_simul, 0, 2);
summary_stats.policies_mean = mean(policies_simul, 2);
summary_stats.policies_std = std(policies_simul, 0, 2);
summary_stats.ModelStats = compute_model_statistics( ...
    simul_block.shocks_simul, idx, policies_ss, n_sectors, endostates_ss, ...
    label, 'shocks_simul');
end
