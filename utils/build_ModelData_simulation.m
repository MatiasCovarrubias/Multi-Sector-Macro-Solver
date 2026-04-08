function [ModelData_simulation, flags] = build_ModelData_simulation( ...
        runtime_results, params, save_label, exp_paths, extra_artifacts)

    if nargin < 5 || isempty(extra_artifacts)
        extra_artifacts = struct();
    end

    ModelData_simulation = struct();
    ModelData_simulation.metadata.save_label = save_label;
    ModelData_simulation.metadata.exp_paths = exp_paths;

    n = params.n_sectors;
    idx = get_variable_indices(n);
    steady_state_aggregates = resolve_steady_state_aggregates(runtime_results);

    if isfield(runtime_results, 'Shocks')
        ModelData_simulation.Shocks = runtime_results.Shocks;
    end
    if isfield(runtime_results, 'rng_state')
        ModelData_simulation.rng_state = runtime_results.rng_state;
    end
    ModelData_simulation.Shared = build_shared_simulation_views(runtime_results);

    %% Process each simulation method
    simul_defs = {
        'SimulFirstOrder',       'FirstOrder',       'ModelStats';
        'SimulSecondOrder',      'SecondOrder',      'ModelStats2nd';
        'SimulPerfectForesight', 'PerfectForesight', 'ModelStatsPF';
        'SimulMITShocks',        'MITShocks',        'ModelStatsMIT';
    };

    flags = struct('has_1storder', false, 'has_2ndorder', false, ...
                   'has_pf', false, 'has_mit', false);
    flag_names = {'has_1storder', 'has_2ndorder', 'has_pf', 'has_mit'};

    for k = 1:size(simul_defs, 1)
        field_src  = simul_defs{k, 1};
        field_dst  = simul_defs{k, 2};
        stats_key  = simul_defs{k, 3};

        if ~isfield(runtime_results, field_src) || isempty(runtime_results.(field_src))
            continue;
        end

        flags.(flag_names{k}) = true;
        ModelData_simulation.(field_dst) = build_simulation_artifact( ...
            runtime_results.(field_src), runtime_results, extra_artifacts, ...
            field_dst, stats_key, idx, steady_state_aggregates);

    end

    ModelData_simulation.metadata.run_flags = flags;
    ModelData_simulation.metadata.has_irfs = false;
end

function simul_block = build_simulation_artifact( ...
        simul_input, runtime_results, extra_artifacts, field_dst, stats_key, idx, steady_state_aggregates)
simul_block = normalize_simulation_block_input(simul_input, idx, field_dst);
simul_block = attach_aggregate_series(simul_block, idx, steady_state_aggregates);
simul_block = attach_summary_stats(simul_block, runtime_results, stats_key, idx, field_dst);
simul_block = merge_extra_artifacts(simul_block, extra_artifacts, field_dst);
end

function simul_block = normalize_simulation_block_input(simul_input, idx, field_dst)
if isstruct(simul_input)
    simul_block = simul_input;
else
    simul_block = struct();
    simul_block.burnin_simul = simul_input(:, 1:0);
    simul_block.shocks_simul = simul_input;
    simul_block.burnout_simul = simul_input(:, 1:0);
    simul_block.burn_in = 0;
    simul_block.burn_out = 0;
    simul_block.T_active = size(simul_input, 2);
    simul_block.T_total = size(simul_input, 2);
end

simul_block.variable_indices = idx;

required_fields = {'burnin_simul', 'shocks_simul', 'burnout_simul', ...
    'variable_indices', 'burn_in', 'burn_out', 'T_active', 'T_total'};
for i = 1:numel(required_fields)
    field_name = required_fields{i};
    assert(isfield(simul_block, field_name), ...
        'build_ModelData_simulation:MissingSimulationField', ...
        '%s is missing required field: %s', field_dst, field_name);
end
end

function simul_block = attach_summary_stats(simul_block, runtime_results, stats_key, idx, field_dst)
if isfield(simul_block, 'summary_stats') && isstruct(simul_block.summary_stats) && ...
        ~isempty(fieldnames(simul_block.summary_stats))
    return;
end

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

if ~isempty(stats_key) && isfield(runtime_results, stats_key) && isstruct(runtime_results.(stats_key))
    summary_stats.ModelStats = runtime_results.(stats_key);
end

assert(isfield(summary_stats, 'ModelStats') && isstruct(summary_stats.ModelStats) && ...
    ~isempty(fieldnames(summary_stats.ModelStats)), ...
    'build_ModelData_simulation:MissingModelStats', ...
    '%s is missing required ModelStats payload: %s', field_dst, stats_key);

simul_block.summary_stats = summary_stats;
end

function simul_block = merge_extra_artifacts(simul_block, extra_artifacts, field_dst)
if ~isstruct(extra_artifacts) || ~isfield(extra_artifacts, field_dst)
    return;
end

extra_fields = fieldnames(extra_artifacts.(field_dst));
for i = 1:numel(extra_fields)
    field_name = extra_fields{i};
    simul_block.(field_name) = extra_artifacts.(field_dst).(field_name);
end
end

function Shared = build_shared_simulation_views(runtime_results)
Shared = struct();

if isfield(runtime_results, 'SolData') && isstruct(runtime_results.SolData)
    Shared.Solution = struct();
    Shared.Solution.StateSpace = struct( ...
        'A', runtime_results.SolData.A, ...
        'B', runtime_results.SolData.B, ...
        'C', runtime_results.SolData.C, ...
        'D', runtime_results.SolData.D);
    Shared.Solution.indices = runtime_results.SolData.indices;

    if isfield(runtime_results, 'steady_state')
        Shared.Solution.steady_state = runtime_results.steady_state;
    end

    Shared.Statistics = struct();
    if isfield(runtime_results.SolData, 'shocks_sd')
        Shared.Statistics.shocks_sd = runtime_results.SolData.shocks_sd;
    end
    if isfield(runtime_results.SolData, 'states_sd')
        Shared.Statistics.states_sd = runtime_results.SolData.states_sd;
    end
    if isfield(runtime_results.SolData, 'policies_sd')
        Shared.Statistics.policies_sd = runtime_results.SolData.policies_sd;
    end
end

if isfield(runtime_results, 'TheoStats') && isstruct(runtime_results.TheoStats) && ...
        ~isempty(fieldnames(runtime_results.TheoStats))
    if ~isfield(Shared, 'Statistics') || isempty(Shared.Statistics)
        Shared.Statistics = struct();
    end
    Shared.Statistics.TheoStats = runtime_results.TheoStats;
end

if isfield(runtime_results, 'steady_state_aggregates') && isstruct(runtime_results.steady_state_aggregates) && ...
        ~isempty(fieldnames(runtime_results.steady_state_aggregates))
    Shared.SteadyStateAggregates = runtime_results.steady_state_aggregates;
end
end

function steady_state_aggregates = resolve_steady_state_aggregates(runtime_results)
if isfield(runtime_results, 'steady_state_aggregates') && ...
        isstruct(runtime_results.steady_state_aggregates) && ...
        ~isempty(fieldnames(runtime_results.steady_state_aggregates))
    steady_state_aggregates = runtime_results.steady_state_aggregates;
else
    required_fields = {'C_ss', 'L_ss', 'GDP_ss', 'I_ss', 'K_ss', 'utility_intratemp_ss'};
    for i = 1:numel(required_fields)
        field_name = required_fields{i};
        assert(isfield(runtime_results, field_name), ...
            'build_ModelData_simulation:MissingSteadyStateAggregate', ...
            'runtime_results is missing required steady-state aggregate: %s', field_name);
    end

    steady_state_aggregates = struct( ...
        'C_ss', runtime_results.C_ss, ...
        'L_ss', runtime_results.L_ss, ...
        'GDP_ss', runtime_results.GDP_ss, ...
        'I_ss', runtime_results.I_ss, ...
        'K_ss', runtime_results.K_ss, ...
        'utility_intratemp_ss', runtime_results.utility_intratemp_ss);
end

required_fields = {'C_ss', 'L_ss', 'GDP_ss', 'I_ss', 'K_ss', 'utility_intratemp_ss'};
for i = 1:numel(required_fields)
    field_name = required_fields{i};
    assert(isfield(steady_state_aggregates, field_name), ...
        'build_ModelData_simulation:MissingSteadyStateAggregate', ...
        'Missing required steady-state aggregate: %s', field_name);
end
end

function simul_block = attach_aggregate_series(simul_block, idx, steady_state_aggregates)
window_defs = {
    'burnin',  'burnin_simul';
    'active',  'shocks_simul';
    'burnout', 'burnout_simul';
};

aggregate_specs = {
    'C',                  idx.c_agg,             steady_state_aggregates.C_ss,                 'log_level';
    'L',                  idx.l_agg,             steady_state_aggregates.L_ss,                 'log_level';
    'GDP',                idx.gdp_agg,           steady_state_aggregates.GDP_ss,               'log_level';
    'I',                  idx.i_agg,             steady_state_aggregates.I_ss,                 'log_level';
    'K',                  idx.k_agg,             steady_state_aggregates.K_ss,                 'log_level';
    'utility_intratemp',  idx.utility_intratemp, steady_state_aggregates.utility_intratemp_ss, 'level';
};

aggregate_series = struct();
for i = 1:size(aggregate_specs, 1)
    field_name = aggregate_specs{i, 1};
    row_idx = aggregate_specs{i, 2};
    steady_state_value = aggregate_specs{i, 3};
    convention = aggregate_specs{i, 4};
    series = struct();

    for j = 1:size(window_defs, 1)
        window_name = window_defs{j, 1};
        block_name = window_defs{j, 2};
        series.(window_name) = normalize_aggregate_series( ...
            simul_block.(block_name)(row_idx, :), steady_state_value, convention, field_name);
    end

    series.full = [series.burnin, series.active, series.burnout];
    aggregate_series.(field_name) = series;
end

simul_block.aggregate_series = aggregate_series;
end

function series = normalize_aggregate_series(raw_series, steady_state_value, convention, field_name)
assert(isnumeric(raw_series) && isrow(raw_series), ...
    'build_ModelData_simulation:InvalidAggregateSeries', ...
    'Aggregate %s must be a numeric row vector.', field_name);

assert(isnumeric(steady_state_value) && isscalar(steady_state_value) && isfinite(steady_state_value), ...
    'build_ModelData_simulation:InvalidSteadyStateAggregate', ...
    'Steady-state value for aggregate %s must be a finite numeric scalar.', field_name);

switch convention
    case 'log_level'
        assert(steady_state_value > 0, ...
            'build_ModelData_simulation:NonpositiveSteadyStateAggregate', ...
            'Steady-state level for aggregate %s must be strictly positive.', field_name);
        series = raw_series - log(steady_state_value);
    case 'level'
        series = raw_series - steady_state_value;
    otherwise
        error('build_ModelData_simulation:UnsupportedAggregateConvention', ...
            'Unsupported aggregate convention for %s: %s', field_name, convention);
end
end
