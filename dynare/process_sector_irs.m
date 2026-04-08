function Results = process_sector_irs(DynareResults, params, ModData, labels, opts)
% PROCESS_SECTOR_IRS Process IRFs for specified sectors
%
% This is the unified function for processing impulse responses. It handles
% both single/few sectors and all-sectors analysis with a single code path.
%
% INPUTS:
%   DynareResults  - Structure from run_dynare_analysis containing IRF data
%   params         - Model parameters structure
%   ModData        - Steady state data structure
%   labels         - Labels struct with fields:
%                    - sector_indices, sector_labels, client_indices, client_labels
%   opts           - Options structure:
%                    - plot_graphs: Generate plots (default: true for <5 sectors)
%                    - save_graphs: Save plots to files (default: false)
%                    - save_intermediate: Save intermediate results (default: false)
%                    - save_interval: How often to save (default: 5)
%                    - exp_paths: Experiment paths structure from setup_experiment_folder
%                    - save_label: Label for files (default: '')
%                    - range_padding: Y-axis padding for plots (default: 0.1)
%                    - ir_plot_length: Number of periods for plots (default: 60)
%                    - shock_description: Description of shock for plot titles (default: '')
%                    - shock_sign: +1 for positive shock, -1 for negative (default: auto-detect)
%
% OUTPUTS:
%   Results - Canonical per-shock IR artifact with:
%             - metadata
%             - entries
%             - summary_stats
%             - labels

%% Input validation
validate_params(params, {'n_sectors'}, 'process_sector_irs');

%% Set defaults
if nargin < 5
    opts = struct();
end

sector_indices = labels.sector_indices;
n_analyzed = numel(sector_indices);
n_sectors = params.n_sectors;

validate_sector_indices(sector_indices, n_sectors, 'process_sector_irs');

opts = normalize_ir_opts(opts, n_analyzed);
context = build_ir_context(ModData, DynareResults, params);
availability = detect_irf_availability(DynareResults);

%% Initialize storage
fprintf('=== PROCESSING IRFs FOR %d SECTORS ===\n', n_analyzed);
if ~isempty(opts.shock_description)
    fprintf('    Shock: %s\n', opts.shock_description);
end

IRFs = cell(n_analyzed, 1);
Stats = initialize_irf_stats(n_analyzed, sector_indices);

%% Process each sector
for idx_pos = 1:n_analyzed
    sector_idx = sector_indices(idx_pos);
    client_idx = labels.client_indices(idx_pos);

    fprintf('Processing sector %d (%d/%d)\n', sector_idx, idx_pos, n_analyzed);

    [IRFs{idx_pos}, sector_stats] = process_single_sector_irf( ...
        idx_pos, sector_idx, client_idx, DynareResults, params, context, availability, opts);
    Stats = assign_sector_stats(Stats, idx_pos, sector_stats);
    print_sector_peak_summary(sector_stats, availability);
    maybe_save_intermediate_results(opts, idx_pos, IRFs, Stats);
end

%% Plot graphs if requested
if opts.plot_graphs
    plot_irf_results(IRFs, labels, sector_indices, opts, availability);
end

%% Print summary
print_irf_statistics_summary(Stats, availability, n_analyzed, opts);

%% Build output structure
Results = struct();
Results.metadata = build_irf_artifact_metadata(opts, labels, availability, IRFs);
Results.entries = finalize_irf_entries(IRFs);
Results.summary_stats = build_irf_summary_stats(Stats);
Results.labels = labels;
Results.shock_description = opts.shock_description;
if isfield(opts, 'shock_config')
    Results.shock_config = opts.shock_config;
    Results.metadata.shock_config = opts.shock_config;
end

end

function opts = normalize_ir_opts(opts, n_analyzed)
default_plot = (n_analyzed <= 5);

opts = set_default(opts, 'plot_graphs', default_plot);
opts = set_default(opts, 'save_graphs', false);
opts = set_default(opts, 'save_intermediate', false);
opts = set_default(opts, 'save_interval', 5);
opts = set_default(opts, 'exp_paths', struct('temp', 'output', 'figures', 'output'));
opts = set_default(opts, 'save_label', '');
opts = set_default(opts, 'range_padding', 0.1);
opts = set_default(opts, 'ir_plot_length', 60);
opts = set_default(opts, 'shock_description', '');
opts = set_default(opts, 'shock_sign', []);
end

function context = build_ir_context(ModData, DynareResults, params)
context = struct();
context.policies_ss = ModData.policies_ss;
context.k_ss = ModData.endostates_ss;
context.n_sectors = params.n_sectors;
context.R = get_ir_row_map();
end

function availability = detect_irf_availability(DynareResults)
availability = struct();
availability.has_1storder = isfield(DynareResults, 'IRSFirstOrder_raw') && ~isempty(DynareResults.IRSFirstOrder_raw);
availability.has_2ndorder = isfield(DynareResults, 'IRSSecondOrder_raw') && ~isempty(DynareResults.IRSSecondOrder_raw);
availability.has_determ = isfield(DynareResults, 'IRSPerfectForesight_raw') && ~isempty(DynareResults.IRSPerfectForesight_raw);
end

function Stats = initialize_irf_stats(n_analyzed, sector_indices)
Stats = struct();
Stats.sector_indices = sector_indices;
Stats.peak_values_firstorder = zeros(n_analyzed, 1);
Stats.peak_values_secondorder = zeros(n_analyzed, 1);
Stats.peak_values_pf = zeros(n_analyzed, 1);
Stats.peak_periods_firstorder = zeros(n_analyzed, 1);
Stats.peak_periods_secondorder = zeros(n_analyzed, 1);
Stats.peak_periods_pf = zeros(n_analyzed, 1);
Stats.half_lives_firstorder = zeros(n_analyzed, 1);
Stats.half_lives_secondorder = zeros(n_analyzed, 1);
Stats.half_lives_pf = zeros(n_analyzed, 1);
Stats.amplifications = zeros(n_analyzed, 1);
Stats.amplifications_2nd = zeros(n_analyzed, 1);
Stats.amplifications_rel = zeros(n_analyzed, 1);
end

function [IRFEntry, sector_stats] = process_single_sector_irf( ...
        idx_pos, sector_idx, client_idx, DynareResults, params, context, availability, opts)

[irs_1st, sectoral_loglin] = maybe_process_irf( ...
    availability.has_1storder, DynareResults, 'IRSFirstOrder_raw', idx_pos, ...
    sector_idx, client_idx, params, context);
[irs_2nd, sectoral_2nd] = maybe_process_irf( ...
    availability.has_2ndorder, DynareResults, 'IRSSecondOrder_raw', idx_pos, ...
    sector_idx, client_idx, params, context);
[irs_pf, sectoral_determ] = maybe_process_irf( ...
    availability.has_determ, DynareResults, 'IRSPerfectForesight_raw', idx_pos, ...
    sector_idx, client_idx, params, context);

shock_sign = resolve_shock_sign(opts.shock_sign, irs_1st, irs_2nd, irs_pf, context.R);
sector_stats = compute_sector_peak_stats(irs_1st, irs_2nd, irs_pf, shock_sign, context.R);

IRFEntry = struct();
IRFEntry.sector_idx = sector_idx;
IRFEntry.client_idx = client_idx;
IRFEntry.first_order = irs_1st;
IRFEntry.second_order = irs_2nd;
IRFEntry.perfect_foresight = irs_pf;
IRFEntry.sectoral_loglin = sectoral_loglin;
IRFEntry.sectoral_secondorder = sectoral_2nd;
IRFEntry.sectoral_determ = sectoral_determ;
IRFEntry.aggregate_first_order = build_aggregate_irf_block(irs_1st, context);
IRFEntry.aggregate_second_order = build_aggregate_irf_block(irs_2nd, context);
IRFEntry.aggregate_perfect_foresight = build_aggregate_irf_block(irs_pf, context);
end

function [IRS, sectoral_data] = maybe_process_irf( ...
        is_available, DynareResults, field_name, idx_pos, sector_idx, client_idx, params, context)
IRS = [];
sectoral_data = [];

if ~is_available
    return;
end

dynare_simul = DynareResults.(field_name){idx_pos};
[IRS, sectoral_data] = process_ir_data(dynare_simul, sector_idx, client_idx, params, ...
    context.n_sectors, context.k_ss, context.policies_ss);
end

function AggregateIRF = build_aggregate_irf_block(IRS, context)
if isempty(IRS)
    AggregateIRF = [];
    return;
end

R = context.R;
AggregateIRF = struct( ...
    'C', IRS(R.C_exp, :), ...
    'I', IRS(R.I_exp, :), ...
    'GDP', IRS(R.GDP_exp, :), ...
    'utility_intratemp', IRS(R.utility_intratemp, :), ...
    'C_exp', IRS(R.C_exp, :), ...
    'I_exp', IRS(R.I_exp, :), ...
    'GDP_exp', IRS(R.GDP_exp, :), ...
    'L', IRS(R.L, :), ...
    'K', IRS(R.K, :));
end

function shock_sign = resolve_shock_sign(explicit_sign, irs_1st, irs_2nd, irs_pf, R)
if ~isempty(explicit_sign)
    shock_sign = explicit_sign;
    return;
end

reference_irf = [];
if ~isempty(irs_1st)
    reference_irf = irs_1st;
elseif ~isempty(irs_2nd)
    reference_irf = irs_2nd;
elseif ~isempty(irs_pf)
    reference_irf = irs_pf;
end

if isempty(reference_irf)
    shock_sign = -1;
    return;
end

if reference_irf(R.A, 1) > 1
    shock_sign = 1;
else
    shock_sign = -1;
end
end

function sector_stats = compute_sector_peak_stats(irs_1st, irs_2nd, irs_pf, shock_sign, R)
sector_stats = struct();
sector_stats.firstorder = compute_peak_stats_for_irf(irs_1st, shock_sign, R);
sector_stats.secondorder = compute_peak_stats_for_irf(irs_2nd, shock_sign, R);
sector_stats.perfect_foresight = compute_peak_stats_for_irf(irs_pf, shock_sign, R);

sector_stats.amplification_abs = sector_stats.perfect_foresight.peak_value - sector_stats.firstorder.peak_value;
sector_stats.amplification_2nd = sector_stats.secondorder.peak_value - sector_stats.firstorder.peak_value;

if sector_stats.firstorder.peak_value > 1e-10
    sector_stats.amplification_rel = ...
        (sector_stats.perfect_foresight.peak_value / sector_stats.firstorder.peak_value - 1) * 100;
else
    sector_stats.amplification_rel = 0;
end

if shock_sign > 0
    sector_stats.sign_str = ' (+)';
else
    sector_stats.sign_str = ' (-)';
end
end

function peak_stats = compute_peak_stats_for_irf(IRS, shock_sign, R)
peak_stats = struct('peak_value', 0, 'peak_period', 0, 'half_life', 0);

if isempty(IRS)
    return;
end

T_stats = min(100, size(IRS, 2));
consumption_series = get_consumption_series(IRS, R, T_stats);
[peak_value, peak_period, half_life] = calculatePeaksAndHalfLives(shock_sign * consumption_series);

peak_stats.peak_value = abs(peak_value);
peak_stats.peak_period = peak_period;
peak_stats.half_life = half_life;
end

function consumption_series = get_consumption_series(IRS, R, T_stats)
if nargin < 3
    T_stats = size(IRS, 2);
end

consumption_series = IRS(R.C_exp, 1:T_stats);
end

function Stats = assign_sector_stats(Stats, idx_pos, sector_stats)
Stats.peak_values_firstorder(idx_pos) = sector_stats.firstorder.peak_value;
Stats.peak_values_secondorder(idx_pos) = sector_stats.secondorder.peak_value;
Stats.peak_values_pf(idx_pos) = sector_stats.perfect_foresight.peak_value;

Stats.peak_periods_firstorder(idx_pos) = sector_stats.firstorder.peak_period;
Stats.peak_periods_secondorder(idx_pos) = sector_stats.secondorder.peak_period;
Stats.peak_periods_pf(idx_pos) = sector_stats.perfect_foresight.peak_period;

Stats.half_lives_firstorder(idx_pos) = sector_stats.firstorder.half_life;
Stats.half_lives_secondorder(idx_pos) = sector_stats.secondorder.half_life;
Stats.half_lives_pf(idx_pos) = sector_stats.perfect_foresight.half_life;

Stats.amplifications(idx_pos) = sector_stats.amplification_abs;
Stats.amplifications_2nd(idx_pos) = sector_stats.amplification_2nd;
Stats.amplifications_rel(idx_pos) = sector_stats.amplification_rel;
end

function print_sector_peak_summary(sector_stats, availability)
if availability.has_2ndorder
    fprintf('  peak: 1st=%.4f, 2nd=%.4f, PF=%.4f%s | amplif: 2nd=%.4f, PF=%.2f%%\n', ...
        sector_stats.firstorder.peak_value, sector_stats.secondorder.peak_value, ...
        sector_stats.perfect_foresight.peak_value, sector_stats.sign_str, ...
        sector_stats.amplification_2nd, sector_stats.amplification_rel);
else
    fprintf('  peak=%.4f (1st), %.4f (PF)%s, amplification=%.2f%%, half-life=%d/%d\n', ...
        sector_stats.firstorder.peak_value, sector_stats.perfect_foresight.peak_value, ...
        sector_stats.sign_str, sector_stats.amplification_rel, ...
        sector_stats.firstorder.half_life, sector_stats.perfect_foresight.half_life);
end
end

function maybe_save_intermediate_results(opts, idx_pos, IRFs, Stats)
if ~(opts.save_intermediate && mod(idx_pos, opts.save_interval) == 0)
    return;
end

intermediate_file = fullfile(opts.exp_paths.temp, ['IRFs_Intermediate_' opts.save_label '.mat']);
save(intermediate_file, 'IRFs', 'Stats');
fprintf('  Saved intermediate results: %s\n', intermediate_file);
end

function plot_irf_results(IRFs, labels, sector_indices, opts, availability)
fprintf('Generating plots...\n');

ax = 0:(opts.ir_plot_length - 1);
n_analyzed = numel(sector_indices);

for idx_pos = 1:n_analyzed
    labels_single = build_single_sector_labels(labels, idx_pos, sector_indices(idx_pos));
    graph_opts = build_graph_opts(opts, sector_indices(idx_pos));

    fprintf('  Plotting sector %d (%d/%d)...\n', sector_indices(idx_pos), idx_pos, n_analyzed);
    call_graphirs_for_sector(IRFs{idx_pos}, labels_single, ax, opts, graph_opts, availability);
end
end

function labels_single = build_single_sector_labels(labels, idx_pos, sector_idx)
labels_single = struct();
labels_single.sector_indices = sector_idx;
labels_single.client_indices = labels.client_indices(idx_pos);
labels_single.sector_labels = labels.sector_labels(idx_pos);
labels_single.client_labels = labels.client_labels(idx_pos);

if isfield(labels, 'sector_labels_latex')
    labels_single.sector_labels_latex = labels.sector_labels_latex(idx_pos);
    labels_single.client_labels_latex = labels.client_labels_latex(idx_pos);
end

if isfield(labels, 'sector_labels_filename')
    labels_single.sector_labels_filename = labels.sector_labels_filename(idx_pos);
end
end

function graph_opts = build_graph_opts(opts, sector_idx)
graph_opts = struct();
graph_opts.figures_folder = opts.exp_paths.figures;
graph_opts.save_figures = opts.save_graphs;
graph_opts.shock_description = opts.shock_description;

if ~isempty(opts.save_label)
    graph_opts.save_label = sprintf('%s_sec%d', opts.save_label, sector_idx);
else
    graph_opts.save_label = sprintf('sec%d', sector_idx);
end
end

function call_graphirs_for_sector(IRFEntry, labels_single, ax, opts, graph_opts, availability)
irs_1st = {IRFEntry.first_order};
irs_2nd = {IRFEntry.second_order};
irs_pf = {IRFEntry.perfect_foresight};
aggregate_irs = struct( ...
    'series_1', [], ...
    'series_2', [], ...
    'series_3', []);

if availability.has_determ && availability.has_1storder && availability.has_2ndorder
    aggregate_irs.series_1 = {IRFEntry.aggregate_perfect_foresight};
    aggregate_irs.series_2 = {IRFEntry.aggregate_first_order};
    aggregate_irs.series_3 = {IRFEntry.aggregate_second_order};
    GraphIRs(irs_pf, irs_1st, irs_2nd, ax, opts.ir_plot_length, ...
        labels_single, opts.range_padding, graph_opts, {}, aggregate_irs);
elseif availability.has_determ && availability.has_1storder
    aggregate_irs.series_1 = {IRFEntry.aggregate_perfect_foresight};
    aggregate_irs.series_2 = {IRFEntry.aggregate_first_order};
    GraphIRs(irs_pf, irs_1st, [], ax, opts.ir_plot_length, ...
        labels_single, opts.range_padding, graph_opts, {}, aggregate_irs);
elseif availability.has_determ && availability.has_2ndorder
    aggregate_irs.series_1 = {IRFEntry.aggregate_perfect_foresight};
    aggregate_irs.series_2 = {IRFEntry.aggregate_second_order};
    GraphIRs(irs_pf, irs_2nd, [], ax, opts.ir_plot_length, ...
        labels_single, opts.range_padding, graph_opts, {'Perfect Foresight', 'Second-Order'}, aggregate_irs);
elseif availability.has_1storder && availability.has_2ndorder
    aggregate_irs.series_1 = {IRFEntry.aggregate_first_order};
    aggregate_irs.series_2 = {IRFEntry.aggregate_second_order};
    GraphIRs(irs_1st, irs_2nd, [], ax, opts.ir_plot_length, ...
        labels_single, opts.range_padding, graph_opts, {'First-Order', 'Second-Order'}, aggregate_irs);
elseif availability.has_determ
    aggregate_irs.series_1 = {IRFEntry.aggregate_perfect_foresight};
    GraphIRs(irs_pf, [], [], ax, opts.ir_plot_length, ...
        labels_single, opts.range_padding, graph_opts, {'Perfect Foresight'}, aggregate_irs);
elseif availability.has_1storder
    aggregate_irs.series_1 = {IRFEntry.aggregate_first_order};
    GraphIRs(irs_1st, [], [], ax, opts.ir_plot_length, ...
        labels_single, opts.range_padding, graph_opts, {'First-Order'}, aggregate_irs);
elseif availability.has_2ndorder
    aggregate_irs.series_1 = {IRFEntry.aggregate_second_order};
    GraphIRs(irs_2nd, [], [], ax, opts.ir_plot_length, ...
        labels_single, opts.range_padding, graph_opts, {'Second-Order'}, aggregate_irs);
end
end

function print_irf_statistics_summary(Stats, availability, n_analyzed, opts)
fprintf('\n=== IRF STATISTICS SUMMARY ===\n');
fprintf('Sectors analyzed: %d\n', n_analyzed);
if ~isempty(opts.shock_description)
    fprintf('Shock: %s\n', opts.shock_description);
end

if availability.has_2ndorder
    fprintf('                    First-Order   Second-Order  Perfect Foresight\n');
    fprintf('Avg peak value:     %.4f        %.4f        %.4f\n', ...
        mean(Stats.peak_values_firstorder), mean(Stats.peak_values_secondorder), mean(Stats.peak_values_pf));
    fprintf('Avg peak period:    %.1f           %.1f           %.1f\n', ...
        mean(Stats.peak_periods_firstorder), mean(Stats.peak_periods_secondorder), mean(Stats.peak_periods_pf));
    fprintf('Avg half-life:      %.1f           %.1f           %.1f\n', ...
        mean(Stats.half_lives_firstorder), mean(Stats.half_lives_secondorder), mean(Stats.half_lives_pf));
    fprintf('Amplification (2nd vs 1st): %.4f\n', mean(Stats.amplifications_2nd));
    fprintf('Amplification (PF vs 1st):  %.4f (%.1f%%)\n', mean(Stats.amplifications), mean(Stats.amplifications_rel));
else
    fprintf('                    First-Order   Perfect Foresight\n');
    fprintf('Avg peak value:     %.4f        %.4f\n', ...
        mean(Stats.peak_values_firstorder), mean(Stats.peak_values_pf));
    fprintf('Avg peak period:    %.1f           %.1f\n', ...
        mean(Stats.peak_periods_firstorder), mean(Stats.peak_periods_pf));
    fprintf('Avg half-life:      %.1f           %.1f\n', ...
        mean(Stats.half_lives_firstorder), mean(Stats.half_lives_pf));
    fprintf('Amplification (PF vs 1st): %.4f (%.1f%%)\n', mean(Stats.amplifications), mean(Stats.amplifications_rel));
end
end

function metadata = build_irf_artifact_metadata(opts, labels, availability, IRFs)
metadata = struct();
metadata.run_flags = struct( ...
    'has_1storder', logical(availability.has_1storder), ...
    'has_2ndorder', logical(availability.has_2ndorder), ...
    'has_pf', logical(availability.has_determ), ...
    'has_mit', false);
metadata.sector_indices = labels.sector_indices;
metadata.ir_horizon = resolve_ir_horizon(opts, IRFs);
metadata.shock_description = opts.shock_description;
metadata.n_sectors = numel(labels.sector_indices);
end

function ir_horizon = resolve_ir_horizon(opts, IRFs)
if isfield(opts, 'ir_horizon') && ~isempty(opts.ir_horizon)
    ir_horizon = opts.ir_horizon;
    return;
end

ir_horizon = 0;
for i = 1:numel(IRFs)
    entry = IRFs{i};
    candidate_fields = {'first_order', 'second_order', 'perfect_foresight'};
    for j = 1:numel(candidate_fields)
        field_name = candidate_fields{j};
        if isfield(entry, field_name) && ~isempty(entry.(field_name))
            ir_horizon = size(entry.(field_name), 2);
            return;
        end
    end
end
end

function entries = finalize_irf_entries(IRFs)
if isempty(IRFs)
    entries = struct('sector_idx', {}, 'client_idx', {}, 'first_order', {}, ...
        'second_order', {}, 'perfect_foresight', {}, 'sectoral_loglin', {}, ...
        'sectoral_secondorder', {}, 'sectoral_determ', {}, ...
        'aggregate_first_order', {}, 'aggregate_second_order', {}, ...
        'aggregate_perfect_foresight', {});
    return;
end

entries = vertcat(IRFs{:});
end

function summary_stats = build_irf_summary_stats(Stats)
summary_stats = struct();
summary_stats.peaks = struct( ...
    'first_order', Stats.peak_values_firstorder(:).', ...
    'second_order', Stats.peak_values_secondorder(:).', ...
    'perfect_foresight', Stats.peak_values_pf(:).');
summary_stats.peak_periods = struct( ...
    'first_order', Stats.peak_periods_firstorder(:).', ...
    'second_order', Stats.peak_periods_secondorder(:).', ...
    'perfect_foresight', Stats.peak_periods_pf(:).');
summary_stats.half_lives = struct( ...
    'first_order', Stats.half_lives_firstorder(:).', ...
    'second_order', Stats.half_lives_secondorder(:).', ...
    'perfect_foresight', Stats.half_lives_pf(:).');
summary_stats.amplifications = struct( ...
    'abs', Stats.amplifications(:).', ...
    'second_order', Stats.amplifications_2nd(:).', ...
    'rel', Stats.amplifications_rel(:).');
end
