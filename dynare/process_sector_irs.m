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
    print_sector_cir_summary(sector_stats, availability);
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

function context = build_ir_context(ModData, ~, params)
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
Stats.cir_values_firstorder = zeros(n_analyzed, 1);
Stats.cir_values_secondorder = zeros(n_analyzed, 1);
Stats.cir_values_pf = zeros(n_analyzed, 1);
Stats.total_effect_signs_firstorder = zeros(n_analyzed, 1);
Stats.total_effect_signs_secondorder = zeros(n_analyzed, 1);
Stats.total_effect_signs_pf = zeros(n_analyzed, 1);
Stats.nonlinear_amp_pf_vs_firstorder = NaN(n_analyzed, 1);
Stats.nonlinear_effect_class_pf_vs_firstorder = NaN(n_analyzed, 1);
end

function [IRFEntry, sector_stats] = process_single_sector_irf( ...
        idx_pos, sector_idx, client_idx, DynareResults, params, context, availability, ~)

[irs_1st, sectoral_loglin] = maybe_process_irf( ...
    availability.has_1storder, DynareResults, 'IRSFirstOrder_raw', idx_pos, ...
    sector_idx, client_idx, params, context);
[irs_2nd, sectoral_2nd] = maybe_process_irf( ...
    availability.has_2ndorder, DynareResults, 'IRSSecondOrder_raw', idx_pos, ...
    sector_idx, client_idx, params, context);
[irs_pf, sectoral_determ] = maybe_process_irf( ...
    availability.has_determ, DynareResults, 'IRSPerfectForesight_raw', idx_pos, ...
    sector_idx, client_idx, params, context);

sector_stats = compute_sector_cir_stats(irs_1st, irs_2nd, irs_pf, context.R);

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
IRFEntry.cir = build_entry_cir_block(sector_stats);
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

function sector_stats = compute_sector_cir_stats(irs_1st, irs_2nd, irs_pf, R)
sector_stats = struct();
sector_stats.firstorder = compute_cir_stats_for_irf(irs_1st, R);
sector_stats.secondorder = compute_cir_stats_for_irf(irs_2nd, R);
sector_stats.perfect_foresight = compute_cir_stats_for_irf(irs_pf, R);
sector_stats.nonlinear_amplification_pf_vs_firstorder = safe_divide_scalar( ...
    sector_stats.perfect_foresight.cumulative_response, ...
    sector_stats.firstorder.cumulative_response);
sector_stats.nonlinear_effect_class_pf_vs_firstorder = classify_nonlinear_effect( ...
    sector_stats.nonlinear_amplification_pf_vs_firstorder);
end

function cir_stats = compute_cir_stats_for_irf(IRS, R)
cir_stats = struct('cumulative_response', 0, 'total_effect_sign', 0);

if isempty(IRS)
    return;
end

consumption_series = get_consumption_series(IRS, R);
cir_stats.cumulative_response = sum(consumption_series);
cir_stats.total_effect_sign = classify_total_effect(cir_stats.cumulative_response);
end

function consumption_series = get_consumption_series(IRS, R)
consumption_series = IRS(R.C_exp, :);
end

function total_effect_sign = classify_total_effect(cumulative_response)
if abs(cumulative_response) < 1e-12
    total_effect_sign = 0;
else
    total_effect_sign = sign(cumulative_response);
end
end

function ratio = safe_divide_scalar(numerator, denominator)
if isfinite(numerator) && isfinite(denominator) && abs(denominator) > 1e-12
    ratio = numerator / denominator;
else
    ratio = NaN;
end
end

function effect_class = classify_nonlinear_effect(ratio)
if ~isfinite(ratio)
    effect_class = NaN;
elseif abs(ratio - 1) < 1e-12
    effect_class = 0;
elseif ratio > 1
    effect_class = 1;
else
    effect_class = -1;
end
end

function Stats = assign_sector_stats(Stats, idx_pos, sector_stats)
Stats.cir_values_firstorder(idx_pos) = sector_stats.firstorder.cumulative_response;
Stats.cir_values_secondorder(idx_pos) = sector_stats.secondorder.cumulative_response;
Stats.cir_values_pf(idx_pos) = sector_stats.perfect_foresight.cumulative_response;

Stats.total_effect_signs_firstorder(idx_pos) = sector_stats.firstorder.total_effect_sign;
Stats.total_effect_signs_secondorder(idx_pos) = sector_stats.secondorder.total_effect_sign;
Stats.total_effect_signs_pf(idx_pos) = sector_stats.perfect_foresight.total_effect_sign;
Stats.nonlinear_amp_pf_vs_firstorder(idx_pos) = sector_stats.nonlinear_amplification_pf_vs_firstorder;
Stats.nonlinear_effect_class_pf_vs_firstorder(idx_pos) = sector_stats.nonlinear_effect_class_pf_vs_firstorder;
end

function cir = build_entry_cir_block(sector_stats)
cir = struct();
cir.response_variable = 'C_exp';
cir.cumulative_responses = struct( ...
    'first_order', sector_stats.firstorder.cumulative_response, ...
    'second_order', sector_stats.secondorder.cumulative_response, ...
    'perfect_foresight', sector_stats.perfect_foresight.cumulative_response);
cir.total_effect_signs = struct( ...
    'first_order', sector_stats.firstorder.total_effect_sign, ...
    'second_order', sector_stats.secondorder.total_effect_sign, ...
    'perfect_foresight', sector_stats.perfect_foresight.total_effect_sign);
cir.nonlinear_amplification = struct( ...
    'pf_vs_first_order', sector_stats.nonlinear_amplification_pf_vs_firstorder);
cir.nonlinear_effect_class = struct( ...
    'pf_vs_first_order', sector_stats.nonlinear_effect_class_pf_vs_firstorder);
end

function print_sector_cir_summary(sector_stats, availability)
if availability.has_2ndorder
    fprintf('  CIR: 1st=%.4f (%s), 2nd=%.4f (%s), PF=%.4f (%s), amp(PF/1st)=%s\n', ...
        sector_stats.firstorder.cumulative_response, format_effect_sign(sector_stats.firstorder.total_effect_sign), ...
        sector_stats.secondorder.cumulative_response, format_effect_sign(sector_stats.secondorder.total_effect_sign), ...
        sector_stats.perfect_foresight.cumulative_response, format_effect_sign(sector_stats.perfect_foresight.total_effect_sign), ...
        format_ratio(sector_stats.nonlinear_amplification_pf_vs_firstorder));
else
    fprintf('  CIR=%.4f (1st, %s), %.4f (PF, %s), amp(PF/1st)=%s\n', ...
        sector_stats.firstorder.cumulative_response, format_effect_sign(sector_stats.firstorder.total_effect_sign), ...
        sector_stats.perfect_foresight.cumulative_response, format_effect_sign(sector_stats.perfect_foresight.total_effect_sign), ...
        format_ratio(sector_stats.nonlinear_amplification_pf_vs_firstorder));
end
end

function label = format_effect_sign(total_effect_sign)
if total_effect_sign > 0
    label = 'positive';
elseif total_effect_sign < 0
    label = 'negative';
else
    label = 'zero';
end
end

function text = format_ratio(value)
if isfinite(value)
    text = sprintf('%.3f', value);
else
    text = 'n/a';
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
    fprintf('Avg CIR:            %.4f        %.4f        %.4f\n', ...
        mean(Stats.cir_values_firstorder), mean(Stats.cir_values_secondorder), mean(Stats.cir_values_pf));
    fprintf('Positive sectors:   %d/%d          %d/%d          %d/%d\n', ...
        sum(Stats.total_effect_signs_firstorder > 0), n_analyzed, ...
        sum(Stats.total_effect_signs_secondorder > 0), n_analyzed, ...
        sum(Stats.total_effect_signs_pf > 0), n_analyzed);
    fprintf('Negative sectors:   %d/%d          %d/%d          %d/%d\n', ...
        sum(Stats.total_effect_signs_firstorder < 0), n_analyzed, ...
        sum(Stats.total_effect_signs_secondorder < 0), n_analyzed, ...
        sum(Stats.total_effect_signs_pf < 0), n_analyzed);
else
    fprintf('                    First-Order   Perfect Foresight\n');
    fprintf('Avg CIR:            %.4f        %.4f\n', ...
        mean(Stats.cir_values_firstorder), mean(Stats.cir_values_pf));
    fprintf('Positive sectors:   %d/%d          %d/%d\n', ...
        sum(Stats.total_effect_signs_firstorder > 0), n_analyzed, ...
        sum(Stats.total_effect_signs_pf > 0), n_analyzed);
    fprintf('Negative sectors:   %d/%d          %d/%d\n', ...
        sum(Stats.total_effect_signs_firstorder < 0), n_analyzed, ...
        sum(Stats.total_effect_signs_pf < 0), n_analyzed);
end
fprintf('Avg nonlinear amplification (PF/1st CIR): %s\n', ...
    format_ratio(finite_mean(Stats.nonlinear_amp_pf_vs_firstorder)));
fprintf('Amplified / attenuated sectors: %d / %d\n', ...
    sum(Stats.nonlinear_effect_class_pf_vs_firstorder > 0), ...
    sum(Stats.nonlinear_effect_class_pf_vs_firstorder < 0));
end

function value = finite_mean(values)
finite_values = values(isfinite(values));
if isempty(finite_values)
    value = NaN;
else
    value = mean(finite_values);
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
        'aggregate_perfect_foresight', {}, 'cir', {});
    return;
end

entries = vertcat(IRFs{:});
end

function summary_stats = build_irf_summary_stats(Stats)
summary_stats = struct();
summary_stats.response_variable = 'C_exp';
summary_stats.cumulative_responses = struct( ...
    'first_order', Stats.cir_values_firstorder(:).', ...
    'second_order', Stats.cir_values_secondorder(:).', ...
    'perfect_foresight', Stats.cir_values_pf(:).');
summary_stats.total_effect_signs = struct( ...
    'first_order', Stats.total_effect_signs_firstorder(:).', ...
    'second_order', Stats.total_effect_signs_secondorder(:).', ...
    'perfect_foresight', Stats.total_effect_signs_pf(:).');
summary_stats.nonlinear_amplification = struct( ...
    'pf_vs_first_order', Stats.nonlinear_amp_pf_vs_firstorder(:).');
summary_stats.nonlinear_effect_class = struct( ...
    'pf_vs_first_order', Stats.nonlinear_effect_class_pf_vs_firstorder(:).');
end
