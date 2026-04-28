function Diagnostics = build_nonlinearity_diagnostics(~, AllShockResults, params, config, ModData)
% BUILD_NONLINEARITY_DIAGNOSTICS Build compact CIR-based IRF diagnostics
%
% The simulation-side diagnostics reconstructed aggregate levels from
% sectoral paths. Those checks are no longer part of the canonical
% postprocessing workflow. The remaining Diagnostics payload exists only to
% support the end-of-run IRF summary table.

Diagnostics = struct();

[first_irf_result, ~] = get_first_valid_irf_result(AllShockResults);
if ~isempty(first_irf_result)
    Diagnostics.has_irfs = true;
    irf_results = get_irf_result_list(AllShockResults);
    n_shocks = numel(irf_results);
    Diagnostics.irf_n_shocks = n_shocks;

    has_2ndorder_irf = has_secondorder_irf_stats(first_irf_result);

    Diagnostics.irf_cir_firstorder = NaN(n_shocks, 1);
    Diagnostics.irf_cir_pf = NaN(n_shocks, 1);
    Diagnostics.irf_nonlinear_amp_pf_vs_firstorder = NaN(n_shocks, 1);
    Diagnostics.irf_nonlinear_amplified_sectors = NaN(n_shocks, 1);
    Diagnostics.irf_nonlinear_attenuated_sectors = NaN(n_shocks, 1);
    Diagnostics.irf_shock_labels = cell(n_shocks, 1);

    if has_2ndorder_irf
        Diagnostics.irf_cir_secondorder = NaN(n_shocks, 1);
    end

    for i = 1:n_shocks
        irf_res = get_irf_result_at(irf_results, i);
        shock_cfg = config.shock_values(i);
        Diagnostics.irf_shock_labels{i} = shock_cfg.label;

        if ~is_valid_irf_result(irf_res)
            continue;
        end
        summary_stats = get_irf_summary_stats(irf_res);
        avg_cir_1st = mean(summary_stats.cumulative_responses.first_order);
        avg_cir_pf = mean(summary_stats.cumulative_responses.perfect_foresight);

        Diagnostics.irf_cir_firstorder(i) = avg_cir_1st * 100;
        Diagnostics.irf_cir_pf(i) = avg_cir_pf * 100;
        Diagnostics.irf_nonlinear_amp_pf_vs_firstorder(i) = finite_mean( ...
            summary_stats.nonlinear_amplification.pf_vs_first_order);
        Diagnostics.irf_nonlinear_amplified_sectors(i) = sum( ...
            summary_stats.nonlinear_effect_class.pf_vs_first_order > 0);
        Diagnostics.irf_nonlinear_attenuated_sectors(i) = sum( ...
            summary_stats.nonlinear_effect_class.pf_vs_first_order < 0);

        if has_2ndorder_irf
            avg_cir_2nd = mean(summary_stats.cumulative_responses.second_order);
            Diagnostics.irf_cir_secondorder(i) = avg_cir_2nd * 100;
        end

    end

    Diagnostics.irf_cir_asymmetry = build_cir_asymmetry_diagnostics( ...
        irf_results, config.shock_values, has_2ndorder_irf);
    Diagnostics.upstreamness = compute_upstreamness(params, ModData);
    Diagnostics.irf_sector_breakdown = build_irf_sector_breakdown( ...
        irf_results, config.shock_values, first_irf_result, Diagnostics.upstreamness);
end

end

function tf = has_secondorder_irf_stats(irf_result)
if ~(isstruct(irf_result) && ~isempty(fieldnames(irf_result)))
    tf = false;
    return;
end

summary_stats = get_irf_summary_stats(irf_result);
tf = isfield(summary_stats, 'cumulative_responses') && ...
    isfield(summary_stats.cumulative_responses, 'second_order') && ...
    ~isempty(summary_stats.cumulative_responses.second_order);
end

function tf = is_valid_irf_result(irf_result)
tf = isstruct(irf_result) && ~isempty(fieldnames(irf_result)) && ...
    isfield(irf_result, 'summary_stats');
end

function summary_stats = get_irf_summary_stats(irf_result)
summary_stats = irf_result.summary_stats;
end

function irf_results = get_irf_result_list(AllShockResults)
irf_results = {};
if ~isstruct(AllShockResults)
    return;
end
if isfield(AllShockResults, 'ShockArtifacts') && iscell(AllShockResults.ShockArtifacts) && ...
        ~isempty(AllShockResults.ShockArtifacts)
    irf_results = AllShockResults.ShockArtifacts;
end
end

function irf_res = get_irf_result_at(irf_results, idx)
if numel(irf_results) < idx
    irf_res = struct();
else
    irf_res = irf_results{idx};
end
end

function asymmetry = build_cir_asymmetry_diagnostics(irf_results, shock_values, has_2ndorder_irf)
pairs = find_positive_negative_shock_pairs(shock_values);
n_pairs = numel(pairs);

asymmetry = struct();
asymmetry.n_pairs = n_pairs;
asymmetry.size_pct = NaN(n_pairs, 1);
asymmetry.negative_label = cell(n_pairs, 1);
asymmetry.positive_label = cell(n_pairs, 1);
asymmetry.first_order = NaN(n_pairs, 1);
asymmetry.perfect_foresight = NaN(n_pairs, 1);
if has_2ndorder_irf
    asymmetry.second_order = NaN(n_pairs, 1);
end

for i = 1:n_pairs
    neg_idx = pairs(i).negative_idx;
    pos_idx = pairs(i).positive_idx;
    neg_res = get_irf_result_at(irf_results, neg_idx);
    pos_res = get_irf_result_at(irf_results, pos_idx);

    asymmetry.size_pct(i) = pairs(i).size_pct;
    asymmetry.negative_label{i} = get_optional_field(shock_values(neg_idx), 'label', '');
    asymmetry.positive_label{i} = get_optional_field(shock_values(pos_idx), 'label', '');

    if ~(is_valid_irf_result(neg_res) && is_valid_irf_result(pos_res))
        continue;
    end

    neg_cirs = neg_res.summary_stats.cumulative_responses;
    pos_cirs = pos_res.summary_stats.cumulative_responses;
    asymmetry.first_order(i) = mean(safe_divide_vectors(neg_cirs.first_order, pos_cirs.first_order));
    asymmetry.perfect_foresight(i) = mean(safe_divide_vectors(neg_cirs.perfect_foresight, pos_cirs.perfect_foresight));
    if has_2ndorder_irf
        asymmetry.second_order(i) = mean(safe_divide_vectors(neg_cirs.second_order, pos_cirs.second_order));
    end
end
end

function SectorBreakdown = build_irf_sector_breakdown(irf_results, shock_values, first_irf_result, Upstreamness)
pairs = find_positive_negative_shock_pairs(shock_values);
sector_indices = resolve_sector_indices(first_irf_result);
sector_labels = resolve_sector_labels(first_irf_result, sector_indices);
n_pairs = numel(pairs);
n_sectors = numel(sector_indices);

SectorBreakdown = struct();
SectorBreakdown.primary_upstreamness_measure = get_optional_field(Upstreamness, 'primary_measure', 'U_M');
SectorBreakdown.sector_indices = sector_indices;
SectorBreakdown.sector_labels = sector_labels;
SectorBreakdown.n_pairs = n_pairs;
SectorBreakdown.rows = initialize_sector_breakdown_rows(n_pairs, n_sectors);

for i = 1:n_pairs
    neg_idx = pairs(i).negative_idx;
    pos_idx = pairs(i).positive_idx;
    neg_res = get_irf_result_at(irf_results, neg_idx);
    pos_res = get_irf_result_at(irf_results, pos_idx);

    row = SectorBreakdown.rows(i);
    row.size_pct = pairs(i).size_pct;
    row.negative_label = get_optional_field(shock_values(neg_idx), 'label', '');
    row.positive_label = get_optional_field(shock_values(pos_idx), 'label', '');
    row.sector_indices = sector_indices;
    row.sector_labels = sector_labels;
    row.upstreamness = get_upstreamness_for_sectors(Upstreamness, sector_indices);

    if is_valid_irf_result(neg_res)
        row.negative_amplification_pf_vs_firstorder = ...
            neg_res.summary_stats.nonlinear_amplification.pf_vs_first_order(:).';
    end

    if is_valid_irf_result(pos_res)
        row.positive_amplification_pf_vs_firstorder = ...
            pos_res.summary_stats.nonlinear_amplification.pf_vs_first_order(:).';
    end

    if is_valid_irf_result(neg_res) && is_valid_irf_result(pos_res)
        neg_pf_cir = neg_res.summary_stats.cumulative_responses.perfect_foresight;
        pos_pf_cir = pos_res.summary_stats.cumulative_responses.perfect_foresight;
        row.pf_asymmetry = safe_divide_vectors(neg_pf_cir, pos_pf_cir);
    end

    row.corr_negative_amplification_upstreamness = safe_corr( ...
        row.negative_amplification_pf_vs_firstorder, row.upstreamness);
    row.corr_positive_amplification_upstreamness = safe_corr( ...
        row.positive_amplification_pf_vs_firstorder, row.upstreamness);
    row.corr_pf_asymmetry_upstreamness = safe_corr(row.pf_asymmetry, row.upstreamness);

    SectorBreakdown.rows(i) = row;
end
end

function rows = initialize_sector_breakdown_rows(n_pairs, n_sectors)
empty_row = struct( ...
    'size_pct', NaN, ...
    'negative_label', '', ...
    'positive_label', '', ...
    'sector_indices', NaN(1, n_sectors), ...
    'sector_labels', {cell(1, n_sectors)}, ...
    'upstreamness', NaN(1, n_sectors), ...
    'negative_amplification_pf_vs_firstorder', NaN(1, n_sectors), ...
    'positive_amplification_pf_vs_firstorder', NaN(1, n_sectors), ...
    'pf_asymmetry', NaN(1, n_sectors), ...
    'corr_negative_amplification_upstreamness', NaN, ...
    'corr_positive_amplification_upstreamness', NaN, ...
    'corr_pf_asymmetry_upstreamness', NaN);

if n_pairs == 0
    rows = repmat(empty_row, 0, 1);
else
    rows = repmat(empty_row, n_pairs, 1);
end
end

function sector_indices = resolve_sector_indices(irf_result)
sector_indices = [];
if isstruct(irf_result) && isfield(irf_result, 'metadata') && ...
        isfield(irf_result.metadata, 'sector_indices')
    sector_indices = irf_result.metadata.sector_indices(:).';
elseif isstruct(irf_result) && isfield(irf_result, 'entries')
    sector_indices = [irf_result.entries.sector_idx];
end
end

function sector_labels = resolve_sector_labels(irf_result, sector_indices)
n_sectors = numel(sector_indices);
sector_labels = arrayfun(@(idx) sprintf('Sector %d', idx), sector_indices, 'UniformOutput', false);

if ~(isstruct(irf_result) && isfield(irf_result, 'labels') && ...
        isfield(irf_result.labels, 'sector_labels'))
    return;
end

candidate_labels = irf_result.labels.sector_labels;
if iscell(candidate_labels) && numel(candidate_labels) == n_sectors
    sector_labels = candidate_labels(:).';
elseif iscell(candidate_labels) && numel(candidate_labels) >= max(sector_indices)
    sector_labels = candidate_labels(sector_indices);
end
end

function values = get_upstreamness_for_sectors(Upstreamness, sector_indices)
values = NaN(1, numel(sector_indices));
if ~(isstruct(Upstreamness) && isfield(Upstreamness, 'has_upstreamness') && ...
        Upstreamness.has_upstreamness && isfield(Upstreamness, 'U_M'))
    return;
end

valid = sector_indices >= 1 & sector_indices <= numel(Upstreamness.U_M);
values(valid) = Upstreamness.U_M(sector_indices(valid));
end

function rho = safe_corr(x, y)
x = x(:);
y = y(:);
valid = isfinite(x) & isfinite(y);
if sum(valid) < 2 || std(x(valid)) == 0 || std(y(valid)) == 0
    rho = NaN;
    return;
end

corr_matrix = corrcoef(x(valid), y(valid));
rho = corr_matrix(1, 2);
end

function pairs = find_positive_negative_shock_pairs(shock_values)
pairs = repmat(struct('size_pct', NaN, 'negative_idx', NaN, 'positive_idx', NaN), 0, 1);

for i = 1:numel(shock_values)
    if get_optional_field(shock_values(i), 'sign', 0) >= 0
        continue;
    end

    size_pct = get_optional_field(shock_values(i), 'size_pct', NaN);
    pos_idx = find_matching_positive_shock(shock_values, size_pct);
    if isnan(pos_idx)
        continue;
    end

    pairs(end + 1) = struct( ... %#ok<AGROW>
        'size_pct', size_pct, ...
        'negative_idx', i, ...
        'positive_idx', pos_idx);
end
end

function pos_idx = find_matching_positive_shock(shock_values, size_pct)
pos_idx = NaN;
if ~isfinite(size_pct)
    return;
end

for j = 1:numel(shock_values)
    candidate_size = get_optional_field(shock_values(j), 'size_pct', NaN);
    candidate_sign = get_optional_field(shock_values(j), 'sign', 0);
    if candidate_sign > 0 && isfinite(candidate_size) && abs(candidate_size - size_pct) < 1e-10
        pos_idx = j;
        return;
    end
end
end

function ratio = safe_divide_vectors(numerator, denominator)
ratio = NaN(size(numerator));
valid = isfinite(numerator) & isfinite(denominator) & abs(denominator) > 1e-12;
ratio(valid) = numerator(valid) ./ denominator(valid);
end

function value = finite_mean(values)
finite_values = values(isfinite(values));
if isempty(finite_values)
    value = NaN;
else
    value = mean(finite_values);
end
end

function value = get_optional_field(s, field_name, default_value)
value = default_value;
if isstruct(s) && isfield(s, field_name)
    value = s.(field_name);
end
end
