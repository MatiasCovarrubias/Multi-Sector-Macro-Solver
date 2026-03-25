function Summary = build_irf_summary_rows(AllShockResults, config, n_shocks)
% BUILD_IRF_SUMMARY_ROWS Build IRF summary rows without printing

Summary = struct();
[first_irf_result, ~] = get_first_valid_irf_result(AllShockResults);
Summary.has_2ndorder_irfs = has_secondorder_irf_stats(first_irf_result);
irf_results = get_irf_result_list(AllShockResults);
Summary.rows = repmat(struct( ...
    'label', '', ...
    'A_level', NaN, ...
    'avg_peak_1st', NaN, ...
    'avg_peak_2nd', NaN, ...
    'avg_peak_pf', NaN, ...
    'avg_amplif_rel', NaN), n_shocks, 1);

for i = 1:n_shocks
    irf_res = get_irf_result_at(irf_results, i);
    shock_cfg = config.shock_values(i);

    Summary.rows(i).label = shock_cfg.label;
    Summary.rows(i).A_level = exp(-shock_cfg.value);
    if ~is_valid_irf_result(irf_res)
        continue;
    end
    summary_stats = get_irf_summary_stats(irf_res);
    Summary.rows(i).avg_peak_1st = mean(summary_stats.peaks.first_order);
    Summary.rows(i).avg_peak_pf = mean(summary_stats.peaks.perfect_foresight);
    Summary.rows(i).avg_amplif_rel = mean(summary_stats.amplifications.rel);

    if Summary.has_2ndorder_irfs
        Summary.rows(i).avg_peak_2nd = mean(summary_stats.peaks.second_order);
    end
end

end

function tf = has_secondorder_irf_stats(irf_result)
if ~(isstruct(irf_result) && ~isempty(fieldnames(irf_result)))
    tf = false;
    return;
end

summary_stats = get_irf_summary_stats(irf_result);
tf = isfield(summary_stats, 'peaks') && isfield(summary_stats.peaks, 'second_order') && ...
    ~isempty(summary_stats.peaks.second_order);
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
