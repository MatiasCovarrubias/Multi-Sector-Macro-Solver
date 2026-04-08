function Diagnostics = build_nonlinearity_diagnostics(ModelData_simulation, AllShockResults, params, config, ModData)
% BUILD_NONLINEARITY_DIAGNOSTICS Build the compact IRF summary diagnostics
%
% The old simulation-side diagnostics reconstructed aggregate levels from
% sectoral paths. Those checks are no longer part of the canonical
% postprocessing workflow. The remaining Diagnostics payload exists only to
% support the end-of-run IRF amplification summary table.

%#ok<INUSD>

Diagnostics = struct();

[first_irf_result, ~] = get_first_valid_irf_result(AllShockResults);
if ~isempty(first_irf_result)
    Diagnostics.has_irfs = true;
    irf_results = get_irf_result_list(AllShockResults);
    n_shocks = numel(irf_results);
    Diagnostics.irf_n_shocks = n_shocks;

    has_2ndorder_irf = has_secondorder_irf_stats(first_irf_result);

    Diagnostics.irf_peak_firstorder = NaN(n_shocks, 1);
    Diagnostics.irf_peak_pf = NaN(n_shocks, 1);
    Diagnostics.irf_amplification_rel = NaN(n_shocks, 1);
    Diagnostics.irf_shock_labels = cell(n_shocks, 1);

    if has_2ndorder_irf
        Diagnostics.irf_peak_secondorder = NaN(n_shocks, 1);
        Diagnostics.irf_amplification_2nd = NaN(n_shocks, 1);
    end

    for i = 1:n_shocks
        irf_res = get_irf_result_at(irf_results, i);
        shock_cfg = config.shock_values(i);
        Diagnostics.irf_shock_labels{i} = shock_cfg.label;

        if ~is_valid_irf_result(irf_res)
            continue;
        end
        summary_stats = get_irf_summary_stats(irf_res);
        avg_peak_1st = mean(summary_stats.peaks.first_order);
        avg_peak_pf = mean(summary_stats.peaks.perfect_foresight);

        Diagnostics.irf_peak_firstorder(i) = avg_peak_1st * 100;
        Diagnostics.irf_peak_pf(i) = avg_peak_pf * 100;

        if has_2ndorder_irf
            avg_peak_2nd = mean(summary_stats.peaks.second_order);
            Diagnostics.irf_peak_secondorder(i) = avg_peak_2nd * 100;
            if avg_peak_1st > 1e-10
                Diagnostics.irf_amplification_2nd(i) = (avg_peak_2nd / avg_peak_1st - 1) * 100;
            else
                Diagnostics.irf_amplification_2nd(i) = 0;
            end
        end

        if isfield(summary_stats, 'amplifications') && isfield(summary_stats.amplifications, 'rel')
            Diagnostics.irf_amplification_rel(i) = mean(summary_stats.amplifications.rel);
        elseif avg_peak_1st > 1e-10
            Diagnostics.irf_amplification_rel(i) = (avg_peak_pf / avg_peak_1st - 1) * 100;
        else
            Diagnostics.irf_amplification_rel(i) = 0;
        end

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
