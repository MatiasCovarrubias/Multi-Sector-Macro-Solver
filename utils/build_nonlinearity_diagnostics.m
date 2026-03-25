function Diagnostics = build_nonlinearity_diagnostics(ModelData_simulation, AllShockResults, params, config, ModData)
% BUILD_NONLINEARITY_DIAGNOSTICS Compute nonlinearity diagnostics without printing

n_sectors = params.n_sectors;
idx = get_variable_indices(n_sectors);

endostates_ss = ModData.endostates_ss;
policies_ss = ModData.policies_ss;
k_ss = endostates_ss(1:n_sectors);

Diagnostics = struct();
Diagnostics.has_firstorder = false;
Diagnostics.has_secondorder = false;
Diagnostics.has_pf = false;
Diagnostics.has_mit = false;
Diagnostics.has_irfs = false;
Diagnostics.simulation_moment_window = 'shocks_simul';

if isfield(ModelData_simulation, 'PerfectForesight') && isfield(ModelData_simulation.PerfectForesight, 'shocks_simul')
    Diagnostics.has_pf = true;
    has_pf = true;

    k_determ = ModelData_simulation.PerfectForesight.shocks_simul(idx.k(1):idx.k(2), :);
    avg_k_logdev = mean(k_determ - k_ss, 2);

    Diagnostics.prealloc_mean_abs_k = mean(abs(avg_k_logdev)) * 100;
    Diagnostics.prealloc_k_mining = avg_k_logdev(1) * 100;

    [max_dev, max_idx] = max(avg_k_logdev);
    [min_dev, min_idx] = min(avg_k_logdev);
    all_labels = SectorLabel(1:n_sectors);
    Diagnostics.prealloc_k_max = max_dev * 100;
    Diagnostics.prealloc_k_max_sector = all_labels.display{max_idx};
    Diagnostics.prealloc_k_min = min_dev * 100;
    Diagnostics.prealloc_k_min_sector = all_labels.display{min_idx};
    Diagnostics.prealloc_k_all = avg_k_logdev * 100;
else
    has_pf = false;
end

if isfield(ModelData_simulation, 'FirstOrder') && isfield(ModelData_simulation.FirstOrder, 'shocks_simul')
    Diagnostics.has_firstorder = true;
    simul_1st = ModelData_simulation.FirstOrder.shocks_simul;
    [C_agg_logdev, I_agg_logdev, GDP_agg_logdev, L_agg_logdev, M_agg_logdev] = ...
        compute_exact_aggregate_logdevs(simul_1st, idx, policies_ss, endostates_ss);

    Diagnostics.firstorder_C_mean_logdev = mean(C_agg_logdev) * 100;
    Diagnostics.firstorder_C_std = std(C_agg_logdev) * 100;
    Diagnostics.firstorder_approx_order = 1;
    Diagnostics.firstorder_Y_mean_logdev = mean(GDP_agg_logdev) * 100;
    Diagnostics.firstorder_Y_std = std(GDP_agg_logdev) * 100;
    Diagnostics.firstorder_L_mean_logdev = mean(L_agg_logdev) * 100;
    Diagnostics.firstorder_I_mean_logdev = mean(I_agg_logdev) * 100;
    Diagnostics.firstorder_M_mean_logdev = mean(M_agg_logdev) * 100;
end

if isfield(ModelData_simulation, 'SecondOrder') && isfield(ModelData_simulation.SecondOrder, 'shocks_simul')
    Diagnostics.has_secondorder = true;
    simul_2nd = ModelData_simulation.SecondOrder.shocks_simul;
    [C_agg_logdev, I_agg_logdev, GDP_agg_logdev, L_agg_logdev, M_agg_logdev] = ...
        compute_exact_aggregate_logdevs(simul_2nd, idx, policies_ss, endostates_ss);

    Diagnostics.secondorder_C_mean_logdev = mean(C_agg_logdev) * 100;
    Diagnostics.secondorder_C_std = std(C_agg_logdev) * 100;
    Diagnostics.secondorder_Y_mean_logdev = mean(GDP_agg_logdev) * 100;
    Diagnostics.secondorder_Y_std = std(GDP_agg_logdev) * 100;
    Diagnostics.secondorder_L_mean_logdev = mean(L_agg_logdev) * 100;
    Diagnostics.secondorder_I_mean_logdev = mean(I_agg_logdev) * 100;
    Diagnostics.secondorder_M_mean_logdev = mean(M_agg_logdev) * 100;
end

if has_pf
    simul_pf = ModelData_simulation.PerfectForesight.shocks_simul;
    [C_agg_logdev, I_agg_logdev, GDP_agg_logdev, L_agg_logdev, M_agg_logdev] = ...
        compute_exact_aggregate_logdevs(simul_pf, idx, policies_ss, endostates_ss);

    Diagnostics.pf_C_mean_logdev = mean(C_agg_logdev) * 100;
    Diagnostics.pf_C_std = std(C_agg_logdev) * 100;
    Diagnostics.pf_Y_mean_logdev = mean(GDP_agg_logdev) * 100;
    Diagnostics.pf_Y_std = std(GDP_agg_logdev) * 100;
    Diagnostics.pf_L_mean_logdev = mean(L_agg_logdev) * 100;
    Diagnostics.pf_I_mean_logdev = mean(I_agg_logdev) * 100;
    Diagnostics.pf_M_mean_logdev = mean(M_agg_logdev) * 100;
end

if isfield(ModelData_simulation, 'MITShocks') && isfield(ModelData_simulation.MITShocks, 'shocks_simul')
    Diagnostics.has_mit = true;
    simul_mit = ModelData_simulation.MITShocks.shocks_simul;
    [C_agg_logdev, I_agg_logdev, GDP_agg_logdev, L_agg_logdev, M_agg_logdev] = ...
        compute_exact_aggregate_logdevs(simul_mit, idx, policies_ss, endostates_ss);

    Diagnostics.mit_C_mean_logdev = mean(C_agg_logdev) * 100;
    Diagnostics.mit_C_std = std(C_agg_logdev) * 100;
    Diagnostics.mit_Y_mean_logdev = mean(GDP_agg_logdev) * 100;
    Diagnostics.mit_Y_std = std(GDP_agg_logdev) * 100;
    Diagnostics.mit_L_mean_logdev = mean(L_agg_logdev) * 100;
    Diagnostics.mit_I_mean_logdev = mean(I_agg_logdev) * 100;
    Diagnostics.mit_M_mean_logdev = mean(M_agg_logdev) * 100;
end

if Diagnostics.has_firstorder && Diagnostics.has_pf
    Diagnostics.precautionary_C = Diagnostics.pf_C_mean_logdev - Diagnostics.firstorder_C_mean_logdev;
    Diagnostics.precautionary_Y = Diagnostics.pf_Y_mean_logdev - Diagnostics.firstorder_Y_mean_logdev;
    Diagnostics.vol_diff_C = Diagnostics.pf_C_std - Diagnostics.firstorder_C_std;
    Diagnostics.vol_diff_Y = Diagnostics.pf_Y_std - Diagnostics.firstorder_Y_std;
end

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
