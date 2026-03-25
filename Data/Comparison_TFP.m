% OFFLINE PREPROCESSING SCRIPT: not called by main.m at runtime.
% This code compares the descriptive statistics of the TFP process
% parameters (the persistence parameter, rho
% and the elements of the variance-covariance matrix) for different ways of
% computing them.

clear all; clc; 
tfp_process = load('TFP_process.mat');

rowNames = {'Mean','St. Deviation','Variance','Min','25th Prctle','Median','75th Prctle','Max'};

rho_field_names = {'modrho', 'modrho_sm', 'modrho_wds', 'modrho_sm_wds', ...
    'modrho_GO', 'modrho_GO_sm', 'modrho_GO_wds', 'modrho_GO_sm_wds', ...
    'modrho_GO_noVA', 'modrho_GO_noVA_sm', 'modrho_GO_noVA_wds', 'modrho_GO_noVA_sm_wds'};
vcv_field_names = {'modvcv', 'modvcv_sm', 'modvcv_wds', 'modvcv_sm_wds', ...
    'modvcv_GO', 'modvcv_GO_sm', 'modvcv_GO_wds', 'modvcv_GO_sm_wds', ...
    'modvcv_GO_noVA', 'modvcv_GO_noVA_sm', 'modvcv_GO_noVA_wds', 'modvcv_GO_noVA_sm_wds'};
resid_field_names = {'ar1resid', 'ar1resid_sm', 'ar1resid_wds', 'ar1resid_sm_wds', ...
    'ar1resid_GO', 'ar1resid_GO_sm', 'ar1resid_GO_wds', 'ar1resid_GO_sm_wds', ...
    'ar1resid_GO_noVA', 'ar1resid_GO_noVA_sm', 'ar1resid_GO_noVA_wds', 'ar1resid_GO_noVA_sm_wds'};
display_names = {'VA', 'VA_sm', 'VA_wds', 'VA_sm_wds', ...
    'GO', 'GO_sm', 'GO_wds', 'GO_sm_wds', ...
    'GO_noVA', 'GO_noVA_sm', 'GO_noVA_wds', 'GO_noVA_sm_wds'};

n_variants = numel(display_names);
rho_matrix = zeros(37, n_variants);
diag_vcv_matrix = zeros(37, n_variants);
shock_std_matrix = zeros(37, n_variants);
avg_shock_correlation = zeros(1, n_variants);
Statistics_OffDiagonal_vcv = cell(1, 36);
offdiag_data = cell(1, n_variants);

for idx = 1:n_variants
    rho_matrix(:, idx) = tfp_process.(rho_field_names{idx});
    vcv_matrix = tfp_process.(vcv_field_names{idx});
    resid_matrix = tfp_process.(resid_field_names{idx});

    diag_vcv_matrix(:, idx) = diag(vcv_matrix);
    shock_std_matrix(:, idx) = std(resid_matrix)';
    avg_shock_correlation(idx) = mean(nonzeros(triu(corr(resid_matrix), 1)));

    for lag = 1:36
        offdiag_data{idx, lag} = diag(vcv_matrix, -lag);
    end
end

Statistics_modrho = build_summary_table(rho_matrix, display_names, rowNames);
Statistics_Diagonal_vcv = build_summary_table(diag_vcv_matrix, display_names, rowNames);
Statistics_shock_std = build_summary_table(shock_std_matrix, display_names, rowNames);
Sectoral_shock_std = array2table(shock_std_matrix, 'VariableNames', display_names);
Avg_shock_correlation = array2table(avg_shock_correlation, 'VariableNames', display_names);

for lag = 1:36
    lag_matrix = zeros(37 - lag, n_variants);
    for idx = 1:n_variants
        lag_matrix(:, idx) = offdiag_data{idx, lag};
    end
    Statistics_OffDiagonal_vcv{1, lag} = build_summary_table(lag_matrix, display_names, rowNames);
end

comparison_pairs = { ...
    'VA', 'modvcv', 'modvcv_wds', 'ar1resid', 'ar1resid_wds'; ...
    'VA_sm', 'modvcv_sm', 'modvcv_sm_wds', 'ar1resid_sm', 'ar1resid_sm_wds'; ...
    'GO', 'modvcv_GO', 'modvcv_GO_wds', 'ar1resid_GO', 'ar1resid_GO_wds'; ...
    'GO_sm', 'modvcv_GO_sm', 'modvcv_GO_sm_wds', 'ar1resid_GO_sm', 'ar1resid_GO_sm_wds'; ...
    'GO_noVA', 'modvcv_GO_noVA', 'modvcv_GO_noVA_wds', 'ar1resid_GO_noVA', 'ar1resid_GO_noVA_wds'; ...
    'GO_noVA_sm', 'modvcv_GO_noVA_sm', 'modvcv_GO_noVA_sm_wds', 'ar1resid_GO_noVA_sm', 'ar1resid_GO_noVA_sm_wds'};

n_pairs = size(comparison_pairs, 1);
pair_labels = comparison_pairs(:, 1);
mean_diag_before = zeros(n_pairs, 1);
mean_diag_after = zeros(n_pairs, 1);
mean_diag_pct_change = zeros(n_pairs, 1);
median_diag_pct_change = zeros(n_pairs, 1);
mean_abs_offdiag_before = zeros(n_pairs, 1);
mean_abs_offdiag_after = zeros(n_pairs, 1);
mean_abs_offdiag_pct_change = zeros(n_pairs, 1);
fro_norm_before = zeros(n_pairs, 1);
fro_norm_after = zeros(n_pairs, 1);
fro_norm_pct_change = zeros(n_pairs, 1);
avg_corr_before = zeros(n_pairs, 1);
avg_corr_after = zeros(n_pairs, 1);
avg_corr_change = zeros(n_pairs, 1);

for pair_idx = 1:n_pairs
    vcv_before = tfp_process.(comparison_pairs{pair_idx, 2});
    vcv_after = tfp_process.(comparison_pairs{pair_idx, 3});
    resid_before = tfp_process.(comparison_pairs{pair_idx, 4});
    resid_after = tfp_process.(comparison_pairs{pair_idx, 5});

    diag_before = diag(vcv_before);
    diag_after = diag(vcv_after);
    offdiag_before = vcv_before(~eye(size(vcv_before)));
    offdiag_after = vcv_after(~eye(size(vcv_after)));

    mean_diag_before(pair_idx) = mean(diag_before);
    mean_diag_after(pair_idx) = mean(diag_after);
    mean_diag_pct_change(pair_idx) = 100 * (mean_diag_after(pair_idx) / mean_diag_before(pair_idx) - 1);
    median_diag_pct_change(pair_idx) = median(100 * (diag_after ./ diag_before - 1));
    mean_abs_offdiag_before(pair_idx) = mean(abs(offdiag_before));
    mean_abs_offdiag_after(pair_idx) = mean(abs(offdiag_after));
    mean_abs_offdiag_pct_change(pair_idx) = 100 * (mean_abs_offdiag_after(pair_idx) / mean_abs_offdiag_before(pair_idx) - 1);
    fro_norm_before(pair_idx) = norm(vcv_before, 'fro');
    fro_norm_after(pair_idx) = norm(vcv_after, 'fro');
    fro_norm_pct_change(pair_idx) = 100 * (fro_norm_after(pair_idx) / fro_norm_before(pair_idx) - 1);
    avg_corr_before(pair_idx) = mean(nonzeros(triu(corr(resid_before), 1)));
    avg_corr_after(pair_idx) = mean(nonzeros(triu(corr(resid_after), 1)));
    avg_corr_change(pair_idx) = avg_corr_after(pair_idx) - avg_corr_before(pair_idx);
end

Winsorization_VCV_Impact = table( ...
    mean_diag_before, mean_diag_after, mean_diag_pct_change, median_diag_pct_change, ...
    mean_abs_offdiag_before, mean_abs_offdiag_after, mean_abs_offdiag_pct_change, ...
    fro_norm_before, fro_norm_after, fro_norm_pct_change, ...
    avg_corr_before, avg_corr_after, avg_corr_change, ...
    'RowNames', pair_labels, ...
    'VariableNames', {'meanDiagVarBefore', 'meanDiagVarAfter', 'meanDiagVarPctChange', ...
    'medianSectorDiagVarPctChange', 'meanAbsOffDiagBefore', 'meanAbsOffDiagAfter', ...
    'meanAbsOffDiagPctChange', 'froNormBefore', 'froNormAfter', 'froNormPctChange', ...
    'avgCorrBefore', 'avgCorrAfter', 'avgCorrChange'});

fprintf('\nWinsorization impact on TFP shock covariance:\n');
disp(Winsorization_VCV_Impact);

clearvars -except Statistics_modrho Statistics_Diagonal_vcv Statistics_OffDiagonal_vcv ...
    Statistics_shock_std Sectoral_shock_std Avg_shock_correlation Winsorization_VCV_Impact;

function summary_table = build_summary_table(data_matrix, variable_names, row_names)
summary_matrix = [mean(data_matrix, 1); std(data_matrix, 0, 1); var(data_matrix, 0, 1); min(data_matrix, [], 1); ...
    prctile(data_matrix, 25, 1); median(data_matrix, 1); prctile(data_matrix, 75, 1); max(data_matrix, [], 1)];
summary_table = array2table(summary_matrix, 'VariableNames', variable_names, 'RowNames', row_names);
end

