function print_summary_table(ModelData)
% PRINT_SUMMARY_TABLE Print the compact end-of-run summary tables.

Summary = build_summary_table_data(ModelData);
ms = Summary.model_stats;
diag = Summary.diagnostics;

fprintf('\n');
fprintf('================================================================================\n');
fprintf('                            SUMMARY TABLE\n');
fprintf('================================================================================\n');
fprintf('Experiment: %s\n', Summary.save_label);

print_model_vs_data_table(ms, Summary.empirical_targets, Summary.config);
print_aggregate_moment_tables(ms);
print_irf_amplification_table(diag);

fprintf('================================================================================\n\n');
end

function print_model_vs_data_table(ms, emp_tgt, config)
fprintf('\n[1] MODEL VS DATA MOMENTS\n');

if ms.n_model_cols == 0
    fprintf('    (No simulation statistics available)\n');
    return;
end

corr_cfg = get_ltfp_summary_config(config.model_type);

fprintf('    %-16s', 'Moment');
print_method_headers(ms, true);
fprintf('\n');
fprintf('    %s\n', repmat('-', 1, 16 + 11 * (ms.n_model_cols + 1)));

rows = {
    'σ(GDP_agg)', 'sigma_VA_agg',   emp_tgt.sigma_VA_agg;
    'σ(C_agg)',   'sigma_C_agg',    emp_tgt.sigma_C_agg;
    'σ(I_agg)',   'sigma_I_agg',    emp_tgt.sigma_I_agg;
    'σ(L_hc_agg)','sigma_L_hc_agg', emp_tgt.sigma_L_agg;
    'corr(L,C)_agg',  'corr_L_C_agg',   get_nested_or_nan(emp_tgt.correlations, 'L_C_agg');
    'corr(I,C)_agg',  'corr_I_C_agg',   get_nested_or_nan(emp_tgt.correlations, 'I_C_agg');
    'σ(GDP) avg', 'sigma_VA_avg',   get_nested_or_nan(emp_tgt, 'sigma_VA_avg');
    'σ(L) avg',   'sigma_L_avg',    get_nested_or_nan(emp_tgt, 'sigma_L_avg');
    'σ(I) avg',   'sigma_I_avg',    get_nested_or_nan(emp_tgt, 'sigma_I_avg');
    'σ(L) emp-wgt', 'sigma_L_avg_empweighted', get_nested_or_nan(emp_tgt, 'sigma_L_avg_empweighted');
    'σ(I) inv-wgt', 'sigma_I_avg_invweighted', get_nested_or_nan(emp_tgt, 'sigma_I_avg_invweighted');
    'σ(Domar)avg', 'sigma_Domar_avg_legacy', get_nested_or_nan(emp_tgt, 'sigma_Domar_avg');
    'avg corr(C,C)', 'avg_pairwise_corr_C', get_nested_or_nan(emp_tgt, 'avg_pairwise_corr_C');
    'avg corr(GDP,GDP)', 'avg_pairwise_corr_VA', get_nested_or_nan(emp_tgt, 'avg_pairwise_corr_VA');
    'avg corr(L,L)', 'avg_pairwise_corr_L', get_nested_or_nan(emp_tgt, 'avg_pairwise_corr_L');
    'avg corr(I,I)', 'avg_pairwise_corr_I', get_nested_or_nan(emp_tgt, 'avg_pairwise_corr_I');
};

for i = 1:size(rows, 1)
    fprintf('    %-16s', rows{i, 1});
    print_method_values(ms, rows{i, 2});
    fprintf('%11s\n', format_stat(rows{i, 3}, 4));
end

if ~isempty(corr_cfg.model_agg_field)
    fprintf('    %-16s', 'corr(L,A)_agg');
    print_method_values(ms, corr_cfg.model_agg_field);
    fprintf('%11s\n', format_stat(get_nested_or_nan(emp_tgt.correlations, corr_cfg.data_agg_field), 4));

    fprintf('    %-16s', 'avg corr(L,A)');
    print_method_values(ms, corr_cfg.model_avg_field);
    fprintf('%11s\n', format_stat(get_nested_or_nan(emp_tgt.correlations, corr_cfg.data_avg_field), 4));
end
end

function print_aggregate_moment_tables(ms)
fprintf('\n[2] AGGREGATE MOMENTS\n');

if ms.n_model_cols == 0
    fprintf('    (No simulation statistics available)\n');
    return;
end

aggregates = {
    'C', 'C';
    'GDP', 'GDP';
    'L', 'L';
    'I', 'I';
    'U_intra', 'utility_intratemp';
};

for i = 1:size(aggregates, 1)
    label = aggregates{i, 1};
    field_name = aggregates{i, 2};

    if ~any_method_has_aggregate_moments(ms, field_name)
        continue;
    end

    fprintf('\n    %s_agg\n', label);
    fprintf('      %-14s%12s%12s%12s%12s\n', 'Method', 'Mean (%)', 'Std (%)', 'Skew', 'Ex. Kurt');
    fprintf('      %s\n', repmat('-', 1, 62));

    print_aggregate_moment_row(ms.has_ms1, '1st', ms.ms1, field_name);
    print_aggregate_moment_row(ms.has_ms2, '2nd', ms.ms2, field_name);
    print_aggregate_moment_row(ms.has_msPF, 'PF',  ms.msPF, field_name);
    print_aggregate_moment_row(ms.has_msMIT, 'MIT', ms.msMIT, field_name);
end
end

function print_irf_amplification_table(diag)
fprintf('\n[3] IRF AMPLIFICATION\n');

if isempty(fieldnames(diag)) || ~isfield(diag, 'has_irfs') || ~diag.has_irfs
    fprintf('    (No IRF results available)\n');
    return;
end

has_2nd = isfield(diag, 'irf_peak_secondorder');

if has_2nd
    fprintf('    %-16s%12s%12s%12s%14s\n', 'Shock', 'Peak(1st)', 'Peak(2nd)', 'Peak(PF)', 'Amplif(PF)');
    fprintf('    %s\n', repmat('-', 1, 66));
    for i = 1:diag.irf_n_shocks
        fprintf('    %-16s%12s%12s%12s%14s\n', ...
            diag.irf_shock_labels{i}, ...
            format_percent(diag.irf_peak_firstorder(i), 3), ...
            format_percent(diag.irf_peak_secondorder(i), 3), ...
            format_percent(diag.irf_peak_pf(i), 3), ...
            format_percent(diag.irf_amplification_rel(i), 1));
    end
else
    fprintf('    %-16s%12s%12s%14s\n', 'Shock', 'Peak(1st)', 'Peak(PF)', 'Amplif(PF)');
    fprintf('    %s\n', repmat('-', 1, 54));
    for i = 1:diag.irf_n_shocks
        fprintf('    %-16s%12s%12s%14s\n', ...
            diag.irf_shock_labels{i}, ...
            format_percent(diag.irf_peak_firstorder(i), 3), ...
            format_percent(diag.irf_peak_pf(i), 3), ...
            format_percent(diag.irf_amplification_rel(i), 1));
    end
end
end

function print_method_headers(ms, include_data)
if ms.has_ms1,   fprintf('%11s', '1st'); end
if ms.has_ms2,   fprintf('%11s', '2nd'); end
if ms.has_msPF,  fprintf('%11s', 'PF'); end
if ms.has_msMIT, fprintf('%11s', 'MIT'); end
if include_data, fprintf('%11s', 'Data'); end
end

function print_method_values(ms, field_name)
if ms.has_ms1,   fprintf('%11s', format_stat(get_nested_or_nan(ms.ms1, field_name), 4)); end
if ms.has_ms2,   fprintf('%11s', format_stat(get_nested_or_nan(ms.ms2, field_name), 4)); end
if ms.has_msPF,  fprintf('%11s', format_stat(get_nested_or_nan(ms.msPF, field_name), 4)); end
if ms.has_msMIT, fprintf('%11s', format_stat(get_nested_or_nan(ms.msMIT, field_name), 4)); end
end

function tf = any_method_has_aggregate_moments(ms, field_name)
tf = has_aggregate_moments(ms.ms1, field_name) || ...
    has_aggregate_moments(ms.ms2, field_name) || ...
    has_aggregate_moments(ms.msPF, field_name) || ...
    has_aggregate_moments(ms.msMIT, field_name);
end

function tf = has_aggregate_moments(method_stats, field_name)
tf = isstruct(method_stats) && isfield(method_stats, 'aggregate_moments') && ...
    isstruct(method_stats.aggregate_moments) && isfield(method_stats.aggregate_moments, field_name);
end

function print_aggregate_moment_row(is_available, method_label, method_stats, field_name)
if ~is_available
    return;
end

moments = struct('mean', NaN, 'std', NaN, 'skewness', NaN, 'kurtosis', NaN);
if has_aggregate_moments(method_stats, field_name)
    moments = method_stats.aggregate_moments.(field_name);
end

fprintf('      %-14s%12s%12s%12s%12s\n', ...
    method_label, ...
    format_percent(100 * moments.mean, 4), ...
    format_percent(100 * moments.std, 4), ...
    format_stat(moments.skewness, 4), ...
    format_stat(moments.kurtosis - 3, 4));
end

function out = format_stat(value, decimals)
if nargin < 2
    decimals = 4;
end

if ~isfinite(value)
    out = 'n/a';
    return;
end

out = sprintf(['%0.', num2str(decimals), 'f'], value);
end

function out = format_percent(value, decimals)
if nargin < 2
    decimals = 4;
end

if ~isfinite(value)
    out = 'n/a';
    return;
end

out = sprintf(['%+0.', num2str(decimals), 'f%%'], value);
end

