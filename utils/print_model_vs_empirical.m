function print_model_vs_empirical(ModelData)
% PRINT_MODEL_VS_EMPIRICAL Display saved model-vs-empirical comparison table

model_stats = get_saved_model_stats(ModelData);
has_theoretical_stats = isfield(ModelData, 'Statistics') && isstruct(ModelData.Statistics) && ...
    isfield(ModelData.Statistics, 'TheoStats') && ~isempty(fieldnames(ModelData.Statistics.TheoStats));

if ~has_theoretical_stats && model_stats.n_model_cols == 0
    return;
end

theo_stats = struct();
if has_theoretical_stats
    theo_stats = ModelData.Statistics.TheoStats;
end

emp_tgt = ModelData.EmpiricalTargets;
corr_cfg = get_ltfp_summary_config(ModelData.metadata.model_type);

fprintf('\n  ┌─ Model vs Empirical Comparison ────────────────────────────────┐\n');
if model_stats.n_model_cols == 0
    fprintf('  │                            Theo       Empirical  Ratio         │\n');
    print_ratio_row('σ(GDP_agg)', theo_stats.sigma_VA_agg, emp_tgt.sigma_VA_agg);
    if isfield(theo_stats, 'sigma_C_agg') && ~isnan(emp_tgt.sigma_C_agg)
        print_ratio_row('σ(C_agg)', theo_stats.sigma_C_agg, emp_tgt.sigma_C_agg);
    end
    print_ratio_row('σ(I_agg)', theo_stats.sigma_I_agg, emp_tgt.sigma_I_agg);
    print_ratio_row('corr(L,C)_agg', get_nested_or_nan(theo_stats, 'corr_L_C_agg'), ...
        get_nested_or_nan(emp_tgt.correlations, 'L_C_agg'));
    print_ratio_row('corr(I,C)_agg', get_nested_or_nan(theo_stats, 'corr_I_C_agg'), ...
        get_nested_or_nan(emp_tgt.correlations, 'I_C_agg'));
    if ~isempty(corr_cfg.model_agg_field)
        print_ratio_row('corr(L,A)_agg', get_nested_or_nan(theo_stats, corr_cfg.model_agg_field), ...
            get_nested_or_nan(emp_tgt.correlations, corr_cfg.data_agg_field));
    end
else
    if isfield(model_stats, 'ms1_source') && strcmp(model_stats.ms1_source, 'theoretical')
        fprintf('  │  1st column uses TheoStats; other columns use saved stats.    │\n');
    else
        fprintf('  │  Using saved simulation stats from ModelData.Statistics.*      │\n');
    end
    fprintf('  │                                                               │\n');
    fprintf('  │                    ');
    if model_stats.has_ms1,   fprintf('   1st   '); end
    if model_stats.has_ms2,   fprintf('   2nd   '); end
    if model_stats.has_msPF,  fprintf('    PF   '); end
    if model_stats.has_msMIT, fprintf('   MIT   '); end
    fprintf('   Data │\n');
    print_multi_row('σ(GDP_agg)', model_stats, 'sigma_VA_agg', emp_tgt.sigma_VA_agg);
    print_multi_row('σ(C_agg)', model_stats, 'sigma_C_agg', emp_tgt.sigma_C_agg);
    print_multi_row('σ(I_agg)', model_stats, 'sigma_I_agg', emp_tgt.sigma_I_agg);
    print_multi_row('σ(L_hc_agg)', model_stats, 'sigma_L_hc_agg', emp_tgt.sigma_L_agg);
    print_multi_row('corr(L,C)_agg', model_stats, 'corr_L_C_agg', get_nested_or_nan(emp_tgt.correlations, 'L_C_agg'));
    print_multi_row('corr(I,C)_agg', model_stats, 'corr_I_C_agg', get_nested_or_nan(emp_tgt.correlations, 'I_C_agg'));
    if ~isempty(corr_cfg.model_agg_field)
        print_multi_row('corr(L,A)_agg', model_stats, corr_cfg.model_agg_field, ...
            get_nested_or_nan(emp_tgt.correlations, corr_cfg.data_agg_field));
    end
    print_multi_row('σ(GDP) avg', model_stats, 'sigma_VA_avg', emp_tgt.sigma_VA_avg);
    print_multi_row('σ(L) avg', model_stats, 'sigma_L_avg', emp_tgt.sigma_L_avg);
    print_multi_row('σ(I) avg', model_stats, 'sigma_I_avg', emp_tgt.sigma_I_avg);
    print_multi_row('avg corr(C,C)', model_stats, 'avg_pairwise_corr_C', get_nested_or_nan(emp_tgt, 'avg_pairwise_corr_C'));
    print_multi_row('avg corr(GDP,GDP)', model_stats, 'avg_pairwise_corr_VA', get_nested_or_nan(emp_tgt, 'avg_pairwise_corr_VA'));
    print_multi_row('avg corr(L,L)', model_stats, 'avg_pairwise_corr_L', get_nested_or_nan(emp_tgt, 'avg_pairwise_corr_L'));
    print_multi_row('avg corr(I,I)', model_stats, 'avg_pairwise_corr_I', get_nested_or_nan(emp_tgt, 'avg_pairwise_corr_I'));
end
fprintf('  └────────────────────────────────────────────────────────────────┘\n');
end

function print_ratio_row(label, model_value, empirical_value)
ratio_value = safe_ratio(model_value, empirical_value);
fprintf('  │  %-18s %10.4f  %10.4f  %10.2f │\n', label, model_value, empirical_value, ratio_value);
end

function print_multi_row(label, model_stats, field_name, empirical_value)
fprintf('  │  %-16s ', label);
if model_stats.has_ms1,   fprintf('%7.4f  ', get_nested_or_nan(model_stats.ms1, field_name)); end
if model_stats.has_ms2,   fprintf('%7.4f  ', get_nested_or_nan(model_stats.ms2, field_name)); end
if model_stats.has_msPF,  fprintf('%7.4f  ', get_nested_or_nan(model_stats.msPF, field_name)); end
if model_stats.has_msMIT, fprintf('%7.4f  ', get_nested_or_nan(model_stats.msMIT, field_name)); end
fprintf('%7.4f │\n', empirical_value);
end

function ratio_value = safe_ratio(numerator, denominator)
if ~isfinite(numerator) || ~isfinite(denominator) || abs(denominator) < 1e-12
    ratio_value = NaN;
else
    ratio_value = numerator / denominator;
end
end
