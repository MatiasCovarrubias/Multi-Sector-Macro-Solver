function print_empirical_targets(emp_tgt)
% PRINT_EMPIRICAL_TARGETS Display empirical target moments
%
% INPUTS:
%   emp_tgt - Empirical targets structure from calib_data.empirical_targets

fprintf('\n  Empirical Targets (HP-filtered, λ=%d, %s aggregation):\n', ...
    emp_tgt.hp_lambda, emp_tgt.aggregation_method);
fprintf('    ── Aggregate volatilities ──\n');
fprintf('    σ(GDP_agg):       %.4f   (aggregate GDP, Törnqvist)\n', emp_tgt.sigma_VA_agg);
if ~isnan(emp_tgt.sigma_C_agg)
    fprintf('    σ(C_agg):         %.4f   (aggregate consumption, NIPA PCE)\n', emp_tgt.sigma_C_agg);
else
    fprintf('    σ(C_agg):         N/A     (Real GDP Components.xls not found)\n');
end
fprintf('    σ(L_agg):         %.4f   (aggregate labor, simple sum)\n', emp_tgt.sigma_L_agg);
fprintf('    σ(I_agg):         %.4f   (aggregate investment, Törnqvist)\n', emp_tgt.sigma_I_agg);
fprintf('    ── Average sectoral volatilities (VA-weighted) ──\n');
fprintf('    σ(VA) avg:        %.4f   (VA-weighted avg of sectoral value-added vol)\n', emp_tgt.sigma_VA_avg);
fprintf('    σ(L) avg:         %.4f   (VA-weighted avg of sectoral labor vol)\n', emp_tgt.sigma_L_avg);
fprintf('    σ(I) avg:         %.4f   (VA-weighted avg of sectoral investment vol)\n', emp_tgt.sigma_I_avg);
fprintf('    ── Average sectoral volatilities (own-variable weighted) ──\n');
fprintf('    σ(VA) own-wgt:    %.4f   (VA-share weighted avg of sectoral value-added vol)\n', emp_tgt.sigma_VA_avg);
fprintf('    σ(L) emp-wgt:     %.4f   (employment-weighted avg of sectoral labor vol)\n', emp_tgt.sigma_L_avg_empweighted);
fprintf('    σ(I) inv-wgt:     %.4f   (investment-weighted avg of sectoral investment vol)\n', emp_tgt.sigma_I_avg_invweighted);
fprintf('    ── Sectoral comovement (avg pairwise corr) ──\n');
if isfield(emp_tgt, 'avg_pairwise_corr_C') && ~isnan(emp_tgt.avg_pairwise_corr_C)
    fprintf('    avg corr(C_j,C_k): %.4f   (sectoral consumption expenditure)\n', emp_tgt.avg_pairwise_corr_C);
end
fprintf('    avg corr(VA_j,VA_k):   %.4f   (sectoral value added)\n', emp_tgt.avg_pairwise_corr_VA);
fprintf('    avg corr(L_j,L_k): %.4f   (sectoral labor)\n', emp_tgt.avg_pairwise_corr_L);
fprintf('    avg corr(I_j,I_k): %.4f   (sectoral investment)\n', emp_tgt.avg_pairwise_corr_I);

if isfield(emp_tgt, 'correlations') && isstruct(emp_tgt.correlations)
    corr_tgt = emp_tgt.correlations;
    fprintf('    ── Aggregate expenditure correlations (data only) ──\n');
    if isfield(corr_tgt, 'L_C_agg') && ~isnan(corr_tgt.L_C_agg)
        fprintf('    corr(L_agg,C_agg):      %.4f   (labor vs NIPA PCE)\n', corr_tgt.L_C_agg);
    end
    if isfield(corr_tgt, 'I_C_agg') && ~isnan(corr_tgt.I_C_agg)
        fprintf('    corr(I_agg,C_agg):      %.4f   (investment vs NIPA PCE)\n', corr_tgt.I_C_agg);
    end
    fprintf('    ── Labor-TFP correlations (data only) ──\n');
    if isfield(corr_tgt, 'L_TFP_agg')
        fprintf('    corr(L_agg,TFP_agg):      %.4f   (VA-based aggregate TFP)\n', corr_tgt.L_TFP_agg);
    end
    if isfield(corr_tgt, 'L_TFP_sectoral_avg_vashare')
        fprintf('    avg corr(L_j,TFP_j):      %.4f   (VA-share weighted sectoral avg)\n', ...
            corr_tgt.L_TFP_sectoral_avg_vashare);
    end
    if isfield(corr_tgt, 'L_TFP_GO_agg')
        fprintf('    corr(L_agg,TFP_GO_agg):   %.4f   (GO-based aggregate TFP)\n', corr_tgt.L_TFP_GO_agg);
    end
    if isfield(corr_tgt, 'L_TFP_GO_sectoral_avg_goshare')
        fprintf('    avg corr(L_j,TFP_GO_j):   %.4f   (GO-share weighted sectoral avg)\n', ...
            corr_tgt.L_TFP_GO_sectoral_avg_goshare);
    end
    if isfield(corr_tgt, 'L_TFP_sectoral')
        fprintf('    sectoral vectors saved in EmpiricalTargets.correlations.*\n');
    end
end

end
