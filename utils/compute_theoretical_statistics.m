function TheoStats = compute_theoretical_statistics(oo_, M_, policies_ss, n_sectors, endostates_ss)
% COMPUTE_THEORETICAL_STATISTICS Extract theoretical moments from the first-order Dynare object.

    %#ok<INUSD>
    TheoStats = struct();

    idx = get_variable_indices(n_sectors);
    n = n_sectors;

    if ~isfield(oo_, 'var') || isempty(oo_.var)
        warning('compute_theoretical_statistics:MissingVarCov', ...
            'oo_.var is not available. Cannot compute theoretical moments.');
        return;
    end

    var_cov = oo_.var;
    if size(var_cov, 1) < idx.n_dynare || size(var_cov, 2) < idx.n_dynare
        warning('compute_theoretical_statistics:VarCovSizeMismatch', ...
            'oo_.var has size %d x %d, but at least %d x %d is required.', ...
            size(var_cov, 1), size(var_cov, 2), idx.n_dynare, idx.n_dynare);
        return;
    end

    ss_of = @(range) policies_ss((range(1):range(2)) - idx.ss_offset);

    p_ss_log = ss_of(idx.p);
    q_ss_log = ss_of(idx.q);
    mout_ss_log = ss_of(idx.mout);
    l_ss_log = ss_of(idx.l);
    i_ss_log = ss_of(idx.i);
    y_ss_log = ss_of(idx.y);
    pk_ss_log = ss_of(idx.pk);

    p_ss = exp(p_ss_log(:));
    q_ss = exp(q_ss_log(:));
    mout_ss = exp(mout_ss_log(:));
    l_ss = exp(l_ss_log(:));
    i_ss = exp(i_ss_log(:));
    y_ss = exp(y_ss_log(:));
    pk_ss = exp(pk_ss_log(:));
    k_ss = exp(endostates_ss(1:n));
    k_ss = k_ss(:);

    cagg_ss_policy = policies_ss(idx.c_agg - idx.ss_offset);
    lagg_ss_policy = policies_ss(idx.l_agg - idx.ss_offset);
    gdp_ss_policy = policies_ss(idx.gdp_agg - idx.ss_offset);
    iagg_ss_policy = policies_ss(idx.i_agg - idx.ss_offset);
    kagg_ss_policy = policies_ss(idx.k_agg - idx.ss_offset);
    utility_ss_policy = policies_ss(idx.utility_intratemp - idx.ss_offset);

    cagg_ss = exp(cagg_ss_policy);
    lagg_ss = exp(lagg_ss_policy);
    gdp_ss = exp(gdp_ss_policy);
    iagg_ss = exp(iagg_ss_policy);
    kagg_ss = exp(kagg_ss_policy);

    va_weights = va_weights_from_y(y_ss);
    go_weights = q_ss / sum(q_ss);
    emp_weights = l_ss / sum(l_ss);
    inv_weights = (i_ss .* pk_ss) / sum(i_ss .* pk_ss);
    omega_Q = (p_ss .* q_ss) / sum(p_ss .* (q_ss - mout_ss));
    omega_Mout = (p_ss .* mout_ss) / sum(p_ss .* (q_ss - mout_ss));
    omega_K = (pk_ss .* k_ss) / sum(pk_ss .* k_ss);

    agg_specs = {
        'C',                  idx.c_agg,             cagg_ss_policy;
        'L',                  idx.l_agg,             lagg_ss_policy;
        'GDP',                idx.gdp_agg,           gdp_ss_policy;
        'I',                  idx.i_agg,             iagg_ss_policy;
        'K',                  idx.k_agg,             kagg_ss_policy;
        'utility_intratemp',  idx.utility_intratemp, utility_ss_policy;
    };

    aggregate_moments = struct();
    for i = 1:size(agg_specs, 1)
        agg_name = agg_specs{i, 1};
        agg_idx = agg_specs{i, 2};
        agg_ss_policy = agg_specs{i, 3};
        aggregate_moments.(agg_name) = build_gaussian_moments( ...
            direct_theoretical_mean(oo_, agg_idx, agg_ss_policy), ...
            safe_diag_std(var_cov, agg_idx));
    end

    sigma_C_agg = aggregate_moments.C.std;
    sigma_L_agg = aggregate_moments.L.std;
    sigma_VA_agg = aggregate_moments.GDP.std;
    sigma_I_agg = aggregate_moments.I.std;
    sigma_K_agg = aggregate_moments.K.std;
    sigma_utility_intratemp_agg = aggregate_moments.utility_intratemp.std;

    Var_a = var_cov(idx.a(1):idx.a(2), idx.a(1):idx.a(2));
    Var_c = var_cov(idx.c(1):idx.c(2), idx.c(1):idx.c(2));
    Var_y = var_cov(idx.y(1):idx.y(2), idx.y(1):idx.y(2));
    Var_i = var_cov(idx.i(1):idx.i(2), idx.i(1):idx.i(2));
    Var_l = var_cov(idx.l(1):idx.l(2), idx.l(1):idx.l(2));
    Var_q = var_cov(idx.q(1):idx.q(2), idx.q(1):idx.q(2));
    Cov_q_gdp = var_cov(idx.q(1):idx.q(2), idx.gdp_agg);
    Cov_l_a = var_cov(idx.l(1):idx.l(2), idx.a(1):idx.a(2));
    Cov_lagg_a = var_cov(idx.l_agg, idx.a(1):idx.a(2));

    corr_CI_agg = safe_corr_from_cov(var_cov(idx.c_agg, idx.i_agg), sigma_C_agg, sigma_I_agg);
    corr_L_C_agg = safe_corr_from_cov(var_cov(idx.l_agg, idx.c_agg), sigma_L_agg, sigma_C_agg);

    sigma_A_VA_agg = sqrt(max(va_weights' * Var_a * va_weights, 0));
    sigma_A_GO_agg = sqrt(max(go_weights' * Var_a * go_weights, 0));
    corr_L_TFP_agg = safe_corr_from_cov(Cov_lagg_a * va_weights, sigma_L_agg, sigma_A_VA_agg);
    corr_L_TFP_GO_agg = safe_corr_from_cov(Cov_lagg_a * go_weights, sigma_L_agg, sigma_A_GO_agg);

    sigma_A_sectoral = sqrt(max(diag(Var_a), 0))';
    sigma_L_sectoral = sqrt(max(diag(Var_l), 0))';
    corr_L_TFP_sectoral = NaN(1, n);
    for j = 1:n
        corr_L_TFP_sectoral(j) = safe_corr_from_cov(var_cov(idx.l(1) + j - 1, idx.a(1) + j - 1), ...
            sigma_L_sectoral(j), sigma_A_sectoral(j));
    end

    sigma_VA_sectoral = sqrt(max(diag(Var_y), 0))';
    sigma_I_sectoral = sqrt(max(diag(Var_i), 0))';
    sigma_VA_avg = weighted_mean_ignore_nan(sigma_VA_sectoral, va_weights);
    sigma_L_avg = weighted_mean_ignore_nan(sigma_L_sectoral, va_weights);
    sigma_I_avg = weighted_mean_ignore_nan(sigma_I_sectoral, va_weights);
    sigma_L_avg_empweighted = weighted_mean_ignore_nan(sigma_L_sectoral, emp_weights);
    sigma_I_avg_invweighted = weighted_mean_ignore_nan(sigma_I_sectoral, inv_weights);

    corr_matrix_C = corr_matrix_from_cov(Var_c);
    corr_matrix_VA = corr_matrix_from_cov(Var_y);
    corr_matrix_L = corr_matrix_from_cov(Var_l);
    corr_matrix_I = corr_matrix_from_cov(Var_i);
    avg_pairwise_corr_C = mean_upper_triangle_ignore_nan(corr_matrix_C);
    avg_pairwise_corr_VA = mean_upper_triangle_ignore_nan(corr_matrix_VA);
    avg_pairwise_corr_L = mean_upper_triangle_ignore_nan(corr_matrix_L);
    avg_pairwise_corr_I = mean_upper_triangle_ignore_nan(corr_matrix_I);

    sigma_Domar_sectoral = NaN(1, n);
    for j = 1:n
        var_domar_j = Var_q(j, j) + var_cov(idx.gdp_agg, idx.gdp_agg) - 2 * Cov_q_gdp(j);
        sigma_Domar_sectoral(j) = sqrt(max(var_domar_j, 0));
    end
    sigma_Domar_avg = weighted_mean_ignore_nan(sigma_Domar_sectoral, go_weights);

    TheoStats.aggregate_definition = 'exact_logdev_to_deterministic_ss';
    TheoStats.variable_convention = 'Y=primary_factors; VA=P_ss*(Q-Mout); GDP=aggregate_VA';
    TheoStats.moment_origin = 'direct_aggregate_policy_tail';
    TheoStats.aggregate_moments = aggregate_moments;
    TheoStats.share_C = cagg_ss / gdp_ss;
    TheoStats.share_I = iagg_ss / gdp_ss;

    TheoStats.sigma_C_agg = sigma_C_agg;
    TheoStats.sigma_L_agg = sigma_L_agg;
    TheoStats.sigma_L_hc_agg = sigma_L_agg;
    TheoStats.sigma_VA_agg = sigma_VA_agg;
    TheoStats.sigma_I_agg = sigma_I_agg;
    TheoStats.sigma_K_agg = sigma_K_agg;
    TheoStats.sigma_utility_intratemp_agg = sigma_utility_intratemp_agg;
    TheoStats.corr_C_I = corr_CI_agg;
    TheoStats.corr_CI_agg = corr_CI_agg;
    TheoStats.corr_L_C_agg = corr_L_C_agg;
    TheoStats.corr_I_C_agg = corr_CI_agg;
    TheoStats.rho_C_agg = direct_theoretical_autocorr(oo_, idx.c_agg);
    TheoStats.rho_L_agg = direct_theoretical_autocorr(oo_, idx.l_agg);
    TheoStats.rho_VA_agg = direct_theoretical_autocorr(oo_, idx.gdp_agg);
    TheoStats.rho_I_agg = direct_theoretical_autocorr(oo_, idx.i_agg);
    TheoStats.rho_K_agg = direct_theoretical_autocorr(oo_, idx.k_agg);
    TheoStats.rho_utility_intratemp_agg = direct_theoretical_autocorr(oo_, idx.utility_intratemp);

    TheoStats.sigma_VA_avg = sigma_VA_avg;
    TheoStats.sigma_VA_sectoral = sigma_VA_sectoral;
    TheoStats.sigma_L_avg = sigma_L_avg;
    TheoStats.sigma_I_avg = sigma_I_avg;
    TheoStats.sigma_L_avg_empweighted = sigma_L_avg_empweighted;
    TheoStats.sigma_I_avg_invweighted = sigma_I_avg_invweighted;
    TheoStats.sigma_Domar_avg = sigma_Domar_avg;
    TheoStats.sigma_Domar_avg_legacy = sigma_Domar_avg;
    TheoStats.sigma_L_sectoral = sigma_L_sectoral;
    TheoStats.sigma_I_sectoral = sigma_I_sectoral;
    TheoStats.sigma_Domar_sectoral = sigma_Domar_sectoral;
    TheoStats.sigma_Domar_sectoral_legacy = sigma_Domar_sectoral;
    TheoStats.domar_definition = 'log_fixed_price_gross_output_share_in_GDP';
    TheoStats.domar_average_weight_definition = 'legacy_normalized_gross_output_weights';
    TheoStats.corr_matrix_C = corr_matrix_C;
    TheoStats.corr_matrix_VA = corr_matrix_VA;
    TheoStats.corr_matrix_L = corr_matrix_L;
    TheoStats.corr_matrix_I = corr_matrix_I;
    TheoStats.avg_pairwise_corr_C = avg_pairwise_corr_C;
    TheoStats.avg_pairwise_corr_VA = avg_pairwise_corr_VA;
    TheoStats.avg_pairwise_corr_L = avg_pairwise_corr_L;
    TheoStats.avg_pairwise_corr_I = avg_pairwise_corr_I;
    TheoStats.corr_L_TFP_agg = corr_L_TFP_agg;
    TheoStats.corr_L_TFP_GO_agg = corr_L_TFP_GO_agg;
    TheoStats.corr_L_TFP_sectoral = corr_L_TFP_sectoral;
    TheoStats.corr_L_TFP_sectoral_avg_vashare = weighted_mean_ignore_nan(corr_L_TFP_sectoral, va_weights);
    TheoStats.corr_L_TFP_sectoral_avg_empshare = weighted_mean_ignore_nan(corr_L_TFP_sectoral, emp_weights);
    TheoStats.va_weights = va_weights;
    TheoStats.go_weights = go_weights;
    TheoStats.emp_weights = emp_weights;
    TheoStats.inv_weights = inv_weights;
    TheoStats.omega_Q = omega_Q;
    TheoStats.omega_Mout = omega_Mout;
    TheoStats.omega_K = omega_K;

    if size(var_cov, 1) >= 5
        TheoStats.sigma_C_legacy = sqrt(max(var_cov(1, 1), 0));
        TheoStats.sigma_L_legacy = sqrt(max(var_cov(2, 2), 0));
        TheoStats.sigma_VA_legacy = sqrt(max(var_cov(3, 3), 0));
        TheoStats.sigma_I_legacy = sqrt(max(var_cov(4, 4), 0));
        TheoStats.sigma_K_legacy = sqrt(max(var_cov(5, 5), 0));

        legacy_std = sqrt(max(diag(var_cov), 0));
        legacy_denom = legacy_std * legacy_std';
        legacy_corr = NaN(size(var_cov));
        valid_mask = legacy_denom > 0;
        legacy_corr(valid_mask) = var_cov(valid_mask) ./ legacy_denom(valid_mask);
        TheoStats.var_cov_agg_legacy = var_cov;
        TheoStats.corr_matrix_agg_legacy = legacy_corr;
    end

    if isfield(oo_, 'autocorr') && ~isempty(oo_.autocorr) && numel(oo_.autocorr) >= 1
        autocorr_lag1 = oo_.autocorr{1};
        if size(autocorr_lag1, 1) >= 5
            TheoStats.rho_C_agg_legacy = autocorr_lag1(1, 1);
            TheoStats.rho_L_agg_legacy = autocorr_lag1(2, 2);
            TheoStats.rho_VA_agg_legacy = autocorr_lag1(3, 3);
            TheoStats.rho_I_agg_legacy = autocorr_lag1(4, 4);
            TheoStats.rho_K_agg_legacy = autocorr_lag1(5, 5);
        end
    end
end

function va_weights = va_weights_from_y(y_ss)
    va_weights = y_ss / sum(y_ss);
end

function moments = build_gaussian_moments(mu, sigma)
    moments = struct( ...
        'mean', mu, ...
        'std', sigma, ...
        'skewness', 0.0, ...
        'kurtosis', 3.0, ...
        'excess_kurtosis', 0.0);
end

function mu = direct_theoretical_mean(oo_, dynare_idx, steady_state_policy_value)
    if isfield(oo_, 'mean') && ~isempty(oo_.mean) && numel(oo_.mean) >= dynare_idx
        mu = oo_.mean(dynare_idx) - steady_state_policy_value;
    else
        mu = 0.0;
    end
end

function rho = direct_theoretical_autocorr(oo_, dynare_idx)
    rho = NaN;
    if ~isfield(oo_, 'autocorr') || isempty(oo_.autocorr) || numel(oo_.autocorr) < 1
        return;
    end

    autocorr_lag1 = oo_.autocorr{1};
    if size(autocorr_lag1, 1) >= dynare_idx && size(autocorr_lag1, 2) >= dynare_idx
        rho = autocorr_lag1(dynare_idx, dynare_idx);
    end
end

function sigma = safe_diag_std(cov_matrix, idx_value)
    sigma = sqrt(max(cov_matrix(idx_value, idx_value), 0));
end

function rho = safe_corr_from_cov(cov_xy, sigma_x, sigma_y)
    if ~isfinite(cov_xy) || ~isfinite(sigma_x) || ~isfinite(sigma_y) || sigma_x <= 0 || sigma_y <= 0
        rho = NaN;
        return;
    end
    rho = cov_xy / (sigma_x * sigma_y);
end

function corr_matrix = corr_matrix_from_cov(cov_matrix)
    n = size(cov_matrix, 1);
    corr_matrix = NaN(n, n);
    sigma = sqrt(max(diag(cov_matrix), 0));
    for i = 1:n
        if sigma(i) > 0
            corr_matrix(i, i) = 1;
        end
        for j = i+1:n
            corr_ij = safe_corr_from_cov(cov_matrix(i, j), sigma(i), sigma(j));
            corr_matrix(i, j) = corr_ij;
            corr_matrix(j, i) = corr_ij;
        end
    end
end

function avg = mean_upper_triangle_ignore_nan(matrix)
    if isempty(matrix)
        avg = NaN;
        return;
    end
    upper_tri = triu(true(size(matrix)), 1);
    values = matrix(upper_tri);
    values = values(isfinite(values));
    if isempty(values)
        avg = NaN;
        return;
    end
    avg = mean(values);
end

function avg = weighted_mean_ignore_nan(values, weights)
    values = values(:);
    weights = weights(:);
    mask = isfinite(values) & isfinite(weights);
    if ~any(mask)
        avg = NaN;
        return;
    end

    values = values(mask);
    weights = weights(mask);
    weights = weights / sum(weights);
    avg = sum(weights .* values);
end

