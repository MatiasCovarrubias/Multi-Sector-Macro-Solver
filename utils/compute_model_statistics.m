function ModelStats = compute_model_statistics(dynare_simul, idx, policies_ss, n_sectors, endostates_ss, sample_label, sample_window)

    if nargin < 6 || isempty(sample_label)
        sample_label = 'Simulation';
    end
    if nargin < 7 || isempty(sample_window)
        sample_window = 'shocks_simul';
    end

    %% Extract simulated series in model-space log levels
    c_simul    = dynare_simul(idx.c(1):idx.c(2), :);
    y_simul    = dynare_simul(idx.y(1):idx.y(2), :);
    l_simul    = dynare_simul(idx.l(1):idx.l(2), :);
    i_simul    = dynare_simul(idx.i(1):idx.i(2), :);
    a_simul    = dynare_simul(idx.a(1):idx.a(2), :);
    iout_simul = dynare_simul(idx.iout(1):idx.iout(2), :);
    mout_simul = dynare_simul(idx.mout(1):idx.mout(2), :);
    p_simul    = dynare_simul(idx.p(1):idx.p(2), :);
    q_simul    = dynare_simul(idx.q(1):idx.q(2), :);
    k_simul    = dynare_simul(idx.k(1):idx.k(2), :);
    cagg_simul = dynare_simul(idx.c_agg, :);
    yagg_simul = dynare_simul(idx.gdp_agg, :);
    lagg_simul = dynare_simul(idx.l_agg, :);
    iagg_simul = dynare_simul(idx.i_agg, :);
    kagg_simul = dynare_simul(idx.k_agg, :);
    magg_simul = NaN(1, size(dynare_simul, 2));
    utility_intratemp_simul = dynare_simul(idx.utility_intratemp, :);

    %% Steady-state levels and weights
    ss_of = @(range) policies_ss((range(1):range(2)) - idx.ss_offset);

    c_ss_log = ss_of(idx.c);
    y_ss_log = ss_of(idx.y);
    l_ss_log = ss_of(idx.l);
    i_ss_log = ss_of(idx.i);
    p_ss_log = ss_of(idx.p);
    q_ss_log = ss_of(idx.q);
    mout_ss_log = ss_of(idx.mout);
    y_ss = exp(ss_of(idx.y));  y_ss = y_ss(:);
    q_ss = exp(ss_of(idx.q));  q_ss = q_ss(:);
    p_ss = exp(ss_of(idx.p));  p_ss = p_ss(:);
    l_ss = exp(ss_of(idx.l));  l_ss = l_ss(:);
    i_ss = exp(ss_of(idx.i));  i_ss = i_ss(:);
    pk_ss = exp(ss_of(idx.pk)); pk_ss = pk_ss(:);

    k_ss = exp(endostates_ss(1:n_sectors));  k_ss = k_ss(:);

    cagg_ss = exp(policies_ss(idx.c_agg - idx.ss_offset));
    yagg_ss = exp(policies_ss(idx.gdp_agg - idx.ss_offset));
    lagg_ss = exp(policies_ss(idx.l_agg - idx.ss_offset));
    iagg_ss = exp(policies_ss(idx.i_agg - idx.ss_offset));
    kagg_ss = exp(policies_ss(idx.k_agg - idx.ss_offset));
    utility_intratemp_ss_log = policies_ss(idx.utility_intratemp - idx.ss_offset);

    va_sector_ss = exp(p_ss_log) .* (exp(q_ss_log) - exp(mout_ss_log));
    va_weights = va_sector_ss' / sum(va_sector_ss);
    go_weights = q_ss' / sum(q_ss);
    emp_weights = l_ss' / sum(l_ss);
    inv_weights = (i_ss .* pk_ss)' / sum(i_ss .* pk_ss);

    c_logdev = c_simul - c_ss_log;
    y_logdev = y_simul - y_ss_log;
    l_logdev = l_simul - l_ss_log;
    i_logdev = i_simul - i_ss_log;

    %% Aggregate log deviations from Dynare aggregate endogenous variables
    C_logdev = cagg_simul - log(cagg_ss);
    I_logdev = iagg_simul - log(iagg_ss);
    GDP_logdev = yagg_simul - log(yagg_ss);
    L_hc_logdev = lagg_simul - log(lagg_ss);
    K_logdev = kagg_simul - log(kagg_ss);
    utility_intratemp_logdev = utility_intratemp_simul - utility_intratemp_ss_log;

    %% Aggregate intermediate use is still reconstructed because it is not a model aggregate endogenous variable
    M_levels = sum(p_ss .* exp(mout_simul), 1);
    M_ss = sum(p_ss .* exp(mout_ss_log));
    M_logdev = log(M_levels) - log(M_ss);

    agg_ss = struct();
    agg_ss.C = cagg_ss;
    agg_ss.I = iagg_ss;
    agg_ss.GDP = yagg_ss;
    agg_ss.L = lagg_ss;
    agg_ss.M = M_ss;
    agg_ss.K = kagg_ss;
    agg_ss.share_C = agg_ss.C / agg_ss.GDP;
    agg_ss.share_I = agg_ss.I / agg_ss.GDP;

    sigma_C_agg  = std(C_logdev);
    sigma_I_agg  = std(I_logdev);
    sigma_VA_agg = std(GDP_logdev);
    sigma_K_agg  = std(K_logdev);
    sigma_L_hc_agg = std(L_hc_logdev);
    sigma_M_agg = std(M_logdev);
    sigma_utility_intratemp_agg = std(utility_intratemp_logdev);

    aggregate_moments = struct();
    aggregate_moments.C = compute_univariate_moments(C_logdev);
    aggregate_moments.I = compute_univariate_moments(I_logdev);
    aggregate_moments.GDP = compute_univariate_moments(GDP_logdev);
    aggregate_moments.L = compute_univariate_moments(L_hc_logdev);
    aggregate_moments.M = compute_univariate_moments(M_logdev);
    aggregate_moments.K = compute_univariate_moments(K_logdev);
    aggregate_moments.utility_intratemp = compute_univariate_moments(utility_intratemp_logdev);

    share_C = agg_ss.share_C;
    share_I = agg_ss.share_I;
    corr_CI = safe_corr(C_logdev, I_logdev);

    %% Legacy aggregates
    sigma_C_pref_agg = std(cagg_simul);
    sigma_I_ces_agg = std(iagg_simul);
    sigma_Y_primaryfactor_agg = std(yagg_simul);
    sigma_L_legacy_agg = std(lagg_simul);
    sigma_M_legacy_agg = std(magg_simul);

    %% Headcount labor aggregate (exact)
    corr_L_C_agg = safe_corr(L_hc_logdev, C_logdev);
    corr_I_C_agg = safe_corr(I_logdev, C_logdev);

    %% Labor-TFP correlations
    A_VA_logdev = va_weights * a_simul;
    omega_Q = (p_ss .* q_ss) / agg_ss.GDP;
    A_GO_logdev = omega_Q' * a_simul;

    corr_L_TFP_sectoral = NaN(1, n_sectors);
    for j = 1:n_sectors
        corr_L_TFP_sectoral(j) = safe_corr(l_simul(j, :), a_simul(j, :));
    end
    corr_L_TFP_agg = safe_corr(L_hc_logdev, A_VA_logdev);
    corr_L_TFP_GO_agg = safe_corr(L_hc_logdev, A_GO_logdev);
    corr_L_TFP_sectoral_avg_vashare = weighted_mean_ignore_nan(corr_L_TFP_sectoral, va_weights);
    corr_L_TFP_sectoral_avg_empshare = weighted_mean_ignore_nan(corr_L_TFP_sectoral, emp_weights);

    %% GDP autocorrelation
    if numel(GDP_logdev) >= 2
        rho_VA_agg = safe_corr(GDP_logdev(1:end-1), GDP_logdev(2:end));
    else
        rho_VA_agg = NaN;
    end

    %% Sectoral statistics
    va_sector_levels = p_ss .* (exp(q_simul) - exp(mout_simul));
    va_sector_logdev = log(va_sector_levels) - log(va_sector_ss);
    [corr_matrix_C, avg_pairwise_corr_C] = safe_corr_matrix_rows(c_logdev);
    [corr_matrix_VA, avg_pairwise_corr_VA] = safe_corr_matrix_rows(va_sector_logdev);
    [corr_matrix_L, avg_pairwise_corr_L] = safe_corr_matrix_rows(l_logdev);
    [corr_matrix_I, avg_pairwise_corr_I] = safe_corr_matrix_rows(i_logdev);

    sigma_VA_sectoral = std(va_sector_logdev, 0, 2)';
    sigma_L_sectoral = std(l_logdev, 0, 2)';
    sigma_I_sectoral = std(i_logdev, 0, 2)';
    sigma_VA_avg = sum(va_weights .* sigma_VA_sectoral);
    sigma_L_avg = sum(va_weights .* sigma_L_sectoral);
    sigma_I_avg = sum(va_weights .* sigma_I_sectoral);

    sigma_L_avg_empweighted = sum(emp_weights .* sigma_L_sectoral);
    sigma_I_avg_invweighted = sum(inv_weights .* sigma_I_sectoral);

    %% Domar weight volatility
    domar_simul = q_simul - repmat(yagg_simul, n_sectors, 1);
    sigma_Domar_sectoral = std(domar_simul, 0, 2)';
    sigma_Domar_avg = sum(go_weights .* sigma_Domar_sectoral);

    %% Pack output
    ModelStats = struct();
    ModelStats.sigma_VA_agg = sigma_VA_agg;
    ModelStats.sigma_C_agg = sigma_C_agg;
    ModelStats.sigma_L_agg = sigma_L_hc_agg;
    ModelStats.sigma_L_hc_agg = sigma_L_hc_agg;
    ModelStats.sigma_I_agg = sigma_I_agg;
    ModelStats.sigma_K_agg = sigma_K_agg;
    ModelStats.sigma_M_agg = sigma_M_agg;
    ModelStats.sigma_utility_intratemp_agg = sigma_utility_intratemp_agg;
    ModelStats.sigma_L_legacy_agg = sigma_L_legacy_agg;
    ModelStats.sigma_M_legacy_agg = sigma_M_legacy_agg;
    ModelStats.sigma_C_pref_agg = sigma_C_pref_agg;
    ModelStats.sigma_I_ces_agg = sigma_I_ces_agg;
    ModelStats.sigma_Y_primaryfactor_agg = sigma_Y_primaryfactor_agg;
    ModelStats.aggregate_definition = 'exact_logdev_to_deterministic_ss';
    ModelStats.sample_window = sample_window;
    ModelStats.aggregate_moments = aggregate_moments;
    ModelStats.share_C = share_C;
    ModelStats.share_I = share_I;
    ModelStats.corr_CI_agg = corr_CI;
    ModelStats.corr_L_C_agg = corr_L_C_agg;
    ModelStats.corr_I_C_agg = corr_I_C_agg;
    ModelStats.rho_VA_agg = rho_VA_agg;
    ModelStats.avg_pairwise_corr_C = avg_pairwise_corr_C;
    ModelStats.avg_pairwise_corr_VA = avg_pairwise_corr_VA;
    ModelStats.avg_pairwise_corr_L = avg_pairwise_corr_L;
    ModelStats.avg_pairwise_corr_I = avg_pairwise_corr_I;
    ModelStats.sigma_VA_avg = sigma_VA_avg;
    ModelStats.sigma_VA_sectoral = sigma_VA_sectoral;
    ModelStats.sigma_L_avg = sigma_L_avg;
    ModelStats.sigma_I_avg = sigma_I_avg;
    ModelStats.sigma_L_avg_empweighted = sigma_L_avg_empweighted;
    ModelStats.sigma_I_avg_invweighted = sigma_I_avg_invweighted;
    ModelStats.sigma_Domar_avg = sigma_Domar_avg;
    ModelStats.sigma_Domar_avg_legacy = sigma_Domar_avg;
    ModelStats.corr_matrix_C = corr_matrix_C;
    ModelStats.sigma_L_sectoral = sigma_L_sectoral;
    ModelStats.sigma_I_sectoral = sigma_I_sectoral;
    ModelStats.sigma_Domar_sectoral = sigma_Domar_sectoral;
    ModelStats.sigma_Domar_sectoral_legacy = sigma_Domar_sectoral;
    ModelStats.domar_definition = 'legacy_go_minus_yagg';
    ModelStats.corr_matrix_VA = corr_matrix_VA;
    ModelStats.corr_matrix_L = corr_matrix_L;
    ModelStats.corr_matrix_I = corr_matrix_I;
    ModelStats.corr_L_TFP_agg = corr_L_TFP_agg;
    ModelStats.corr_L_TFP_GO_agg = corr_L_TFP_GO_agg;
    ModelStats.corr_L_TFP_sectoral = corr_L_TFP_sectoral;
    ModelStats.corr_L_TFP_sectoral_avg_vashare = corr_L_TFP_sectoral_avg_vashare;
    ModelStats.corr_L_TFP_sectoral_avg_empshare = corr_L_TFP_sectoral_avg_empshare;
    ModelStats.va_weights = va_weights;
    ModelStats.go_weights = go_weights;
    ModelStats.emp_weights = emp_weights;
    ModelStats.inv_weights = inv_weights;
end

function moments = compute_univariate_moments(x)
    x = x(:);
    x = x(isfinite(x));
    moments = struct('mean', NaN, 'std', NaN, 'skewness', NaN, 'kurtosis', NaN);
    if isempty(x)
        return;
    end

    mu = mean(x);
    centered = x - mu;
    sigma = std(x);

    moments.mean = mu;
    moments.std = sigma;

    if sigma == 0
        return;
    end

    z = centered / sigma;
    moments.skewness = mean(z .^ 3);
    moments.kurtosis = mean(z .^ 4);
end

function rho = safe_corr(x, y)
    x = x(:);
    y = y(:);
    mask = isfinite(x) & isfinite(y);
    if nnz(mask) < 2
        rho = NaN;
        return;
    end

    x = x(mask);
    y = y(mask);
    if std(x) == 0 || std(y) == 0
        rho = NaN;
        return;
    end

    rho_mat = corrcoef(x, y);
    rho = rho_mat(1, 2);
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

function [corr_matrix, avg] = safe_corr_matrix_rows(data)
    n = size(data, 1);
    corr_matrix = NaN(n, n);
    for i = 1:n
        corr_matrix(i, i) = 1;
        for j = i+1:n
            rho_ij = safe_corr(data(i, :), data(j, :));
            corr_matrix(i, j) = rho_ij;
            corr_matrix(j, i) = rho_ij;
        end
    end

    upper_tri = triu(true(n), 1);
    values = corr_matrix(upper_tri);
    values = values(isfinite(values));
    if isempty(values)
        avg = NaN;
    else
        avg = mean(values);
    end
end
