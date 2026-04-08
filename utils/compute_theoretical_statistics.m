function TheoStats = compute_theoretical_statistics(oo_, M_, policies_ss, n_sectors, endostates_ss)
% COMPUTE_THEORETICAL_STATISTICS Build first-order TheoStats from Dynare moments.
%
% This routine intentionally keeps the theoretical workflow minimal:
% - std moments are read directly from oo_.var
% - lag-1 autocorrelation is read from oo_.autocorr{1}
% - skewness/kurtosis are Gaussian-by-construction for first-order objects
%
% All required aggregate fields are always returned so downstream code can
% rely on a stable TheoStats schema.

    %#ok<INUSD>
    idx = get_variable_indices(n_sectors);

    aggregate_specs = {
        'C',                 'c_agg',              idx.c_agg;
        'L',                 'l_agg',              idx.l_agg;
        'GDP',               'gdp_agg',            idx.gdp_agg;
        'I',                 'i_agg',              idx.i_agg;
        'K',                 'k_agg',              idx.k_agg;
        'utility_intratemp', 'utility_intratemp',  idx.utility_intratemp;
    };

    if isfield(oo_, 'var') && ~isempty(oo_.var)
        var_cov = real(oo_.var);
    else
        var_cov = [];
        warning('compute_theoretical_statistics:MissingVarCov', ...
            'oo_.var is not available. TheoStats moments are returned as NaN.');
    end

    aggregate_moments = struct();
    sigma_values = struct();
    rho_values = struct();
    dynare_indices = struct();
    missing_sigma_names = {};

    for i = 1:size(aggregate_specs, 1)
        agg_name = aggregate_specs{i, 1};
        dynare_name = aggregate_specs{i, 2};
        fallback_idx = aggregate_specs{i, 3};

        dynare_idx = resolve_endo_index(M_, dynare_name, fallback_idx);
        dynare_indices.(agg_name) = dynare_idx;

        sigma_val = extract_std_from_cov(var_cov, dynare_idx);
        rho_val = extract_lag1_autocorr(oo_, dynare_idx);

        aggregate_moments.(agg_name) = build_gaussian_moments(0.0, sigma_val);
        sigma_values.(agg_name) = sigma_val;
        rho_values.(agg_name) = rho_val;
        if ~isfinite(sigma_val)
            missing_sigma_names{end + 1} = agg_name; %#ok<AGROW>
        end
    end

    if ~isempty(missing_sigma_names)
        warning('compute_theoretical_statistics:MissingAggregateStd', ...
            'Could not read theoretical std from oo_.var for aggregates: %s', ...
            strjoin(missing_sigma_names, ', '));
    end

    corr_CI_agg = extract_corr_from_cov(var_cov, dynare_indices.C, dynare_indices.I);
    corr_L_C_agg = extract_corr_from_cov(var_cov, dynare_indices.L, dynare_indices.C);

    cagg_ss = exp(policies_ss(idx.c_agg - idx.ss_offset));
    iagg_ss = exp(policies_ss(idx.i_agg - idx.ss_offset));
    gdp_ss = exp(policies_ss(idx.gdp_agg - idx.ss_offset));

    TheoStats = struct();
    TheoStats.aggregate_definition = 'exact_logdev_to_deterministic_ss';
    TheoStats.variable_convention = 'Y=primary_factors; VA=P_ss*(Q-Mout); GDP=aggregate_VA';
    TheoStats.moment_origin = 'dynare_oo_var';
    TheoStats.aggregate_moments = aggregate_moments;
    TheoStats.share_C = safe_ratio(cagg_ss, gdp_ss);
    TheoStats.share_I = safe_ratio(iagg_ss, gdp_ss);

    TheoStats.sigma_C_agg = sigma_values.C;
    TheoStats.sigma_L_agg = sigma_values.L;
    TheoStats.sigma_L_hc_agg = sigma_values.L;
    TheoStats.sigma_VA_agg = sigma_values.GDP;
    TheoStats.sigma_I_agg = sigma_values.I;
    TheoStats.sigma_K_agg = sigma_values.K;
    TheoStats.sigma_utility_intratemp_agg = sigma_values.utility_intratemp;

    TheoStats.corr_C_I = corr_CI_agg;
    TheoStats.corr_CI_agg = corr_CI_agg;
    TheoStats.corr_L_C_agg = corr_L_C_agg;
    TheoStats.corr_I_C_agg = corr_CI_agg;

    TheoStats.rho_C_agg = rho_values.C;
    TheoStats.rho_L_agg = rho_values.L;
    TheoStats.rho_VA_agg = rho_values.GDP;
    TheoStats.rho_I_agg = rho_values.I;
    TheoStats.rho_K_agg = rho_values.K;
    TheoStats.rho_utility_intratemp_agg = rho_values.utility_intratemp;

    TheoStats.endo_index_lookup = dynare_indices;
    TheoStats = attach_compatibility_placeholders(TheoStats, n_sectors, var_cov);
end

function idx = resolve_endo_index(M_, target_name, fallback_idx)
idx = fallback_idx;

if ~isfield(M_, 'endo_names') || isempty(M_.endo_names)
    return;
end

raw_names = M_.endo_names;
if ischar(raw_names)
    endo_names = cellstr(raw_names);
elseif isstring(raw_names)
    endo_names = cellstr(raw_names);
elseif iscell(raw_names)
    endo_names = raw_names;
else
    return;
end

for i = 1:numel(endo_names)
    endo_names{i} = strtrim(char(endo_names{i}));
end

hit = find(strcmp(endo_names, target_name), 1, 'first');
if ~isempty(hit)
    idx = hit;
end
end

function sigma = extract_std_from_cov(var_cov, dynare_idx)
sigma = NaN;
if isempty(var_cov) || ~isfinite(dynare_idx) || dynare_idx < 1
    return;
end
if size(var_cov, 1) < dynare_idx || size(var_cov, 2) < dynare_idx
    return;
end
sigma = sqrt(max(real(var_cov(dynare_idx, dynare_idx)), 0));
end

function rho = extract_lag1_autocorr(oo_, dynare_idx)
rho = NaN;
if ~isfinite(dynare_idx) || dynare_idx < 1
    return;
end
if ~isfield(oo_, 'autocorr') || isempty(oo_.autocorr) || numel(oo_.autocorr) < 1
    return;
end

autocorr_lag1 = oo_.autocorr{1};
if size(autocorr_lag1, 1) < dynare_idx || size(autocorr_lag1, 2) < dynare_idx
    return;
end
rho = real(autocorr_lag1(dynare_idx, dynare_idx));
end

function rho = extract_corr_from_cov(var_cov, idx_x, idx_y)
rho = NaN;
if isempty(var_cov)
    return;
end
if ~isfinite(idx_x) || ~isfinite(idx_y) || idx_x < 1 || idx_y < 1
    return;
end
if size(var_cov, 1) < idx_x || size(var_cov, 2) < idx_y
    return;
end

sigma_x = extract_std_from_cov(var_cov, idx_x);
sigma_y = extract_std_from_cov(var_cov, idx_y);
if ~isfinite(sigma_x) || ~isfinite(sigma_y) || sigma_x <= 0 || sigma_y <= 0
    return;
end

cov_xy = real(var_cov(idx_x, idx_y));
rho = cov_xy / (sigma_x * sigma_y);
end

function ratio = safe_ratio(num, den)
if ~isfinite(num) || ~isfinite(den) || abs(den) < 1e-12
    ratio = NaN;
else
    ratio = num / den;
end
end

function moments = build_gaussian_moments(mu, sigma)
moments = struct( ...
    'mean', mu, ...
    'std', sigma, ...
    'skewness', 0.0, ...
    'kurtosis', 3.0, ...
    'excess_kurtosis', 0.0);
end

function TheoStats = attach_compatibility_placeholders(TheoStats, n_sectors, var_cov)
TheoStats.sigma_VA_avg = NaN;
TheoStats.sigma_VA_sectoral = NaN(1, n_sectors);
TheoStats.sigma_L_avg = NaN;
TheoStats.sigma_I_avg = NaN;
TheoStats.sigma_L_avg_empweighted = NaN;
TheoStats.sigma_I_avg_invweighted = NaN;
TheoStats.sigma_Domar_avg = NaN;
TheoStats.sigma_Domar_avg_legacy = NaN;
TheoStats.sigma_L_sectoral = NaN(1, n_sectors);
TheoStats.sigma_I_sectoral = NaN(1, n_sectors);
TheoStats.sigma_Domar_sectoral = NaN(1, n_sectors);
TheoStats.sigma_Domar_sectoral_legacy = NaN(1, n_sectors);
TheoStats.domar_definition = 'log_fixed_price_gross_output_share_in_GDP';
TheoStats.domar_average_weight_definition = 'legacy_normalized_gross_output_weights';

TheoStats.corr_matrix_C = [];
TheoStats.corr_matrix_VA = [];
TheoStats.corr_matrix_L = [];
TheoStats.corr_matrix_I = [];
TheoStats.avg_pairwise_corr_C = NaN;
TheoStats.avg_pairwise_corr_VA = NaN;
TheoStats.avg_pairwise_corr_L = NaN;
TheoStats.avg_pairwise_corr_I = NaN;
TheoStats.corr_L_TFP_agg = NaN;
TheoStats.corr_L_TFP_GO_agg = NaN;
TheoStats.corr_L_TFP_sectoral = NaN(1, n_sectors);
TheoStats.corr_L_TFP_sectoral_avg_vashare = NaN;
TheoStats.corr_L_TFP_sectoral_avg_empshare = NaN;

TheoStats.va_weights = NaN(1, n_sectors);
TheoStats.go_weights = NaN(1, n_sectors);
TheoStats.emp_weights = NaN(1, n_sectors);
TheoStats.inv_weights = NaN(1, n_sectors);
TheoStats.omega_Q = NaN(1, n_sectors);
TheoStats.omega_Mout = NaN(1, n_sectors);
TheoStats.omega_K = NaN(1, n_sectors);

TheoStats.sigma_C_legacy = TheoStats.sigma_C_agg;
TheoStats.sigma_L_legacy = TheoStats.sigma_L_agg;
TheoStats.sigma_VA_legacy = TheoStats.sigma_VA_agg;
TheoStats.sigma_I_legacy = TheoStats.sigma_I_agg;
TheoStats.sigma_K_legacy = TheoStats.sigma_K_agg;
TheoStats.rho_C_agg_legacy = TheoStats.rho_C_agg;
TheoStats.rho_L_agg_legacy = TheoStats.rho_L_agg;
TheoStats.rho_VA_agg_legacy = TheoStats.rho_VA_agg;
TheoStats.rho_I_agg_legacy = TheoStats.rho_I_agg;
TheoStats.rho_K_agg_legacy = TheoStats.rho_K_agg;

if isempty(var_cov)
    TheoStats.var_cov_agg_legacy = [];
    TheoStats.corr_matrix_agg_legacy = [];
else
    TheoStats.var_cov_agg_legacy = var_cov;
    TheoStats.corr_matrix_agg_legacy = corr_matrix_from_cov(var_cov);
end
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
        if sigma(i) > 0 && sigma(j) > 0
            corr_ij = real(cov_matrix(i, j)) / (sigma(i) * sigma(j));
            corr_matrix(i, j) = corr_ij;
            corr_matrix(j, i) = corr_ij;
        end
    end
end
end

