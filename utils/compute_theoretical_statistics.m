function TheoStats = compute_theoretical_statistics(oo_, M_, policies_ss, n_sectors, endostates_ss)
% COMPUTE_THEORETICAL_STATISTICS Extract theoretical moments from Dynare solution
%
% Computes theoretical variances for EXPENDITURE-BASED aggregates using
% the state-space representation. This replaces the old approach that used
% Dynare's legacy aggregates (cagg, yagg, iagg, magg).
%
% NEW AGGREGATES (expenditure-based, matching model specification):
%   C_constP = Σ_j P̄_j × C_j        (consumption expenditure)
%   I_constP = Σ_j P̄_j × I_j^out    (investment expenditure)  
%   GDP_constP = Σ_j P̄_j × (Q_j - M_j^out) = C_constP + I_constP
%
% The variance of log-deviations is computed using first-order approximation:
%   var(log X_constP) ≈ ω' × Var(x̃) × ω
% where ω is the vector of steady-state expenditure shares.
%
% INPUTS:
%   oo_         - Dynare output structure (contains oo_.dr for state-space)
%   M_          - Dynare model structure
%   policies_ss - Steady state policies
%   n_sectors   - Number of sectors
%
% OUTPUTS:
%   TheoStats - Structure with theoretical moments for expenditure-based aggregates

    TheoStats = struct();
    
    % Get variable indices
    idx = get_variable_indices(n_sectors);
    n = n_sectors;
    
    %% ===== COMPUTE THEORETICAL VARIANCE FROM STATE-SPACE =====
    % The state-space is: x_{t+1} = A x_t + B ε_t, y_t = C x_t + D ε_t
    % Unconditional variance: Σ_x = A Σ_x A' + B Σ_ε B' (solve Lyapunov)
    %                         Σ_y = C Σ_x C' + D Σ_ε D'
    
    if ~isfield(oo_, 'dr') || ~isfield(oo_.dr, 'ghx')
        warning('compute_theoretical_statistics:NoDr', ...
            'oo_.dr not available. Cannot compute theoretical moments.');
        return;
    end
    
    % Get shock covariance matrix
    Sigma_eps = M_.Sigma_e;  % n_shocks x n_shocks
    
    % Get state transition matrices from Dynare
    % ghx: policy function w.r.t. states (n_vars x n_states)
    % ghu: policy function w.r.t. shocks (n_vars x n_shocks)
    ghx = oo_.dr.ghx;
    ghu = oo_.dr.ghu;
    
    % Get indices for states (k, a) and convert to DR order
    n_states = 2 * n;  % k and a
    
    % State indices in decision rule order
    k_dr = oo_.dr.inv_order_var(1:n);
    a_dr = oo_.dr.inv_order_var(n+1:2*n);
    state_dr_idx = [k_dr; a_dr];
    
    % Extract state transition for states only
    A_full = ghx(state_dr_idx, :);
    B_full = ghu(state_dr_idx, :);
    
    % Solve Lyapunov equation for state variance: Σ_x = A Σ_x A' + B Σ_ε B'
    BQB = B_full * Sigma_eps * B_full';
    try
        Sigma_x = dlyap(A_full, BQB);
    catch
        % Fallback: iterative solution
        Sigma_x = BQB;
        for iter = 1:500
            Sigma_x_new = A_full * Sigma_x * A_full' + BQB;
            if max(abs(Sigma_x_new(:) - Sigma_x(:))) < 1e-10
                break;
            end
            Sigma_x = Sigma_x_new;
        end
    end
    
    % Extract indices for sectoral variables in DR order
    c_dr = oo_.dr.inv_order_var(idx.c(1):idx.c(2));
    y_dr = oo_.dr.inv_order_var(idx.y(1):idx.y(2));
    i_dr = oo_.dr.inv_order_var(idx.i(1):idx.i(2));
    iout_dr = oo_.dr.inv_order_var(idx.iout(1):idx.iout(2));
    q_dr = oo_.dr.inv_order_var(idx.q(1):idx.q(2));
    mout_dr = oo_.dr.inv_order_var(idx.mout(1):idx.mout(2));
    l_dr = oo_.dr.inv_order_var(idx.l(1):idx.l(2));
    
    % Get policy matrices for these variables
    C_c = ghx(c_dr, :);      D_c = ghu(c_dr, :);
    C_y = ghx(y_dr, :);      D_y = ghu(y_dr, :);
    C_i = ghx(i_dr, :);      D_i = ghu(i_dr, :);
    C_iout = ghx(iout_dr, :); D_iout = ghu(iout_dr, :);
    C_q = ghx(q_dr, :);       D_q = ghu(q_dr, :);
    C_mout = ghx(mout_dr, :); D_mout = ghu(mout_dr, :);
    C_l = ghx(l_dr, :);       D_l = ghu(l_dr, :);
    
    % Compute variance of sectoral variables
    % Var(y) = C Σ_x C' + D Σ_ε D'
    Var_c = C_c * Sigma_x * C_c' + D_c * Sigma_eps * D_c';
    Var_y = C_y * Sigma_x * C_y' + D_y * Sigma_eps * D_y';
    Var_i = C_i * Sigma_x * C_i' + D_i * Sigma_eps * D_i';
    Var_iout = C_iout * Sigma_x * C_iout' + D_iout * Sigma_eps * D_iout';
    Var_q = C_q * Sigma_x * C_q' + D_q * Sigma_eps * D_q';
    Var_mout = C_mout * Sigma_x * C_mout' + D_mout * Sigma_eps * D_mout';
    Var_l = C_l * Sigma_x * C_l' + D_l * Sigma_eps * D_l';
    
    % Covariance between q and mout (needed for GDP = Q - M^out)
    Cov_q_mout = C_q * Sigma_x * C_mout' + D_q * Sigma_eps * D_mout';
    Cov_c_l = C_c * Sigma_x * C_l' + D_c * Sigma_eps * D_l';
    
    %% ===== COMPUTE EXPENDITURE-BASED AGGREGATE VARIANCES =====
    % For log-linearized variables, var(Σ ω_j x̃_j) = ω' Var(x̃) ω
    % where ω_j = (P̄_j × X_j^ss) / (Σ_k P̄_k × X_k^ss) is expenditure share
    
    % Get steady-state levels and prices
    p_ss_idx = (idx.p(1):idx.p(2)) - idx.ss_offset;
    c_ss_idx = (idx.c(1):idx.c(2)) - idx.ss_offset;
    iout_ss_idx = (idx.iout(1):idx.iout(2)) - idx.ss_offset;
    q_ss_idx = (idx.q(1):idx.q(2)) - idx.ss_offset;
    mout_ss_idx = (idx.mout(1):idx.mout(2)) - idx.ss_offset;
    l_ss_idx = (idx.l(1):idx.l(2)) - idx.ss_offset;
    y_ss_idx = (idx.y(1):idx.y(2)) - idx.ss_offset;
    pk_ss_idx = (idx.pk(1):idx.pk(2)) - idx.ss_offset;
    
    p_ss = exp(policies_ss(p_ss_idx));  p_ss = p_ss(:);
    c_ss = exp(policies_ss(c_ss_idx));  c_ss = c_ss(:);
    iout_ss = exp(policies_ss(iout_ss_idx));  iout_ss = iout_ss(:);
    q_ss = exp(policies_ss(q_ss_idx));  q_ss = q_ss(:);
    mout_ss = exp(policies_ss(mout_ss_idx));  mout_ss = mout_ss(:);
    l_ss = exp(policies_ss(l_ss_idx));  l_ss = l_ss(:);
    y_ss = exp(policies_ss(y_ss_idx));  y_ss = y_ss(:);
    pk_ss = exp(policies_ss(pk_ss_idx));  pk_ss = pk_ss(:);
    k_ss = exp(endostates_ss(1:n));  k_ss = k_ss(:);
    va_weights = va_weights_from_y(y_ss);
    
    % Expenditure weights for consumption
    % log(C_constP) ≈ Σ_j ω_j^C × c̃_j where ω_j^C = (P̄_j × C_j^ss) / C_constP_ss
    C_constP_ss = sum(p_ss .* c_ss);
    omega_C = (p_ss .* c_ss) / C_constP_ss;
    
    % Expenditure weights for investment
    I_constP_ss = sum(p_ss .* iout_ss);
    omega_I = (p_ss .* iout_ss) / I_constP_ss;
    
    % Expenditure weights for GDP = Q - M^out
    GDP_constP_ss = sum(p_ss .* (q_ss - mout_ss));
    omega_Q = (p_ss .* q_ss) / GDP_constP_ss;
    omega_Mout = (p_ss .* mout_ss) / GDP_constP_ss;
    
    % Labor weights (simple sum for headcount)
    L_hc_ss = sum(l_ss);
    omega_L = l_ss / L_hc_ss;

    % Domar weights for TFP aggregation
    domar_weights = omega_Q;
    
    % Capital weights (value-weighted: pk_ss × K_ss)
    K_val_ss = sum(pk_ss .* k_ss);
    omega_K = (pk_ss .* k_ss) / K_val_ss;
    
    % Theoretical variances of expenditure-based aggregates
    % var(log C_constP) = ω_C' × Var(c̃) × ω_C
    sigma_C_agg = sqrt(omega_C' * Var_c * omega_C);
    sigma_I_agg = sqrt(omega_I' * Var_iout * omega_I);
    
    % var(log GDP) = ω_Q' Var(q̃) ω_Q + ω_M' Var(m̃) ω_M - 2 ω_Q' Cov(q̃,m̃) ω_M
    % Note: GDP = Q - M^out, so covariance term is subtracted
    var_GDP = omega_Q' * Var_q * omega_Q + omega_Mout' * Var_mout * omega_Mout ...
              - 2 * omega_Q' * Cov_q_mout * omega_Mout;
    sigma_VA_agg = sqrt(max(var_GDP, 0));  % Ensure non-negative
    
    % Labor aggregate (headcount)
    sigma_L_agg = sqrt(omega_L' * Var_l * omega_L);
    cov_L_C = omega_C' * Cov_c_l * omega_L;
    corr_L_C_agg = safe_corr_from_cov(cov_L_C, sigma_C_agg, sigma_L_agg);

    % TFP state covariance and labor-TFP correlations
    Var_a = Sigma_x(n+1:2*n, n+1:2*n);
    Cov_l_a = C_l * Sigma_x(:, n+1:2*n);

    var_A_VA = va_weights' * Var_a * va_weights;
    sigma_A_VA_agg = sqrt(max(var_A_VA, 0));
    cov_L_A_VA = omega_L' * Cov_l_a * va_weights;
    corr_L_TFP_agg = safe_corr_from_cov(cov_L_A_VA, sigma_L_agg, sigma_A_VA_agg);

    var_A_GO = domar_weights' * Var_a * domar_weights;
    sigma_A_GO_agg = sqrt(max(var_A_GO, 0));
    cov_L_A_GO = omega_L' * Cov_l_a * domar_weights;
    corr_L_TFP_GO_agg = safe_corr_from_cov(cov_L_A_GO, sigma_L_agg, sigma_A_GO_agg);

    corr_L_TFP_sectoral = NaN(1, n);
    sigma_A_sectoral = sqrt(max(diag(Var_a), 0))';
    for j = 1:n
        corr_L_TFP_sectoral(j) = safe_corr_from_cov(Cov_l_a(j, j), sqrt(max(Var_l(j, j), 0)), sigma_A_sectoral(j));
    end
    corr_L_TFP_sectoral_avg_vashare = weighted_mean_ignore_nan(corr_L_TFP_sectoral, va_weights);
    corr_L_TFP_sectoral_avg_empshare = weighted_mean_ignore_nan(corr_L_TFP_sectoral, omega_L');

    sigma_VA_sectoral = sqrt(max(diag(Var_y), 0))';
    sigma_VA_avg = weighted_mean_ignore_nan(sigma_VA_sectoral, va_weights);

    corr_matrix_C = corr_matrix_from_cov(Var_c);
    corr_matrix_VA = corr_matrix_from_cov(Var_y);
    corr_matrix_L = corr_matrix_from_cov(Var_l);
    corr_matrix_I = corr_matrix_from_cov(Var_i);

    avg_pairwise_corr_C = mean_upper_triangle_ignore_nan(corr_matrix_C);
    avg_pairwise_corr_VA = mean_upper_triangle_ignore_nan(corr_matrix_VA);
    avg_pairwise_corr_L = mean_upper_triangle_ignore_nan(corr_matrix_L);
    avg_pairwise_corr_I = mean_upper_triangle_ignore_nan(corr_matrix_I);
    
    % Capital aggregate (k is the first n states in Sigma_x)
    Var_k = Sigma_x(1:n, 1:n);
    sigma_K_agg = sqrt(omega_K' * Var_k * omega_K);
    
    % Store in output structure
    TheoStats.sigma_C_agg = sigma_C_agg;
    TheoStats.sigma_L_agg = sigma_L_agg;
    TheoStats.sigma_L_hc_agg = sigma_L_agg;
    TheoStats.sigma_VA_agg = sigma_VA_agg;
    TheoStats.sigma_I_agg = sigma_I_agg;
    TheoStats.sigma_K_agg = sigma_K_agg;
    TheoStats.sigma_VA_avg = sigma_VA_avg;
    TheoStats.sigma_VA_sectoral = sigma_VA_sectoral;
    TheoStats.avg_pairwise_corr_C = avg_pairwise_corr_C;
    TheoStats.avg_pairwise_corr_VA = avg_pairwise_corr_VA;
    TheoStats.avg_pairwise_corr_L = avg_pairwise_corr_L;
    TheoStats.avg_pairwise_corr_I = avg_pairwise_corr_I;
    TheoStats.corr_matrix_C = corr_matrix_C;
    TheoStats.corr_matrix_VA = corr_matrix_VA;
    TheoStats.corr_matrix_L = corr_matrix_L;
    TheoStats.corr_matrix_I = corr_matrix_I;
    TheoStats.corr_L_TFP_agg = corr_L_TFP_agg;
    TheoStats.corr_L_TFP_GO_agg = corr_L_TFP_GO_agg;
    TheoStats.corr_L_TFP_sectoral = corr_L_TFP_sectoral;
    TheoStats.corr_L_TFP_sectoral_avg_vashare = corr_L_TFP_sectoral_avg_vashare;
    TheoStats.corr_L_TFP_sectoral_avg_empshare = corr_L_TFP_sectoral_avg_empshare;
    
    % Aggregates reported directly by Dynare's stoch_simul call
    if isfield(oo_, 'var') && ~isempty(oo_.var) && size(oo_.var, 1) >= 5
        var_cov = oo_.var;
        TheoStats.sigma_C_legacy = sqrt(var_cov(1, 1));  % c_agg
        TheoStats.sigma_L_legacy = sqrt(var_cov(2, 2));  % l_agg
        TheoStats.sigma_VA_legacy = sqrt(var_cov(3, 3)); % gdp_agg
        TheoStats.sigma_I_legacy = sqrt(var_cov(4, 4));  % i_agg
        TheoStats.sigma_K_legacy = sqrt(var_cov(5, 5));  % k_agg
        
        % Full variance-covariance matrix of legacy aggregates
        TheoStats.var_cov_agg_legacy = var_cov;
        
        % Correlations for the reported stoch_simul aggregate block
        std_vec = sqrt(diag(var_cov));
        corr_matrix = var_cov ./ (std_vec * std_vec');
        TheoStats.corr_matrix_agg_legacy = corr_matrix;
    end
    
    % Compute expenditure-based aggregate covariances for correlation
    % Cov(C, I) for the new aggregates
    Cov_c_iout = C_c * Sigma_x * C_iout' + D_c * Sigma_eps * D_iout';
    cov_CI = omega_C' * Cov_c_iout * omega_I;
    TheoStats.corr_C_I = safe_corr_from_cov(cov_CI, sigma_C_agg, sigma_I_agg);
    TheoStats.corr_L_C_agg = corr_L_C_agg;
    TheoStats.corr_I_C_agg = TheoStats.corr_C_I;
    
    % Autocorrelations (legacy, from Dynare)
    if isfield(oo_, 'autocorr') && ~isempty(oo_.autocorr) && numel(oo_.autocorr) >= 1
        autocorr_lag1 = oo_.autocorr{1};
        
        % These are autocorrelations of the reported stoch_simul aggregate block
        if size(autocorr_lag1, 1) >= 5
            TheoStats.rho_C_agg_legacy = autocorr_lag1(1, 1);
            TheoStats.rho_L_agg_legacy = autocorr_lag1(2, 2);
            TheoStats.rho_VA_agg_legacy = autocorr_lag1(3, 3);
            TheoStats.rho_I_agg_legacy = autocorr_lag1(4, 4);
            TheoStats.rho_K_agg_legacy = autocorr_lag1(5, 5);
        end
    end
    
    % Store expenditure weights
    TheoStats.omega_C = omega_C;
    TheoStats.omega_I = omega_I;
    TheoStats.omega_Q = omega_Q;
    TheoStats.omega_Mout = omega_Mout;
    TheoStats.omega_K = omega_K;
    
    % Store VA weights from steady state
    TheoStats.va_weights = va_weights;
    
end

function va_weights = va_weights_from_y(y_ss)
    va_weights = y_ss / sum(y_ss);
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

