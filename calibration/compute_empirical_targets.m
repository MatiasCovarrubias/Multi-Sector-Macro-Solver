function empirical_targets = compute_empirical_targets(VA_raw, EMP_raw, InvRaw, VAn, Invn, GO_raw, VA_for_domar, Cons_agg, TFP_raw, TFP_GO_raw, GOn_raw, Cons_sectoral)

    if nargin < 9
        TFP_raw = [];
    end
    if nargin < 10
        TFP_GO_raw = [];
    end
    if nargin < 11
        GOn_raw = [];
    end
    if nargin < 12
        Cons_sectoral = [];
    end

    hp_lambda = 100;
    epsilon = 1e-10;

    VA_clean = max(VA_raw, epsilon);
    EMP_clean = max(EMP_raw, epsilon);
    Inv_clean = max(InvRaw, epsilon);
    VAn_clean = max(VAn, epsilon);
    Invn_clean = max(Invn, epsilon);
    if ~isempty(Cons_sectoral)
        Cons_sectoral_clean = max(Cons_sectoral, epsilon);
    else
        Cons_sectoral_clean = [];
    end

    [yrnum, n_sectors] = size(VA_clean);

    %% Weights from nominal shares
    VAsh = VAn_clean ./ repmat(sum(VAn_clean, 2), 1, n_sectors);
    Invsh = Invn_clean ./ repmat(sum(Invn_clean, 2), 1, n_sectors);

    va_weights = mean(VAsh, 1);
    va_weights = va_weights / sum(va_weights);

    EMPsh = EMP_clean ./ repmat(sum(EMP_clean, 2), 1, n_sectors);
    emp_weights = mean(EMPsh, 1);
    emp_weights = emp_weights / sum(emp_weights);

    inv_weights = mean(Invsh, 1);
    inv_weights = inv_weights / sum(inv_weights);

    %% Tornqvist aggregation
    aggVA = ones(yrnum, 1);
    for t = 1:yrnum-1
        gr = 0;
        for j = 1:n_sectors
            gr = gr + 0.5 * (VAsh(t,j) + VAsh(t+1,j)) * (log(VA_clean(t+1,j)) - log(VA_clean(t,j)));
        end
        aggVA(t+1) = aggVA(t) * exp(gr);
    end

    aggInv = ones(yrnum, 1);
    for t = 1:yrnum-1
        gr = 0;
        for j = 1:n_sectors
            gr = gr + 0.5 * (Invsh(t,j) + Invsh(t+1,j)) * (log(Inv_clean(t+1,j)) - log(Inv_clean(t,j)));
        end
        aggInv(t+1) = aggInv(t) * exp(gr);
    end

    aggEMP = sum(EMP_clean, 2);

    %% Aggregate volatilities (HP-filtered)
    cyc_VA_agg = hp_cycle_from_levels(aggVA, hp_lambda, epsilon);
    sigma_VA_agg = std(cyc_VA_agg);

    cyc_Inv_agg = hp_cycle_from_levels(aggInv, hp_lambda, epsilon);
    sigma_I_agg = std(cyc_Inv_agg);

    cyc_EMP_agg = hp_cycle_from_levels(aggEMP, hp_lambda, epsilon);
    sigma_L_agg = std(cyc_EMP_agg);

    if ~isempty(Cons_agg)
        cyc_Cons_agg = hp_cycle_from_levels(Cons_agg, hp_lambda, epsilon);
        sigma_C_agg = std(cyc_Cons_agg);
    else
        sigma_C_agg = NaN;
    end

    %% Sectoral volatilities and comovement
    cyc_VA_sectoral = hp_cycles_from_levels_matrix(VA_clean, hp_lambda, epsilon);
    sigma_VA_sectoral = std(cyc_VA_sectoral, 0, 1);
    sigma_VA_avg = sum(va_weights .* sigma_VA_sectoral);
    [corr_matrix_VA, avg_pairwise_corr_VA] = safe_corr_matrix_columns(cyc_VA_sectoral);

    cyc_L_sectoral = hp_cycles_from_levels_matrix(EMP_clean, hp_lambda, epsilon);
    sigma_L_sectoral = std(cyc_L_sectoral, 0, 1);
    sigma_L_avg = sum(va_weights .* sigma_L_sectoral);
    sigma_L_avg_empweighted = sum(emp_weights .* sigma_L_sectoral);
    [corr_matrix_L, avg_pairwise_corr_L] = safe_corr_matrix_columns(cyc_L_sectoral);

    cyc_I_sectoral = hp_cycles_from_levels_matrix(Inv_clean, hp_lambda, epsilon);
    sigma_I_sectoral = std(cyc_I_sectoral, 0, 1);
    sigma_I_avg = sum(va_weights .* sigma_I_sectoral);
    sigma_I_avg_invweighted = sum(inv_weights .* sigma_I_sectoral);
    [corr_matrix_I, avg_pairwise_corr_I] = safe_corr_matrix_columns(cyc_I_sectoral);

    if ~isempty(Cons_sectoral_clean)
        cyc_C_sectoral = hp_cycles_from_levels_matrix(Cons_sectoral_clean, hp_lambda, epsilon);
        [corr_matrix_C, avg_pairwise_corr_C] = safe_corr_matrix_columns(cyc_C_sectoral);
    else
        corr_matrix_C = [];
        avg_pairwise_corr_C = NaN;
    end

    %% Domar weight volatility
    GO_clean = max(GO_raw, epsilon);
    VA_domar_clean = max(VA_for_domar, epsilon);
    [~, n_sectors_domar] = size(GO_clean);

    go_weights = mean(GO_clean, 1) / sum(mean(GO_clean, 1));
    Domar = GO_clean ./ repmat(sum(VA_domar_clean, 2), 1, n_sectors_domar);

    sigma_Domar_sectoral = zeros(1, n_sectors_domar);
    for i = 1:n_sectors_domar
        sigma_Domar_sectoral(i) = std(hp_cycle_from_levels(Domar(:, i), hp_lambda, epsilon));
    end
    sigma_Domar_avg = sum(go_weights .* sigma_Domar_sectoral);

    %% Labor-TFP correlations
    correlations = struct();
    correlations.method = 'corr(HP-filter(log x), HP-filter(log y))';
    correlations.hp_lambda = hp_lambda;
    correlations.L_C_agg = NaN;
    correlations.I_C_agg = NaN;

    if ~isempty(Cons_agg)
        correlations.L_C_agg = safe_corr(cyc_EMP_agg, cyc_Cons_agg);
        correlations.I_C_agg = safe_corr(cyc_Inv_agg, cyc_Cons_agg);
    end

    if ~isempty(TFP_raw)
        TFP_clean = max(TFP_raw, epsilon);
        aggTFP = tornqvist_chain_index(TFP_clean, VAsh, epsilon);
        L_TFP_sectoral = compute_sector_correlations(EMP_clean, TFP_clean, hp_lambda, epsilon);

        correlations.L_TFP_sectoral = L_TFP_sectoral;
        correlations.L_TFP_agg = safe_corr(cyc_EMP_agg, hp_cycle_from_levels(aggTFP, hp_lambda, epsilon));
        correlations.L_TFP_sectoral_avg_vashare = weighted_mean_ignore_nan(L_TFP_sectoral, va_weights);
        correlations.L_TFP_sectoral_avg_empshare = weighted_mean_ignore_nan(L_TFP_sectoral, emp_weights);
    end

    if ~isempty(TFP_GO_raw)
        TFP_GO_clean = max(TFP_GO_raw, epsilon);
        if isempty(GOn_raw)
            GOn_clean = GO_clean;
        else
            GOn_clean = max(GOn_raw, epsilon);
        end
        go_corr_weights = mean(GOn_clean ./ repmat(sum(GOn_clean, 2), 1, n_sectors), 1);
        go_corr_weights = go_corr_weights / sum(go_corr_weights);
        domar_growth_weights = GOn_clean ./ repmat(sum(VAn_clean, 2), 1, n_sectors);
        aggTFP_GO = tornqvist_chain_index(TFP_GO_clean, domar_growth_weights, epsilon);
        L_TFP_GO_sectoral = compute_sector_correlations(EMP_clean, TFP_GO_clean, hp_lambda, epsilon);

        correlations.L_TFP_GO_sectoral = L_TFP_GO_sectoral;
        correlations.L_TFP_GO_agg = safe_corr(cyc_EMP_agg, hp_cycle_from_levels(aggTFP_GO, hp_lambda, epsilon));
        correlations.L_TFP_GO_sectoral_avg_goshare = weighted_mean_ignore_nan(L_TFP_GO_sectoral, go_corr_weights);
        correlations.L_TFP_GO_sectoral_avg_empshare = weighted_mean_ignore_nan(L_TFP_GO_sectoral, emp_weights);
    end

    %% Pack output
    empirical_targets = struct();
    empirical_targets.hp_lambda = hp_lambda;
    empirical_targets.va_weights = va_weights;
    empirical_targets.go_weights = go_weights;
    empirical_targets.emp_weights = emp_weights;
    empirical_targets.inv_weights = inv_weights;
    empirical_targets.aggregation_method = 'tornqvist';

    empirical_targets.sigma_VA_agg = sigma_VA_agg;
    empirical_targets.sigma_C_agg = sigma_C_agg;
    empirical_targets.sigma_L_agg = sigma_L_agg;
    empirical_targets.sigma_I_agg = sigma_I_agg;

    empirical_targets.sigma_VA_avg = sigma_VA_avg;
    empirical_targets.sigma_L_avg = sigma_L_avg;
    empirical_targets.sigma_I_avg = sigma_I_avg;
    empirical_targets.sigma_L_avg_empweighted = sigma_L_avg_empweighted;
    empirical_targets.sigma_I_avg_invweighted = sigma_I_avg_invweighted;
    empirical_targets.sigma_Domar_avg = sigma_Domar_avg;

    empirical_targets.avg_pairwise_corr_C = avg_pairwise_corr_C;
    empirical_targets.avg_pairwise_corr_VA = avg_pairwise_corr_VA;
    empirical_targets.avg_pairwise_corr_L = avg_pairwise_corr_L;
    empirical_targets.avg_pairwise_corr_I = avg_pairwise_corr_I;

    empirical_targets.sigma_VA_sectoral = sigma_VA_sectoral;
    empirical_targets.sigma_L_sectoral = sigma_L_sectoral;
    empirical_targets.sigma_I_sectoral = sigma_I_sectoral;
    empirical_targets.sigma_Domar_sectoral = sigma_Domar_sectoral;
    empirical_targets.corr_matrix_C = corr_matrix_C;
    empirical_targets.corr_matrix_VA = corr_matrix_VA;
    empirical_targets.corr_matrix_L = corr_matrix_L;
    empirical_targets.corr_matrix_I = corr_matrix_I;
    empirical_targets.correlations = correlations;
end

function cyc = hp_cycle_from_levels(series, hp_lambda, epsilon)
    log_series = log(max(series(:), epsilon));
    [~, cyc] = hpfilter(log_series, 'Smoothing', hp_lambda);
end

function cyc_mat = hp_cycles_from_levels_matrix(levels, hp_lambda, epsilon)
    [T, n] = size(levels);
    cyc_mat = NaN(T, n);
    for i = 1:n
        cyc_mat(:, i) = hp_cycle_from_levels(levels(:, i), hp_lambda, epsilon);
    end
end

function agg = tornqvist_chain_index(levels, weights, epsilon)
    [T, n] = size(levels);
    agg = ones(T, 1);

    for t = 1:T-1
        growth_t = 0;
        for j = 1:n
            growth_t = growth_t + 0.5 * (weights(t, j) + weights(t + 1, j)) * ...
                (log(max(levels(t + 1, j), epsilon)) - log(max(levels(t, j), epsilon)));
        end
        agg(t + 1) = agg(t) * exp(growth_t);
    end
end

function rho_vec = compute_sector_correlations(lhs_levels, rhs_levels, hp_lambda, epsilon)
    n_sectors = size(lhs_levels, 2);
    rho_vec = NaN(1, n_sectors);

    for i = 1:n_sectors
        lhs_cycle = hp_cycle_from_levels(lhs_levels(:, i), hp_lambda, epsilon);
        rhs_cycle = hp_cycle_from_levels(rhs_levels(:, i), hp_lambda, epsilon);
        rho_vec(i) = safe_corr(lhs_cycle, rhs_cycle);
    end
end

function rho = safe_corr(x, y)
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

function [corr_matrix, avg] = safe_corr_matrix_columns(data)
    n = size(data, 2);
    corr_matrix = NaN(n, n);
    for i = 1:n
        corr_matrix(i, i) = 1;
        for j = i+1:n
            rho_ij = safe_corr(data(:, i), data(:, j));
            corr_matrix(i, j) = rho_ij;
            corr_matrix(j, i) = rho_ij;
        end
    end
    avg = mean_upper_triangle_ignore_nan(corr_matrix);
end

function avg = weighted_mean_ignore_nan(values, weights)
    mask = isfinite(values) & isfinite(weights);
    if ~any(mask)
        avg = NaN;
        return;
    end

    weights = weights(mask);
    values = values(mask);
    weight_sum = sum(weights);
    if weight_sum == 0
        avg = NaN;
        return;
    end

    weights = weights / weight_sum;
    avg = sum(weights .* values);
end
