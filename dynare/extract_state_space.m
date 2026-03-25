function SolData = extract_state_space(oo_, M_, n_sectors, policies_ss)
    dim = n_sectors;
    inv = oo_.dr.inv_order_var;
    idx = get_variable_indices(n_sectors);

    %% Variable index ranges in Dynare ordering
    k_ind      = idx.k;
    a_ind      = idx.a;
    c_ind      = idx.c;
    l_ind      = idx.l;
    pk_ind     = idx.pk;
    pm_ind     = idx.pm;
    m_ind      = idx.m;
    mout_ind   = idx.mout;
    i_ind      = idx.i;
    iout_ind   = idx.iout;
    p_ind      = idx.p;
    q_ind      = idx.q;
    y_ind      = idx.y;
    c_util_ind = idx.c_util;
    l_util_ind = idx.l_util;
    cagg_ind   = idx.c_agg;
    lagg_ind   = idx.l_agg;
    gdpagg_ind = idx.gdp_agg;
    iagg_ind   = idx.i_agg;
    kagg_ind   = idx.k_agg;
    utility_intratemp_ind = idx.utility_intratemp;

    %% Map to Dynare internal ordering
    k_inv      = [inv(k_ind(1)),    inv(k_ind(2))];
    a_inv      = [inv(a_ind(1)),    inv(a_ind(2))];
    c_inv      = [inv(c_ind(1)),    inv(c_ind(2))];
    l_inv      = [inv(l_ind(1)),    inv(l_ind(2))];
    pk_inv     = [inv(pk_ind(1)),   inv(pk_ind(2))];
    pm_inv     = [inv(pm_ind(1)),   inv(pm_ind(2))];
    m_inv      = [inv(m_ind(1)),    inv(m_ind(2))];
    mout_inv   = [inv(mout_ind(1)), inv(mout_ind(2))];
    i_inv      = [inv(i_ind(1)),    inv(i_ind(2))];
    iout_inv   = [inv(iout_ind(1)), inv(iout_ind(2))];
    p_inv      = [inv(p_ind(1)),    inv(p_ind(2))];
    q_inv      = [inv(q_ind(1)),    inv(q_ind(2))];
    y_inv      = [inv(y_ind(1)),    inv(y_ind(2))];
    c_util_inv = inv(c_util_ind);
    l_util_inv = inv(l_util_ind);
    cagg_inv   = inv(cagg_ind);
    lagg_inv   = inv(lagg_ind);
    gdpagg_inv = inv(gdpagg_ind);
    iagg_inv   = inv(iagg_ind);
    kagg_inv   = inv(kagg_ind);
    utility_intratemp_inv = inv(utility_intratemp_ind);

    %% State transition: S(t) = A*S(t-1) + B*e(t)
    ghx = oo_.dr.ghx;
    ghu = oo_.dr.ghu;

    A = [ghx(k_inv(1):k_inv(2),:); ghx(a_inv(1):a_inv(2),:)];
    B = [ghu(k_inv(1):k_inv(2),:); ghu(a_inv(1):a_inv(2),:)];

    %% Policy functions: X(t) = C*S(t-1) + D*e(t)
    pol_rows_x = [c_inv; l_inv; pk_inv; pm_inv; m_inv; mout_inv; ...
                  i_inv; iout_inv; p_inv; q_inv; y_inv];
    agg_rows = [
        c_util_inv;
        l_util_inv;
        cagg_inv;
        lagg_inv;
        gdpagg_inv;
        iagg_inv;
        kagg_inv;
        utility_intratemp_inv
    ];

    C_blocks = cell(size(pol_rows_x, 1) + numel(agg_rows), 1);
    D_blocks = cell(size(pol_rows_x, 1) + numel(agg_rows), 1);
    bi = 0;
    for r = 1:size(pol_rows_x, 1)
        bi = bi + 1;
        C_blocks{bi} = ghx(pol_rows_x(r,1):pol_rows_x(r,2), :);
        D_blocks{bi} = ghu(pol_rows_x(r,1):pol_rows_x(r,2), :);
    end
    for r = 1:numel(agg_rows)
        bi = bi + 1;
        C_blocks{bi} = ghx(agg_rows(r), :);
        D_blocks{bi} = ghu(agg_rows(r), :);
    end

    C = vertcat(C_blocks{:});
    D = vertcat(D_blocks{:});

    assert(size(C, 1) == numel(policies_ss), ...
        'extract_state_space:PolicyDimensionMismatch', ...
        'StateSpace.C has %d rows but policies_ss has %d elements.', ...
        size(C, 1), numel(policies_ss));

    %% Pack output
    SolData = struct();
    SolData.A = A;
    SolData.B = B;
    SolData.C = C;
    SolData.D = D;
    SolData.k_ss = oo_.steady_state(1:dim);
    SolData.policies_ss = policies_ss;
    SolData.indices = struct('k', k_ind, 'a', a_ind, 'c', c_ind, 'l', l_ind, ...
        'pk', pk_ind, 'pm', pm_ind, 'm', m_ind, 'mout', mout_ind, ...
        'i', i_ind, 'iout', iout_ind, 'p', p_ind, 'q', q_ind, 'y', y_ind, ...
        'c_util', c_util_ind, 'l_util', l_util_ind, 'c_agg', cagg_ind, ...
        'l_agg', lagg_ind, 'gdp_agg', gdpagg_ind, 'i_agg', iagg_ind, ...
        'k_agg', kagg_ind, 'utility_intratemp', utility_intratemp_ind);
end
