function SolData = extract_state_space(oo_, M_, n_sectors, policies_ss)
    dim = n_sectors;
    inv = oo_.dr.inv_order_var;

    %% Variable index ranges in Dynare ordering
    k_ind      = [1,        dim];
    a_ind      = [dim+1,    2*dim];
    c_ind      = [2*dim+1,  3*dim];
    l_ind      = [3*dim+1,  4*dim];
    pk_ind     = [4*dim+1,  5*dim];
    pm_ind     = [5*dim+1,  6*dim];
    m_ind      = [6*dim+1,  7*dim];
    mout_ind   = [7*dim+1,  8*dim];
    i_ind      = [8*dim+1,  9*dim];
    iout_ind   = [9*dim+1,  10*dim];
    p_ind      = [10*dim+1, 11*dim];
    q_ind      = [11*dim+1, 12*dim];
    y_ind      = [12*dim+1, 13*dim];
    cagg_ind   = 13*dim+1;
    lagg_ind   = 13*dim+2;
    yagg_ind   = 13*dim+3;
    iagg_ind   = 13*dim+4;
    magg_ind   = 13*dim+5;

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
    cagg_inv   = inv(cagg_ind);
    lagg_inv   = inv(lagg_ind);
    yagg_inv   = inv(yagg_ind);
    iagg_inv   = inv(iagg_ind);
    magg_inv   = inv(magg_ind);

    %% State transition: S(t) = A*S(t-1) + B*e(t)
    ghx = oo_.dr.ghx;
    ghu = oo_.dr.ghu;

    A = [ghx(k_inv(1):k_inv(2),:); ghx(a_inv(1):a_inv(2),:)];
    B = [ghu(k_inv(1):k_inv(2),:); ghu(a_inv(1):a_inv(2),:)];

    %% Policy functions: X(t) = C*S(t-1) + D*e(t)
    pol_rows_x = [c_inv; l_inv; pk_inv; pm_inv; m_inv; mout_inv; ...
                  i_inv; iout_inv; p_inv; q_inv; y_inv];
    agg_rows = [cagg_inv; lagg_inv; yagg_inv; iagg_inv; magg_inv];

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
        'cagg', cagg_ind, 'lagg', lagg_ind, 'yagg', yagg_ind, ...
        'iagg', iagg_ind, 'magg', magg_ind);
end
