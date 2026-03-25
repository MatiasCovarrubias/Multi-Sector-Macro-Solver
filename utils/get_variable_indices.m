function idx = get_variable_indices(n)
% GET_VARIABLE_INDICES Returns variable index ranges for the model
%
% Single source of truth for variable indexing. Uses Dynare ordering
% (states first, then policies) as the canonical form.
%
% INPUTS:
%   n - Number of sectors (typically 37)
%
% OUTPUTS:
%   idx - Structure with index ranges:
%
%   STATES:
%     idx.k     - Capital: [1, n]
%     idx.a     - TFP: [n+1, 2*n]
%
%   POLICIES (11 sectoral blocks + 8 aggregates):
%     idx.c     - Consumption: [2*n+1, 3*n]
%     idx.l     - Labor: [3*n+1, 4*n]
%     idx.pk    - Capital price: [4*n+1, 5*n]
%     idx.pm    - Intermediate price: [5*n+1, 6*n]
%     idx.m     - Intermediate input: [6*n+1, 7*n]
%     idx.mout  - Intermediate output: [7*n+1, 8*n]
%     idx.i     - Investment: [8*n+1, 9*n]
%     idx.iout  - Investment output: [9*n+1, 10*n]
%     idx.p     - Price: [10*n+1, 11*n]
%     idx.q     - Gross output: [11*n+1, 12*n]
%     idx.y     - Value added: [12*n+1, 13*n]
%
%   AGGREGATES:
%     idx.c_util            - Utility consumption aggregate: 13*n+1
%     idx.l_util            - Utility labor aggregate: 13*n+2
%     idx.c_agg             - Expenditure consumption aggregate: 13*n+3
%     idx.l_agg             - Headcount labor aggregate: 13*n+4
%     idx.gdp_agg           - Expenditure GDP aggregate: 13*n+5
%     idx.i_agg             - Expenditure investment aggregate: 13*n+6
%     idx.k_agg             - Deterministic-price capital aggregate: 13*n+7
%     idx.utility_intratemp - Intratemporal utility: 13*n+8
%
%   OFFSET for policies_ss:
%     idx.ss_offset = 2*n (subtract from Dynare index to get policies_ss index)
%
%   Example: To get y from policies_ss:
%     policies_ss(idx.y(1) - idx.ss_offset : idx.y(2) - idx.ss_offset)

    idx = struct();
    
    % States (Dynare indices 1 to 2n)
    idx.k = [1, n];
    idx.a = [n+1, 2*n];
    
    % Policies (Dynare indices 2n+1 onwards)
    idx.c = [2*n+1, 3*n];
    idx.l = [3*n+1, 4*n];
    idx.pk = [4*n+1, 5*n];
    idx.pm = [5*n+1, 6*n];
    idx.m = [6*n+1, 7*n];
    idx.mout = [7*n+1, 8*n];
    idx.i = [8*n+1, 9*n];
    idx.iout = [9*n+1, 10*n];
    idx.p = [10*n+1, 11*n];
    idx.q = [11*n+1, 12*n];
    idx.y = [12*n+1, 13*n];
    
    % Aggregates
    idx.c_util = 13*n + 1;
    idx.l_util = 13*n + 2;
    idx.c_agg = 13*n + 3;
    idx.l_agg = 13*n + 4;
    idx.gdp_agg = 13*n + 5;
    idx.i_agg = 13*n + 6;
    idx.k_agg = 13*n + 7;
    idx.utility_intratemp = 13*n + 8;
    
    % Offset: subtract this from Dynare indices to get policies_ss indices
    idx.ss_offset = 2*n;
    
    % Number of sectors (for convenience)
    idx.n = n;
    
    % Total dimensions
    idx.n_states = 2*n;
    idx.n_policies = 11*n + 8;  % policies_ss length
    idx.n_dynare = 13*n + 8;    % Dynare simulation length
end
