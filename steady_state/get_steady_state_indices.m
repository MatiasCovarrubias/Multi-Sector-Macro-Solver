function idx = get_steady_state_indices(n)
% GET_STEADY_STATE_INDICES Index map for the steady-state solver vector.

idx = struct();

idx.c = [1, n];
idx.l = [n + 1, 2 * n];
idx.pk = [2 * n + 1, 3 * n];
idx.pm = [3 * n + 1, 4 * n];
idx.m = [4 * n + 1, 5 * n];
idx.mout = [5 * n + 1, 6 * n];
idx.i = [6 * n + 1, 7 * n];
idx.iout = [7 * n + 1, 8 * n];
idx.p = [8 * n + 1, 9 * n];
idx.q = [9 * n + 1, 10 * n];
idx.y = [10 * n + 1, 11 * n];

idx.c_util = 11 * n + 1;
idx.l_util = 11 * n + 2;

idx.n = n;
idx.n_base = 11 * n + 2;
end
