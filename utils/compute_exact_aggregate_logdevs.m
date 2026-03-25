function [C_agg_logdev, I_agg_logdev, GDP_agg_logdev, L_agg_logdev, M_agg_logdev, K_agg_logdev, ss] = ...
        compute_exact_aggregate_logdevs(simul, idx, policies_ss, endostates_ss)
% COMPUTE_EXACT_AGGREGATE_LOGDEVS Build exact aggregate log deviations.
%
% Input simulation paths are assumed to be in model-space log levels, as
% returned by Dynare. This helper reconstructs exact aggregate levels using
% deterministic steady-state levels and then renormalizes them to
% log-deviations from deterministic steady state.

epsilon = 1e-12;

ss_of = @(range) policies_ss((range(1):range(2)) - idx.ss_offset);

p_ss = exp(ss_of(idx.p));         p_ss = p_ss(:);
c_ss = exp(ss_of(idx.c));         c_ss = c_ss(:);
iout_ss = exp(ss_of(idx.iout));   iout_ss = iout_ss(:);
q_ss = exp(ss_of(idx.q));         q_ss = q_ss(:);
mout_ss = exp(ss_of(idx.mout));   mout_ss = mout_ss(:);
l_ss = exp(ss_of(idx.l));         l_ss = l_ss(:);
pk_ss = exp(ss_of(idx.pk));       pk_ss = pk_ss(:);
k_ss = exp(endostates_ss(1:numel(pk_ss))); k_ss = k_ss(:);

c_logdev = simul(idx.c(1):idx.c(2), :) - ss_of(idx.c);
iout_logdev = simul(idx.iout(1):idx.iout(2), :) - ss_of(idx.iout);
q_logdev = simul(idx.q(1):idx.q(2), :) - ss_of(idx.q);
mout_logdev = simul(idx.mout(1):idx.mout(2), :) - ss_of(idx.mout);
l_logdev = simul(idx.l(1):idx.l(2), :) - ss_of(idx.l);
k_logdev = simul(idx.k(1):idx.k(2), :) - endostates_ss(1:numel(pk_ss));

ss = struct();
ss.C = sum(p_ss .* c_ss);
ss.I = sum(p_ss .* iout_ss);
ss.GDP = sum(p_ss .* (q_ss - mout_ss));
ss.L = sum(l_ss);
ss.M = sum(p_ss .* mout_ss);
ss.K = sum(pk_ss .* k_ss);
ss.share_C = ss.C / ss.GDP;
ss.share_I = ss.I / ss.GDP;

C_levels = sum(p_ss .* (c_ss .* exp(c_logdev)), 1);
I_levels = sum(p_ss .* (iout_ss .* exp(iout_logdev)), 1);
GDP_levels = sum(p_ss .* (q_ss .* exp(q_logdev) - mout_ss .* exp(mout_logdev)), 1);
L_levels = sum(l_ss .* exp(l_logdev), 1);
M_levels = sum(p_ss .* (mout_ss .* exp(mout_logdev)), 1);
K_levels = sum(pk_ss .* (k_ss .* exp(k_logdev)), 1);

C_agg_logdev = logdev_from_levels(C_levels, ss.C, 'C', epsilon);
I_agg_logdev = logdev_from_levels(I_levels, ss.I, 'I', epsilon);
GDP_agg_logdev = logdev_from_levels(GDP_levels, ss.GDP, 'GDP', epsilon);
L_agg_logdev = logdev_from_levels(L_levels, ss.L, 'L', epsilon);
M_agg_logdev = logdev_from_levels(M_levels, ss.M, 'M', epsilon);
K_agg_logdev = logdev_from_levels(K_levels, ss.K, 'K', epsilon);
end

function x_logdev = logdev_from_levels(levels, steady_state_level, label, epsilon)
assert(isfinite(steady_state_level) && steady_state_level > 0, ...
    'compute_exact_aggregate_logdevs:InvalidSteadyStateLevel', ...
    'Steady-state level for %s must be finite and strictly positive.', label);
assert(all(isfinite(levels)), ...
    'compute_exact_aggregate_logdevs:NonFiniteAggregateLevel', ...
    'Aggregate %s contains non-finite simulated levels.', label);

tolerance = max(epsilon, 1e-10 * max(1, abs(steady_state_level)));
min_level = min(levels);
if min_level <= 0
    assert(min_level >= -tolerance, ...
        'compute_exact_aggregate_logdevs:NonpositiveAggregateLevel', ...
        'Aggregate %s reached a nonpositive level below tolerance (min=%.4e).', ...
        label, min_level);
    levels = max(levels, tolerance);
end

x_logdev = log(levels) - log(steady_state_level);
end
