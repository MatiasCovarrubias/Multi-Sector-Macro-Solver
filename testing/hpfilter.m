function [trend, cycle] = hpfilter(y, lambda)
T = length(y);
col = size(y, 2) > size(y, 1);
if col
    y = y(:);
end
e = ones(T, 1);
D = spdiags([e, -2 * e, e], 0:2, T - 2, T);
trend = (speye(T) + lambda * (D' * D)) \ y;
cycle = y - trend;
if col
    trend = trend';
    cycle = cycle';
end
end
