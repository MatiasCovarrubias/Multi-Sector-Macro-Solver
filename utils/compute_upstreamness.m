function Upstreamness = compute_upstreamness(params, ModData)
% COMPUTE_UPSTREAMNESS Calculate steady-state sector upstreamness measures

Upstreamness = initialize_upstreamness_output(params, ModData);
if ~Upstreamness.has_upstreamness
    return;
end

n = Upstreamness.n_sectors;
idx = get_variable_indices(n);
policies_ss = exp(ModData.policies_ss(:));

Gamma_M = get_matrix_param(params, ModData, 'Gamma_M', 'parGamma_M');
Gamma_I = get_matrix_param(params, ModData, 'Gamma_I', 'parGamma_I');
sigma_m = get_scalar_param(params, ModData, 'sigma_m', 'parsigma_m');
sigma_I = get_scalar_param(params, ModData, 'sigma_I', 'parsigma_I');

if isempty(Gamma_M) || isempty(Gamma_I) || ~isfinite(sigma_m) || ~isfinite(sigma_I) || ...
        numel(policies_ss) < idx.n_policies
    Upstreamness.has_upstreamness = false;
    return;
end

Pk = get_policy_block(policies_ss, idx.pk, idx.ss_offset);
Pm = get_policy_block(policies_ss, idx.pm, idx.ss_offset);
M = get_policy_block(policies_ss, idx.m, idx.ss_offset);
Mout = get_policy_block(policies_ss, idx.mout, idx.ss_offset);
Inv = get_policy_block(policies_ss, idx.i, idx.ss_offset);
P = get_policy_block(policies_ss, idx.p, idx.ss_offset);
Q = get_policy_block(policies_ss, idx.q, idx.ss_offset);

Delta_M = Gamma_M .* ((P .^ (-sigma_m)) * (Pm .^ sigma_m)') .* ((1 ./ Q) * M');
Delta_I = Gamma_I .* ((P .^ (-sigma_I)) * (Pk .^ sigma_I)') .* ((1 ./ Q) * Inv');

identity = eye(n);
ones_vec = ones(n, 1);

Upstreamness.U_M = (identity - Delta_M) \ ones_vec;
Upstreamness.U_I = (identity - Delta_I) \ ones_vec;
Upstreamness.U_simple = Mout ./ Q;
Upstreamness.Delta_M = Delta_M;
Upstreamness.Delta_I = Delta_I;
Upstreamness.primary_measure = 'U_M';
Upstreamness.has_upstreamness = true;
end

function Upstreamness = initialize_upstreamness_output(params, ModData)
n = resolve_n_sectors(params, ModData);
Upstreamness = struct( ...
    'has_upstreamness', false, ...
    'primary_measure', 'U_M', ...
    'n_sectors', n, ...
    'U_M', NaN(n, 1), ...
    'U_I', NaN(n, 1), ...
    'U_simple', NaN(n, 1), ...
    'Delta_M', NaN(n, n), ...
    'Delta_I', NaN(n, n));

if n > 0 && isstruct(ModData) && isfield(ModData, 'policies_ss') && ~isempty(ModData.policies_ss)
    Upstreamness.has_upstreamness = true;
end
end

function n = resolve_n_sectors(params, ModData)
n = 0;
if isstruct(params) && isfield(params, 'n_sectors') && isfinite(params.n_sectors)
    n = double(params.n_sectors);
elseif isstruct(ModData) && isfield(ModData, 'parameters') && ...
        isfield(ModData.parameters, 'parn_sectors') && isfinite(ModData.parameters.parn_sectors)
    n = double(ModData.parameters.parn_sectors);
end
end

function value = get_scalar_param(params, ModData, params_field, moddata_field)
value = NaN;
if isstruct(params) && isfield(params, params_field)
    value = params.(params_field);
elseif isstruct(ModData) && isfield(ModData, 'parameters') && isfield(ModData.parameters, moddata_field)
    value = ModData.parameters.(moddata_field);
end
end

function value = get_matrix_param(params, ModData, params_field, moddata_field)
value = [];
if isstruct(params) && isfield(params, params_field)
    value = params.(params_field);
elseif isstruct(ModData) && isfield(ModData, 'parameters') && isfield(ModData.parameters, moddata_field)
    value = ModData.parameters.(moddata_field);
end
end

function values = get_policy_block(policies_ss, dynare_range, ss_offset)
values = policies_ss(dynare_range(1) - ss_offset:dynare_range(2) - ss_offset);
values = values(:);
end
