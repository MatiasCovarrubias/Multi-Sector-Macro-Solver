function [ModData, params, cache_info] = load_or_build_cached_steady_state(params, runtime_config, exp_paths, cache_opts)
% LOAD_OR_BUILD_CACHED_STEADY_STATE Shared steady-state cache policy

if nargin < 4 || isempty(cache_opts)
    cache_opts = struct();
end

cache_opts = set_default(cache_opts, 'project_root', '');
cache_opts = set_default(cache_opts, 'allow_project_cache', false);
cache_opts = set_default(cache_opts, 'verbose', false);
cache_opts = set_default(cache_opts, 'save_local_cache', true);

current_params = params;

local_ss_file = fullfile(exp_paths.experiment, 'steady_state.mat');
project_ss_file = '';
if cache_opts.allow_project_cache && ~isempty(cache_opts.project_root)
    project_ss_file = fullfile(cache_opts.project_root, 'calibrated_steady_state.mat');
end

cache_info = struct( ...
    'source', 'recomputed', ...
    'path', local_ss_file, ...
    'saved', false, ...
    'elapsed_seconds', NaN);

if ~runtime_config.force_recalibrate
    if exist(local_ss_file, 'file') == 2
        loaded = load(local_ss_file, 'ModData', 'params');
        if is_steady_state_cache_compatible(loaded.params, current_params)
            [ModData, params] = reconcile_cached_steady_state( ...
                loaded.ModData, loaded.params, current_params);
            cache_info.source = 'experiment_cache';
            cache_info.path = local_ss_file;
            return;
        end
    end

    if ~isempty(project_ss_file) && exist(project_ss_file, 'file') == 2
        loaded = load(project_ss_file, 'ModData', 'params');
        if is_steady_state_cache_compatible(loaded.params, current_params)
            [ModData, params] = reconcile_cached_steady_state( ...
                loaded.ModData, loaded.params, current_params);
            cache_info.source = 'project_cache';
            cache_info.path = project_ss_file;
            return;
        end
    end
end

calib_opts = struct();
calib_opts.gridpoints = runtime_config.gridpoints;
calib_opts.verbose = cache_opts.verbose;
calib_opts.sol_guess_file = runtime_config.sol_guess_file;
calib_opts.fsolve_options = runtime_config.fsolve_options;

tic;
[ModData, params] = calibrate_steady_state(params, calib_opts);
cache_info.elapsed_seconds = toc;

if cache_opts.save_local_cache
    save(local_ss_file, 'ModData', 'params');
    cache_info.saved = true;
end
end

function [ModData, params] = reconcile_cached_steady_state(ModData, cached_params, current_params)
params = cached_params;

runtime_fields = { ...
    'delta', 'smooth', 'wds', 'covariance_scale', 'tfp_suffix', 'rho', ...
    'Sigma_A_full', 'Sigma_A', 'tfp_process_fields', 'shock_scaling'};

for i = 1:numel(runtime_fields)
    field_name = runtime_fields{i};
    if isfield(current_params, field_name)
        params.(field_name) = current_params.(field_name);
    end
end

params.GHH = resolve_preference_flag(current_params);

if isfield(ModData, 'parameters')
    if isfield(current_params, 'rho')
        ModData.parameters.parrho = current_params.rho;
    end
    if isfield(current_params, 'Sigma_A')
        ModData.parameters.parSigma_A = current_params.Sigma_A;
    end
    if isfield(current_params, 'delta')
        ModData.parameters.pardelta = current_params.delta;
    end
    ModData.parameters.parGHH = double(params.GHH);
end
end

function tf = is_steady_state_cache_compatible(cached_params, current_params)
cached_signature = build_steady_state_signature(cached_params);
current_signature = build_steady_state_signature(current_params);
tf = are_steady_state_signatures_equal(cached_signature, current_signature);
end

function GHH = resolve_preference_flag(params)
GHH = true;
if isfield(params, 'GHH')
    GHH = logical(params.GHH);
end
end

function signature = build_steady_state_signature(params)
signature = struct();

signature.n_sectors = get_signature_field(params, 'n_sectors');
signature.model_type = get_signature_field(params, 'model_type');

signature.beta = get_signature_field(params, 'beta');
signature.theta = get_signature_field(params, 'theta');
signature.eps_c = get_signature_field(params, 'eps_c');
signature.eps_l = get_signature_field(params, 'eps_l');
signature.GHH = resolve_preference_flag(params);

signature.sigma_c = get_signature_field(params, 'sigma_c');
signature.sigma_m = get_signature_field(params, 'sigma_m');
signature.sigma_q = get_signature_field(params, 'sigma_q');
signature.sigma_y = get_signature_field(params, 'sigma_y');
signature.sigma_I = get_signature_field(params, 'sigma_I');
signature.sigma_l = get_signature_field(params, 'sigma_l');

signature.conssh_data = get_signature_field(params, 'conssh_data');
signature.capsh_data = get_signature_field(params, 'capsh_data');
signature.vash_data = get_signature_field(params, 'vash_data');
signature.ionet_data = get_signature_field(params, 'ionet_data');
signature.invnet_data = get_signature_field(params, 'invnet_data');
end

function value = get_signature_field(params, field_name)
value = [];
if isfield(params, field_name)
    value = params.(field_name);
end
end

function tf = are_steady_state_signatures_equal(signature_a, signature_b)
tf = isequaln(signature_a, signature_b);
end
