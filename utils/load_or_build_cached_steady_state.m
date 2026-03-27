function [ModData, params, cache_info] = load_or_build_cached_steady_state(params, runtime_config, exp_paths, cache_opts)
% LOAD_OR_BUILD_CACHED_STEADY_STATE Shared steady-state cache policy

if nargin < 4 || isempty(cache_opts)
    cache_opts = struct();
end

cache_opts = set_default(cache_opts, 'project_root', '');
cache_opts = set_default(cache_opts, 'allow_project_cache', false);
cache_opts = set_default(cache_opts, 'verbose', false);

current_params = params;

local_ss_file = fullfile(exp_paths.experiment, 'steady_state.mat');
project_ss_file = '';
if cache_opts.allow_project_cache && ~isempty(cache_opts.project_root)
    project_ss_file = fullfile(cache_opts.project_root, 'calibrated_steady_state.mat');
end

cache_info = struct('source', 'recomputed', 'path', local_ss_file);

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

[ModData, params] = calibrate_steady_state(params, calib_opts);
save(local_ss_file, 'ModData', 'params');
end

function [ModData, params] = reconcile_cached_steady_state(ModData, cached_params, current_params)
params = cached_params;

runtime_fields = { ...
    'delta', 'n_sectors', 'model_type', 'smooth', 'wds', 'covariance_scale', ...
    'tfp_suffix', 'rho', 'Sigma_A_full', 'Sigma_A', 'tfp_process_fields', ...
    'shock_scaling', 'conssh_data', 'capsh_data', 'vash_data', ...
    'ionet_data', 'invnet_data'};

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
tf = resolve_preference_flag(cached_params) == resolve_preference_flag(current_params);
end

function GHH = resolve_preference_flag(params)
GHH = true;
if isfield(params, 'GHH')
    GHH = logical(params.GHH);
end
end
