function params = params_config()
% PARAMS_CONFIG User-facing parameter configuration for main.m

% Set to false to use the values below.
use_defaults = false;

params = struct();
params.GHH = false;
params.beta = 0.96;
params.eps_l = 0.5;
params.eps_c = 0.33;
params.theta = 1;
params.phi = 4.0;
params.sigma_c = 0.5;
params.sigma_m = 0.001;
params.sigma_q = 0.5;
params.sigma_y = 0.8;
params.sigma_I = 0.5;
params.sigma_l = 0.065;

if use_defaults
    params = params_config_defaults();
end
end
