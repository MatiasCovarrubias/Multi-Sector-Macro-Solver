function params = params_config_defaults()
% PARAMS_CONFIG_DEFAULTS Canonical parameter defaults for main.m

params = struct();
params.GHH = true;
params.beta = 0.96;
params.eps_l = 0.05;
params.eps_c = 0.5;
params.theta = 1;
params.phi = 3.5;
params.sigma_c = 0.5;
params.sigma_m = 0.001;
params.sigma_q = 0.5;
params.sigma_y = 0.6;
params.sigma_I = 0.5;
params.sigma_l = 0.5;
end
