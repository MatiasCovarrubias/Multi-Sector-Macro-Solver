function ModelData = build_ModelData(config, save_label, sector_indices, n_shocks, ...
                                     calib_data, labels, params, ModData, exp_paths)
    ModelData = struct();

    %% Metadata
    ModelData.metadata.date = config.date;
    ModelData.metadata.exp_label = config.exp_label;
    ModelData.metadata.save_label = save_label;
    ModelData.metadata.model_type = config.model_type;
    ModelData.metadata.smooth = calib_data.smooth;
    ModelData.metadata.wds = calib_data.wds;
    ModelData.metadata.covariance_scale = calib_data.covariance_scale;
    ModelData.metadata.tfp_suffix = calib_data.tfp_suffix;
    ModelData.metadata.sector_indices = sector_indices;
    ModelData.metadata.sector_labels = labels.sector_labels;
    ModelData.metadata.config = config;
    ModelData.metadata.exp_paths = exp_paths;
    ModelData.metadata.n_shocks = n_shocks;
    ModelData.metadata.run_flags = struct('has_1storder', false, 'has_2ndorder', false, ...
        'has_pf', false, 'has_mit', false);
    ModelData.metadata.has_irfs = false;
    ModelData.metadata.has_diagnostics = false;

    %% Calibration
    ModelData.calibration = calib_data;
    ModelData.params = params;
    ModelData.EmpiricalTargets = calib_data.empirical_targets;

    %% Steady state
    ModelData.SteadyState.parameters = ModData.parameters;
    ModelData.SteadyState.policies_ss = ModData.policies_ss;
    ModelData.SteadyState.endostates_ss = ModData.endostates_ss;
    ModelData.SteadyState.C_ss = ModData.C_ss;
    ModelData.SteadyState.L_ss = ModData.L_ss;
    ModelData.SteadyState.GDP_ss = ModData.GDP_ss;
    ModelData.SteadyState.I_ss = ModData.I_ss;
    ModelData.SteadyState.K_ss = ModData.K_ss;
    ModelData.SteadyState.utility_intratemp_ss = ModData.utility_intratemp_ss;

end
