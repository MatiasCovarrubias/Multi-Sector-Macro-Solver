function save_experiment_results(config, params_config_snapshot, exp_paths, ModelData, ModelData_simulation, ModelData_IRs, flags, has_irfs, source_files)
    fprintf('\n--- Save ---\n');

    if ~config.save_results
        fprintf('Saving disabled\n');
        return;
    end

    if nargin < 9 || isempty(source_files)
        source_files = struct('runtime', '', 'params', '');
    end

    validate_pipeline_outputs(ModelData, ModelData_simulation, ModelData_IRs, ...
        flags, has_irfs, 'save_experiment_results');

    runtime_config_snapshot = config;
    filename_runtime_config = fullfile(exp_paths.experiment, 'runtime_config.mat');
    save(filename_runtime_config, 'runtime_config_snapshot');
    fprintf('runtime_config: %s\n', filename_runtime_config);

    filename_params_config = fullfile(exp_paths.experiment, 'params_config.mat');
    save(filename_params_config, 'params_config_snapshot');
    fprintf('params_config: %s\n', filename_params_config);

    if isfield(source_files, 'runtime') && exist(source_files.runtime, 'file') == 2
        runtime_config_copy = fullfile(exp_paths.experiment, 'runtime_config.m');
        copyfile(source_files.runtime, runtime_config_copy);
        fprintf('runtime_config source: %s\n', runtime_config_copy);
    end

    if isfield(source_files, 'params') && exist(source_files.params, 'file') == 2
        params_config_copy = fullfile(exp_paths.experiment, 'params_config.m');
        copyfile(source_files.params, params_config_copy);
        fprintf('params_config source: %s\n', params_config_copy);
    end

    filename_model = fullfile(exp_paths.experiment, 'ModelData.mat');
    save(filename_model, 'ModelData');
    fprintf('ModelData: %s\n', filename_model);

    if flags.has_1storder || flags.has_2ndorder || flags.has_pf || flags.has_mit
        filename_simul = fullfile(exp_paths.experiment, 'ModelData_simulation.mat');
        save(filename_simul, 'ModelData_simulation');
        fprintf('ModelData_simulation: %s\n', filename_simul);
    end

    if has_irfs
        filename_irs = fullfile(exp_paths.experiment, 'ModelData_IRs.mat');
        save(filename_irs, 'ModelData_IRs');
        fprintf('ModelData_IRs: %s\n', filename_irs);
    end
end
