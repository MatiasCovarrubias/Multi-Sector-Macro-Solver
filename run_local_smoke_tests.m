function Results = run_local_smoke_tests(cfg)
% RUN_LOCAL_SMOKE_TESTS Fast local checks for the MATLAB pipeline

current_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder, 'testing'), '-begin');

if nargin < 1 || isempty(cfg)
    cfg = build_test_defaults();
end

old_dir = pwd;
cleanup_dir = onCleanup(@() cd(old_dir));
cd(cfg.project_root);

env = setup_test_environment(struct('require_dynare', cfg.run_dynare_smoke));

log_file = start_log(cfg.log_dir, 'local_smoke');
cleanup_log = onCleanup(@() stop_log());

fprintf('\nRunning local smoke tests in %s\n', cfg.project_root);

Results = struct();
Results.ok = true;
Results.env = env;
Results.log_file = log_file;
Results.stages = struct();

shared = struct();

run_stage('setup', @stage_setup);

if stage_passed('setup')
    run_stage('load_calibration', @stage_load_calibration);
else
    record_skip('load_calibration', 'setup failed');
end

if stage_passed('load_calibration')
    run_stage('steady_state', @stage_steady_state);
else
    record_skip('steady_state', 'load_calibration failed');
end

if cfg.run_dynare_smoke && stage_passed('steady_state')
    run_stage('dynare_smoke', @stage_dynare_smoke);
elseif cfg.run_dynare_smoke
    record_skip('dynare_smoke', 'steady_state failed');
else
    record_skip('dynare_smoke', 'cfg.run_dynare_smoke = false');
end

print_stage_summary();

    function stage_setup()
        required_files = {'calibration_data.mat', 'TFP_process.mat', 'beadat_37sec.mat'};
        missing_files = {};
        for i = 1:numel(required_files)
            if exist(required_files{i}, 'file') ~= 2
                missing_files{end+1} = required_files{i}; %#ok<AGROW>
            end
        end
        assert(isempty(missing_files), ...
            'run_local_smoke_tests:MissingFile', ...
            'Missing required local data files: %s', strjoin(missing_files, ', '));

        validate_sector_indices(cfg.sector_indices, 37, 'run_local_smoke_tests');
        shared.runtime_config = build_test_runtime_config(struct( ...
            'model_type', cfg.model_type, ...
            'covariance_scale', cfg.covariance_scale, ...
            'sector_indices', cfg.sector_indices, ...
            'force_recalibrate', cfg.force_recalibrate, ...
            'date', cfg.smoke.date, ...
            'exp_label', cfg.smoke.exp_label, ...
            'gridpoints', cfg.smoke.gridpoints, ...
            'ir_horizon', cfg.smoke.ir_horizon, ...
            'ir_plot_length', cfg.smoke.ir_plot_length, ...
            'simul_T', cfg.smoke.simul_T, ...
            'simul_burn_in', cfg.smoke.simul_burn_in, ...
            'simul_burn_out', cfg.smoke.simul_burn_out, ...
            'shock_sizes_pct', cfg.smoke.shock_sizes_pct, ...
            'run_firstorder_simul', true, ...
            'run_secondorder_simul', false, ...
            'run_pf_simul', false, ...
            'run_mit_shocks_simul', false, ...
            'run_firstorder_irs', false, ...
            'run_secondorder_irs', false, ...
            'run_pf_irs', false));
        shared.save_label = strcat(shared.runtime_config.date, shared.runtime_config.exp_label);
    end

    function stage_load_calibration()
        params = build_default_params();
        [shared.calib_data, shared.params] = load_calibration_data( ...
            params, cfg.sector_indices, cfg.model_type, shared.runtime_config.shock_scaling, ...
            shared.runtime_config.smooth, shared.runtime_config.wds, ...
            shared.runtime_config.covariance_scale);
        shared.labels = shared.calib_data.labels;

        assert(isfield(shared.calib_data, 'empirical_targets'), ...
            'run_local_smoke_tests:MissingEmpiricalTargets', ...
            'Calibration data must include empirical targets.');
        assert(isfield(shared.params, 'n_sectors') && shared.params.n_sectors == 37, ...
            'run_local_smoke_tests:InvalidSectorCount', ...
            'Expected 37 sectors after loading calibration data.');
    end

    function stage_steady_state()
        shared.exp_paths = setup_experiment_folder(shared.save_label);
        [shared.ModData, shared.params] = load_or_build_cached_steady_state( ...
            shared.params, shared.runtime_config, shared.exp_paths, struct( ...
            'project_root', cfg.project_root, ...
            'allow_project_cache', cfg.allow_project_cache, ...
            'verbose', false));

        assert(isfield(shared.ModData, 'policies_ss') && isfield(shared.ModData, 'endostates_ss'), ...
            'run_local_smoke_tests:InvalidSteadyState', ...
            'Steady-state output is missing required fields.');
    end

    function stage_dynare_smoke()
        opts = build_dynare_opts(shared.runtime_config, cfg.sector_indices, 'base');
        opts.verbose = false;
        opts.run_firstorder_simul = true;
        opts.run_secondorder_simul = false;
        opts.run_firstorder_irs = false;
        opts.run_secondorder_irs = false;
        opts.run_pf_irs = false;
        opts.run_pf_simul = false;
        opts.run_mit_shocks_simul = false;

        params_for_dynare = shared.params;
        params_for_dynare.IRshock = shared.runtime_config.shock_values(1).value;
        shared.BaseResults = run_dynare_analysis(shared.ModData, params_for_dynare, opts);

        assert(isfield(shared.BaseResults, 'SimulFirstOrder'), ...
            'run_local_smoke_tests:MissingFirstOrderSimulation', ...
            'Dynare smoke test did not produce SimulFirstOrder.');
        assert(isfield(shared.BaseResults.SimulFirstOrder, 'summary_stats'), ...
            'run_local_smoke_tests:MissingModelStats', ...
            'Dynare smoke test did not attach summary_stats to SimulFirstOrder.');
        assert(isfield(shared.BaseResults.SimulFirstOrder.summary_stats, 'ModelStats'), ...
            'run_local_smoke_tests:MissingModelStats', ...
            'Dynare smoke test did not produce first-order ModelStats.');
        assert(size(shared.BaseResults.SimulFirstOrder.burnin_simul, 2) == shared.runtime_config.simul_burn_in, ...
            'run_local_smoke_tests:InvalidBurnInWindow', ...
            'First-order burn-in window length is inconsistent.');
        assert(size(shared.BaseResults.SimulFirstOrder.shocks_simul, 2) == shared.runtime_config.simul_T, ...
            'run_local_smoke_tests:InvalidShockWindow', ...
            'First-order active shock window length is inconsistent.');
        assert(size(shared.BaseResults.SimulFirstOrder.burnout_simul, 2) == shared.runtime_config.simul_burn_out, ...
            'run_local_smoke_tests:InvalidBurnOutWindow', ...
            'First-order burn-out window length is inconsistent.');
        assert(strcmp(shared.BaseResults.SimulFirstOrder.summary_stats.ModelStats.sample_window, 'shocks_simul'), ...
            'run_local_smoke_tests:InvalidMomentSample', ...
            'First-order moments must use shocks_simul.');
    end

    function run_stage(stage_name, fn)
        started_at = tic;
        try
            fn();
            record_stage(stage_name, true, toc(started_at), '');
        catch ME
            Results.ok = false;
            record_stage(stage_name, false, toc(started_at), ME.message, ME);
            if ~cfg.continue_on_failure
                rethrow(ME);
            end
        end
    end

    function record_skip(stage_name, reason)
        stage_key = make_stage_key(stage_name);
        Results.stages.(stage_key) = struct( ...
            'passed', false, ...
            'skipped', true, ...
            'elapsed_seconds', 0, ...
            'message', reason);
    end

    function record_stage(stage_name, passed, elapsed_seconds, message, ME)
        if nargin < 5
            ME = [];
        end

        stage_key = make_stage_key(stage_name);
        stage_info = struct();
        stage_info.passed = passed;
        stage_info.skipped = false;
        stage_info.elapsed_seconds = elapsed_seconds;
        stage_info.message = message;
        if ~isempty(ME)
            stage_info.identifier = ME.identifier;
            stage_info.stack = format_stack(ME);
        end
        Results.stages.(stage_key) = stage_info;
    end

    function passed = stage_passed(stage_name)
        stage_key = make_stage_key(stage_name);
        passed = isfield(Results.stages, stage_key) && ...
            isfield(Results.stages.(stage_key), 'passed') && ...
            Results.stages.(stage_key).passed;
    end

    function print_stage_summary()
        stage_names = fieldnames(Results.stages);
        fprintf('\nSmoke test stage summary:\n');
        for i = 1:numel(stage_names)
            stage_info = Results.stages.(stage_names{i});
            if stage_info.skipped
                fprintf('  %s: SKIPPED (%s)\n', stage_names{i}, stage_info.message);
            elseif stage_info.passed
                fprintf('  %s: PASS (%.2fs)\n', stage_names{i}, stage_info.elapsed_seconds);
            else
                fprintf('  %s: FAIL (%s)\n', stage_names{i}, stage_info.message);
            end
        end
        fprintf('Log file: %s\n', Results.log_file);
    end
end

function log_file = start_log(log_dir, prefix)
if exist(log_dir, 'dir') ~= 7
    mkdir(log_dir);
end

timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
log_file = fullfile(log_dir, sprintf('%s_%s.log', prefix, timestamp));
diary(log_file);
end

function stop_log()
try
    diary off;
catch
end
end

function stack_lines = format_stack(ME)
stack_lines = cell(numel(ME.stack), 1);
for i = 1:numel(ME.stack)
    frame = ME.stack(i);
    stack_lines{i} = sprintf('%s:%d', frame.name, frame.line);
end
end

function stage_key = make_stage_key(stage_name)
stage_key = regexprep(stage_name, '[^A-Za-z0-9_]', '_');
if ~isempty(stage_key) && isstrprop(stage_key(1), 'digit')
    stage_key = ['x' stage_key];
end
end
