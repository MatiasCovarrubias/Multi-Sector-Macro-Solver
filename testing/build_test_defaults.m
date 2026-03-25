function cfg = build_test_defaults()
% BUILD_TEST_DEFAULTS Default settings for local MATLAB test workflows

testing_root = fileparts(mfilename('fullpath'));
project_root = fileparts(testing_root);
experiments_root = fullfile(project_root, 'experiments');

cfg = struct();
cfg.project_root = project_root;
cfg.testing_root = testing_root;
cfg.experiments_root = experiments_root;

cfg.sector_indices = 1;
cfg.model_type = 'VA';
cfg.covariance_scale = 1.0;
cfg.force_recalibrate = false;
cfg.allow_project_cache = true;
cfg.continue_on_failure = false;
cfg.run_dynare_smoke = false;

cfg.log_dir = fullfile(experiments_root, 'test_logs');
cfg.fixture_dir = fullfile(experiments_root, 'test_fixtures');
cfg.fixture_bundle = fullfile(cfg.fixture_dir, 'packaging_fixture.mat');
cfg.stage_inputs_fixture = fullfile(cfg.fixture_dir, 'stage_inputs_fixture.mat');
cfg.steady_state_fixture = fullfile(cfg.fixture_dir, 'steady_state_fixture.mat');
cfg.replay_root = fullfile(experiments_root, 'fixture_replays');

cfg.smoke = struct();
cfg.smoke.date = "_Mar_2026";
cfg.smoke.exp_label = "_local_smoke";
cfg.smoke.gridpoints = 4;
cfg.smoke.ir_horizon = 12;
cfg.smoke.ir_plot_length = 8;
cfg.smoke.simul_T = 8;
cfg.smoke.simul_burn_in = 2;
cfg.smoke.simul_burn_out = 2;
cfg.smoke.shock_sizes_pct = 5;

cfg.fixture = struct();
cfg.fixture.date = "_Mar_2026";
cfg.fixture.exp_label = "_fixture";
cfg.fixture.gridpoints = 4;
cfg.fixture.ir_horizon = 20;
cfg.fixture.ir_plot_length = 12;
cfg.fixture.simul_T = 16;
cfg.fixture.simul_burn_in = 0;
cfg.fixture.simul_burn_out = 0;
cfg.fixture.shock_sizes_pct = 5;
cfg.fixture.run_firstorder_simul = true;
cfg.fixture.run_firstorder_irs = true;
cfg.fixture.run_pf_irs = true;
end
