function Summary = run_debug_loop(cfg)
% RUN_DEBUG_LOOP Run the local MATLAB test workflow with persistent logs

current_folder = fileparts(mfilename('fullpath'));
addpath(fullfile(current_folder, 'testing'), '-begin');

if nargin < 1 || isempty(cfg)
    cfg = build_test_defaults();
end

old_dir = pwd;
cleanup_dir = onCleanup(@() cd(old_dir));
cd(cfg.project_root);

if exist(cfg.log_dir, 'dir') ~= 7
    mkdir(cfg.log_dir);
end

timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
suite_log = fullfile(cfg.log_dir, sprintf('debug_loop_%s.log', timestamp));
diary(suite_log);
cleanup_log = onCleanup(@() stop_log());

fprintf('\nStarting debug loop\n');

Summary = struct();
Summary.suite_log = suite_log;
Summary.smoke = [];
Summary.fixtures = [];
Summary.ok = true;

try
    Summary.smoke = run_local_smoke_tests(cfg);
catch ME
    Summary.ok = false;
    Summary.smoke = struct('ok', false, 'message', ME.message);
    print_error('smoke', ME);
end

if exist(cfg.fixture_bundle, 'file') == 2
    try
        Summary.fixtures = run_fixture_tests(cfg);
    catch ME
        Summary.ok = false;
        Summary.fixtures = struct('ok', false, 'message', ME.message);
        print_error('fixtures', ME);
    end
else
    Summary.fixtures = struct( ...
        'ok', false, ...
        'skipped', true, ...
        'message', sprintf('Fixture bundle not found: %s', cfg.fixture_bundle));
    fprintf('Fixture replay skipped: %s\n', Summary.fixtures.message);
end

Summary.ok = Summary.ok && Summary.smoke.ok && ...
    (~isfield(Summary.fixtures, 'skipped') || Summary.fixtures.skipped || Summary.fixtures.ok);

fprintf('Suite log: %s\n', suite_log);
fprintf('Overall status: %s\n', bool_str(Summary.ok));
end

function print_error(stage_name, ME)
fprintf('\n%s failed: %s\n', stage_name, ME.message);
for i = 1:numel(ME.stack)
    fprintf('  at %s:%d\n', ME.stack(i).name, ME.stack(i).line);
end
end

function stop_log()
try
    diary off;
catch
end
end

function value = bool_str(flag)
if flag
    value = 'PASS';
else
    value = 'FAIL';
end
end
