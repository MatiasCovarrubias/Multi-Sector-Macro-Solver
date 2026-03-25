function env = setup_test_environment(opts)
% SETUP_TEST_ENVIRONMENT Prepare paths for MATLAB-side test workflows

if nargin < 1
    opts = struct();
end

testing_root = fileparts(mfilename('fullpath'));
project_root = fileparts(testing_root);

addpath(testing_root, '-begin');
addpath(fullfile(project_root, 'utils'), '-begin');
setup_runtime_paths(project_root);

opts = set_default(opts, 'require_dynare', false);

env = struct();
env.project_root = project_root;
env.dynare_path = '';
env.has_dynare = false;

dynare_location = which('dynare');
env.has_dynare = ~isempty(dynare_location);

if env.has_dynare
    env.dynare_path = dynare_location;
elseif opts.require_dynare
    error('setup_test_environment:DynareMissing', ...
        'Dynare must already be on the MATLAB path before running this test.');
end
end
