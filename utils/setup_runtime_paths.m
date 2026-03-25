function setup_runtime_paths(current_folder)
% SETUP_RUNTIME_PATHS Add only active runtime folders to the MATLAB path
%
% Deliberately excludes offline preprocessing (`Data`) and legacy (`OLD`)
% folders so runtime scripts only resolve against supported pipeline code.

active_dirs = {'plotting', 'dynare', 'steady_state', 'calibration'};

for i = 1:numel(active_dirs)
    addpath(fullfile(current_folder, active_dirs{i}), '-begin');
end
end
