function assign_base_workspace_vars(vars_struct)
% ASSIGN_BASE_WORKSPACE_VARS Assign a struct of variables into MATLAB base workspace

var_names = fieldnames(vars_struct);
for i = 1:numel(var_names)
    assignin('base', var_names{i}, vars_struct.(var_names{i}));
end
end
