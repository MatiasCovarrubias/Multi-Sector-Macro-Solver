function value = read_base_workspace_var(var_name)
% READ_BASE_WORKSPACE_VAR Read a variable from MATLAB base workspace

value = evalin('base', var_name);
end
