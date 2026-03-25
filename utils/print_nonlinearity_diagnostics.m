function Diagnostics = print_nonlinearity_diagnostics(ModelData_simulation, AllShockResults, params, config, ModData)
% PRINT_NONLINEARITY_DIAGNOSTICS Legacy wrapper kept silent on purpose.
%
% Nonlinearity diagnostics are still computed for packaging, but they are no
% longer printed during normal runs.

Diagnostics = build_nonlinearity_diagnostics(ModelData_simulation, AllShockResults, params, config, ModData);

end

