function Diagnostics = print_nonlinearity_diagnostics(ModelData_simulation, AllShockResults, params, config, ModData)
% PRINT_NONLINEARITY_DIAGNOSTICS Legacy wrapper kept silent on purpose.
%
% Only the compact CIR-based IRF summary is still computed for
% packaging, and it is not printed during normal runs.

Diagnostics = build_nonlinearity_diagnostics(ModelData_simulation, AllShockResults, params, config, ModData);

end

