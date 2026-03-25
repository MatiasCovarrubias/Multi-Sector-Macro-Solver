function endo_simul = run_generated_dynare_job(session, mod_name, workspace_vars)
% RUN_GENERATED_DYNARE_JOB Execute a generated Dynare job and return oo_.endo_simul

if nargin >= 3 && ~isempty(workspace_vars)
    assign_base_workspace_vars(workspace_vars);
end

clear_generated_dynare_outputs();
cleanup = onCleanup(@() clear_generated_dynare_outputs()); %#ok<NASGU>
run_dynare_mod(session.dynare_folder, mod_name);
assert(evalin('base', 'exist(''oo_'', ''var'')') == 1, ...
    'run_generated_dynare_job:MissingDynareOutput', ...
    'Dynare job %s did not populate oo_ in the base workspace.', mod_name);

oo_results = read_base_workspace_var('oo_');
assert(isstruct(oo_results) && isfield(oo_results, 'endo_simul') && ~isempty(oo_results.endo_simul), ...
    'run_generated_dynare_job:MissingSimulationOutput', ...
    'Dynare job %s did not populate oo_.endo_simul.', mod_name);

endo_simul = oo_results.endo_simul;
end

function clear_generated_dynare_outputs()
evalin('base', ['if exist(''oo_'', ''var''), clear(''oo_''); end; ' ...
                'if exist(''Simulated_time_series'', ''var''), clear(''Simulated_time_series''); end']);
end
