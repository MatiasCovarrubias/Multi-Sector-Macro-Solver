function setup_dynare_runtime(session, ModData, params, opts)
% SETUP_DYNARE_RUNTIME Prepare the Dynare-specific runtime boundary
%
% All runtime side effects needed by Dynare live here: generated config
% files, temp `.mat` payloads, and base workspace variables.

assign_model_parameters_to_base(ModData.parameters);
write_dynare_runtime_files(session, ModData, opts);

base_vars = struct();
base_vars.params = params;
base_vars.policies_ss = session.policies_ss;
base_vars.k_ss = session.endostates_ss;
base_vars.C_ss = session.steady_state_aggregates.C_ss;
base_vars.L_ss = session.steady_state_aggregates.L_ss;
base_vars.GDP_ss = session.steady_state_aggregates.GDP_ss;
base_vars.I_ss = session.steady_state_aggregates.I_ss;
base_vars.K_ss = session.steady_state_aggregates.K_ss;
base_vars.utility_intratemp_ss = session.steady_state_aggregates.utility_intratemp_ss;
base_vars.parn_sectors = session.n_sectors;
assign_base_workspace_vars(base_vars);
end

function assign_model_parameters_to_base(parameters)
params_vars = struct2cell(parameters);
params_names = fieldnames(parameters);
for i = 1:numel(params_vars)
    assignin('base', params_names{i}, params_vars{i});
end
end

function write_dynare_runtime_files(session, ModData, opts)
N = opts.ir_horizon;
ax = 0:N-1;
policies_ss = ModData.policies_ss;
k_ss = ModData.endostates_ss;
idx = session.idx;
p_ss_log = policies_ss((idx.p(1) - idx.ss_offset):(idx.p(2) - idx.ss_offset));
pk_ss_log = policies_ss((idx.pk(1) - idx.ss_offset):(idx.pk(2) - idx.ss_offset));
modstruct_path = fullfile(session.dynare_folder, 'ModStruct_temp.mat');
runtime_payload = collect_runtime_payload( ...
    ModData.parameters, policies_ss, k_ss, p_ss_log, pk_ss_log, N, ax, ...
    session.steady_state_aggregates);
save(modstruct_path, '-struct', 'runtime_payload');

model_config_path = fullfile(session.dynare_folder, 'model_config.mod');
fid = fopen(model_config_path, 'w');
fprintf(fid, '@#define n_sectors = %d\n', session.n_sectors);
switch opts.model_type
    case 'VA'
        fprintf(fid, '@#define MODEL_TYPE = 1\n');
    case {'GO', 'GO_noVA'}
        fprintf(fid, '@#define MODEL_TYPE = 2\n');
    otherwise
        fclose(fid);
        error('setup_dynare_runtime:UnsupportedModelType', ...
            'Unsupported model_type: %s', char(opts.model_type));
end
fclose(fid);

base_model_include_path = fullfile(session.dynare_folder, 'base_model_include.mod');
fid = fopen(base_model_include_path, 'w');
fprintf(fid, '@#include "%s"\n', get_base_model_filename(opts.model_type));
fclose(fid);
end

function runtime_payload = collect_runtime_payload(parameters, policies_ss, k_ss, p_ss_log, pk_ss_log, N, ax, ...
        steady_state_aggregates)
runtime_payload = struct();
runtime_payload.policies_ss = policies_ss;
runtime_payload.k_ss = k_ss;
runtime_payload.p_ss_log = p_ss_log;
runtime_payload.pk_ss_log = pk_ss_log;
runtime_payload.N = N;
runtime_payload.ax = ax;
runtime_payload.C_ss = steady_state_aggregates.C_ss;
runtime_payload.L_ss = steady_state_aggregates.L_ss;
runtime_payload.GDP_ss = steady_state_aggregates.GDP_ss;
runtime_payload.I_ss = steady_state_aggregates.I_ss;
runtime_payload.K_ss = steady_state_aggregates.K_ss;
runtime_payload.utility_intratemp_ss = steady_state_aggregates.utility_intratemp_ss;

parameter_names = fieldnames(parameters);
for i = 1:numel(parameter_names)
    param_name = parameter_names{i};
    if numel(param_name) >= 3 && strcmp(param_name(1:3), 'par')
        runtime_payload.(param_name) = parameters.(param_name);
    end
end

end

function base_model_filename = get_base_model_filename(model_type)
switch char(model_type)
    case 'VA'
        base_model_filename = 'ProdNetRbc_base.mod';
    case {'GO', 'GO_noVA'}
        base_model_filename = 'ProdNetRbc_base_GO.mod';
    otherwise
        error('setup_dynare_runtime:UnsupportedModelType', ...
            'Unsupported model_type: %s', char(model_type));
end
end
