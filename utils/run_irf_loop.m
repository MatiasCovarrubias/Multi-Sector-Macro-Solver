function [AllShockResults, BaseResults] = run_irf_loop(config, sector_indices, ...
        ModData, params, BaseResults, AllShockResults, exp_paths, labels)

    n_shocks = numel(config.shock_values);
    fprintf('\nIRF Analysis: %d shocks, sectors [%s]\n', n_shocks, num2str(sector_indices));
    shared_solution_cache = extract_solution_cache(BaseResults);
    announce_cache_status(shared_solution_cache);

    for shock_idx = 1:n_shocks
        shock_config = config.shock_values(shock_idx);
        params.IRshock = shock_config.value;

        fprintf('\n[IRF] Shock %d/%d | %s | value=%.4f\n', ...
            shock_idx, n_shocks, shock_config.label, shock_config.value);
        fprintf('      Description: %s\n', shock_config.description);

        dynare_opts = build_dynare_opts(config, sector_indices, 'irf');
        dynare_opts.solution_cache = shared_solution_cache;
        dynare_opts.runtime_context = struct( ...
            'stage', 'IRF', ...
            'shock_idx', shock_idx, ...
            'n_shocks', n_shocks, ...
            'shock_label', shock_config.label, ...
            'shock_description', shock_config.description);

        tic;
        DynareResults = run_dynare_analysis(ModData, params, dynare_opts);
        elapsed_irf = toc;
        shared_solution_cache = merge_solution_cache(shared_solution_cache, DynareResults);

        AllShockResults.DynareResults{shock_idx} = DynareResults;
        BaseResults = update_base_irf_fields(BaseResults, DynareResults, shock_idx);

        ir_opts = build_ir_processing_opts(config, exp_paths, shock_config);

        IRShockArtifact = process_sector_irs(DynareResults, params, ModData, labels, ir_opts);
        IRShockArtifact.shock_config = shock_config;
        if ~isfield(AllShockResults, 'ShockArtifacts') || isempty(AllShockResults.ShockArtifacts)
            AllShockResults.ShockArtifacts = cell(n_shocks, 1);
        end
        print_ir_artifact_debug(shock_idx, n_shocks, IRShockArtifact, AllShockResults);
        AllShockResults.ShockArtifacts{shock_idx} = IRShockArtifact;

        fprintf('  Done (%.1fs)\n', elapsed_irf);
    end

    fprintf('All %d shocks processed\n', n_shocks);
end

function BaseResults = update_base_irf_fields(BaseResults, DynareResults, shock_idx)
if shock_idx ~= 1
    return;
end

if ~isfield(BaseResults, 'TheoStats') && isfield(DynareResults, 'TheoStats')
    BaseResults.TheoStats = DynareResults.TheoStats;
end

if ~isfield(BaseResults, 'steady_state_aggregates') && isfield(DynareResults, 'steady_state_aggregates')
    BaseResults.steady_state_aggregates = DynareResults.steady_state_aggregates;
    BaseResults.C_ss = DynareResults.C_ss;
    BaseResults.L_ss = DynareResults.L_ss;
    BaseResults.GDP_ss = DynareResults.GDP_ss;
    BaseResults.I_ss = DynareResults.I_ss;
    BaseResults.K_ss = DynareResults.K_ss;
    BaseResults.utility_intratemp_ss = DynareResults.utility_intratemp_ss;
end

if ~isfield(BaseResults, 'steady_state') && isfield(DynareResults, 'steady_state')
    BaseResults.steady_state = DynareResults.steady_state;
end
end

function ir_opts = build_ir_processing_opts(config, exp_paths, shock_config)
ir_opts = struct();
ir_opts.plot_graphs = config.plot_irs;
ir_opts.save_graphs = config.save_results;
ir_opts.save_intermediate = config.save_results;
ir_opts.save_interval = 5;
ir_opts.exp_paths = exp_paths;
ir_opts.save_label = shock_config.label;
ir_opts.ir_plot_length = config.ir_plot_length;
ir_opts.ir_horizon = config.ir_horizon;
ir_opts.shock_description = shock_config.description;
ir_opts.shock_config = shock_config;
end

function cache = extract_solution_cache(RuntimeResults)
cache = struct();

if isempty(RuntimeResults) || ~isstruct(RuntimeResults)
    return;
end

cache = merge_solution_cache(cache, RuntimeResults);
end

function cache = merge_solution_cache(cache, RuntimeResults)
if isempty(RuntimeResults) || ~isstruct(RuntimeResults)
    return;
end

cache = maybe_copy_solution(cache, RuntimeResults, ...
    {'oo_1st', 'M_1st', 'options_1st'}, 'first_order');
cache = maybe_copy_solution(cache, RuntimeResults, ...
    {'oo_2nd', 'M_2nd', 'options_2nd'}, 'second_order');
cache = maybe_copy_pf_irf_cache(cache, RuntimeResults);
end

function cache = maybe_copy_solution(cache, RuntimeResults, source_fields, target_field)
has_all_fields = all(isfield(RuntimeResults, source_fields));
if ~has_all_fields
    return;
end

values = cellfun(@(field_name) RuntimeResults.(field_name), source_fields, 'UniformOutput', false);
if any(cellfun(@isempty, values))
    return;
end

cache.(target_field) = struct( ...
    'oo', values{1}, ...
    'M', values{2}, ...
    'options', values{3});
end

function cache = maybe_copy_pf_irf_cache(cache, RuntimeResults)
if ~isfield(RuntimeResults, 'pf_irf_cache') || isempty(RuntimeResults.pf_irf_cache)
    return;
end

candidate = RuntimeResults.pf_irf_cache;
required_fields = {'signature', 'session'};
if ~(isstruct(candidate) && all(isfield(candidate, required_fields)))
    return;
end

if isfield(candidate, 'was_reused')
    candidate = rmfield(candidate, 'was_reused');
end

cache.pf_irf = candidate;
end

function announce_cache_status(cache)
has_first = isfield(cache, 'first_order');
has_second = isfield(cache, 'second_order');
has_pf = isfield(cache, 'pf_irf');

if has_first || has_second || has_pf
    fprintf('[IRF] Reusing cached solutions: 1st=%s, 2nd=%s, PF=%s\n', ...
        logical_str(has_first), logical_str(has_second), logical_str(has_pf));
else
    fprintf('[IRF] No cached IR solutions available; first requested shock will solve them once.\n');
end
end

function s = logical_str(flag)
if flag
    s = 'yes';
else
    s = 'no';
end
end

function print_ir_artifact_debug(shock_idx, n_shocks, artifact, AllShockResults)
fprintf('[run_irf_loop] shock %d/%d | artifact_fields=%s\n', ...
    shock_idx, n_shocks, join_fieldnames(fieldnames(artifact)));
if isfield(artifact, 'metadata')
    fprintf('[run_irf_loop] shock %d/%d | metadata_fields=%s\n', ...
        shock_idx, n_shocks, join_fieldnames(fieldnames(artifact.metadata)));
end
if isfield(artifact, 'entries') && ~isempty(artifact.entries)
    fprintf('[run_irf_loop] shock %d/%d | first_entry_fields=%s\n', ...
        shock_idx, n_shocks, join_fieldnames(fieldnames(artifact.entries(1))));
end
if shock_idx > 1 && isfield(AllShockResults, 'ShockArtifacts') && ...
        numel(AllShockResults.ShockArtifacts) >= 1 && ...
        ~isempty(AllShockResults.ShockArtifacts{1})
    reference_artifact = AllShockResults.ShockArtifacts{1};
    if isstruct(reference_artifact)
        reference_fields = fieldnames(reference_artifact);
        candidate_fields = fieldnames(artifact);
        if ~isequal(reference_fields, candidate_fields)
            fprintf(2, '[run_irf_loop] shock %d top-level field mismatch vs shock 1\n', shock_idx);
            fprintf(2, '  shock 1 fields: %s\n', join_fieldnames(reference_fields));
            fprintf(2, '  shock %d fields: %s\n', shock_idx, join_fieldnames(candidate_fields));
        end
    end
end
end

function joined = join_fieldnames(fields)
if isempty(fields)
    joined = '(none)';
    return;
end
joined = strjoin(fields(:).', ', ');
end
