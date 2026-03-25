function ModelData_IRs = build_ModelData_IRs(AllShockResults, config, save_label, sector_indices, n_shocks)
    ModelData_IRs = struct();
    ModelData_IRs.save_label = save_label;

    if ~has_irf_artifacts(AllShockResults)
        return;
    end

    n_sectors_analyzed = numel(sector_indices);

    ModelData_IRs.sector_indices = sector_indices;
    ModelData_IRs.ir_horizon = config.ir_horizon;
    ModelData_IRs.metadata = struct();
    ModelData_IRs.metadata.run_flags = build_ir_run_flags(config);
    ModelData_IRs.metadata.has_irfs = true;
    ModelData_IRs.shocks = [];

    for i = 1:n_shocks
        shock_artifact = get_runtime_irf_artifact(AllShockResults, i);
        assert(is_valid_irf_artifact(shock_artifact), ...
            'build_ModelData_IRs:MissingShockArtifact', ...
            'Missing canonical IR artifact for shock %d.', i);
        packaged_shock = package_shock_artifact( ...
            shock_artifact, config.shock_values(i), sector_indices, ...
            ModelData_IRs.metadata.run_flags, config.ir_horizon, n_sectors_analyzed);
        if i == 1
            ModelData_IRs.shocks = repmat(packaged_shock, n_shocks, 1);
        else
            assert_matching_struct_fields(ModelData_IRs.shocks(1), packaged_shock, i);
        end
        try
            ModelData_IRs.shocks(i) = packaged_shock;
        catch ME
            print_assignment_failure_debug(ModelData_IRs.shocks, packaged_shock, i, ME);
            rethrow(ME);
        end
    end
end

function run_flags = build_ir_run_flags(config)
run_flags = struct( ...
    'has_1storder', logical(config.run_firstorder_irs), ...
    'has_2ndorder', logical(config.run_secondorder_irs), ...
    'has_pf', logical(config.run_pf_irs), ...
    'has_mit', false);
end

function tf = has_irf_artifacts(AllShockResults)
tf = false;
for i = 1:numel(get_irf_artifact_list(AllShockResults))
    artifact = get_runtime_irf_artifact(AllShockResults, i);
    if is_valid_irf_artifact(artifact)
        tf = true;
        return;
    end
end
end

function artifacts = get_irf_artifact_list(AllShockResults)
artifacts = {};
if isstruct(AllShockResults)
    if isfield(AllShockResults, 'ShockArtifacts') && iscell(AllShockResults.ShockArtifacts)
        artifacts = AllShockResults.ShockArtifacts;
    end
end
end

function artifact = get_runtime_irf_artifact(AllShockResults, idx)
artifacts = get_irf_artifact_list(AllShockResults);
if numel(artifacts) < idx
    artifact = struct();
else
    artifact = artifacts{idx};
end
end

function tf = is_valid_irf_artifact(artifact)
tf = isstruct(artifact) && ~isempty(fieldnames(artifact)) && ...
    isfield(artifact, 'metadata') && isfield(artifact, 'entries') && ...
    isfield(artifact, 'summary_stats');
end

function shock_artifact = package_shock_artifact(artifact, shock_config, sector_indices, run_flags, ir_horizon, n_sectors_analyzed)
assert(numel(artifact.entries) == n_sectors_analyzed, ...
    'build_ModelData_IRs:UnexpectedEntryCount', ...
    'Shock artifact entry count does not match analyzed sector count.');
assert(numel(artifact.summary_stats.peaks.first_order) == n_sectors_analyzed, ...
    'build_ModelData_IRs:UnexpectedSummarySize', ...
    'Shock artifact summary size does not match analyzed sector count.');

shock_artifact = struct( ...
    'label', get_optional_field(shock_config, 'label', ''), ...
    'value', get_optional_field(shock_config, 'value', []), ...
    'size_pct', get_optional_field(shock_config, 'size_pct', []), ...
    'sign', get_optional_field(shock_config, 'sign', []), ...
    'A_level', get_optional_field(shock_config, 'A_level', []), ...
    'description', get_optional_field(shock_config, 'description', ''), ...
    'sector_indices', sector_indices, ...
    'run_flags', run_flags, ...
    'metadata', build_shock_metadata(artifact.metadata, run_flags, sector_indices, ir_horizon), ...
    'entries', artifact.entries, ...
    'summary_stats', artifact.summary_stats);

shock_artifact.metadata.n_entries = numel(artifact.entries);
shock_artifact.metadata.shock_config = shock_config;
if isempty(shock_artifact.metadata.shock_description)
    shock_artifact.metadata.shock_description = shock_artifact.description;
end
end

function metadata = build_shock_metadata(source_metadata, run_flags, sector_indices, ir_horizon)
metadata = struct();
if isstruct(source_metadata) && ~isempty(fieldnames(source_metadata))
    metadata_fields = fieldnames(source_metadata);
    for i = 1:numel(metadata_fields)
        field_name = metadata_fields{i};
        metadata.(field_name) = source_metadata.(field_name);
    end
end
metadata.run_flags = run_flags;
metadata.sector_indices = sector_indices;
metadata.ir_horizon = ir_horizon;
end

function value = get_optional_field(s, field_name, default_value)
value = default_value;
if isstruct(s) && isfield(s, field_name)
    value = s.(field_name);
end
end

function assert_matching_struct_fields(reference_struct, candidate_struct, idx)
reference_fields = fieldnames(reference_struct);
candidate_fields = fieldnames(candidate_struct);
if isequal(reference_fields, candidate_fields)
    return;
end

fprintf(2, '[build_ModelData_IRs] field mismatch before assigning shock %d\n', idx);
fprintf(2, '  reference fields: %s\n', join_fieldnames(reference_fields));
fprintf(2, '  candidate fields: %s\n', join_fieldnames(candidate_fields));
fprintf(2, '  only in reference: %s\n', join_fieldnames(setdiff(reference_fields, candidate_fields, 'stable')));
fprintf(2, '  only in candidate: %s\n', join_fieldnames(setdiff(candidate_fields, reference_fields, 'stable')));
error('build_ModelData_IRs:FieldMismatch', ...
    'Packaged shock %d fields do not match the first packaged shock.', idx);
end

function print_assignment_failure_debug(existing_struct, packaged_shock, idx, ME)
fprintf(2, '[build_ModelData_IRs] assignment failed for shock %d: %s\n', idx, ME.message);
if isempty(existing_struct)
    fprintf(2, '  existing struct array is empty before assignment.\n');
    return;
end

fprintf(2, '  existing fields: %s\n', join_fieldnames(fieldnames(existing_struct(1))));
fprintf(2, '  packaged fields: %s\n', join_fieldnames(fieldnames(packaged_shock)));
end

function joined = join_fieldnames(fields)
if isempty(fields)
    joined = '(none)';
    return;
end
joined = strjoin(fields(:).', ', ');
end
