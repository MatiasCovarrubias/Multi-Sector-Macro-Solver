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
    ModelData_IRs.metadata.run_flags = resolve_ir_run_flags(AllShockResults, config);
    ModelData_IRs.metadata.has_irfs = true;
    packaged_shocks = cell(n_shocks, 1);

    for i = 1:n_shocks
        shock_artifact = get_runtime_irf_artifact(AllShockResults, i);
        assert(is_valid_irf_artifact(shock_artifact), ...
            'build_ModelData_IRs:MissingShockArtifact', ...
            'Missing canonical IR artifact for shock %d.', i);
        packaged_shocks{i} = package_shock_artifact( ...
            shock_artifact, config.shock_values(i), sector_indices, config.ir_horizon, n_sectors_analyzed);
    end

    ModelData_IRs.shocks = vertcat(packaged_shocks{:});
end

function run_flags = build_ir_run_flags(config)
run_flags = struct( ...
    'has_1storder', logical(config.run_firstorder_irs), ...
    'has_2ndorder', logical(config.run_secondorder_irs), ...
    'has_pf', logical(config.run_pf_irs), ...
    'has_mit', false);
end

function run_flags = resolve_ir_run_flags(AllShockResults, config)
run_flags = struct();

artifacts = get_irf_artifact_list(AllShockResults);
for i = 1:numel(artifacts)
    artifact = get_runtime_irf_artifact(AllShockResults, i);
    if is_valid_irf_artifact(artifact) && isfield(artifact, 'metadata') && ...
            isstruct(artifact.metadata) && isfield(artifact.metadata, 'run_flags')
        run_flags = artifact.metadata.run_flags;
        return;
    end
end

run_flags = build_ir_run_flags(config);
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

function shock_artifact = package_shock_artifact(artifact, shock_config, sector_indices, ir_horizon, n_sectors_analyzed)
assert(numel(artifact.entries) == n_sectors_analyzed, ...
    'build_ModelData_IRs:UnexpectedEntryCount', ...
    'Shock artifact entry count does not match analyzed sector count.');
assert(numel(artifact.summary_stats.peaks.first_order) == n_sectors_analyzed, ...
    'build_ModelData_IRs:UnexpectedSummarySize', ...
    'Shock artifact summary size does not match analyzed sector count.');

metadata = build_shock_metadata(artifact.metadata, sector_indices, ir_horizon);

shock_artifact = struct( ...
    'label', get_optional_field(shock_config, 'label', ''), ...
    'value', get_optional_field(shock_config, 'value', []), ...
    'size_pct', get_optional_field(shock_config, 'size_pct', []), ...
    'sign', get_optional_field(shock_config, 'sign', []), ...
    'A_level', get_optional_field(shock_config, 'A_level', []), ...
    'description', get_optional_field(shock_config, 'description', ''), ...
    'sector_indices', sector_indices, ...
    'run_flags', metadata.run_flags, ...
    'metadata', metadata, ...
    'entries', artifact.entries, ...
    'summary_stats', artifact.summary_stats);

shock_artifact.metadata.n_entries = numel(artifact.entries);
shock_artifact.metadata.shock_config = shock_config;
if isempty(shock_artifact.metadata.shock_description)
    shock_artifact.metadata.shock_description = shock_artifact.description;
end
end

function metadata = build_shock_metadata(source_metadata, sector_indices, ir_horizon)
metadata = struct();
if isstruct(source_metadata) && ~isempty(fieldnames(source_metadata))
    metadata_fields = fieldnames(source_metadata);
    for i = 1:numel(metadata_fields)
        field_name = metadata_fields{i};
        metadata.(field_name) = source_metadata.(field_name);
    end
end

if ~isfield(metadata, 'run_flags') || isempty(metadata.run_flags)
    metadata.run_flags = struct('has_1storder', false, 'has_2ndorder', false, 'has_pf', false, 'has_mit', false);
end
metadata.sector_indices = sector_indices;
metadata.ir_horizon = ir_horizon;
end

function value = get_optional_field(s, field_name, default_value)
value = default_value;
if isstruct(s) && isfield(s, field_name)
    value = s.(field_name);
end
end
