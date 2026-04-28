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
    ModelData_IRs.cir_asymmetry = build_cir_asymmetry( ...
        ModelData_IRs.shocks, config.shock_values, n_sectors_analyzed);
end

function run_flags = build_ir_run_flags(config)
run_flags = struct( ...
    'has_1storder', logical(config.run_firstorder_irs), ...
    'has_2ndorder', logical(config.run_secondorder_irs), ...
    'has_pf', logical(config.run_pf_irs), ...
    'has_mit', false);
end

function run_flags = resolve_ir_run_flags(AllShockResults, config)
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
assert(numel(artifact.summary_stats.cumulative_responses.first_order) == n_sectors_analyzed, ...
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

function asymmetry = build_cir_asymmetry(shocks, shock_values, n_sectors_analyzed)
pairs = find_positive_negative_shock_pairs(shock_values);

asymmetry = struct();
asymmetry.description = 'CIR asymmetry: negative-shock CIR divided by positive-shock CIR.';
asymmetry.interpretation = 'Values below -1 indicate negative asymmetry.';
asymmetry.n_pairs = numel(pairs);
asymmetry.rows = initialize_asymmetry_rows(numel(pairs), n_sectors_analyzed);

for i = 1:numel(pairs)
    neg_idx = pairs(i).negative_idx;
    pos_idx = pairs(i).positive_idx;

    row = asymmetry.rows(i);
    row.size_pct = pairs(i).size_pct;
    row.negative_shock_idx = neg_idx;
    row.positive_shock_idx = pos_idx;
    row.negative_label = get_optional_field(shock_values(neg_idx), 'label', '');
    row.positive_label = get_optional_field(shock_values(pos_idx), 'label', '');

    row.ratio = compute_asymmetry_ratio( ...
        shocks(neg_idx).summary_stats.cumulative_responses, ...
        shocks(pos_idx).summary_stats.cumulative_responses, ...
        n_sectors_analyzed);
    row.negative_asymmetry = compute_negative_asymmetry_flags(row.ratio);

    asymmetry.rows(i) = row;
end
end

function rows = initialize_asymmetry_rows(n_pairs, n_sectors_analyzed)
empty_method_vectors = build_method_vector_struct(NaN(1, n_sectors_analyzed));
empty_flags = build_method_vector_struct(false(1, n_sectors_analyzed));
empty_row = struct( ...
    'size_pct', NaN, ...
    'negative_shock_idx', NaN, ...
    'positive_shock_idx', NaN, ...
    'negative_label', '', ...
    'positive_label', '', ...
    'ratio', empty_method_vectors, ...
    'negative_asymmetry', empty_flags);

if n_pairs == 0
    rows = repmat(empty_row, 0, 1);
else
    rows = repmat(empty_row, n_pairs, 1);
end
end

function values = build_method_vector_struct(default_vector)
values = struct( ...
    'first_order', default_vector, ...
    'second_order', default_vector, ...
    'perfect_foresight', default_vector);
end

function pairs = find_positive_negative_shock_pairs(shock_values)
pairs = repmat(struct('size_pct', NaN, 'negative_idx', NaN, 'positive_idx', NaN), 0, 1);

for i = 1:numel(shock_values)
    if get_optional_field(shock_values(i), 'sign', 0) >= 0
        continue;
    end

    size_pct = get_optional_field(shock_values(i), 'size_pct', NaN);
    pos_idx = find_matching_positive_shock(shock_values, size_pct);
    if isnan(pos_idx)
        continue;
    end

    pairs(end + 1) = struct( ... %#ok<AGROW>
        'size_pct', size_pct, ...
        'negative_idx', i, ...
        'positive_idx', pos_idx);
end
end

function pos_idx = find_matching_positive_shock(shock_values, size_pct)
pos_idx = NaN;
if ~isfinite(size_pct)
    return;
end

for j = 1:numel(shock_values)
    candidate_size = get_optional_field(shock_values(j), 'size_pct', NaN);
    candidate_sign = get_optional_field(shock_values(j), 'sign', 0);
    if candidate_sign > 0 && isfinite(candidate_size) && abs(candidate_size - size_pct) < 1e-10
        pos_idx = j;
        return;
    end
end
end

function ratio = compute_asymmetry_ratio(negative_cirs, positive_cirs, n_sectors_analyzed)
ratio = struct();
method_fields = {'first_order', 'second_order', 'perfect_foresight'};

for i = 1:numel(method_fields)
    method = method_fields{i};
    numerator = get_method_vector(negative_cirs, method, n_sectors_analyzed);
    denominator = get_method_vector(positive_cirs, method, n_sectors_analyzed);
    ratio.(method) = safe_divide_vectors(numerator, denominator);
end
end

function flags = compute_negative_asymmetry_flags(ratio)
flags = struct( ...
    'first_order', ratio.first_order < -1, ...
    'second_order', ratio.second_order < -1, ...
    'perfect_foresight', ratio.perfect_foresight < -1);
end

function values = get_method_vector(method_struct, method, n_sectors_analyzed)
values = NaN(1, n_sectors_analyzed);
if isstruct(method_struct) && isfield(method_struct, method)
    candidate = method_struct.(method);
    if isnumeric(candidate) && numel(candidate) == n_sectors_analyzed
        values = candidate(:).';
    end
end
end

function ratio = safe_divide_vectors(numerator, denominator)
ratio = NaN(size(numerator));
valid = isfinite(numerator) & isfinite(denominator) & abs(denominator) > 1e-12;
ratio(valid) = numerator(valid) ./ denominator(valid);
end
