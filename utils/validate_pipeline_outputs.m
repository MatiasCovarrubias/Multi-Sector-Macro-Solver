function validate_pipeline_outputs(ModelData, ModelData_simulation, ModelData_IRs, flags, has_irfs, context)
% VALIDATE_PIPELINE_OUTPUTS Validate packaged outputs against runtime flags

if nargin < 6
    context = 'validate_pipeline_outputs';
end

expects_diagnostics = has_irfs;
if expects_diagnostics
    has_diagnostics = isfield(ModelData, 'Diagnostics') && isstruct(ModelData.Diagnostics) && ...
        ~isempty(fieldnames(ModelData.Diagnostics));
    if ~has_diagnostics
        error('validate_pipeline_outputs:MissingDiagnostics', ...
            '[%s] ModelData.Diagnostics must be attached before validation.', context);
    end
end

validate_ModelData(ModelData, context);
assert_metadata_flags(ModelData.metadata, flags, has_irfs, context, 'ModelData.metadata');

if flags.has_1storder || flags.has_2ndorder || flags.has_pf || flags.has_mit
    validate_ModelData_simulation(ModelData_simulation, flags, context);
end

if has_irfs
    irf_flags = build_irf_flags_from_config(ModelData.metadata.config);
    if ~isfield(ModelData_IRs, 'metadata') || isempty(ModelData_IRs.metadata)
        ModelData_IRs.metadata = struct();
    end
    if ~isfield(ModelData_IRs.metadata, 'run_flags') || isempty(ModelData_IRs.metadata.run_flags)
        ModelData_IRs.metadata.run_flags = irf_flags;
    end
    if ~isfield(ModelData_IRs.metadata, 'has_irfs')
        ModelData_IRs.metadata.has_irfs = has_irfs;
    end
    validate_ModelData_IRs(ModelData_IRs, context);
    assert_metadata_flags(ModelData_IRs.metadata, irf_flags, has_irfs, context, 'ModelData_IRs.metadata');
end
end

function assert_metadata_flags(metadata, flags, has_irfs, context, metadata_name)
assert(islogical(metadata.has_irfs) || isnumeric(metadata.has_irfs), ...
    'validate_pipeline_outputs:InvalidHasIrfsFlag', ...
    '[%s] %s.has_irfs must be logical-like.', context, metadata_name);
assert(isequal(logical(metadata.has_irfs), logical(has_irfs)), ...
    'validate_pipeline_outputs:HasIrfsMismatch', ...
    '[%s] %s.has_irfs does not match pipeline has_irfs.', context, metadata_name);

expected_fields = fieldnames(flags);
for i = 1:numel(expected_fields)
    field_name = expected_fields{i};
    assert(isfield(metadata.run_flags, field_name) && ...
        isequal(metadata.run_flags.(field_name), flags.(field_name)), ...
        'validate_pipeline_outputs:RunFlagMismatch', ...
        '[%s] %s.run_flags.%s does not match pipeline flags.', ...
        context, metadata_name, field_name);
end
end

function irf_flags = build_irf_flags_from_config(config)
required_fields = {'run_firstorder_irs', 'run_secondorder_irs', 'run_pf_irs'};
for i = 1:numel(required_fields)
    assert(isfield(config, required_fields{i}), ...
        'validate_pipeline_outputs:MissingIRFConfigFlag', ...
        'ModelData.metadata.config is missing required IRF flag: %s', required_fields{i});
end

irf_flags = struct( ...
    'has_1storder', logical(config.run_firstorder_irs), ...
    'has_2ndorder', logical(config.run_secondorder_irs), ...
    'has_pf', logical(config.run_pf_irs), ...
    'has_mit', false);
end
