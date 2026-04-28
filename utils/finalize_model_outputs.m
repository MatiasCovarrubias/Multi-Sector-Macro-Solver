function [ModelData, ModelData_simulation, ModelData_IRs] = finalize_model_outputs( ...
        ModelData, ModelData_simulation, ModelData_IRs, flags, has_irfs, Diagnostics)
% FINALIZE_MODEL_OUTPUTS Attach run metadata and IRF summary diagnostics before saving

if nargin < 6
    Diagnostics = [];
end

ModelData.metadata.run_flags = flags;
ModelData.metadata.has_irfs = has_irfs;

ModelData_simulation.metadata.run_flags = flags;
ModelData_simulation.metadata.has_irfs = has_irfs;

if ~isempty(fieldnames(ModelData_IRs))
    if ~isfield(ModelData_IRs, 'metadata') || isempty(ModelData_IRs.metadata)
        ModelData_IRs.metadata = struct();
    end
    if ~isfield(ModelData_IRs.metadata, 'run_flags') || isempty(ModelData_IRs.metadata.run_flags)
        ModelData_IRs.metadata.run_flags = build_irf_flags_from_config(ModelData.metadata.config);
    end
    ModelData_IRs.metadata.has_irfs = has_irfs;
end

if ~isempty(Diagnostics)
    ModelData.Diagnostics = Diagnostics;
end

has_diagnostics = isfield(ModelData, 'Diagnostics') && isstruct(ModelData.Diagnostics) && ...
    ~isempty(fieldnames(ModelData.Diagnostics));
ModelData.metadata.has_diagnostics = has_diagnostics;

end

function irf_flags = build_irf_flags_from_config(config)
required_fields = {'run_firstorder_irs', 'run_secondorder_irs', 'run_pf_irs'};
for i = 1:numel(required_fields)
    assert(isfield(config, required_fields{i}), ...
        'finalize_model_outputs:MissingIRFConfigFlag', ...
        'ModelData.metadata.config is missing required IRF flag: %s', required_fields{i});
end

irf_flags = struct( ...
    'has_1storder', logical(config.run_firstorder_irs), ...
    'has_2ndorder', logical(config.run_secondorder_irs), ...
    'has_pf', logical(config.run_pf_irs), ...
    'has_mit', false);
end
