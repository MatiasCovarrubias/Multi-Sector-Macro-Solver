function ModelData = attach_simulation_statistics(ModelData, ModelData_simulation)
% ATTACH_SIMULATION_STATISTICS Copy packaged simulation views into ModelData

if ~isstruct(ModelData_simulation) || isempty(fieldnames(ModelData_simulation))
    return;
end

if ~isfield(ModelData, 'Statistics') || isempty(ModelData.Statistics)
    ModelData.Statistics = struct();
end

if isfield(ModelData_simulation, 'Shared') && isstruct(ModelData_simulation.Shared)
    shared = ModelData_simulation.Shared;
    if isfield(shared, 'Solution') && isstruct(shared.Solution) && ...
            ~isempty(fieldnames(shared.Solution))
        ModelData.Solution = shared.Solution;
    end
    if isfield(shared, 'Statistics') && isstruct(shared.Statistics) && ...
            ~isempty(fieldnames(shared.Statistics))
        shared_stat_fields = fieldnames(shared.Statistics);
        for i = 1:numel(shared_stat_fields)
            field_name = shared_stat_fields{i};
            ModelData.Statistics.(field_name) = shared.Statistics.(field_name);
        end
    end
end

method_fields = {'FirstOrder', 'SecondOrder', 'PerfectForesight', 'MITShocks'};
for i = 1:numel(method_fields)
    field_name = method_fields{i};
    if ~isfield(ModelData_simulation, field_name) || ~isstruct(ModelData_simulation.(field_name))
        continue;
    end
    if ~isfield(ModelData_simulation.(field_name), 'summary_stats') || ...
            ~isstruct(ModelData_simulation.(field_name).summary_stats) || ...
            isempty(fieldnames(ModelData_simulation.(field_name).summary_stats))
        continue;
    end

    ModelData.Statistics.(field_name) = ModelData_simulation.(field_name).summary_stats;
end

end
