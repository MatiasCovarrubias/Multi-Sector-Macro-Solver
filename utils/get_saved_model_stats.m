function model_stats = get_saved_model_stats(ModelData)
% GET_SAVED_MODEL_STATS Collect saved model statistics by simulation method

stats_root = struct();
if isstruct(ModelData) && isfield(ModelData, 'Statistics') && isstruct(ModelData.Statistics)
    stats_root = ModelData.Statistics;
end

model_stats = struct();
model_stats.has_ms1 = has_method_stats(stats_root, 'FirstOrder');
model_stats.has_ms2 = has_method_stats(stats_root, 'SecondOrder');
model_stats.has_msPF = has_method_stats(stats_root, 'PerfectForesight');
model_stats.has_msMIT = has_method_stats(stats_root, 'MITShocks');
model_stats.n_model_cols = model_stats.has_ms1 + model_stats.has_ms2 + ...
    model_stats.has_msPF + model_stats.has_msMIT;

if model_stats.has_ms1, model_stats.ms1 = stats_root.FirstOrder.ModelStats; else, model_stats.ms1 = struct(); end
if model_stats.has_ms2, model_stats.ms2 = stats_root.SecondOrder.ModelStats; else, model_stats.ms2 = struct(); end
if model_stats.has_msPF, model_stats.msPF = stats_root.PerfectForesight.ModelStats; else, model_stats.msPF = struct(); end
if model_stats.has_msMIT, model_stats.msMIT = stats_root.MITShocks.ModelStats; else, model_stats.msMIT = struct(); end
end

function tf = has_method_stats(stats_root, field_name)
tf = isfield(stats_root, field_name) && isstruct(stats_root.(field_name)) && ...
    isfield(stats_root.(field_name), 'ModelStats') && isstruct(stats_root.(field_name).ModelStats) && ...
    ~isempty(fieldnames(stats_root.(field_name).ModelStats));
end
