function Summary = build_summary_table_data(ModelData)
% BUILD_SUMMARY_TABLE_DATA Gather summary-table inputs from packaged outputs

Summary = struct();
Summary.config = ModelData.metadata.config;
Summary.params = ModelData.params;
Summary.empirical_targets = ModelData.EmpiricalTargets;
Summary.save_label = ModelData.metadata.save_label;
Summary.n_shocks = ModelData.metadata.n_shocks;
Summary.run_status = struct( ...
    'has_1storder_simul', ModelData.metadata.run_flags.has_1storder, ...
    'has_2ndorder_simul', ModelData.metadata.run_flags.has_2ndorder, ...
    'has_pf_simul', ModelData.metadata.run_flags.has_pf, ...
    'has_mit_simul', ModelData.metadata.run_flags.has_mit, ...
    'has_irfs', ModelData.metadata.has_irfs, ...
    'save_results', ModelData.metadata.config.save_results);

Summary.has_theoretical_stats = isfield(ModelData, 'Statistics') && isstruct(ModelData.Statistics) && ...
    isfield(ModelData.Statistics, 'TheoStats') && ~isempty(fieldnames(ModelData.Statistics.TheoStats));
if Summary.has_theoretical_stats
    Summary.theoretical_stats = ModelData.Statistics.TheoStats;
else
    Summary.theoretical_stats = struct();
end

Summary.model_stats = get_saved_model_stats(ModelData);

Summary.has_diagnostics = isfield(ModelData, 'Diagnostics') && ~isempty(fieldnames(ModelData.Diagnostics));
if Summary.has_diagnostics
    Summary.diagnostics = ModelData.Diagnostics;
else
    Summary.diagnostics = struct();
end

end
