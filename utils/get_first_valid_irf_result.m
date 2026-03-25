function [irf_result, irf_idx] = get_first_valid_irf_result(AllShockResults)
% GET_FIRST_VALID_IRF_RESULT Return the first non-empty IRF result

irf_result = [];
irf_idx = [];

irf_results = get_irf_result_list(AllShockResults);
if isempty(irf_results)
    return;
end

for i = 1:numel(irf_results)
    candidate = irf_results{i};
    if is_valid_irf_result(candidate)
        irf_result = candidate;
        irf_idx = i;
        return;
    end
end
end

function irf_results = get_irf_result_list(AllShockResults)
irf_results = {};
if ~isstruct(AllShockResults)
    return;
end
if isfield(AllShockResults, 'ShockArtifacts') && iscell(AllShockResults.ShockArtifacts) && ...
        ~isempty(AllShockResults.ShockArtifacts)
    irf_results = AllShockResults.ShockArtifacts;
end
end

function tf = is_valid_irf_result(candidate)
tf = isstruct(candidate) && ~isempty(fieldnames(candidate)) && ...
    isfield(candidate, 'summary_stats');
end
