% OFFLINE PREPROCESSING SCRIPT: not called by main.m at runtime.
% This code loads measured TFP levels (constructed separately in Stata),
% detrends them, and estimates sector-by-sector AR(1) persistence together
% with full and diagonal shock covariance variants.
clc; clear all;

detorder = 4;
model_types = {'VA', 'GO', 'GO_noVA'};
series_variants = {'baseline', 'smooth', 'wds', 'smooth_wds'};
tfp_process = struct();

for model_idx = 1:numel(model_types)
    model_type = model_types{model_idx};
    for variant_idx = 1:numel(series_variants)
        series_variant = series_variants{variant_idx};
        csv_file = get_tfp_csv_filename(model_type, series_variant);
        tfp_levels_table = readtable(csv_file);
        tfp_levels = table2array(tfp_levels_table(:, 2:end));

        [ar1coeff, ar1resid] = estimate_tfp_process(tfp_levels, detorder);
        [rho_field, resid_field, vcv_field] = get_tfp_process_field_names(model_type, series_variant);

        tfp_process.(rho_field) = ar1coeff;
        tfp_process.(resid_field) = ar1resid;
        tfp_process.(vcv_field) = cov(ar1resid);
    end
end

save('TFP_process.mat', '-struct', 'tfp_process');

function [ar1coeff, ar1resid] = estimate_tfp_process(tfp_levels, detorder)
dim = size(tfp_levels, 2);
yrnum = size(tfp_levels, 1);
timetrend = (1:yrnum)';
detTFP = zeros(size(tfp_levels));

for sector_idx = 1:dim
    polycoeff = polyfit(timetrend, log(tfp_levels(:, sector_idx)), detorder);
    detTFP(:, sector_idx) = log(tfp_levels(:, sector_idx)) - polyval(polycoeff, timetrend);
end

ar1coeff = zeros(dim, 1);
ar1resid = zeros(yrnum - 1, dim);
for sector_idx = 1:dim
    ar1coeff(sector_idx) = arma_mlear(detTFP(:, sector_idx), 1, 0);
    ar1resid(:, sector_idx) = detTFP(2:end, sector_idx) - ...
        ar1coeff(sector_idx) * detTFP(1:end-1, sector_idx);
end
end

function csv_file = get_tfp_csv_filename(model_type, series_variant)
switch model_type
    case 'VA'
        base_file = 'TFP_37';
    case 'GO'
        base_file = 'TFP_GO_37';
    case 'GO_noVA'
        base_file = 'TFP_GO_noVA_37';
    otherwise
        error('TFPprep:InvalidModelType', 'Invalid model type: %s', model_type);
end

switch series_variant
    case 'baseline'
        csv_file = [base_file, '.csv'];
    case 'smooth'
        csv_file = [base_file, '_sm.csv'];
    case 'wds'
        csv_file = [base_file, '_wds.csv'];
    case 'smooth_wds'
        csv_file = [base_file, '_sm_wds.csv'];
    otherwise
        error('TFPprep:InvalidSeriesVariant', ...
            'Invalid TFP series variant: %s', series_variant);
end
end

function [rho_field, resid_field, vcv_field] = get_tfp_process_field_names(model_type, series_variant)
switch model_type
    case 'VA'
        model_suffix = '';
    case 'GO'
        model_suffix = '_GO';
    case 'GO_noVA'
        model_suffix = '_GO_noVA';
    otherwise
        error('TFPprep:InvalidModelType', 'Invalid model type: %s', model_type);
end

switch series_variant
    case 'baseline'
        variant_suffix = '';
    case 'smooth'
        variant_suffix = '_sm';
    case 'wds'
        variant_suffix = '_wds';
    case 'smooth_wds'
        variant_suffix = '_sm_wds';
    otherwise
        error('TFPprep:InvalidSeriesVariant', ...
            'Invalid TFP series variant: %s', series_variant);
end

field_suffix = [model_suffix, variant_suffix];
rho_field = ['modrho', field_suffix];
resid_field = ['ar1resid', field_suffix];
vcv_field = ['modvcv', field_suffix];
end
