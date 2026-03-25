function [calib_data, params] = load_calibration_data(params, sector_indices, model_type, shock_scaling, smooth, wds, covariance_scale)

%% Defaults
if nargin < 3, model_type = 'VA'; end
if nargin < 4, shock_scaling = struct('sectors', [], 'factor', 1.0); end
if nargin < 5 || isempty(smooth), smooth = false; end
if nargin < 6 || isempty(wds), wds = false; end
if nargin < 7 || isempty(covariance_scale), covariance_scale = 1.0; end

validate_sector_indices(sector_indices, 37, 'load_calibration_data');
assert(ismember(model_type, {'VA', 'GO', 'GO_noVA'}), 'Invalid model_type: %s', model_type);
[smooth, wds, covariance_scale, tfp_suffix] = normalize_tfp_switches( ...
    smooth, wds, covariance_scale);

%% Load raw data
load('calibration_data.mat');
tfp_process = load('TFP_process.mat');

if ~exist('beadat_37sec.mat', 'file')
    error('load_calibration_data:MissingBEAData', 'beadat_37sec.mat not found.');
end
load('beadat_37sec.mat');

n_sectors = 37;

%% Process shares and networks
conssh_data = mean(Cons47bea ./ repmat(sum(Cons47bea, 2), 1, n_sectors))';

capshbea(capshbea < 0.05) = 0.05;
capsh_data = mean(capshbea, 2);

vash_data = (mean(VA47bea ./ GO47bea))';

invnet_data = mean(invmat, 3);
invnet_data(invnet_data < 0.001) = 0.001;
invnet_data = invnet_data ./ sum(invnet_data, 1);

ionet_data = mean(IOmat47bea ./ repmat(sum(IOmat47bea), [n_sectors 1 1]), 3);
ionet_data(ionet_data < 0.001) = 0.001;
ionet_data = ionet_data ./ sum(ionet_data, 1);

%% TFP process
params.delta = mean(depratbea, 2);
params.n_sectors = n_sectors;
params.model_type = model_type;
params.smooth = smooth;
params.wds = wds;
params.covariance_scale = covariance_scale;
params.tfp_suffix = tfp_suffix;

[rho_field, sigma_field] = get_tfp_process_field_names(model_type, smooth, wds);
if ~isfield(tfp_process, rho_field) || ~isfield(tfp_process, sigma_field)
    error('load_calibration_data:MissingTFPProcessFields', ...
        'TFP_process.mat is missing required fields `%s` and/or `%s`.', rho_field, sigma_field);
end
params.rho = tfp_process.(rho_field);
params.Sigma_A_full = tfp_process.(sigma_field);
params.Sigma_A = apply_covariance_scale(params.Sigma_A_full, covariance_scale);
params.tfp_process_fields = struct('rho', rho_field, 'Sigma_A_full', sigma_field);

%% Shock scaling
if ~isempty(shock_scaling.sectors) && shock_scaling.factor ~= 1.0
    validate_sector_indices(shock_scaling.sectors, n_sectors, 'shock_scaling');
    D = eye(n_sectors);
    for i = 1:numel(shock_scaling.sectors)
        D(shock_scaling.sectors(i), shock_scaling.sectors(i)) = shock_scaling.factor;
    end
    params.Sigma_A = D * params.Sigma_A * D;
    params.shock_scaling = shock_scaling;
end

%% Store in params
params.conssh_data = conssh_data;
params.capsh_data = capsh_data;
params.vash_data = vash_data;
params.ionet_data = ionet_data;
params.invnet_data = invnet_data;

%% Aggregate consumption (NIPA)
Cons_agg = [];
if exist('Real GDP Components.xls', 'file')
    Cons_agg = readmatrix('Real GDP Components.xls', 'Sheet', 1, 'Range', 'C9:BU9')';
elseif exist('Data/Real GDP Components.xls', 'file')
    Cons_agg = readmatrix('Data/Real GDP Components.xls', 'Sheet', 1, 'Range', 'C9:BU9')';
end

%% Empirical targets
empirical_targets = compute_empirical_targets( ...
    VA_raw, EMP_raw, InvRaw, VAn, Invn, GO47bea, VA47bea, Cons_agg, TFP, TFP_GO, GOn, Cons47bea);

%% Client rankings and labels
[client_indices, ranking] = compute_client_rankings(ionet_data, sector_indices, n_sectors);

sector_label_struct = SectorLabel(sector_indices);
client_label_struct = SectorLabel(client_indices);

labels = struct();
labels.sector_indices = sector_indices;
labels.sector_labels = sector_label_struct.display;
labels.sector_labels_latex = sector_label_struct.latex;
labels.sector_labels_filename = sector_label_struct.filename;
labels.client_indices = client_indices;
labels.client_labels = client_label_struct.display;
labels.client_labels_latex = client_label_struct.latex;
labels.client_labels_filename = client_label_struct.filename;
labels.ranking = ranking;

%% Output
calib_data = struct();
calib_data.conssh_data = conssh_data;
calib_data.capsh_data = capsh_data;
calib_data.vash_data = vash_data;
calib_data.ionet_data = ionet_data;
calib_data.invnet_data = invnet_data;
calib_data.labels = labels;
calib_data.empirical_targets = empirical_targets;
calib_data.model_type = model_type;
calib_data.smooth = smooth;
calib_data.wds = wds;
calib_data.covariance_scale = covariance_scale;
calib_data.tfp_suffix = tfp_suffix;
calib_data.shock_scaling = shock_scaling;

end

function [smooth, wds, covariance_scale, tfp_suffix] = normalize_tfp_switches(smooth, wds, covariance_scale)
if ischar(smooth) || isstring(smooth)
    legacy_variant = char(string(smooth));
    assert(ismember(legacy_variant, {'baseline', 'smooth', 'wds', 'smooth_wds'}), ...
        'Invalid legacy TFP series variant: %s', legacy_variant);
    smooth = ismember(legacy_variant, {'smooth', 'smooth_wds'});
    wds = ismember(legacy_variant, {'wds', 'smooth_wds'});
end

assert(isscalar(smooth) && (islogical(smooth) || isnumeric(smooth)), ...
    'Invalid smooth flag: must be a logical scalar.');
assert(isscalar(wds) && (islogical(wds) || isnumeric(wds)), ...
    'Invalid wds flag: must be a logical scalar.');
assert(isscalar(covariance_scale) && isnumeric(covariance_scale) && isfinite(covariance_scale), ...
    'Invalid covariance_scale: must be a finite numeric scalar.');
assert(covariance_scale >= 0 && covariance_scale <= 1, ...
    'Invalid covariance_scale: must lie in [0, 1].');

smooth = logical(smooth);
wds = logical(wds);
covariance_scale = double(covariance_scale);

tfp_suffix = '';
if smooth
    tfp_suffix = [tfp_suffix, '_sm'];
end
if wds
    tfp_suffix = [tfp_suffix, '_wds'];
end
end

function [rho_field, sigma_field] = get_tfp_process_field_names(model_type, smooth, wds)
switch model_type
    case 'VA'
        base_name = '';
    case 'GO'
        base_name = '_GO';
    case 'GO_noVA'
        base_name = '_GO_noVA';
    otherwise
        error('load_calibration_data:InvalidModelType', 'Invalid model_type: %s', model_type);
end

series_suffix = '';
if smooth
    series_suffix = [series_suffix, '_sm'];
end
if wds
    series_suffix = [series_suffix, '_wds'];
end

rho_field = ['modrho', base_name, series_suffix];
sigma_field = ['modvcv', base_name, series_suffix];
end

function Sigma_scaled = apply_covariance_scale(Sigma_full, covariance_scale)
Sigma_diag = diag(diag(Sigma_full));
Sigma_scaled = Sigma_diag + covariance_scale * (Sigma_full - Sigma_diag);
end
