function cfg = get_ltfp_summary_config(model_type)
% GET_LTFP_SUMMARY_CONFIG Shared labor-TFP reporting field map

cfg = struct('label', '', 'model_agg_field', '', 'model_avg_field', '', ...
    'model_sectoral_field', '', 'data_agg_field', '', 'data_avg_field', '', ...
    'data_sectoral_field', '');

switch model_type
    case 'VA'
        cfg.label = 'VA';
        cfg.model_agg_field = 'corr_L_TFP_agg';
        cfg.model_avg_field = 'corr_L_TFP_sectoral_avg_vashare';
        cfg.model_sectoral_field = 'corr_L_TFP_sectoral';
        cfg.data_agg_field = 'L_TFP_agg';
        cfg.data_avg_field = 'L_TFP_sectoral_avg_vashare';
        cfg.data_sectoral_field = 'L_TFP_sectoral';
    case {'GO', 'GO_noVA'}
        cfg.label = 'GO';
        cfg.model_agg_field = 'corr_L_TFP_GO_agg';
        cfg.model_avg_field = 'corr_L_TFP_sectoral_avg_vashare';
        cfg.model_sectoral_field = 'corr_L_TFP_sectoral';
        cfg.data_agg_field = 'L_TFP_GO_agg';
        cfg.data_avg_field = 'L_TFP_sectoral_avg_vashare';
        cfg.data_sectoral_field = 'L_TFP_GO_sectoral';
end
end
