function runtime_config = build_test_runtime_config(overrides)
% BUILD_TEST_RUNTIME_CONFIG Shared runtime-config builder for tests/fixtures

if nargin < 1 || isempty(overrides)
    overrides = struct();
end

runtime_config = build_default_config();
runtime_config.plot_irs = false;
runtime_config.save_results = false;
runtime_config.fsolve_options = build_quiet_fsolve_options();

override_fields = fieldnames(overrides);
for i = 1:numel(override_fields)
    field_name = override_fields{i};
    runtime_config.(field_name) = overrides.(field_name);
end

runtime_config.shock_values = build_shock_values(runtime_config.shock_sizes_pct);
end
