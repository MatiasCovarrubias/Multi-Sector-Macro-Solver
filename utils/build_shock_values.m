function shock_values = build_shock_values(shock_sizes_pct)
    validate_shock_sizes_pct(shock_sizes_pct);
    shock_idx = 0;
    for s = 1:numel(shock_sizes_pct)
        pct = shock_sizes_pct(s);
        A_neg = 1 - pct/100;
        A_pos = 1/A_neg;

        shock_idx = shock_idx + 1;
        shock_values(shock_idx) = struct(... %#ok<SAGROW>
            'value', -log(A_neg), 'label', sprintf('neg%dpct', pct), ...
            'size_pct', pct, 'sign', -1, 'A_level', A_neg, ...
            'description', sprintf('-%d%% TFP (A=%.2f)', pct, A_neg));

        shock_idx = shock_idx + 1;
        shock_values(shock_idx) = struct(... %#ok<SAGROW>
            'value', log(A_neg), 'label', sprintf('pos%dpct', pct), ...
            'size_pct', pct, 'sign', +1, 'A_level', A_pos, ...
            'description', sprintf('+%.1f%% TFP (A=%.4f)', (A_pos-1)*100, A_pos));
    end

    fprintf('\nShocks (%d sizes -> %d total):\n', numel(shock_sizes_pct), numel(shock_values));
    for i = 1:numel(shock_values)
        fprintf('  %s: IR=%.4f, A=%.4f\n', shock_values(i).label, shock_values(i).value, shock_values(i).A_level);
    end
end

function validate_shock_sizes_pct(shock_sizes_pct)
    assert(isnumeric(shock_sizes_pct) && isvector(shock_sizes_pct) && ~isempty(shock_sizes_pct), ...
        'build_shock_values:InvalidShockSizes', ...
        'shock_sizes_pct must be a non-empty numeric vector.');
    assert(all(isfinite(shock_sizes_pct)), ...
        'build_shock_values:NonFiniteShockSize', ...
        'shock_sizes_pct must contain only finite values.');
    assert(all(shock_sizes_pct > 0 & shock_sizes_pct < 100), ...
        'build_shock_values:OutOfRangeShockSize', ...
        'All shock_sizes_pct entries must lie strictly between 0 and 100.');
end
