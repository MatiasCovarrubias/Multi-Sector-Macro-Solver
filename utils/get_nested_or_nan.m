function value = get_nested_or_nan(s, field_name)
% GET_NESTED_OR_NAN Return a field value or NaN when absent

if isstruct(s) && isfield(s, field_name)
    value = s.(field_name);
else
    value = NaN;
end
end
