function data = make_simple_meter_table(timestamps)
%MAKE_SIMPLE_METER_TABLE Creates a single-house meter data table for testing.
%   All rows belong to house_id=1, one row per timestamp.

    n = numel(timestamps);
    data = table(ones(n, 1), timestamps(:), repmat(70.0, n, 1), ...
        'VariableNames', {'house_id', 'timestamp', 'T_supply_C'});
end