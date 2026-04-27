function data = make_two_house_meter_table(timestamps)
%MAKE_TWO_HOUSE_METER_TABLE Creates a two-house meter data table for testing.
%   Both houses have one row per timestamp.

    n = numel(timestamps);
    house_ids = [ones(n, 1); 2*ones(n, 1)];
    ts = [timestamps(:); timestamps(:)];
    T_supply = repmat(70.0, 2*n, 1);

    data = table(house_ids, ts, T_supply, ...
        'VariableNames', {'house_id', 'timestamp', 'T_supply_C'});
end