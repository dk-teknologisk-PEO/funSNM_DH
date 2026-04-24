function data = make_test_meter_table(timestamps)
%MAKE_TEST_METER_TABLE Creates a minimal meter data table for testing.
    n = numel(timestamps);
    n_per_house = floor(n/2);
    
    house_ids = [repmat(1, n_per_house, 1); repmat(2, n - n_per_house, 1)];
    ts = [timestamps(1:n_per_house)'; timestamps(1:n-n_per_house)'];
    T_supply = repmat(70.0, numel(house_ids), 1);
    
    data = table(house_ids, ts, T_supply, ...
        'VariableNames', {'house_id', 'timestamp', 'T_supply_C'});
end