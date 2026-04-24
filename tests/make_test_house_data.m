function house_data = make_test_house_data()
%MAKE_TEST_HOUSE_DATA Creates a single test house data table row.
    house_data = table();
    house_data.house_id = 1;
    house_data.T_supply_C = 68.5;
    house_data.T_main_ukf_C = 70.0;
    house_data.flow_kg_h = 200;
    house_data.length_service_m = 10;
    house_data.x_pos_m = 5;
    house_data.timestamp = datetime(2019,11,15,12,0,0);
end