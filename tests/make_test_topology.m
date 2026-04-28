function topology = make_test_topology()
%MAKE_TEST_TOPOLOGY Creates a minimal topology struct for testing.

    topology.pipe_parameters.main_pipe.insulation_W_m_K = 0.21;
    topology.pipe_parameters.csac_pipe.insulation_W_m_K = 0.16;

    topology.cul_de_sacs = struct('id', {0, 1, 2, 3}, ...
        'dist_on_main_m', {24.04, 53.59, 126.09, 165.83}, ...
        'length_m', {164.64, 174.18, 155.16, 149.63});
end