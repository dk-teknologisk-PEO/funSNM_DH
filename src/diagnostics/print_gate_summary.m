function print_gate_summary(csac, time, house_id, can_update, innovation, dynamic_gate, startup_gate, filter_name)
    effective_gate = max(dynamic_gate, startup_gate);
    fprintf('%s | CSAC %d | %s | house %d | can_update=%d | innovation=%.3f | dynamic_gate=%.3f | effective_gate=%.3f\n', ...
        filter_name, csac, string(time), house_id, can_update, innovation, dynamic_gate, effective_gate);
end