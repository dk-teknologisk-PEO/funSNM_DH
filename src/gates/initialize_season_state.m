function season_state = initialize_season_state()
%INITIALIZE_SEASON_STATE Creates the initial heating season gate state struct.

    season_state.active = false;
    season_state.count = 0;
    season_state.start_date = NaT;
    season_state.start_timestep = NaN;
    season_state.end_timestep = NaN;
    season_state.last_gate_check_date = NaT;
    season_state.last_active_date = NaT;
    season_state.last_real_season_end_date = NaT;
    season_state.last_season_duration_days = 0;
    season_state.cooldown_active = false;
    season_state.active_days = 0;
    season_state.inactive_days = 0;
    season_state.snapshot_timestep = NaN;
end