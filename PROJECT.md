# funSNM_DH

**State estimation of district heating (DH) service pipes and utility meters using an Unscented Kalman Filter (UKF) and a Particle Filter (PF).**

The project estimates, per consumer (house), two unobservable quantities from standard utility-meter readings:

1. The temperature **offset** of the utility-meter supply-temperature sensor (sensor bias / drift).
2. The **U-value** (insulation conductance, W/m/K) of the service pipe between the cul-de-sac (CSAC) junction and the house.

In parallel it estimates shared U-values for the CSAC pipes and the main pipe, and tracks the heating season. Estimation runs network-by-network and CSAC-by-CSAC across multiple simulated DH networks (input data is produced by a separate "DH simulation" project).

---

## Repository layout

```
funSNM_DH/
├── README.md
├── .gitignore                       (ignores *.png *.mat *.csv *.json *.fig)
├── structure.txt                    (snapshot of folder tree)
│
├── main.m                           (current production entry point – UKF only)
├── main_UKF_and_PF.m                (UKF + PF side-by-side comparison)
├── main_documented.m                (older heavily-commented UKF+PF script)
├── test_drift_scenarios.m           (sensor-drift sweep: linear / step)
├── test_reference_sensor.m          (reference-sensor placement sweep)
├── test_reference_sensor_direct.m   (diagnostic: noise-free ref. sensor
│                                     bypassing T_inlet_fitted propagation;
│                                     compares 3 scenarios across all networks)
├── inspect_data.m                   (export physical inputs + diagnostic plot)
├── inspect_csac_summary.m           (per-CSAC console summary of flows / ΔT)
│
├── config/
│   └── config.json                  (single source of truth for all tuning)
│
├── data/
│   └── weather/                     (Aalborg.mat, Aalborg_air.mat,
│                                     Aalborg_soil.mat)
│
├── src/
│   ├── CSACs/                       (cul-de-sac state container + timestep loop)
│   ├── data_handling/               (import, weather, drift injection)
│   ├── diagnostics/                 (KPIs, plots, loggers, summaries)
│   ├── gates/                       (flow / season / rejection gating)
│   ├── kalman_filter/               (UKF: sigma points, predict, update, master offset)
│   ├── network_model/               (pipe physics: main / service temperature & U)
│   └── particle_filter/             (PF predict + update + init)
│
├── tests/                           (unit tests + helpers; run_all_tests.m)
│
└── results/                         (auto-generated, .gitignored)
    ├── YYYY-MM-DD_HHMM/{ukf,pf}/    (one run per timestamp)
    ├── drift_test/<scenario>/
    ├── sensor_test/, sensor_test_direct/
    └── inspection_data/
```

> External data path: `config.project.paths` points to
> `C:/Users/PEO/Documents/GitHub/DH simulation/{network_topologies,results}` —
> the `funSNM_DH` repo expects those to be present on disk; they are produced
> by the companion **DH simulation** project.

---

## Configuration — `config/config.json`

All tuning lives in one file under `project.*`. Highlights:

- **`time`** – simulation window (default `2018-01-05 → 2020-12-31`).
- **`location`** – `Aalborg` (drives the weather files under `data/weather/`).
- **`datasets`** – list of network IDs to process (default `[1,2,3,4]`).
- **`paths`** – external topology / data folders, weather folder.
- **`cutoff`** – gating thresholds: `flow_cutoff` (kg/h),
  `delta_T_gate_threshold`, `U_min/max`, `offset_min/max`, `alpha_min`.
- **`initialization.ukf`** – measurement / process noise, initial state
  uncertainty for `[offset, U]`.
- **`initialization.pf`** – `num_particles` (default 1000).
- **`initialization`** (misc) – `max_air_temperature`, `min_active_houses`,
  `max_innovation_C`, `max_nis`, `max_delta_T_change_rate`,
  `innovation_gate_N_sigma`, `hibernation_reset_threshold_hours`.
- **`csac_U_estimation`** / **`main_pipe_U_estimation`** – gain-based slow
  adaptation of shared CSAC and main-pipe U-values with warmup, correlation
  and convergence checks.
- **`main_pipe_coupling`** – blends CSAC-local main-pipe estimate with the
  global one.
- **`heating_season`** + **`heating_season_gate`** – NIS-driven season
  start/stop, lookback windows, P-inflation between seasons.
- **`master_offset`** – mean-centring step (`gamma`, `deadzone`) keeping the
  CSAC-wide mean offset near zero.
- **`simulation`** – initial-guess distributions for offsets and U-values
  (drives `rng` seeding in init code).
- **`consecutive_rejection`** – inflate `P` after repeated gate rejections.
- **`debug`** – flags to disable stability gate / mean-centring and to print
  per-update info.

---

## Entry points

### `main.m` (recommended)

Production UKF orchestration. For each network:

1. `importData` → meter data, network data, topology.
2. `initialize_all_csacs` → cell array of CSAC state structs + true trajectories.
3. Shared `U_csac` and `U_main` initialised once and kept consistent across CSACs.
4. Time-step loop: `process_csac_timestep` per CSAC, then slow adaptation
   of `U_csac` and `U_main`, season management and master-offset application.
5. Post-processing: `post_process_all_csacs`, plots, KPI summaries (CSV).

Outputs go to `results/YYYY-MM-DD_HHMM/ukf/`.

### `main_UKF_and_PF.m`

Runs the same pipeline with both UKF and PF for direct comparison. Outputs
to `results/<timestamp>/ukf/` and `results/<timestamp>/pf/`.

### `main_documented.m`

Older, more heavily commented variant kept as reference / documentation.

### Experiment scripts

- **`test_drift_scenarios.m`** – injects a known offset drift (linear or
  step) on one house and compares KPIs across scenarios. Output:
  `results/drift_test/<scenario>/`.
- **`test_reference_sensor.m`** – tests adding noisy reference temperature
  sensors at main-pipe / CSAC junctions (none / single / dual combinations).
- **`test_reference_sensor_direct.m`** – diagnostic variant that uses a
  noise-free reference and bypasses `T_inlet_fitted` propagation along the
  main pipe. Runs every network in `datasets` under three scenarios:
  (a) `no_sensor` – no reference sensors, CSACs run independently;
  (b) `reference_direct` – reference sensors active from day one;
  (c) `reference_after_convergence` – reference sensors switched on only
  after `CONVERGENCE_DAYS` (default 30) of elapsed simulation time, so each
  CSAC's UKF is allowed to converge first. Network topologies are pre-scanned
  directly from `network_<NN>.json` (no meter-data import) to determine the
  true maximum number of CSACs, so reporting works regardless of network
  size. Outputs per scenario/network: KPI plots, CSV summary, and per-CSAC
  time-series plots of meter-offset error and service-pipe U-value error
  over the full simulation span under
  `results/sensor_test_direct/net<N>/<scenario>/`.
  **Key finding (May 2026):** using the reference as a hard inlet input
  (scenario b) degrades the estimates — network-averaged offset TW-MAE grows
  from 0.110 K (baseline) to 0.378 K, with a systematic offset bias of
  −0.26 K and a doubled innovation-gate rejection rate. Scenario (c) is
  statistically indistinguishable from the baseline. The reference's value
  is therefore metrological, not statistical: it should be used as an
  SI-traceability anchor (uncertainty-budget input, offline calibration of
  the inlet model, augmented slow-bias state, or shadow audit channel)
  rather than as a UKF input. See Section 5.3 of the A4.1.5 report for the
  full discussion.
- **`inspect_data.m`** – exports physical inputs (flow, T_supply, T_air,
  T_soil, T_main estimate, ΔT) per CSAC to CSV and plots them.
- **`inspect_csac_summary.m`** – prints per-CSAC tables of pipe geometry,
  flow statistics and monthly activity.

---

## Source modules (`src/`)

### `CSACs/`

- `initialize_all_csacs.m` – builds the per-CSAC state container (meter
  data slice, ground-truth tables, UKF state per house, season state,
  drift schedule, …).
- `initialize_csac_state.m` – one-CSAC init helper.
- `process_csac_timestep.m` – the inner per-timestep step: gather active
  houses, run flow / season / innovation gates, predict & update UKF per
  house, log results.

### `kalman_filter/`

UKF for the 2-D state `x = [offset; U]`:

- `initialize_kalman_filter.m`, `generate_sigma_points.m`,
  `predict_measurement_ukf.m`, `update_house_ukf_gated.m`,
  `update_ukf_house.m`, `update_snapshot.m`.
- `apply_master_offset.m` – soft mean-centring of all offsets in a CSAC.

### `particle_filter/`

`initialize_pf_state.m`, `predict_measurement_pf.m`, `update_pf_house.m`.

### `network_model/` — pipe physics / forward model

- `get_supply_temp.m`, `get_main_temp.m` – analytical service-pipe heat
  loss (supply ↔ main-junction temperature).
- `calculate_main_pipe_temp.m`, `estimate_main_pipe_temp.m`,
  `main_temp_estimator{,_global,_global_weighted}.m` – main-pipe
  temperature profile estimators.
- `estimate_csac_U.m`, `estimate_main_pipe_U.m`,
  `compute_csac_slope_signal.m` – slow adaptation of pipe U-values.

### `gates/`

- `flow_validity_gate.m`, `check_consecutive_rejections.m`.
- `initialize_season_state.m`, `manage_heating_season.m`,
  `check_heating_season_gate.m`, `apply_season_actions.m`.

### `data_handling/`

- `importData.m`, `importAalborgData.m` – load topology + meter + network
  CSVs/MATs produced by the DH simulation project.
- `get_air_temperature.m`, `get_soil_temperature.m`, `city_position.m`,
  `get_10km_grid.m` – weather access (Aalborg).
- `generate_true_trajectories.m`, `apply_offset_drift_to_data.m` – inject
  ground-truth drift / step scenarios.
- `build_daily_T_air_max_table.m` – daily max-air lookup for season gate.
- `load_reference_sensor_data.m` – read `J_Main_<id>_s` columns as
  reference-sensor inputs.

### `diagnostics/`

- `initialize_logger.m`, `update_snapshot.m` (in kalman_filter) – per-house
  per-timestep logs.
- `compute_house_kpis.m`, `compute_and_save_network_kpis.m` – KPI table
  (TW-MAE, final error, bias, std, rejection %, U_csac error, …).
  **TW-MAE** = _Time-Weighted Mean Absolute Error_: each absolute error is
  weighted by the elapsed time Δt to the next valid sample (capped at 48 h
  to prevent long off-season gaps from dominating), then divided by the
  total weighted time. Equivalent to ordinary MAE for regularly-sampled,
  continuously-active series, but fairer in the presence of seasonal gaps
  and rejected samples.
- `plot_diagnostics.m`, `plot_csac_U_diagnostics.m`,
  `plot_kpi_bar_charts.m`, `post_process_all_csacs.m`,
  `print_csac_summary.m`, `print_full_gate_summary.m`.

---

## Tests (`tests/`)

Lightweight MATLAB unit tests using `assert_true`, `assert_false`,
`assert_near` plus `make_*` fixture builders (test config, topology,
meter tables, CSAC arrays, gate params, UKF states, …).

Run them all with:

```matlab
addpath('tests'); run_all_tests
```

Covered modules include: `apply_master_offset`, `apply_offset_drift_to_data`,
`apply_season_actions`, `build_daily_T_air_max_table`,
`check_consecutive_rejections`, `compute_csac_slope_signal`,
`compute_house_kpis`, `estimate_main_pipe_U`, `estimate_main_pipe_temp`,
`generate_true_trajectories`, `initialize_season_state`,
`manage_heating_season`, `print_csac_summary`, `update_house_ukf_gated`,
`update_shared_U_csac`, `update_snapshot`.

---

## Outputs (`results/`)

Auto-created, git-ignored. Each main run uses a `YYYY-MM-DD_HHMM`
timestamp folder; UKF and PF outputs land in `ukf/` and `pf/` subfolders:

- `Diagnostics_Network_CSAC_<net>_<csac>.png` – per-CSAC diagnostic plots.
- `<filter>_csac<idx>_filter_outputs.csv` – time series of estimated
  offset, U, P-diagonals, gate decisions per house.
- Experiment scripts add `drift_test/<scenario>/`,
  `sensor_test{,_direct}/<scenario>/`, `inspection_data/`.

---

## Dependencies

- **MATLAB** (developed against a recent release; uses `jsondecode`,
  `timetable`/`table`, `splitapply`, `datestr`, `waitbar`).
- No MATLAB toolboxes are strictly required by the core math, but
  Statistics-style helpers are used in places.
- External data from the companion **DH simulation** project (topology
  JSON/MAT files + simulated meter/network CSVs/MATs) referenced via
  `config.project.paths`.

---

## Typical workflow

1. Make sure the **DH simulation** outputs exist at the paths declared in
   `config.json`.
2. Edit `config.json` (datasets to run, season parameters, debug flags).
3. Run `main.m` (or `main_UKF_and_PF.m` for filter comparison).
4. Inspect generated plots and `*_filter_outputs.csv` under
   `results/<timestamp>/`.
5. For sensor-bias / sensor-placement studies, run
   `test_drift_scenarios.m`, `test_reference_sensor.m` or
   `test_reference_sensor_direct.m`.
6. Run `tests/run_all_tests.m` after refactors.
