# funSNM_DH

**Service-pipe & utility-meter state estimation for district heating networks**, using an Unscented Kalman Filter (UKF) and a Particle Filter (PF).

The toolbox jointly estimates, from measured supply/return temperatures, flows and energies at the substation and at each consumer:

- per-meter **temperature offsets** (drift / calibration error of the utility meters),
- per-service-pipe **U-values** (heat loss coefficients), and
- shared parameters: **CSAC** (consumer-side ambient/soil correction) and the **main-pipe U-value**.

It is intended for diagnostics of meter drift and pipe degradation in small district-heating branches where a few high-quality reference signals are available.

---

## Quick start

**Requirements**

- MATLAB R2021a or newer (uses `arguments` blocks, `string` arrays, `timetable`).
- Statistics and Machine Learning Toolbox (for the particle filter / sampling).
- No external dependencies; all code is under `src/`.

**Run the default experiment**

1. Clone the repo and open the folder in MATLAB.
2. Edit `config/config.json` if you want to point at a different dataset or change filter settings (see [PROJECT.md](PROJECT.md) for the full schema).
3. From the repo root, run:

   ```matlab
   main
   ```

   This loads the data referenced in `config.json`, runs the UKF + PF, and writes figures and `.mat` files to `results/`.

**Other entry points**

- `main_UKF_and_PF.m` — side-by-side comparison of UKF and PF on the same data.
- `main_documented.m` — heavily commented walk-through, good for first-time readers.
- `test_drift_scenarios.m` — synthetic drift injection to check estimator behaviour.
- `test_reference_sensor.m` — sensitivity studies on adding noisy reference temperature sensors at main-pipe / CSAC junctions.
- `test_reference_sensor_direct.m` — diagnostic study comparing three reference-sensor strategies across all networks: (a) no reference, (b) reference used as the service-pipe inlet from day one, (c) reference activated only after 30 days of UKF convergence. Key finding: using the reference as a hard inlet input degrades meter-offset and U-value estimates (scenario b is ~3–4× worse than baseline; scenario c is statistically indistinguishable from baseline). The reference is best used as an **SI-traceability anchor** — for uncertainty budgeting, offline calibration of the inlet model, or a shadow-audit channel — rather than as an estimator input. See Section 5.3 of the A4.1.5 report.
- `inspect_data.m`, `inspect_csac_summary.m` — quick data / result inspection helpers.

**Run the unit tests**

```matlab
cd tests
run_all_tests
```

---

## Repository layout (one-liner)

```
config/   JSON configuration (single source of truth for a run)
data/     Input measurements and weather data (read-only)
src/      Library code: CSACs, kalman_filter, particle_filter,
          network_model, gates, data_handling, diagnostics
tests/    MATLAB unit tests (run via run_all_tests)
results/  Generated figures and .mat outputs (git-ignored)
*.m       Entry-point scripts at the repo root
```

For a module-by-module description and the full `config.json` schema, see **[PROJECT.md](PROJECT.md)**.

---

## Typical workflow

1. Drop a new measurement file into `data/` and update the path in `config/config.json`.
2. Choose a filter and tune its noise covariances in the `ukf` / `pf` sections of `config.json`.
3. Run `main` (or `main_UKF_and_PF` to compare both filters).
4. Inspect `results/` — offset and U-value trajectories, innovation diagnostics, gate-rejection statistics.
5. Re-run with adjusted gates / priors as needed.

---

## Contact

Peter Friis Østergaard — peo@teknologisk.dk
Danish Technological Institute, Metrology, Aarhus
