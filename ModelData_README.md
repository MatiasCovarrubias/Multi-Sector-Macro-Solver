# ModelData Structures

Saved objects:

- `ModelData.mat`
- `ModelData_simulation.mat`
- `ModelData_IRs.mat`

## `ModelData`

Core model object: metadata, calibration, steady state, solution views, summary statistics, and diagnostics.

| Field | Type | Description |
|-------|------|-------------|
| `metadata.date` | string | Experiment date label |
| `metadata.exp_label` | string | Experiment label |
| `metadata.save_label` | string | Combined experiment label |
| `metadata.model_type` | string | `'VA'`, `'GO'`, or `'GO_noVA'` |
| `metadata.smooth` | logical | `true` when the smoothed TFP series is selected |
| `metadata.wds` | logical | `true` when winsorized TFP growth is selected |
| `metadata.diagVCV` | logical | `true` when the shock covariance is forced diagonal |
| `metadata.tfp_suffix` | string | TFP data suffix implied by `smooth` and `wds` |
| `metadata.sector_indices` | vector | Analyzed sector indices |
| `metadata.sector_labels` | struct | Sector label data |
| `metadata.config` | struct | Full runtime configuration |
| `metadata.exp_paths` | struct | Experiment paths |
| `metadata.n_shocks` | int | Number of shock configurations |
| `metadata.run_flags` | struct | Final run-status flags: `{has_1storder, has_2ndorder, has_pf, has_mit}` |
| `metadata.has_irfs` | logical | Whether IRFs were packaged for this run |
| `metadata.has_diagnostics` | logical | Whether the compact IRF summary diagnostics are attached in `ModelData.Diagnostics` |
| `calibration` | struct | Full calibration data |
| `params` | struct | Model parameters |
| `EmpiricalTargets` | struct | Empirical target moments |
| `SteadyState.parameters` | struct | Structural parameters at SS |
| `SteadyState.policies_ss` | vector | Log steady-state policy values |
| `SteadyState.endostates_ss` | vector | Log steady-state states |
| `SteadyState.C_ss` | scalar | Aggregate consumption expenditure at SS |
| `SteadyState.L_ss` | scalar | Aggregate labor headcount at SS |
| `SteadyState.GDP_ss` | scalar | Aggregate GDP expenditure at SS |
| `SteadyState.I_ss` | scalar | Aggregate investment expenditure at SS |
| `SteadyState.K_ss` | scalar | Aggregate capital at steady state |
| `SteadyState.utility_intratemp_ss` | scalar | Intratemporal utility at SS |
| `Solution.StateSpace` | struct | State-space matrices (`A`, `B`, `C`, `D`) |
| `Solution.indices` | struct | Variable indices in state space |
| `Solution.steady_state` | vector | Dynare steady-state vector |
| `Statistics.TheoStats` | struct | First-order theoretical moments, read directly from aggregate endogenous variables when available |
| `Statistics.shocks_sd` | vector | Shock standard deviations |
| `Statistics.states_sd` | vector | State variable std devs |
| `Statistics.policies_sd` | vector | Policy variable std devs |
| `Statistics.FirstOrder` | struct | First-order simulation summary stats |
| `Statistics.SecondOrder` | struct | Second-order simulation summary stats |
| `Statistics.PerfectForesight` | struct | Perfect-foresight simulation summary stats |
| `Statistics.MITShocks` | struct | MIT-shocks simulation summary stats |
| `Diagnostics` | struct | Compact CIR-based IRF diagnostics used by the end-of-run summary table |

## `ModelData_simulation`

Full simulation object. Contains shared solution views, shock draws, and method-specific simulation paths.

| Field | Type | Description |
|-------|------|-------------|
| `metadata.save_label` | string | Experiment label |
| `metadata.exp_paths` | struct | Output folder paths |
| `metadata.run_flags` | struct | Final run-status flags |
| `metadata.has_irfs` | logical | Whether IRFs were also packaged |
| `rng_state` | struct | RNG state for reproducibility |
| `Shared.Solution` | struct | Shared state-space solution |
| `Shared.Statistics` | struct | Shared theoretical statistics |

### Shared shocks

| Field | Type | Description |
|-------|------|-------------|
| `Shocks.data` | matrix `[simul_T × n_sectors]` | Master active shock draws from `mvnrnd` |
| `Shocks.T` | int | Number of active shock periods |
| `Shocks.rng_state` | struct | RNG state at time of generation |
| `Shocks.Sigma_A` | matrix `[n × n]` | Shock covariance matrix |
| `Shocks.usage.<Method>` | struct | `{start, end}` row range used by each method |

### Method blocks

Each of `FirstOrder`, `SecondOrder`, `PerfectForesight`, and `MITShocks` uses the same schema:

| Field | Type | Description |
|-------|------|-------------|
| `<Method>.burnin_simul` | matrix `[n_vars × burn_in]` | Stored burn-in window |
| `<Method>.shocks_simul` | matrix `[n_vars × T_active]` | Active-shock window |
| `<Method>.burnout_simul` | matrix `[n_vars × burn_out]` | Stored burn-out window |
| `<Method>.variable_indices` | struct | Variable index mapping |
| `<Method>.burn_in` | int | Burn-in periods |
| `<Method>.burn_out` | int | Burn-out periods |
| `<Method>.T_active` | int | Active-shock periods |
| `<Method>.T_total` | int | `burn_in + T_active + burn_out` |
| `<Method>.summary_stats` | struct | Simulation summary stats |
| `<Method>.aggregate_series` | struct | Canonical aggregate simulation series stored as deviations from deterministic steady state |

### Window conventions

- `burnin_simul` stores the zero-shock lead-in window.
- `shocks_simul` is the active-shock window.
- `burnout_simul` stores the zero-shock transition window after the active shocks.

The common timing fields are:

- `config.simul_T`
- `config.simul_burn_in`
- `config.simul_burn_out`

### MIT special case

- `MITShocks.burn_in` is always `0`
- `MITShocks.burnin_simul` is empty
- `MITShocks.T_total = T_active + burn_out`

### Aggregate moment convention

Simulation-side aggregate moments are computed from:

1. Read aggregate `C`, `I`, `GDP`, `L`, `K`, and `utility_intratemp` directly from the Dynare aggregate endogenous-variable rows.
2. Convert `C`, `I`, `GDP`, `L`, and `K` to `log(X_t) - log(X_ss_det)`.
3. Convert `utility_intratemp` to deviation from its deterministic steady-state policy value.
4. Reconstruct aggregate `M` separately, because it is not stored as its own aggregate endogenous variable.
5. Compute moments on `shocks_simul` only.

The stored `aggregate_series` blocks follow the same canonical convention used by the aggregate moment code:

- `C`, `I`, `GDP`, `L`, and `K` are stored as `log(X_t) - log(X_ss_det)`
- `utility_intratemp` is stored as deviation from its deterministic steady-state policy value
- each series is stored by window (`burnin`, `active`, `burnout`) and as a concatenated `full` path

Sectoral value added in the moment code is also intended to use fixed steady-state prices:

- `VA_j(t) = \bar P_j (Q_j(t) - M^{out}_j(t))`

This is the convention expected by the April 2026 Python analysis layer when it recomputes benchmark moments for cross-checking.

In `ModelData.Statistics.<Method>.ModelStats` this appears as:

- `sample_window = 'shocks_simul'`
- `aggregate_definition = 'exact_logdev_to_deterministic_ss'`
- `aggregate_moments.C`
- `aggregate_moments.I`
- `aggregate_moments.GDP`
- `aggregate_moments.L`
- `aggregate_moments.K`

Legacy comparison fields are labeled explicitly, for example:

- `sigma_L_legacy_agg`
- `sigma_M_legacy_agg`
- `sigma_Domar_avg_legacy`

Variable naming inside `ModelStats` now follows:

- `Y` = sectoral primary factors
- `VA` = fixed-price sectoral value added, `VA_j(t) = \bar P_j (Q_j(t) - M^{out}_j(t))`
- `GDP` = aggregate fixed-price value added read from the aggregate policy variable

For reporting targets, the relevant sectoral volatility object is `sigma_VA_avg`.
There is no separate reported `sigma(Y)` object in the MATLAB-to-Python contract.

Additional metadata fields document this explicitly:

- `variable_convention = 'Y=primary_factors; VA=P_ss*(Q-Mout); GDP=aggregate_VA'`
- `domar_definition = 'log_fixed_price_gross_output_share_in_GDP'`
- `domar_average_weight_definition = 'legacy_normalized_gross_output_weights'`

Each `aggregate_moments.<X>` block stores:

- `mean`
- `std`
- `skewness`
- `kurtosis`

`PerfectForesight.CapitalStats` is attached when the PF simulation runs.

## `ModelData_IRs`

Impulse-response object. Contains one canonical artifact per shock.

### Top-level fields

| Field | Type | Description |
|-------|------|-------------|
| `save_label` | string | Experiment label |
| `sector_indices` | vector | Analyzed sector indices |
| `ir_horizon` | int | IRF horizon (e.g., 200) |
| `shocks` | struct array `(n_shocks)` | Per-shock IR artifacts |

### Shock artifact: `shocks(i)`

| Field | Type | Description |
|-------|------|-------------|
| `label` | string | Short label, e.g. `'neg20pct'` |
| `value` | double | IRshock value (log deviation from SS) |
| `size_pct` | numeric scalar | Shock size in percent (e.g., 12.5 or 50) |
| `sign` | int | −1 (negative) or +1 (positive) |
| `A_level` | double | TFP level: A = exp(−value) |
| `description` | string | Human-readable description |
| `sector_indices` | vector | Analyzed sectors for this shock artifact |
| `run_flags` | struct | IR methods available for this shock artifact |
| `metadata` | struct | Shock metadata |
| `entries` | struct array | Per-sector IR entries |
| `summary_stats` | struct | Per-shock summary stats |

### Per-sector entry: `shocks(i).entries(j)`

Each entry is a struct with:

| Field | Type | Description |
|-------|------|-------------|
| `sector_idx` | int | Sector number |
| `first_order` | matrix `[29 × T]` | First-order IRF |
| `second_order` | matrix `[29 × T]` | Second-order IRF |
| `perfect_foresight` | matrix `[29 × T]` | Perfect-foresight IRF |
| `sectoral_loglin` | struct | Full-sector first-order log-deviation blocks retained for compatibility/debugging |
| `sectoral_secondorder` | struct | Full-sector second-order log-deviation blocks retained for compatibility/debugging |
| `sectoral_determ` | struct | Full-sector perfect-foresight log-deviation blocks retained for compatibility/debugging |
| `aggregate_first_order` | struct | Stored aggregate first-order IR series |
| `aggregate_second_order` | struct | Stored aggregate second-order IR series |
| `aggregate_perfect_foresight` | struct | Stored aggregate perfect-foresight IR series |
| `cir` | struct | Per-sector CIR statistics for this shock |

### Per-sector CIR block: `shocks(i).entries(j).cir`

Each sector entry stores the CIR values directly, so callers do not need to recompute them from the full IR matrices.

| Field | Type | Description |
|-------|------|-------------|
| `response_variable` | string | Currently `'C_exp'`, row 2 of the IR matrix |
| `cumulative_responses.first_order` | scalar | Sum over all periods of the first-order `C_exp` IR |
| `cumulative_responses.second_order` | scalar | Sum over all periods of the second-order `C_exp` IR, or zero when unavailable |
| `cumulative_responses.perfect_foresight` | scalar | Sum over all periods of the perfect-foresight `C_exp` IR |
| `total_effect_signs.first_order` | scalar | Sign of the first-order CIR: `-1`, `0`, or `+1` |
| `total_effect_signs.second_order` | scalar | Sign of the second-order CIR: `-1`, `0`, or `+1` |
| `total_effect_signs.perfect_foresight` | scalar | Sign of the perfect-foresight CIR: `-1`, `0`, or `+1` |
| `nonlinear_amplification.pf_vs_first_order` | scalar | `CIR(PF) / CIR(first order)` for this sector and shock |
| `nonlinear_effect_class.pf_vs_first_order` | scalar | `+1` amplification, `-1` attenuation, `0` neutral, `NaN` unavailable |

`sectoral_*` blocks can contain:

- `C_all`
- `L_all`
- `Iout_all`
- `Q_all`
- `Mout_all`
- `K_all`

Each `aggregate_*` block stores:

- `C`
- `I`
- `GDP`
- `utility_intratemp`
- `C_exp`
- `I_exp`
- `GDP_exp`
- `L`
- `K`

Here `C`, `I`, `GDP`, `L`, and `K` are read directly from the Dynare aggregate endogenous variables, stored as log deviations from deterministic steady state. The `*_exp` fields are compatibility aliases for the same direct aggregate IR rows. The `sectoral_*` payloads are not used to rebuild these aggregate IRs in the current packaging path.

### Per-shock summary stats: `shocks(i).summary_stats`

`summary_stats` collects the same per-sector CIR quantities into vectors ordered like `shocks(i).sector_indices`.

| Field | Type | Description |
|-------|------|-------------|
| `response_variable` | string | Currently `'C_exp'` |
| `cumulative_responses.first_order` | row vector `[1 × n_analyzed]` | First-order CIR by analyzed sector |
| `cumulative_responses.second_order` | row vector `[1 × n_analyzed]` | Second-order CIR by analyzed sector |
| `cumulative_responses.perfect_foresight` | row vector `[1 × n_analyzed]` | Perfect-foresight CIR by analyzed sector |
| `total_effect_signs.first_order` | row vector `[1 × n_analyzed]` | Sign of first-order CIR by analyzed sector |
| `total_effect_signs.second_order` | row vector `[1 × n_analyzed]` | Sign of second-order CIR by analyzed sector |
| `total_effect_signs.perfect_foresight` | row vector `[1 × n_analyzed]` | Sign of perfect-foresight CIR by analyzed sector |
| `nonlinear_amplification.pf_vs_first_order` | row vector `[1 × n_analyzed]` | PF-vs-first-order nonlinear amplification by analyzed sector |
| `nonlinear_effect_class.pf_vs_first_order` | row vector `[1 × n_analyzed]` | Amplification/attenuation class by analyzed sector |

### Cross-shock CIR asymmetry: `ModelData_IRs.cir_asymmetry`

Positive and negative shocks are paired by `size_pct`.

| Field | Type | Description |
|-------|------|-------------|
| `description` | string | Definition of the ratio |
| `interpretation` | string | Rule of thumb for negative asymmetry |
| `n_pairs` | int | Number of matched positive/negative shock-size pairs |
| `rows(k).size_pct` | scalar | Matched shock size |
| `rows(k).negative_label` | string | Negative-shock label |
| `rows(k).positive_label` | string | Positive-shock label |
| `rows(k).ratio.first_order` | row vector `[1 × n_analyzed]` | Negative first-order CIR divided by positive first-order CIR |
| `rows(k).ratio.second_order` | row vector `[1 × n_analyzed]` | Negative second-order CIR divided by positive second-order CIR |
| `rows(k).ratio.perfect_foresight` | row vector `[1 × n_analyzed]` | Negative PF CIR divided by positive PF CIR |
| `rows(k).negative_asymmetry.<method>` | logical row vector | `true` when the corresponding ratio is below `-1` |

### Diagnostics used by the printed summary

`ModelData.Diagnostics` stores a compact view of the IR analysis:

- `upstreamness`: steady-state upstreamness computed in MATLAB from log steady-state policies converted to levels. It stores `U_M`, `U_I`, `U_simple`, `Delta_M`, and `Delta_I`; `U_M` is the primary measure used for correlations.
- `irf_sector_breakdown.rows(k)`: one row per matched shock size, with sector vectors for negative-shock amplification, positive-shock amplification, PF asymmetry, upstreamness, and the three upstreamness correlations printed in the final summary.

### IRF row map

| Row | Variable | Description |
|-----|----------|-------------|
| 1 | `A_ir` | TFP level (not deviation) |
| 2 | `C_ir` | Aggregate consumption from direct aggregate endogenous variable (log dev from SS) |
| 3 | `I_ir` | Aggregate investment from direct aggregate endogenous variable (log dev from SS) |
| 4 | `Cj_ir` | Sectoral consumption |
| 5 | `Pj_ir` | Sectoral price |
| 6 | `Ioutj_ir` | Sectoral investment output |
| 7 | `Moutj_ir` | Sectoral intermediate output |
| 8 | `Lj_ir` | Sectoral labor |
| 9 | `Ij_ir` | Sectoral investment input |
| 10 | `Mj_ir` | Sectoral intermediate input |
| 11 | `Yj_ir` | Sectoral output |
| 12 | `Qj_ir` | Sectoral Tobin's Q |
| 13 | `A_client_ir` | Client TFP level (not deviation) |
| 14 | `Cj_client_ir` | Client consumption |
| 15 | `Pj_client_ir` | Client price |
| 16 | `Ioutj_client_ir` | Client investment output |
| 17 | `Moutj_client_ir` | Client intermediate output |
| 18 | `Lj_client_ir` | Client labor |
| 19 | `Ij_client_ir` | Client investment input |
| 20 | `Mj_client_ir` | Client intermediate input |
| 21 | `Yj_client_ir` | Client output |
| 22 | `Qj_client_ir` | Client Tobin's Q |
| 23 | `Kj_ir` | Sectoral capital (log dev from SS) |
| 24 | `GDP_ir` | Aggregate GDP from direct aggregate endogenous variable (log dev from SS) |
| 25 | `Pmj_client_ir` | Client intermediate price (log dev from SS) |
| 26 | `gammaij_client_ir` | Client expenditure share deviation |
| 27 | `L_ir` | Aggregate labor headcount (log dev from SS) |
| 28 | `K_ir` | Aggregate capital with deterministic price weights (log dev from SS) |
| 29 | `utility_intratemp_ir` | Intratemporal utility (level deviation from SS) |

Rows 2, 3, and 24 are expenditure-based aggregates. Rows 27 and 28 are the model-implied aggregate `L` and `K` endogenouses, and row 29 is `utility_intratemp`.
