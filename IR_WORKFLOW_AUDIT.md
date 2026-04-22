# IR Workflow Audit

This note traces the impulse-response workflow across the two codebases:

- MATLAB/Dynare in `Multi-Sector-Macro-Solver`
- Python/JAX in `jaxecon`

The focus is the post-solution IR workflow:

- what exact experiment is being run
- whether positive and negative IRs are computed correctly
- whether any normalization or treatment could distort the final responses
- whether there are likely bugs or misleading conventions

## Executive takeaways

1. I did **not** find a clear arithmetic sign bug in either pipeline. Both explicitly compute negative and positive TFP-state shocks.
2. The biggest source of confusion is **shock semantics**: a configured shock size `p` means
   - negative leg: `A = 1 - p`
   - positive leg: `A = 1 / (1 - p)`
   so the positive leg is the log-symmetric counterpart, **not** a literal `A = 1 + p`.
3. The biggest source of cross-codebase mismatch is **conditioning**:
   - Dynare IRs are deterministic-steady-state responses.
   - Python defaults to a **GIR averaged over ergodic draws** (`use_gir = True`, `long_simulation = True`), not a deterministic-steady-state IR.
4. I did **not** find evidence that `covariance_scale` or `shock_scaling` leaks into the deterministic IR size. Those affect the stochastic-simulation side, not the IR amplitude.
5. I found two smaller issues:
   - MATLAB shock labels use `%d`, so non-integer configured sizes can be mislabeled in text output.
   - Python `analysis_hooks.py` prints sectoral row references that are off by one relative to the canonical MATLAB row map, though the actual MATLAB-to-Python extraction map appears correct.

## Shock-size semantics to keep in mind

The configured shock size is easiest to think of as the **negative leg size**. The positive leg is built to be symmetric in log space.

| Configured size | Negative TFP level | Positive TFP level actually used |
| --- | ---: | ---: |
| `12.5` | `0.8750` | `1.142857...` |
| `25` | `0.7500` | `1.333333...` |

So a configuration that looks like "plus/minus 25%" is **not** `0.75` versus `1.25`. It is `0.75` versus `1.3333`.

This convention is used in both codebases.

## Dynare workflow

### Exact experiment

The MATLAB/Dynare IRs are **not** Dynare's built-in `irf` objects. Dynare is used to solve the model, and the IR paths are then generated manually in MATLAB.

Core path:

- `runtime_config.m`
- `utils/build_shock_values.m`
- `main.m`
- `utils/run_irf_loop.m`
- `dynare/run_dynare_analysis.m`
- `dynare/process_sector_irs.m`
- `utils/process_ir_data.m`

The exact experiment is:

1. Pick a sector `j` from `config.sector_indices`.
2. Pick a shock size from `config.shock_sizes_pct`.
3. Build two experiments in `build_shock_values.m`:
   - negative: `A_neg = 1 - pct/100`
   - positive: `A_pos = 1 / A_neg`
4. Convert the experiment into `params.IRshock`.
5. Shock the **initial TFP state** of sector `j`.
6. Simulate forward with **zero future exogenous shocks**.

For perturbation IRs, `compute_perturbation_irfs()` does:

```matlab
ss_shocked = oo_.steady_state;
ss_shocked(n_sectors + sector_idx) = -params.IRshock;
shockssim_ir = zeros([opts.ir_horizon, n_sectors]);
IRS_all{ii} = simult_(M_, options_, ss_shocked, oo_.dr, shockssim_ir, order);
```

For perfect-foresight IRs, `solve_pf_irf_for_sector()` does:

```matlab
a_init = zeros(pf_irf_session.n_sectors, 1);
a_init(sector_idx) = -ir_shock;
oo_local.steady_state(idx.a(1):idx.a(2)) = a_init;
oo_local.exo_simul = pf_irf_session.zero_exo;
perfect_foresight_solver(...)
```

So the Dynare-side IR experiment is:

- a one-time perturbation to the **initial log TFP state**
- followed by a deterministic transition with **zero future shocks**

It is **not** an innovation IR to `e_t` with nonzero `shockssim_ir`.

### Positive versus negative IRs

The sign logic is internally consistent.

- `utils/build_shock_values.m` creates both legs explicitly.
- `utils/run_irf_loop.m` loops over all shock configs and assigns `params.IRshock = shock_config.value`.
- `run_dynare_analysis.m` then applies the shock through `-params.IRshock`.

This means:

- negative TFP shock: `params.IRshock > 0`, so `a_j(0) = -params.IRshock < 0`
- positive TFP shock: `params.IRshock < 0`, so `a_j(0) = -params.IRshock > 0`

I did not find a sign-flip bug in the stored IR series.

`process_sector_irs.m` does auto-detect shock sign from the TFP path, but only for summary statistics:

- if `A_ir(1) > 1`, treat as positive
- else treat as negative

That sign is then used to report positive peak magnitudes. It does **not** flip the stored IR paths themselves.

### Normalization and treatment

I did not find a hidden normalization by shock size.

Important details:

- `normalize_simulation_matrix()` only fixes orientation (`n_vars x T` versus `T x n_vars`).
- `process_ir_data.m` mixes **levels** and **log deviations**:
  - row 1 and row 13 (`A_ir`, `A_client_ir`) are converted with `exp(...)`, so they are **TFP levels**
  - most other rows are `dynare_simul - steady_state`, i.e. deviations in the model's logged representation
  - `gammaij_client_ir` is a **log share deviation**
- `covariance_scale` and `shock_scaling` operate on the stochastic shock covariance side, not on `params.IRshock`.

So there is no obvious amplitude leak from covariance scaling into the deterministic IR experiment.

### Likely issues or confusing conventions

1. **The positive-shock label is semantically misleading.**  
   A configured `25` means the positive leg is `A = 1.3333`, not `A = 1.25`. This is intentional, but easy to misread when looking at labels or filenames.

2. **Non-integer shock sizes can be mislabeled in MATLAB text output.**  
   `build_shock_values.m` uses `%d` in `sprintf(...)`, so a size like `12.5` is rendered as an integer in labels such as `neg%dpct` / `pos%dpct`.

3. **Second-order IRs are run without pruning.**  
   That is not necessarily wrong, but it is a place to be careful if long-horizon second-order responses ever look unstable.

4. **The experiment is an initial-state IR, not an innovation IR.**  
   If anyone expected a Dynare "one-period shock innovation" experiment, that is not what this code is currently doing.

## Global-solution workflow

### Exact experiment

Core path:

- `DEQN/analysis.py`
- `DEQN/analysis/GIR.py`
- `DEQN/econ_models/RbcProdNet_April2026/analysis_hooks.py`
- `DEQN/econ_models/RbcProdNet_April2026/matlab_irs.py`
- `DEQN/econ_models/RbcProdNet_April2026/plot_helpers.py`

The Python IR path does the following:

1. Discover IR shock sizes from the MATLAB IR object when available.
2. Choose the TFP states to shock (`n_sectors + sector_idx`).
3. For each sector, each shock size, and each sign:
   - build a counterfactual initial state
   - simulate forward with **zero future shocks**
   - subtract the baseline path from the counterfactual path
4. Either:
   - average those differences over many draws from a reference sample (`GIR`), or
   - compute the response from a single stochastic steady state (`IR_stoch_ss`)

The exact counterfactual-state logic in `GIR.py` is:

- negative shock: add `log(1 - shock_size)` to the log TFP level
- positive shock, default mode: add `-log(1 - shock_size)`

So the Python side uses the same log-symmetric convention as MATLAB by default.

### What the default Python IR actually is

This is the most important conceptual point for comparison with Dynare.

Current defaults in `DEQN/analysis.py` are:

- `long_simulation = True`
- `use_gir = True`
- `ergodic_price_aggregation = False`

So the default Python line is **not** a deterministic-steady-state IR. It is a **generalized impulse response averaged over draws from the nonlinear ergodic sample**.

If `use_gir = False`, the Python line becomes a **stochastic-steady-state IR**, which is still not the deterministic-steady-state object used by Dynare.

That means the standard overlay compares:

- dashed MATLAB benchmark: deterministic-steady-state IR
- solid Python line: GIR or stochastic-steady-state IR

This is not a bug, but it is a real source of mismatch.

### Positive versus negative IRs

The Python sign logic also looks correct.

`GIR_fn()` explicitly loops over:

- each requested shock size
- `shock_sign in ["pos", "neg"]`

and stores separate keys such as:

- `pos_12.5`
- `neg_12.5`
- `pos_12.5_stochss`
- `neg_12.5_stochss`

I did not find evidence that positive and negative responses are being inferred indirectly or flipped after the fact.

### Normalization and treatment

I did not find a hidden normalization by shock size on the Python side either.

Important details:

- The IR is always `counterfactual - baseline`.
- It is **not** divided by the shock magnitude.
- For aggregate variables, the model returns log-deviation objects such as `Agg. Consumption`, `Agg. GDP`, etc.
- Figures multiply both Python and MATLAB series by `100` for display as percent-like units.
- When `ergodic_price_aggregation = True`, Python re-aggregates both the DEQN series and the MATLAB overlays using fixed ergodic prices. That is intentional, not a leak.

There is also an explicit timing alignment choice:

- for most MATLAB overlays, Python skips the initial MATLAB point because MATLAB period 0 is the shocked initial state with no policy response yet
- for `Kj`, Python does **not** skip the initial point

That capital exception looks intentional and economically defensible because capital is predetermined, but it is a special case that should be remembered when checking alignment.

### Likely issues or confusing conventions

1. **The default Python object is not the same object as the Dynare benchmark.**  
   This is the main conceptual mismatch in the comparison.

2. **`analysis.py` still contains a legacy-looking `shock_size = 0.2` entry.**  
   The actual IR computation uses `ir_shock_sizes`, usually discovered from MATLAB, so that config entry is mostly noise and can mislead a reader.

3. **`analysis_hooks.py` prints row references that are off by one relative to the canonical MATLAB row map.**  
   The actual extraction logic in `matlab_irs.py` appears correct, so this looks like a documentation/console issue rather than a data bug.

## Cross-codebase comparison checklist

If the goal is an apples-to-apples comparison, the main checks are:

1. **Same shock meaning**  
   Confirm that "25" means `0.75` versus `1.3333`, not `0.75` versus `1.25`.

2. **Same conditioning object**  
   Dynare benchmark = deterministic steady state.  
   Python default = GIR over ergodic draws.

3. **Same timing alignment**  
   MATLAB overlays usually use `skip_initial = True`; capital uses a different rule.

4. **Same aggregation convention**  
   Direct aggregate rows versus fixed-price re-aggregation should not be mixed accidentally.

5. **Same response unit**  
   Raw stored IRs are mostly log deviations; plotting multiplies by `100` for display.

## Bottom line

My current read is:

- I do **not** see a clear bug in how positive versus negative shocks are computed.
- I do **not** see a hidden normalization leaking into deterministic IR size.
- The two biggest reasons results can look puzzling are:
  - the **log-symmetric shock convention**, especially on the positive side
  - the fact that the Python default IR is **not anchored at the deterministic steady state**

If you want the next step, I would isolate the comparison in this order:

1. hold fixed one sector and one shock size
2. verify the impact TFP level in both codebases
3. compare Dynare to Python with `use_gir = False`
4. if needed, add a Python deterministic-steady-state IR so the comparison is fully like-for-like
